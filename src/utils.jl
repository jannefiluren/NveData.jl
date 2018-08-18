
"""
    read_metadata()

Read metadata stored in excel sheet (Serier til Avrenningskart 2016.xlsx)
"""
function read_metadata()

    file_name = joinpath(dirname(@__FILE__), "../raw/metadata.csv")

    df_meta = CSV.read(file_name; delim = ';', missingstring = "NULL")

    filter!(row -> !ismissing(row[:drainage_basin_key]), df_meta)

    regine_area = df_meta[:regine_area]
    main_no = df_meta[:main_no]

    df_meta[:regine_main] = ["$(regine_area[i]).$(main_no[i])" for i in eachindex(regine_area)]

    return df_meta
    
end


"""
    subset_metadata(df_meta, stat_list)

Select metadata using station selection
"""
function subset_metadata(df_meta, stat_list)

    df_meta = filter(row -> row[:regine_main] in stat_list, df_meta)
    
    return df_meta

end


"""
    read_dbk_ind()

Read file with information about drainage basin key and indices indicating
the location of the watersheds in the senorge grid (SeNorge.txt)
"""
function read_dbk_ind()

    file_name = joinpath(dirname(@__FILE__), "../raw/SeNorge.txt")

    lines = readlines(file_name)

    dbk_ind = OrderedDict{Float64, Array{Int64, 1}}()

    for line in lines        
        str_splitted = split(line)
        key = parse(Float64, str_splitted[2])
        values = map(x -> parse(Int64,x) + 1, str_splitted[3:end])
        dbk_ind[key] = values
    end

    return dbk_ind

end


"""
    subset_dbk_ind(dict_in, df_meta)

Select dictionary linking senorge indices with drainage basin key using metadata table
"""
function subset_dbk_ind(dict_in, df_meta)

    dbk_sel = df_meta[:drainage_basin_key]

    dict_out = OrderedDict{Float64, Array{Int, 1}}()

    for dbk in dbk_sel

        dict_out[dbk] = dict_in[dbk]

    end

    return dict_out

end


"""
    watershed_info(dbk_ind)

Collect information (elevation, landuse, ...) about gridcells in a watershed.
"""
function watershed_info(dbk_ind, elev_breaks = 0:200:4000)

    # Initilize dictonary

    df_all = OrderedDict{Float64, DataFrame}()

    for (dbk, ind) in dbk_ind
        df_all[dbk] = DataFrame(senorge_ind = ind, area = 1.0)
    end

    # Add landuse and elevation information to dataframes

    raster_names = ["lus_unclassified", "lus_glacier", "lus_agriculture", "lus_bog",
    "lus_lake", "lus_forest", "lus_bare_mountain", "lus_urban", "elevation"]

    for raster in raster_names

        file_name = joinpath(dirname(@__FILE__), "../raw/$(raster).asc")
        
        data_raster = read_esri_raster(file_name)
        
        data_vec = raster2vec(data_raster)

        for (dbk, senorge_ind) in dbk_ind

            df_all[dbk][Meta.parse(raster)] = data_vec[senorge_ind]

        end

    end

    # Remove watershed if elevation data is missing

    for (dbk, df_wsh) in df_all
        
        if any(isnan.(df_wsh[:elevation]))
            delete!(df_all, dbk)
            info("Removed watershed with drainage basin key $dbk because elevation data is missing!")
        end

    end

    # Add indicies for elevation bands to dataframes

    for (dbk, df_wsh) in df_all

        elevation = df_wsh[:elevation]

        elev_ind = similar(elevation)

        elev_band = 1
        
        for ielev in 1:length(elev_breaks)-1
            
            ind = findall(elev_breaks[ielev] .<= elevation .< elev_breaks[ielev+1])
            
            if !isempty(ind)
                elev_ind[ind] .= elev_band
                elev_band += 1
            end
            
        end
        
        df_all[dbk][:elev_ind] = elev_ind
        
    end

    return df_all

end


"""
    aggregate_elevations!(metdata_elev, itime, df_geo, senorge_raw, scale_func; agg_func = mean)

Average data over elevation bands.
"""
function aggregate_elevations!(metdata_elev, itime, df_geo, senorge_raw, scale_func; agg_func = mean) 
    
    for dbk in keys(df_geo)

        senorge_ind = convert(Array{Int64, 1}, df_geo[dbk][:senorge_ind])

        elev_ind = convert(Array{Int64, 1}, df_geo[dbk][:elev_ind])

        elev_iter = sort(unique(elev_ind))

        for ielev in elev_iter

            iagg = senorge_ind[elev_ind .== ielev]

            metdata_elev[dbk][itime, ielev] = agg_func(scale_func(senorge_raw[iagg]))

        end

    end

    return nothing

end


"""
    read_senorge_data(time_vec, met_var, df_geo; senorge = "V2.0", scale_func = x -> x)

Read senorge meteorological data.
"""
function read_senorge_data(time_vec, met_var, df_geo; senorge = "V2.0", scale_func = x -> x)

    # Allocate memory for outputs

    metdata_elev = OrderedDict{Float64, Array{Float64, 2}}()
    
    for (dbk, df_tmp) in df_geo

        nrows = length(time_vec)
        ncols = length(unique(df_tmp[:elev_ind]))

        metdata_elev[dbk] = zeros(nrows, ncols)

    end

    senorge_raw = zeros(Float64, 1550*1195)

    # Loop over time

    for itime = 1:length(time_vec)

        date_str1 = Dates.format(time_vec[itime], "yyyy")
        date_str2 = Dates.format(time_vec[itime], "yyyy_mm_dd")

        if senorge == "V2.0"

            path_met = "//hdata/grid/metdata/met_obs_v2.0"
            
            file_name = joinpath(path_met, "$met_var", "$date_str1", "$(met_var)_$date_str2.bil")    

        end

        print("Processing data for $date_str2 \n")
        
        try

            read_bil!(senorge_raw, file_name)

            aggregate_elevations!(metdata_elev, itime, df_geo, senorge_raw, scale_func)

        catch

            error("No $met_var data available from $(time_vec[i]) and onwards")

        end

    end

    return metdata_elev, time_vec

end


"""
    write_met_data(save_folder, file_desc, time_vec, met_data, df_meta)

Write meteorological data to files
"""
function write_met_data(save_folder, file_desc, time_vec, met_data, df_meta)

    for row in eachrow(df_meta)

        try

            dbk = row[:drainage_basin_key]

            regine_main = row[:regine_main]

            file_path = joinpath(save_folder, "$(replace(regine_main, "." => "_"))_data")

            mkpath(file_path)

            file_name = joinpath(file_path, "$file_desc.txt")
            data = hcat(Dates.format.(time_vec, "yyyy-mm-dd HH:MM"), round.(met_data[dbk], digits = 2))

            f = open(file_name, "w")
            writedlm(f, data, ";", quotes=false)
            close(f)

        catch

            warn("Unable to write $file_desc data for station $stat_name")

        end

    end

end


"""
    write_update(save_folder, file_desc, time_old, met_old, time_new, met_new, df_meta)

Combine old and new senorge data and write to file
"""
function write_update(save_folder, file_desc, time_old, met_old, time_new, met_new, df_meta)

    for row in eachrow(df_meta)

        try

            regine_main = row[:regine_main]

            dbk = row[:drainage_basin_key]

            file_path = joinpath(save_folder, "$(replace(regine_main, "." => "_"))_data")

            mkpath(file_path)

            met_data = vcat(met_old[regine_main], met_new[dbk])
            time_vec = vcat(time_old, time_new)

            file_name = joinpath(file_path, "$file_desc.txt")
            data = hcat(Dates.format(time_vec, "yyyy-mm-dd HH:MM"), round.(met_data, 2))

            f = open(file_name, "w")
            writedlm(f, data, ";", quotes=false)
            close(f)

        catch

            warn("Unable to write $file_desc data for station $regine_main")

        end

    end

end


"""
    metadata_elevbands(save_folder, df_geo, df_meta)

Compute average coverage of the different landuse classes, average elevation
and total area for the individual elevation zones
"""
function metadata_elevbands(save_folder, df_geo, df_meta)

    regine_main = convert(Array, df_meta[:regine_main])

    for (dbk, df_tmp) in df_geo

        # Aggregate data for elevation bands
        
        df_agg = aggregate(df_tmp, :elev_ind, [mean, sum])

        colnames = [Symbol(replace(string(name), "Statistics." => "")) for name in names(df_agg)]

        names!(df_agg, colnames)
        
        keep = [:elev_ind, :lus_unclassified_mean, :lus_glacier_mean, :lus_agriculture_mean,
                :lus_bog_mean, :lus_lake_mean, :lus_forest_mean, :lus_bare_mountain_mean,
                :lus_urban_mean, :elevation_mean, :area_sum]

        df_save = df_agg[keep]
        
        # Save data to files
        
        irow = findall(df_meta[:drainage_basin_key] .== dbk)
        
        regine_main = df_meta[:regine_main][irow][1]

        file_path = joinpath(save_folder, "$(replace(regine_main, "." => "_"))_data")

        mkpath(file_path)
        
        file_name = joinpath(file_path, "metadata.txt")
                    
        CSV.write(file_name, df_save; delim = ';')

    end

end


"""
    process_runoff_data(stat_list, time_vec, df_meta, path_runoff, save_folder)

Process runoff data obtained using shell scripts.
"""
function process_runoff_data(stat_list, time_vec, df_meta, path_runoff, save_folder)

    # File names for all available stations

    regine_area = convert(Array{Int64}, df_meta[:regine_area])
    main_no = convert(Array{Int64}, df_meta[:main_no])
    point_no = convert(Array{Int64}, df_meta[:point_no])
    param_key = convert(Array{Int64}, df_meta[:param_key])
    version_no_end = convert(Array{Int64}, df_meta[:version_no_end])

    stat_all = ["$(regine_area[i]).$(main_no[i])" for i in eachindex(regine_area)]

    # Catchment area

    area_total = df_meta[:area_total]

    # Loop over stations in input list

    for stat_name in stat_list

        @info "Save runoff data for station: $(stat_name)"

        # Load data for one station

        istat = findall(stat_all .== stat_name)
        istat = istat[1]

        file_name = "$(regine_area[istat]).$(main_no[istat]).$(point_no[istat]).$(param_key[istat]).$(version_no_end[istat])"

        data = readdlm(joinpath(path_runoff, file_name), skipstart = 2) 

        # Match desired time period using data frames and joins

        df_runoff = DataFrame(Time = map(DateTime, data[:,1]), Runoff = data[:,2])

        df_ref = DataFrame(Time = time_vec)

        df_final = join(df_ref, df_runoff, on = :Time, kind = :left)

	    sort!(df_final, [:Time])

        # Save data to file

        try

            time_save = df_final[:Time]

            runoff_save = df_final[:Runoff]

            area = tryparse(Float64, replace(area_total[istat], "," => "."))

            runoff_save = (runoff_save * 86400 * 1000) / (area * 1e6)  # Convert from m3/s to mm/day

            runoff_save = convert(Array{Union{Missing, Float64}, 1}, runoff_save)

            runoff_save[ismissing.(runoff_save)] .= -999.0
            
            data = hcat(Dates.format.(time_save, "yyyy-mm-dd HH:MM"), round.(runoff_save; digits = 2))

            file_path = joinpath(save_folder, "$(replace(stat_name, "." => "_"))_data")
            file_name = joinpath(file_path, "runoff.txt")

            mkpath(file_path)

            f = open(file_name, "w")
            writedlm(f, data, ";", quotes=false)
            close(f)
            
        catch

            @warn "Unable to write runoff data for station $stat_name"

        end

    end

end


"""
     read_old_data(save_folder, met_var; crop = 10)
    
Read existing data.
"""
function read_old_data(save_folder, met_var; crop = 10)

    stat_list = readdir(save_folder)
    stat_list = [replace(x, "_data" => "") for x in stat_list]
    stat_list = [replace(x, "_" => ".") for x in stat_list]

    data = OrderedDict()
    
    time = []

    for stat_name in stat_list

        file_path = joinpath(save_folder, "$(replace(stat_name, "." => "_"))_data", "$met_var.txt")
        
        str  = readline(file_path)
        nsep = length(collect((m.match for m = eachmatch(r";", str))))
        tmp  = CSV.read(file_path, delim = ";", header = false,
                        dateformat="yyyy-mm-dd HH:MM", allowmissing=:none, 
                        types = vcat(DateTime, repeat([Float64], nsep)))

        time = tmp[1:end-crop, 1]
        tmp  = convert(Array{Float64}, tmp[1:end-crop, 2:end])

        data[stat_name] = tmp

    end

    return data, time

end


"""
    find_last_senorge(senorge_version)

Find last available senorge data.
"""
function find_last_senorge(senorge_version)

    time_start = floor(now(), Dates.Day)
    time_stop = time_start + Dates.Day(10)

    time_final = []

    for time_check = time_start:time_stop

        date_str1 = Dates.format(time_check, "yyyy")
        date_str2 = Dates.format(time_check, "yyyy_mm_dd")

        if senorge_version == "V2.0"
            path_met = "//hdata/grid/metdata/met_obs_v2.0"
            file_name = joinpath(path_met, "rr", "$date_str1", "rr_$date_str2.bil")
        end

        if isfile(file_name)
            time_final = time_check
        else
            break
        end
        
    end

    return time_final

end
