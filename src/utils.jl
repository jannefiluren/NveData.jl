
"""
    read_metadata()

Read metadata stored in excel sheet (Serier til Avrenningskart 2016.xlsx).
"""
function read_metadata()

    file_name = joinpath(dirname(@__FILE__), "../raw/Serier til Avrenningskart 2016.xlsx")

    df_meta = readxlsheet(DataFrame, file_name, "Ark1")

end


"""
    read_dbk_ind()

Read file with information about drainage basin key and indices indicating
the location of the watersheds in the senorge grid (SeNorge.txt).
"""
function read_dbk_ind()

    file_name = joinpath(dirname(@__FILE__), "../raw/SeNorge.txt")

    lines = readlines(file_name)

    dbk_ind = OrderedDict()

    for line in lines
        str_splitted = split(line)
        dbk_ind[str_splitted[2]] = map(x -> parse(Int64,x) + 1, str_splitted[3:end])
    end

    return dbk_ind

end


"""
    link_elev_ind(dbk_ind, dbk_elev, elev_breaks = 0:200:4000)

For each catchment, identify which grid cells belong to the different
elevation bands.
"""
function link_elev_ind(dbk_ind, dbk_elev, elev_breaks = 0:200:4000)

    elev_ind = OrderedDict()

    for (dbk, ind) in dbk_ind

        tmp = OrderedDict()

        for ielev in 1:length(elev_breaks)-1

            ikeep = elev_breaks[ielev] .<= dbk_elev[dbk] .< elev_breaks[ielev+1]

            if any(ikeep)
                tmp[ielev] = ind[ikeep]
            end

        end

        elev_ind[dbk] = tmp

    end

    return elev_ind

end


"""
    link_stat_dbk(df_meta)

Link station name with drainage_basin_key.
"""
function link_stat_dbk(df_meta)

    regine_area = [Int(regine_area) for regine_area in df_meta[:regine_area]]
    main_no = [Int(main_no) for main_no in df_meta[:main_no]]
    drainage_basin_key = [Int(drainage_basin_key) for drainage_basin_key in df_meta[:drainage_basin_key]]

    stat_dbk = OrderedDict()

    for irow = 1:size(df_meta, 1)

        stat_dbk[ "$(regine_area[irow]).$(main_no[irow])"] = "$(drainage_basin_key[irow])"

    end

    return stat_dbk

end


"""
    link_dbk_elev(dbk_ind, dem_vec)

Link drainage basin key and elevations. Remove watersheds lacking elevation data.
"""
function link_dbk_elev(dbk_ind, dem_vec)

    dbk_elev = OrderedDict()

    for (dbk, ind) in dbk_ind

        if any(isnan.(dem_vec[ind]))
            delete!(dbk_ind, dbk)
        else
            dbk_elev[dbk] = dem_vec[ind]
        end

    end

    return dbk_elev, dbk_ind

end


"""
    read_esri_raster(file_name)

Read esri raster file.
"""
function read_esri_raster(file_name)

    lines = readlines(file_name)

    res = Dict()

    for irow in 1:6
        str_splitted = split(lines[irow])
        res[str_splitted[1]] = parse(Int64, str_splitted[2])
    end

    raster = readdlm(file_name, ' ', skipstart=6)

    res["data"] = convert(Array{Float64,2}, raster[:, 1:end-1])

    return res

end


"""
    raster2vec(raster)

Convert a esri raster to a vector and replace missing data with not-a-number.
"""
function raster2vec(raster)

    data = raster["data"]
    data[data .== raster["NODATA_value"]] = NaN
    data = transpose(data)
    vec = data[:]

    return vec

end


"""
    read_bil(file_name)

Read binary file with senorge data.
"""
function read_bil(file_name)

    data = zeros(Float64, 1550*1195)
    
    read_bil!(data, file_name)

    return data

end

function read_bil!(data, file_name)

    @assert length(data) == 1550*1195

    f = open(file_name, "r") 
    i = 1
    while ~eof(f)
        data[i] = read(f, Int16)
        i += 1
    end
    close(f)

end


"""
    aggregate_elevations(elev_ind, data_vec, scale_func; agg_func = mean)

Average data over elevation bands.
"""
function aggregate_elevations(elev_ind, data_vec, scale_func; agg_func = mean)

    res = OrderedDict()

    for reg_nr in keys(elev_ind)

        res_single = zeros(length(elev_ind[reg_nr]))

        i = 1

        for ind in values(elev_ind[reg_nr])
            res_single[i] = agg_func(scale_func(data_vec[ind]))
            i += 1
        end

        res[reg_nr] = transpose(res_single)

    end

    return res

end

function aggregate_elevations!(res, elev_ind, data_vec, scale_func; agg_func = mean)

    for reg_nr in keys(elev_ind)

        res_single = zeros(length(elev_ind[reg_nr]))

        i = 1

        for ind in values(elev_ind[reg_nr])
            res_single[i] = agg_func(scale_func(data_vec[ind]))
            i += 1
        end

        res[reg_nr] = vcat(res[reg_nr], transpose(res_single))

    end

end


"""
    read_senorge_data(time_vec, met_var, elev_ind; senorge = "V2.0", scale_func = x -> x)

Read senorge meteorological data.
"""
function read_senorge_data(time_vec, met_var, elev_ind; senorge = "V2.0", scale_func = x -> x)

    data_out = []
    data_grid = []
    time_out = DateTime[]

    for i = 1:length(time_vec)

        date_str1 = Dates.format(time_vec[i], "yyyy")
        date_str2 = Dates.format(time_vec[i], "yyyy_mm_dd")

        if senorge == "V2.0"

            path_met = "//hdata/grid/metdata/met_obs_v2.0"
            
            file_name = joinpath(path_met, "$met_var", "$date_str1", "$(met_var)_$date_str2.bil")    

        end

        print("Processing data for $date_str2 \n")
        
        try
            
            if i == 1

                data_grid = read_bil(file_name)
                data_out = aggregate_elevations(elev_ind, data_grid, scale_func)

            else
                                
                read_bil!(data_grid, file_name)
                aggregate_elevations!(data_out, elev_ind, data_grid, scale_func)
                
            end

            push!(time_out, time_vec[i])

        catch

            warn("No $met_var data available from $(time_vec[i]) and onwards")
            break

        end

    end

    return data_out, time_out

end


"""
    write_met_data(stat_list, save_folder, file_desc, time_vec, met_data, stat_dbk)

Write meteorological data to files
"""
function write_met_data(stat_list, save_folder, file_desc, time_vec, met_data, stat_dbk)

    for stat_name in stat_list

        try

            dbk = stat_dbk[stat_name]

            file_path = joinpath(save_folder, "$(replace(stat_name, ".", "_"))_data")

            mkpath(file_path)

            file_name = joinpath(file_path, "$file_desc.txt")
            data = hcat(Dates.format(time_vec, "yyyy-mm-dd HH:MM"), round(met_data[dbk], 2))

            f = open(file_name, "w")
            writedlm(f, data, ";", quotes=false)
            close(f)

        catch

            warn("Unable to write $file_desc data for station $stat_name")

        end

    end

end


"""
    write_update(stat_list, save_folder, file_desc, time_old, met_old, time_new, met_new, stat_dbk)

Combine old and new senorge data and write to file.
"""
function write_update(stat_list, save_folder, file_desc, time_old, met_old, time_new, met_new, stat_dbk)

    for stat_name in stat_list

        try

            dbk = stat_dbk[stat_name]

            file_path = joinpath(save_folder, "$(replace(stat_name, ".", "_"))_data")

            mkpath(file_path)

            met_data = vcat(met_old[stat_name], met_new[dbk])
            time_vec = vcat(time_old, time_new)

            file_name = joinpath(file_path, "$file_desc.txt")
            data = hcat(Dates.format(time_vec, "yyyy-mm-dd HH:MM"), round(met_data, 2))

            f = open(file_name, "w")
            writedlm(f, data, ";", quotes=false)
            close(f)

        catch

            warn("Unable to write $file_desc data for station $stat_name")

        end

    end

end


"""
    metadata_elevbands(stat_list, save_folder, elev_ind, stat_dbk)

Compute average coverage of the different landuse classes, average elevation
and total area for the individual elevation zones.
"""
function metadata_elevbands(stat_list, save_folder, elev_ind, stat_dbk)

    # Process grids

    raster_names = ["lus_unclassified", "lus_glacier", "lus_agriculture", "lus_bog",
                   "lus_lake", "lus_forest", "lus_bare_mountain", "lus_urban",
                   "elevation", "area"]
    
    agg_type = [mean, mean, mean, mean, mean, mean, mean, mean, mean, sum]

    res_dict = Dict(key => DataFrame() for key in keys(elev_ind))

    for i in eachindex(raster_names)

        file_name = joinpath(Pkg.dir("NveData"), "raw/$(raster_names[i]).asc")

        data_raster = read_esri_raster(file_name)

        data_vec = raster2vec(data_raster)

        data_out = aggregate_elevations(elev_ind, data_vec, x -> x; agg_func = agg_type[i])

        for (dbk, values) in data_out

            res_dict[dbk][parse(raster_names[i])] = values[:, 1]

        end

    end

    # Write to files

    for stat_name in stat_list

        try

            dbk = stat_dbk[stat_name]

            file_path = joinpath(save_folder, "$(replace(stat_name, ".", "_"))_data")

            mkpath(file_path)

            file_name = joinpath(file_path, "metadata.txt")
            
            writetable(file_name, res_dict[dbk], separator = ';')
            
        catch

            warn("Unable to write metadata for elevation bands for station $stat_name")

        end

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

        # Load data for one station

        istat = find(stat_all .== stat_name)
        istat = istat[1]

        file_name = "$(regine_area[istat]).$(main_no[istat]).$(point_no[istat]).$(param_key[istat]).$(version_no_end[istat])"

        data = readdlm(joinpath(path_runoff, file_name), skipstart = 2) 

        # Match desired time period using data frames and joins

        df_all = DataFrame(Time = map(DateTime, data[:,1]), Runoff = data[:,2])

        df_subset = DataFrame(Time = time_vec)

        df_final = join(df_subset, df_all, on = :Time, kind = :left)

        # Save data to file

        try

            time_save = df_final[:Time]

            runoff_save = df_final[:Runoff]

            runoff_save = (runoff_save * 86400 * 1000) / (area_total[istat] * 1e6)  # Convert from m3/s to mm/day

            runoff_save[isna(runoff_save)] = -999.0
            runoff_save[runoff_save .< 0.0] = -999.0

            data = hcat(Dates.format(time_save, "yyyy-mm-dd HH:MM"), round(runoff_save, 2))

            file_path = joinpath(save_folder, "$(replace(stat_name, ".", "_"))_data")
            file_name = joinpath(file_path, "Q_obs.txt")

            mkpath(file_path)

            f = open(file_name, "w")
            writedlm(f, data, ";", quotes=false)
            close(f)
            
        catch

            warn("Unable to write runoff data for station $stat_name")

        end

    end

end


"""
     read_old_data(stat_list, save_folder, met_var; crop = 10)
    
Read existing data.
"""
function read_old_data(stat_list, save_folder, met_var; crop = 10)

    data = OrderedDict()
    
    time = []

    for stat_name in stat_list

        file_path = joinpath(save_folder, "$(replace(stat_name, ".", "_"))_data", "$met_var.txt")

        str  = readline(file_path)
        nsep = length(matchall(r";", str))
        tmp  = CSV.read(file_path, delim = ";", header = false,
                        dateformat="yyyy-mm-dd HH:MM", nullable = false, types = vcat(DateTime, repmat([Float64], nsep)))

        time = Array(tmp[1:end-crop, 1])
        tmp  = Array(tmp[1:end-crop, 2:end])

        data[stat_name] = tmp

    end

    return data, time

end


"""
    find_last_senorge(senorge)

Find last available senorge data.
"""
function find_last_senorge(senorge)

    time_start = floor(now(), Dates.Day)
    time_stop = time_start + Dates.Day(10)

    time_final = []

    for time_check = time_start:time_stop

        date_str1 = Dates.format(time_check, "yyyy")
        date_str2 = Dates.format(time_check, "yyyy_mm_dd")

        if senorge == "V2.0"
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
