
"""
    update_dataset(stat_list, save_folder)

Update an existing initial dataset.
"""
function update_dataset(stat_list, save_folder)

    # Path to runoff data

    path_runoff = "//hdata/fou/Avrenningskart/Data/runoff_all/"
    
    # Read metadata from excel sheet

    df_meta = read_metadata()

    # Link drainage basin key with watershed indices

    dbk_ind = read_dbk_ind()

    # Read digital elevation model

    file_name = joinpath(Pkg.dir("NveData"), "raw/elevation.asc")

    dem_raster = read_esri_raster(file_name)

    dem_vec = raster2vec(dem_raster)

    # Link drainage basin key with elevations

    dbk_elev, dbk_ind = link_dbk_elev(dbk_ind, dem_vec)

    # For each catchment, link indicies of grid cells to individual elevation bands

    elev_ind = link_elev_ind(dbk_ind, dbk_elev)

    # Link station name and drainage basin key

    stat_dbk = link_stat_dbk(df_meta)

    # Read existing data

    tair_old, time_old = read_old_data(stat_list, save_folder, "Tair")

    prec_old, time_old = read_old_data(stat_list, save_folder, "Prec")

    # Find last available senorge data

    time_first = time_old[end] + Dates.Day(1)
    time_last = find_last_senorge("V2.0")

    time_vec = time_first:time_last

    # Read new meteorological data

    met_var = "tm"

    tair_new, time_new = read_senorge_data(time_vec, met_var, elev_ind; scale_func = x -> (x-2732.0)/10.0)

    met_var = "rr"

    prec_new, time_new = read_senorge_data(time_vec, met_var, elev_ind; scale_func = x -> x/10.0)

    # Write data to files

    write_update(stat_list, save_folder, "Tair", time_old, tair_old, time_new, tair_new, stat_dbk)

    write_update(stat_list, save_folder, "Prec", time_old, prec_old, time_new, prec_new, stat_dbk)

    # Process runoff data

    time_first = time_old[1]
    
    time_vec = time_first:time_last

    process_runoff_data(stat_list, time_vec, df_meta, path_runoff, save_folder)

end

