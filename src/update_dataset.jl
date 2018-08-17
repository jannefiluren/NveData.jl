
"""
    update_dataset(stat_list, save_folder; path_runoff)

Update an existing initial dataset.
"""
function update_dataset(stat_list, save_folder;
    path_runoff = "/home/jmg/flood_forecasting/runoff_observed")

    # Read metadata from excel sheet

    df_meta = read_metadata()

    df_meta = subset_metadata(df_meta, stat_list)

    # Dictonary that links drainage basin keys to senorge indices

    dbk_ind = read_dbk_ind()

    dbk_ind = subset_dbk_ind(dbk_ind, df_meta)

    # Dictonary that links drainage basin keys to dataframes with
    # geographic information (landuse, elevation...)

    df_geo = watershed_info(dbk_ind)

    # Read existing data

    tair_old, time_old = read_old_data(save_folder, "tair")

    prec_old, time_old = read_old_data(save_folder, "prec")

    # Find last available senorge data

    time_first = time_old[end] + Dates.Day(1)
    time_last = find_last_senorge("V2.0")

    time_vec = time_first:time_last

    # Read new meteorological data

    met_var = "tm"

    tair_new, time_new = read_senorge_data(time_vec, met_var, df_geo; scale_func = x -> (x.-2732.0)./10.0)

    met_var = "rr"

    prec_new, time_new = read_senorge_data(time_vec, met_var, df_geo; scale_func = x -> x./10.0)

    # Write data to files

    write_update(save_folder, "tair", time_old, tair_old, time_new, tair_new, df_meta)

    write_update(save_folder, "prec", time_old, prec_old, time_new, prec_new, df_meta)

    # Process runoff data

    time_first = time_old[1]
    
    time_vec = time_first:time_last

    process_runoff_data(stat_list, time_vec, df_meta, path_runoff, save_folder)

end

