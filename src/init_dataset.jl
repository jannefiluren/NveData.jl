
"""
    init_dataset(stat_list, time_vec, save_folder; path_runoff)

Create an initial dataset.
"""
function init_dataset(stat_list, time_vec, save_folder;
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

    # Process metadata for the individual elevation bands

    metadata_elevbands(save_folder, df_geo, df_meta)

    # Read meteorological data

    met_var = "tm"

    tair_data, tair_time = read_senorge_data(time_vec, met_var, df_geo; scale_func = x -> (x .- 2732.0)./10.0)

    met_var = "rr"

    prec_data, prec_time = read_senorge_data(time_vec, met_var, df_geo; scale_func = x -> x/10.0)

    # Write data to files

    write_met_data(save_folder, "tair", time_vec, tair_data, df_meta)

    write_met_data(save_folder, "prec", time_vec, prec_data, df_meta)

    # Process runoff data

    process_runoff_data(stat_list, time_vec, df_meta, path_runoff, save_folder)

end

