
"""
    init_dataset(stat_list, time_vec, save_folder)

Create an initial dataset.
"""
function init_dataset(stat_list, time_vec, save_folder)

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

    # Process metadata for the individual elevation bands

    metadata_elevbands(stat_list, save_folder, elev_ind, stat_dbk)

    # Read meteorological data

    met_var = "tm"

    tair_data, tair_time = read_senorge_data(time_vec, met_var, elev_ind; scale_func = x -> (x-2732.0)/10.0)

    met_var = "rr"

    prec_data, prec_time = read_senorge_data(time_vec, met_var, elev_ind; scale_func = x -> x/10.0)

    # Write data to files

    write_met_data(stat_list, save_folder, "Tair", time_vec, tair_data, stat_dbk)

    write_met_data(stat_list, save_folder, "Prec", time_vec, prec_data, stat_dbk)

    # Process runoff data

    process_runoff_data(stat_list, time_vec, df_meta, path_runoff, save_folder)

end

