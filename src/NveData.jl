module NveData

using DataStructures
using DataFrames
using CSV

export watershed_info

export read_esri_raster
export raster2vec
export read_bil, read_bil!

export read_dbk_ind, subset_dbk_ind
export read_metadata, subset_metadata

export read_senorge_data
export aggregate_elevations, aggregate_elevations!
export senorge_info

export write_met_data
export write_update
export metadata_elevbands
export process_runoff_data

export init_dataset, update_dataset

export find_last_senorge

export read_old_data

include("utils.jl")
include("utils_raster.jl")
include("init_dataset.jl")
include("update_dataset.jl")

end
