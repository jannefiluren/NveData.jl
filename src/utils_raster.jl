
"""
senorge_info()

Get coordinate information about senorge grids
"""
function senorge_info()

nrows, ncols = 1550, 1195

xllcorner, yllcorner = -75000, 6450000

senorge_ind = reshape(collect(1:nrows*ncols), ncols, nrows)'

xvec = collect(xllcorner:1000:(xllcorner+(ncols-1)*1000))'

xcoord = repmat(xvec, nrows, 1)

yvec = collect(yllcorner+(nrows-1)*1000:-1000:yllcorner)

ycoord = repmat(yvec, 1, ncols)

return senorge_ind, xcoord, ycoord

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
data[data .== raster["NODATA_value"]] .= NaN
data = permutedims(data)
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
