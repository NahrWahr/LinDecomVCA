using ArchGDAL
using JLD2

function LoadTIF(path::String="~/Pixxel/project/testdata/Hyperion_Canada.tif")::Array{UInt8}
    # Requires ArchGDAL
    Data = ArchGDAL.readraster(path)

    return Data
end

# Select bands which are to be discarded or selected
# SelectedBands = (hyperion_bands[set]) * 1e-3

# Uncomment line below to load a different TIF image.
# Pass the path of the image to $LoadTIF function

#Index, X = Discard0(IndexFlat(LoadTIF()[:,100:end-100,:]))

# Save the proccessed Index and X from above to reduce times if using same file often
#@save "/home/rnarwar/Pixxel/project/testdata/IndARD.jld2" Index X

# IndARD.jld2 contains variables :Index and :X.
# :Index is cartesian position of pixel in image.
# :X is the actual spectral data of each pixel, i.e. the spectral vector.
# @load "/home/rnarwar/Pixxel/project/testdata/IndARD.jld2"
# X = permutedims(X)[set,:]
