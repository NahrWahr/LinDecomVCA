"""
        * Ind::CartesianIndex(x,y) where x,y are pixel position
        * Vec::Vector{pType <: UInt8} where Vec is spectral data of pixel
        """
struct IPixel{pType <: UInt8}
    Ind::CartesianIndex
    Vec::Vector{pType}
end

"""
    * (Removes pixels which are 0 at all bands.)
         * (Returns the data in the format of Tuple)
    """
function Discard0(X::Matrix{IPixel})::Tuple{Vector{CartesianIndex}, Matrix{UInt8}}
    n = size(X,1)
    Out = Vector{Vector{UInt8}}(undef,0)
    Index = Vector{CartesianIndex{2}}(undef,0)
    for i in 1:n
	if sum(X[i].Vec) != 0
	    push!(Out, X[i].Vec)
	    push!(Index, X[i].Ind)
	end
    end
    Out = permutedims(hcat(Out...))
    return Index,Out
end

"""
    * (Reshapes(flattens) the 3D array into a 2D array of dimensions L and M×N)
        * (Stores the Cartesian Indices of the pixels position)
        * (Returns the data in the format of IPixel (IndexedPixel))
    """
@inline function IndexFlat(Data::Array{UInt8, 3})::Matrix{IPixel}
    a,b,c = size(Data)

    A = Matrix{IPixel}(undef,a*b,1)
    @inbounds for i=1:a
	@inbounds for j=1:b
	    A[b*(i-1)+j] = IPixel(CartesianIndex(i,j), Data[i,j,:])
	end
    end
    return A
end
