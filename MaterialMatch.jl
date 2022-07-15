function SelectSpectra(sp)
    n = length(sp)
    Spectra = Vector{Any}(undef, n)
    for sno=1:n
        x,y, = Vector.(map(x->reinterpret(Float32, base64decode(x)),(get(sp[sno],"XData",1),get(sp[sno],"YData",1))))

        CBand = BandMatch(x)
        Spectra[sno] = [x[CBand],y[CBand]]
        #Spectra[sno] = [x,y]
    end
    return Spectra
end

function BandMatch(x::Vector{Float32})::Vector{Int32}
    CorresBand = Vector{Int32}(undef, length(SelectedBands))
    for i in eachindex(SelectedBands)
        Band = findmin(abs.(x .- SelectedBands[i]))
        if first(Band) < 0.002
            CorresBand[i] = last(Band)
        else
            return ones(length(SelectedBands))
        end
    end
    return CorresBand
end

function FindMaterial(ObservedSpectra,
                      SpectraLibrary)

    materialmatch = Vector{Int32}(undef,NoEM2Extract)

    for j in eachindex(ObservedSpectra)

        diff = Vector{Float32}(undef,length(SpectraLibrary))
        for i in eachindex(SpectraLibrary)

            # Calculate difference in materials and VCA extracted spectra.
            diff[i] = norm(SpectraLibrary[i] - ObservedSpectra[j])

        end

        # Find the closest match and store its index
        index = findmin(diff)
        if first(index) < 1
            materialmatch[j] = last(index)
        else
            materialmatch[j] = 3
        end

    end
    return materialmatch
end

# Below we load the entire spectral library as a dictionary
#sm = JSON.parsefile("/EcoStress/SamplesInfo.json")
#sp = JSON.parsefile("/EcoStress/SpectralData.json")

# Read the spectral library dictionary/load processed library

#SpectralLibrary = last.(SelectSpectra(sp))

# Scale and center the spectra
#SpectralLibrary = map(x -> x/norm(x), map(x -> x .- minimum(x), SpectralLibrary))

# Scale and center the observed spectra
#ObservedSpectra = map(x -> x/norm(x), map(x -> x .- minimum(x), ObservedSpectra))

#materialmatch = FindMaterial(ObservedSpectra, SpectralLibrary)
