function IdentifySpectra(ObsSpec, SpecLibrary)

    materialmatch = Vector{Int32}(undef,NoEM2Extract)

    for j in eachindex(ObsSpec)

        diff = Vector{Float32}(undef,length(SpecLibrary))

        for i in eachindex(SpecLibrary)
            # Calculate difference in materials and VCA extracted spectra.
            diff[i] = norm(SpecLibrary[i] - ObsSpec[j])

        end

        # Find the closest match and store its index
        index = findmin(diff)

        if first(index) < 1 # Match must be this close
            materialmatch[j] = last(index)

        else
            materialmatch[j] = 1 # Set material to unknown
        end
    end

    return materialmatch
end
