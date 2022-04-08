#=
# This file contains the implementations for finding expectation values
=#
# TODO: UPDATE THIS FOR MULTIPLE PARTICLES
#   Probably should split the functionality into separate functions - 
#   a single particle set and a multi-particle set
#   Some should probably only be available if there are multiple particles
function x1_expectation(beads::Matrix{Float64}, nTsl::Int64, nPar::Int64)
    x1_av = 0.0

    for j = 1:nTsl
        for i = 1:nPar
            x1_av += beads[j,i]
        end
    end
    
    x1_av = x1_av / nTsl
    return x1_av
end

function x2_expectation(beads::Matrix{Float64}, nTsl::Int64, nPar::Int64)
    x2_av = 0.0
    
    for j = 1:nTsl
        for i = 1:nPar
            x2_av += beads[j,i] * beads[j,i]
        end
    end

    x2_av = x2_av / nTsl
    return x2_av
end

#function Px_expectation(beads::Matrix{Float64}, nTsl::Int64, nPar::Int64)
