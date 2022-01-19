using LoopVectorization
using VectorizationBase
using Random

@inline function ExtPot(lam::Float64, R::VectorizationBase.Vec{4, Float64})
    constant = 1/(4 * lam)
    RR = 0.0
    @turbo for i in eachindex(R,R)
        RR += R[i] * R[i]
    end
end

function main()
    nTsl = 2^11
    nPar = 2^10

    lam = 0.5

    potentials1 = zeros(nTsl, nPar)
    potentials2 = zeros(nTsl, nPar)
    pos = randn(nTsl, nPar)
    
    println("Turbo:")
    #=
    @time begin
        @turbo for tSlice = 1:nTsl
            for ptcl = 1:nPar
                potentials1[tSlice, ptcl] = ExtPot(lam, pos[tSlice, ptcl])
            end
        end
    end
    =#
    k = 1/(4 * lam)
    @time begin
        @turbo for tSlice = 1:nTsl
            for ptcl = 1:nPar
                potentials1[tSlice, ptcl] = k * pos[tSlice, ptcl] * pos[tSlice, ptcl]
            end
        end
    end

    println("\nNon-turbo:")
    #=
    @time begin
        @turbo for tSlice = 1:nTsl
            for ptcl = 1:nPar
                potentials2[tSlice, ptcl] = ExtPot(lam, pos[tSlice, ptcl])
            end
        end
    end
    =#
    @time begin
        @inbounds for tSlice = 1:nTsl
            @simd for ptcl = 1:nPar
                potentials2[tSlice, ptcl] = k * pos[tSlice, ptcl] * pos[tSlice, ptcl]
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
