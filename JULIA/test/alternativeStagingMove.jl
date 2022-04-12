@inbounds function StagingMove(Manent::Function, Param::Params, Path::Paths,
                                ptcl::Int64, rng::MersenneTwister)
    ### Create, initialize variables, vectors, etc.
    oldAction       = 0.0::Float64
    newAction       = 0.0::Float64
    edgeExclusion   = 1::Int64
    m               = 16
    tau             = Param.tau
    oldBeads        = zeros(Float64, m-1)
    oldPotentials   = zeros(Float64, m-1)
    oldDeterminants = zeros(Float64, m-1)

    # Choose the start and end of the stage
    alpha_start     = rand(rng, 1:Param.nTsl)
    alpha_end       = ModTslice((alpha_start + m), Param.nTsl)

    ### Compute old action values, store old values in respective objects
    for a = 1:m-1
        tSlice              = ModTslice(alpha_start + a, Param.nTsl)
        oldBeads[a]         = Path.beads[tSlice, ptcl]
        oldPotentials[a]    = Path.potentials[tSlice, ptcl]
        oldDeterminants[a]  = Path.determinants[tSlice]
        oldAction           += ComputeAction(Param, Path, tSlice)
    end

    ### Perform the staging shift and calculate the new action
    for a = 1:m-1
        tSlice                  = ModTslice((alpha_start + a), Param.nTsl)
        temp                    = ModTslice((tSlice - 1), Param.nTsl)
        temp == 0 ? tSlicem1 = Param.nTsl : tSlicem1 = temp
        tau1                    = (m - a) * tau
        avex                    = (tau1 * Path.beads[tSlicem1, ptcl] + 
                                   tau * Path.beads[alpha_end, ptcl]) / (tau + tau1)
        sigma2                  = 2.0 * Param.lam / (1.0 / tau + 1.0 / tau1)
        Path.beads[tSlice, ptcl]= avex + sqrt(sigma2) * randn(rng)
        UpdatePotential(Param, Path, tSlice, ptcl)
        if Param.nPar != 1
            BuildDeterminantMatrix(Param, Path, tSlice)
            UpdateManent(Manent, Param, Path, tSlice, ptcl)
        end

        newAction               += ComputeAction(Param, Path, tSlice)
    end

    ### Metropolis, accept update or restore old values
    if rand(rng) < exp(-(newAction - oldAction))
        Path.numAcceptStaging += 1
    else
        for a = 1:m-1
            tSlice = ModTslice((alpha_start + a), Param.nTsl)
            Path.beads[tSlice, ptcl]        = oldBeads[a]
            Path.potentials[tSlice, ptcl]   = oldPotentials[a]
            Path.determinants[tSlice]       = oldDeterminants[a]
        end
    end
end
