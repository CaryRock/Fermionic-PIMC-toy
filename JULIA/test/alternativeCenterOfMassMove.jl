@inbounds function CenterOfMassMove(Manent::Function, Param::Params, 
                                Path::Paths, ptcl::Int64, rng::MersenneTwister)
    ### Create basic variables, parameters ####################################
    delta = Param.delta
    shift = Shift(rng, delta)
    oldAction = 0.0::Float64
    newAction = 0.0::Float64

    ### Initialize containers of old values, parameters, arrays, etc. #########
    # Compute the old action
    for tSlice = 1:Param.nTsl
        oldAction += ComputeAction(Param, Path, tSlice)
    end
    
    # Save old dets, potentials, and determinants arrays to avoid recalculating them
    oldPotentials = zeros(Float64, Param.nTsl, Param.nPar)
    oldDeterminants = zeros(Float64, Param.nTsl)
    oldDets = zeros(Float64, Param.nPar, Param.nPar, Param.nTsl)

    for tSlice = 1:Param.nTsl
        oldDeterminants[tSlice] = Path.determinants[tSlice]

        for i = 1:Param.nPar
            oldPotentials[tSlice, i] = Path.potentials[tSlice, i]

            for j = 1:Param.nPar
                oldDets[i, j, tSlice] = Path.dets[i, j, tSlice]
            end
        end
    end

    ### Add shift to chosen worldline, shift beads, recalculate relavant values
    for tSlice = 1:Param.nTsl
        Path.beads[tSlice, ptcl] += shift
        BuildDeterminantMatrix(Param, Path, tSlice)
        UpdatePotential(Param, Path, tSlice, ptcl)
    end

    if Param.nPar != 1
        # Update determinants for multiple particles - leave alone for N = 1
        for tSlice = 1:Param.nTsl
            UpdateManent(Manent, Param, Path, tSlice, ptcl)
        end
    end

    # Compute value for new action
    for tSlice = 1:Param.nTsl
        newAction += ComputeAction(Param, Path, tSlice)
    end

    ### Metropolis, accept update or restore old values
    if rand(rng) < exp(-abs(newAction - oldAction))
        Path.numAcceptCOM += 1
    else
        for tSlice = 1:Param.nTsl
            Path.beads[tSlice, ptcl] -= shift
            Path.determinants[tSlice] = oldDeterminants[tSlice]

            for i = 1:Param.nPar
                Path.potentials[tSlice, i] = oldPotentials[tSlice, i]

                for j = 1:Param.nPar
                    Path.dets[i, j, tSlice] = oldDets[i, j, tSlice]
                end
            end
        end
    end
end
