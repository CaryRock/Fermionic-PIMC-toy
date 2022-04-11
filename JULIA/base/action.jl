# Chooses the function that interact with the action of the particle(s)
function WhichManent(Manent::Function, Determinant::Function, 
        Permanant::Function, Boltzmannant::Function, bosons::Bool, boltzmannons::Bool)
    if bosons && !boltzmannons
        return Permanant
    elseif boltzmannons && !bosons
        return Boltzmannant
    elseif bosons && boltzmannons
        println("Please choose one of \'--bosons\' or \'--boltzmannons\', not both.")
        exit()
    else
        return Determinant
    end
end

# For recursion reasons, this should probably be changed to
# @inbounds function Determinant(Param.tau::Float64, Param.nTsl::Int64, Param.nPar::Int64, beads::Array{Float64,}, tSlice::Int64)
@inbounds function Determinant(Param::Params, Path::Paths,tSlice::Int64)
    Neg1o2tau = -1.0 / (2.0 * Param.tau)
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    if (Param.nPar == 1)
        #return 1.0
        return MathConstants.e
    elseif (Param.nPar == 2)
        return ( 
                exp(Neg1o2tau * ( (Path.beads[tSlice, 1] - Path.beads[tModPlus, 1] )^2 + 
                                 (Path.beads[tSlice, 2] - Path.beads[tModPlus, 2] )^2 ) ) - 
                    exp(Neg1o2tau * ( (Path.beads[tSlice, 1] - Path.beads[tModPlus, 2] )^2 + 
                                     (Path.beads[tSlice, 2] - Path.beads[tModPlus, 1] )^2 ) )
               )
    else
        println("This part isn't done yet!")
        exit()
    end
end

# For recursion reasons, this should probably be changed to
# @inbounds function Permanent(Param.tau::Float64, Param.nTsl::Int64, Param.nPar::Int64, beads::Array{Float64,}, tSlice::Int64)
@inbounds function Permanent(Param::Params, Path::Paths, tSlice::Int64)
    Neg1o2tau = -1.0 / (2.0 * Param.tau)
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    if (Param.nPar == 1)
        #return 1.0
        return MathConstants.e
    elseif (Param.nPar == 2)
        return ( 
                exp(Neg1o2tau * ( (Path.beads[tSlice, 1] - Path.beads[tModPlus, 1] )^2 + 
                                 (Path.beads[tSlice, 2] - Path.beads[tModPlus, 2] )^2 ) ) + 
                    exp(Neg1o2tau * ( (Path.beads[tSlice, 1] - Path.beads[tModPlus, 2] )^2 + 
                                     (Path.beads[tSlice, 2] - Path.beads[tModPlus, 1] )^2 ) )
               )
    else
        println("This part isn't done yet!")
        exit()
    end
end

@inline function Boltzmannant(Param::Params, Path::Paths, tSlice::Int64)
    return MathConstants.e
end

@inbounds function InstantiateManents(Manent::Function, Param::Params, Path::Paths)
    for ptcl = 1:Param.nPar
        for tSlice = 1:Param.nTsl
            Path.determinants[tSlice, ptcl] = Manent(Param, Path, tSlice)
        end
    end
end

@inbounds function InstantiatePotentials(Param::Params, Path::Paths)
    for tSlice = 1:Param.nTsl
        for ptcl = 1:Param.nPar
            tModPlus = ModTslice(tSlice+1,Param.nTsl)
            vextT = ExtPotential(Param.lam, Path.beads[tSlice,ptcl])
            vextTPlus = ExtPotential(Param.lam, Path.beads[tModPlus,ptcl])

            if (CutOff(Path.beads[tSlice,ptcl],Path.beads[tModPlus,ptcl]))
                Path.potentials[tSlice,ptcl] = 0.0
            else
                #Path.potentials[tSlice,ptcl] = Param.tau / 2.0 * vextT
                Path.potentials[tSlice,ptcl] = Param.tau / 2.0 * (vextT + vextTPlus) 
                #Path.potentials[tSlice,ptcl] = ExtPotential(Param.lam,Path.beads[tSlice,ptcl]) 
            end
        end
    end
end

# This should be something like exp(ComputeAction)*Path.determinants[] for 
# computing the density
@inbounds function ComputeAction(Param::Params, Path::Paths, tSlice::Int64)
    # Computes the potential action of a particle along its worldline
    action = 0.0 #+ 0.0im
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    for ptcl = 1:Param.nPar
        if (CutOff(Path.beads[tSlice,ptcl],Path.beads[tModPlus,ptcl]))
        else
            action += (Path.potentials[tSlice, ptcl] + Path.potentials[tModPlus, ptcl]) + 
                abs(log(Complex(Path.determinants[tSlice, ptcl])))
#            action += ( Path.potentials[tSlice,ptcl] + Path.potentials[tModPlus,ptcl] ) * 
#                abs(log(Complex((Path.determinants[tSlice,ptcl]))))
        end
    end
    return action
end

@inbounds function UpdatePotential(Param::Params, Path::Paths, tSlice::Int64, ptcl::Int64)
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    vextT = ExtPotential(Param.lam, Path.beads[tSlice, ptcl])
    vexTPlus = ExtPotential(Param.lam, Path.beads[tModPlus, ptcl])

    Path.potentials[tSlice,ptcl] = Param.tau / 2.0 * (vextT + vexTPlus)
end

@inbounds function UpdateManent(Manent::Function, Param::Params, Path::Paths, tSlice::Int64, ptcl::Int64)
    tModMinus = ModTslice(tSlice - 1, Param.nTsl)
    
    Path.determinants[tModMinus, ptcl] = Manent(Param, Path, tModMinus)
    Path.determinants[tSlice, ptcl] = Manent(Param, Path, tSlice)
end
