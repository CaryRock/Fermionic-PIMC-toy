# Contains the functions that interact with the action of the particle(s)
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

@inbounds function Determinant(Param::Params, Path::Paths,tSlice::Int64)
    # Just short-circuit the whole thing for Boltzmannons
    return 1.0
    tau = Param.tau
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    if (Param.nPar == 1)
        return exp(-1/(2*tau) * 
                   (Path.beads[tSlice,1] - Path.beads[tModPlus,1])^2 )
    end

    for ptclRow = 1:Param.nTsl
        for ptclCol = 1:Param.nPar
            # Iterate over the beads - recall, Julia is Column-major
            if (CutOff(Path.beads[tSlice,ptclRow], Path.beads[tModPlus,ptclCol]))
                Path.determinants[ptclRow,ptclCol] = 1.0
            else
                Path.determinants[ptclRow,ptclCol] = exp(-1/(2*tau) * 
                    (Path.beads[tSlice,ptclRow] - Path.beads[tModPlus,ptclCol])^2)
            end
        end
    end
    return det(Path.determinants)
end

@inbounds function Permanent(Param::Params, Path::Paths, tSlice::Int64)
    Neg1o2tau = -1/(2 * Param.tau)
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    if (Param.nPar == 1)
        return exp(Neg1o2tau * 
                   (Path.beads[tSlice,1] - Path.beads[tModPlus,1])^2 )
    elseif (Param.nPar == 2)
        return ( 
                exp(Neg1o2tau * (Path.beads[tSlice, 1] + Path.beads[tModPlus, 2] )^2) + 
                exp(Neg1o2tau * (Path.beads[tSlice, 2] + Path.beads[tModPlus, 1] )^2)
               )
    else
        println("This part isn't done yet!")
        exit()
    end
end

@inline function Boltzmannant(Param::Params, Path::Paths, tSlice::Int64)
    return 1.0
end
#=
@inbounds function InitializeDeterminants(Param::Params, Path::Paths)
    for tSlice = 1:Param.nTsl
        for ptcl = 1:Param.nPar
            Path.determinants[tSlice,ptcl] = Determinant(Param, Path, tSlice)
        end
    end
end
# Note: Potential "simplification": add as an argument "Manent::Function"to a 
# more general form of either of these two and combine both into the same. 
@inbounds function InitializePermanents(Param::Params, Path::Paths)
    for tSlice = 1:Param.nTsl
        for ptcl = 1:Param.nPar
            Path.determinants[tSlice, ptcl] = Permanent(Param, Path, tSlice) #TODO: Is this a todo? It's convenient, though
        end
    end
end
=#
@inbounds function UpdateManents(Manent::Function, Param::Params, Path::Paths)
    for ptcl = 1:Param.nPar
        for tSlice = 1:Param.nTsl
            Path.determinants[tSlice, ptcl] = Manent(Param, Path, tSlice)
        end
    end
end

@inbounds function InitializeBoltzmannant(Param::Params, Path::Paths)
    for tSlice = 1:Param.nTsl
        for ptcl = 1:Param.nPar
            Path.determinants[tSlice, ptcl] = 1.0
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
                Path.potentials[tSlice,ptcl] = ExtPotential(Param.lam,Path.beads[tSlice,ptcl]) 
            end
        end
    end
end

# This should be something like exp(ComputeAction)*Path.determinants[] for 
# computing the density
function ComputeAction(Param::Params, Path::Paths, tSlice::Int64)
    # Computes the potential action of a particle along its worldline
    action = 0.0
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    @inbounds for ptcl = 1:Param.nPar
        if (CutOff(Path.beads[tSlice,ptcl],Path.beads[tModPlus,ptcl]))
            action += 0.0
        else
            action += ( Path.potentials[tSlice,ptcl] + Path.potentials[tModPlus,ptcl] ) * 
                        log(Path.determinants[tSlice,ptcl])
        end
    end
    return action
end

Base.@propagate_inbounds function AltComputeAction(Param::Params, Path::Paths, action::Float64)
    # Computes the potential action of a particle along its worldline
    #action = 0.0
    result = zeros(Float64, Param.nPar * Param.nTsl)

    for ptcl = 1:Param.nPar
        for tSlice = 1:Param.nTsl
            tModPlus = ModTslice(tSlice + 1, Param.nTsl)
            if ( CutOff(Path.beads[tSlice, ptcl], Path.beads[tModPlus, ptcl]) )
                #action += 0.0
                result[(ptcl - 1) * Param.nTsl + tSlice ] = 0.0
            else
                #action += ( Path.potentials[tSlice, ptcl] + Path.potentials[tModPlus, ptcl] ) * 
                #            Path.determinants[tSlice, ptcl]
                #
                result[(ptcl - 1) * Param.nTsl + tSlice ] = ( Path.potentials[tSlice, ptcl] + Path.potentials[tModPlus, ptcl] ) * 
                            Path.determinants[tSlice, ptcl]
            end
        end
    end
    
    for ptcl = 1:Param.nPar
        for tSlice = 1:Param.nTsl
            action += result[(ptcl - 1) * Param.nTsl + tSlice ]
        end
    end

    #return action
end

# TODO: FIGURE OUT HOW TO MAKE THIS FUNCTION WORK WITH ExtPotential() AND @turbo
@inline function UpdatePotential(Path::Paths, tSlice::Int64, ptcl::Int64, lam::Float64)
    @inbounds Path.potentials[tSlice,ptcl] = ExtPotential(lam,Path.beads[tSlice,ptcl])
end

function UpdateDeterminant(Param::Params, Path::Paths, tSlice::Int64, ptcl::Int64)
    tModMinus = ModTslice(tSlice - 1, Param.nTsl)

    Path.determinants[tModMinus,ptcl] = Determinant(Param, Path, tModMinus)
    Path.determinants[tSlice,ptcl] = Determinant(Param, Path, tSlice)
end
# Could these be combined into a "UpdateManent"?
function UpdatePermanent(Param::Params, Path::Paths, tSlice::Int64, ptcl::Int64)
    tModMinus = ModTslice(tSlice - 1, Patam.nTsl)

    Path.determinants[tModMinus, ptcl] = Permanent(Param, Path, tModMinus)
    Path.determinants[tSlice, ptcl] = Permanent(Param, Path, tSlice)
end

function UpdateManent(Manent::Function, Param::Params, Path::Paths, tSlice::Int64, ptcl::Int64)
    tModMinus = ModTslice(tSlice - 1, Param.nTsl)

    Path.determinants[tModMinus, ptcl] = Manent(Param, Path, tModMinus)
    Path.determinants[tSlice, ptcl] = Manent(Param, Path, tSlice)
end
