# Contains the functions that interact with the action of the particle(s)
function WhichManent(Manent::Function, Determinant::Function, Permanant::Function, boson::Bool, boltzmannon::Bool)
    if boson
        return Permanant
    else
        return Determinant
    end
end

# TODO: SEE IF THIS FUNCTION CAN BE COMPILED/WRAPPED IN C++ CODE USING BOOST
@inbounds function Determinant(Param::Params, Path::Paths,tSlice::Int64)
    # Just short-circuit the whole thing for Boltzmannons
    Neg1o2tau = -1.0/(2.0 * Param.tau)
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    if (Param.nPar == 1)    # Scalar value
        return exp( Neg1o2tau * 
                   (Path.beads[tSlice,1] - Path.beads[tModPlus,1])^2 )
    
    elseif (Param.nPar == 2)    # Simple 2x2 matrix, return simple value
        return (
    exp(Neg1o2tau * ((Path.beads[tSlice, 1] - Path.beads[tModPlus, 1])^2 + (Path.beads[tSlice, 2] - Path.beads[tModPlus, 2] )^2)) -
    exp(Neg1o2tau * ((Path.beads[tSlice, 2] - Path.beads[tModPlus, 1])^2 + (Path.beads[tModPlus, 1] - Path.beads[tSlice, 2] )^2))
               )
    else
        println("This part isn't done yet!")
        exit()
    end
end

# TODO: SEE IF THIS FUNCTION CAN BE COMPILED/WRAPPED IN C++ CODE USING BOOST
@inbounds function Permanent(Param::Params, Path::Paths, tSlice::Int64)
    Neg1o2tau = -1.0/(2.0 * Param.tau)
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    if (Param.nPar == 1)
        return exp(Neg1o2tau * 
                   (Path.beads[tSlice,1] - Path.beads[tModPlus,1])^2 )
    elseif (Param.nPar == 2)
        return ( 
    exp(Neg1o2tau * ((Path.beads[tSlice, 1] - Path.beads[tModPlus, 1])^2 + (Path.beads[tSlice, 2] - Path.beads[tModPlus, 2] )^2)) +
    exp(Neg1o2tau * ((Path.beads[tSlice, 2] - Path.beads[tModPlus, 1])^2 + (Path.beads[tModPlus, 1] - Path.beads[tSlice, 2] )^2))
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
    for ptcl = 1:Param.nPar
        for tSlice = 1:Param.nTsl
            Path.determinants[tSlice,ptcl] = Determinant(Param, Path, tSlice)
        end
    end
end
# Note: Potential "simplification": add as an argument "Manent::Function"to a 
# more general form of either of these two and combine both into the same. 
@inbounds function InitializePermanents(Param::Params, Path::Paths)
    for ptcl = 1:Param.nPar
        for tSlice = 1:Param.nTsl
            Path.determinants[tSlice, ptcl] = Permanent(Param, Path, tSlice)
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
                        Path.determinants[tSlice,ptcl]
        end
    end
    return action
end

Base.@propagate_inbounds function AltComputeAction(Param::Params, Path::Paths, action::Float64)
    # Computes the potential action of a particle along its worldline
    result = zeros(Float64, Param.nPar * Param.nTsl)

    for ptcl = 1:Param.nPar
        for tSlice = 1:Param.nTsl
            tModPlus = ModTslice(tSlice + 1, Param.nTsl)
            if ( CutOff(Path.beads[tSlice, ptcl], Path.beads[tModPlus, ptcl]) )
                #action += 0.0
                result[(ptcl - 1) * Param.nTsl + tSlice ] = 0.0
            else
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
end

@inline function UpdatePotential(Path::Paths, tSlice::Int64, ptcl::Int64, lam::Float64)
    @inbounds Path.potentials[tSlice,ptcl] = ExtPotential(lam,Path.beads[tSlice,ptcl])
end

### These functions are meant to update a all the determiinants as shifting
# a worldline leads to all the beads changing distance. Changed bead distance
# means new determinants.


### These functions are meant to update single determinants as opposed to a 
# whole worldline of them. 
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
    Path.determinants[tSlice, ptcl] = Manent(Param, Path, tModMinus)
end
