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

function BuildDeterminantTensor(Param::Params, beads::Array{Float64,}, dets::Array{Float64,3})
    Neg1o2tau = -1.0 / (2.0 * Param.tau)

    for tSlice = 1:Param.nTsl
        tModPlus = ModTslice(tSlice + 1, Param.nTsl)
        for i = 1:Param.nPar
            for j = 1:Param.nPar
                dets[i, j, tSlice] = exp(Neg1o2tau * (beads[tSlice, i] - beads[tModPlus, j])^2)
            end
        end
    end
end

# Functionally identical to BuildDeterminantTensor, but instead of building the whole tensor, just a single
# matrix slice of it - useful for iteration
function BuildDeterminantMatrix(Param::Params, beads::Array{Float64,}, dets::Array{Float64,3}, tSlice::Int64)
    Net1o2tau = -1.0 / (2.0 * Param.tau)
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    for i = 1:Param.nPar
        for j = 1:Param.nPar
            dets[tSlice,i,j] = exp(Neg1o2tau * (Path.beads[tSlice, i] - Path.beads[tModPlus, j])^2 )
        end
    end
end

# For recursion reasons, this should probably be changed to
@inbounds function Determinant(Param::Params, beads::Array{Float64,}, tSlice::Int64)
    #@inbounds function Determinant(Param::Params, Path::Paths, tSlice::Int64)
    Neg1o2tau = -1.0 / (2.0 * Param.tau)
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    if (Param.nPar == 1)
        #return 1.0
        return MathConstants.e
    elseif (Param.nPar == 2)
        return det(beads)   # Technically, 'beads' here is the slice of the 
        # array that has the exponential terms (like immediately below)
    elseif (Param.nPar == -1)
        return ( 
                exp(Neg1o2tau * ( (beads[tSlice, 1] - beads[tModPlus, 1] )^2 + 
                                 (beads[tSlice, 2] - beads[tModPlus, 2] )^2 ) ) - 
                    exp(Neg1o2tau * ( (beads[tSlice, 1] - beads[tModPlus, 2] )^2 + 
                                     (beads[tSlice, 2] - beads[tModPlus, 1] )^2 ) )
               )
    elseif (Param.nPar == 3)
        return (

               )
    else
        println("This part isn't done yet!")
        exit()
    end
end

# For recursion reasons, this should probably be changed to
# @inbounds function Permanant(Param.tau::Float64, Param.nTsl::Int64, Param.nPar::Int64, beads::Array{Float64,}, tSlice::Int64)
@inbounds function Permanant(Param::Params, beads::Array{Float64,}, tSlice::Int64)
    Neg1o2tau = -1.0 / (2.0 * Param.tau)
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    if (Param.nPar == 1)
        #return 1.0
        return MathConstants.e
    elseif (Param.nPar == 2)
        return ( 
                exp(Neg1o2tau * ( (beads[tSlice, 1] - beads[tModPlus, 1] )^2 + 
                                 (beads[tSlice, 2] - beads[tModPlus, 2] )^2 ) ) + 
                    exp(Neg1o2tau * ( (beads[tSlice, 1] - beads[tModPlus, 2] )^2 + 
                                     (beads[tSlice, 2] - beads[tModPlus, 1] )^2 ) )
               )
    else
        println("This part isn't done yet!")
        exit()
    end
end

@inline function Boltzmannant(Param::Params, Path::Paths, tSlice::Int64)
    return MathConstants.e
end

@inbounds function InstantiateManents(Manent::Function, Param::Params, beads::Array{Float64,})
    for tSlice = 1:Param.nTsl
        for ptcl = 1:Param.nPar
            Path.determinants[tSlice, ptcl] = Manent(Param, beads, tSlice)
        end
    end
end

@inbounds function InstantiateManentsDets(Manent::Function, Param::Params, dets::Array{Float64, 3})
    for tSlice = 1:Param.nTsl
        for i = 1:Param.nPar
            for j = 1:Param.nPar
                Path.determinants[tSlice, i] = Manent(Param, dets, tSlice)
            end
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
@inbounds function ComputeAction(Param::Params, Path::Paths, tSlice::Int64)
    # Computes the potential action of a particle along its worldline
    action = 0.0 + 0.0im
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    for ptcl = 1:Param.nPar
        if (CutOff(Path.beads[tSlice,ptcl],Path.beads[tModPlus,ptcl]))
        else
            action += ( Path.potentials[tSlice,ptcl] + Path.potentials[tModPlus,ptcl] ) * 
                abs(log(Complex((Path.determinants[tSlice,ptcl]))))
        end
    end
    return action
end

@inline function UpdatePotential(Path::Paths, tSlice::Int64, ptcl::Int64, lam::Float64)
    @inbounds Path.potentials[tSlice,ptcl] = ExtPotential(lam,Path.beads[tSlice,ptcl])
end

@inbounds function UpdateManent(Manent::Function, Param::Params, beads::Array{Float64,}, tSlice::Int64, ptcl::Int64)
    tModMinus = ModTslice(tSlice - 1, Param.nTsl)
    
    Path.determinants[tModMinus, ptcl] = Manent(Param, beads, tModMinus)
    Path.determinants[tSlice, ptcl] = Manent(Param, beads, tSlice)
end
