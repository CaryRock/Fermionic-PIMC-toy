# Chooses the function that interact with the action of the particle(s)
function WhichManent(Manent::Function, Determinant::Function, 
        Permanant::Function, Boltzmannant::Function, bosons::Bool, boltzmannons::Bool)
    if bosons && !boltzmannons
        return Permanant
    elseif  boltzmannons && !bosons
        return Boltzmannant
    elseif
        println("Please choose one of \'--bosons\' or \'--boltzmannons\', not both.")
        exit()
    else
        return Determinant
end

function BuildDeterminantTensor(Param::Params, beads::Array{Float64,}, dets::Array{Float64,3})
    Neg1o2tau = -1.0 / (2.0 * Param.tau)

    for tSlice = 1:Param.nTsl
        tModPlus = ModTslice(tSlice + 1, Param.nTsl)
        for i = 1:Param.nPar
            for j = 1:Param.nPar
                dets[tSlice, i, j] = exp(Neg1o2tau * (beads[tSlice, i] - beads[tModPlus, j])^2 )
            end
        end
    end
end

# Functionally identical to BuildDeterminantTensor, but instead of building the whole tensor, just a single
# matrix slice of it - useful for iteration
function BuildDeterminantMatrix(Param::Params, beads::Array{Float64,}, dets::Array{Float64,3}, tSlice::Int64)
    Neg1o2tau = -1.0 / (2.0 * Param.tau)
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    for i = 1:Param.nPar
        for j = 1:Param.nPar
            dets[tSlice, i, j] = exp(Neg1o2tau * (beads[tSlice, i] - beads[tModPlus, j])^2)
        end
    end
end

# For recursion reasons, this should probably be changed to
@inbounds function Determinant(Param::Params, dets::Array{Float64,3}, tSlice::Int64)
    if (Param.nPar != 1)
        #tModPlus = ModTslice(tSlice + 1, Param.nTsl)
        return det(dets[tSlice])
    elseif (Param.nPar == 1)
        return MathConstants.e
    else
        return println("Something has gone wrong here! Exiting...")
        exit()
    end
end

# For recursion reasons, this should probably be changed to
# @inbounds function Permanant(Param.tau::Float64, Param.nTsl::Int64, Param.nPar::Int64, beads::Array{Float64,}, tSlice::Int64)
@inbounds function Permanant(Param::Params, beads::Array{Float64,}, tSlice::Int64)
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)
end

@inline function Boltzmannant(Param::Params, Path::Paths, tSlice::Int64)
    return MathConstants.e
end

@inbounds function InstantiateManents(Manent::Function, Param::Params, beads::Array{Float64,})
end

@inbounds function InstantiatePotentials(Param::Params, Path::Paths)
end

# This should be something like exp(ComputeAction)*Path.determinants[] for 
# computing the density
@inbounds function ComputeAction(Param::Params, Path::Paths, tSlice::Int64)
end

@inline function UpdatePotential(Path::Paths, tSlice::Int64, ptcl::Int64, lam::Float64)
    @inbounds Path.potentials[tSlice,ptcl] = ExtPotential(lam,Path.beads[tSlice,ptcl])
end

@inbounds function UpdateManent(Manent::Function, Param::Params, beads::Array{Float64,}, tSlice::Int64, ptcl::Int64)
    tModMinus = ModTslice(tSlice - 1, Param.nTsl)
    
    Path.determinants[tModMinus, ptcl] = Manent(Param, beads, tModMinus)
    Path.determinants[tSlice, ptcl] = Manent(Param, beads, tSlice)
end
