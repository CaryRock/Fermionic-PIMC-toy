#=
Chooses the function that interact with the action of the particle(s)
Input: Function handle, Determinant/Permanant/Boltzmannant function handles,
        bools for commandline options for which particles being simulated
=#
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

#=
Initializes/rebuilds the "Determinant" tensor, which stores all the (uncomputed)
values of the determinants. Meant to be fed to the Mannet() function that was
determined above.
Inputs: Parameters, particle beads matrix, uncomputed/empty determinant tensor
=#
function BuildDeterminantTensor(Param::Params, beads::Array{Float64,}, dets::Array{Float64,3})
    Neg1o2tau = -1.0 / (2.0 * Param.tau)

    for tSlice = 1:Param.nTsl
        tModPlus = ModTslice(tSlice + 1, Param.nTsl)    # Update every iteration
        for i = 1:Param.nPar
            for j = 1:Param.nPar
                dets[i, j, tSlice] = exp(Neg1o2tau * (beads[tSlice, i] - beads[tModPlus, j])^2)
            end
        end
    end
end

#=
Functionally identical to BuildDeterminantTensor, but instead of building the 
whole tensor, just a single matrix slice of it - useful for iteration or for
rebuilding a single element/determinant
Inputs: Parameters, particle beads matrix, determinant tensor, specific time 
        slice to rebuild
=#
function BuildDeterminantMatrix(Param::Params, beads::Array{Float64,}, dets::Array{Float64,3}, tSlice::Int64)
    Net1o2tau = -1.0 / (2.0 * Param.tau)
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    for i = 1:Param.nPar
        for j = 1:Param.nPar
            dets[i,j,tSlice] = exp(Neg1o2tau * (Path.beads[tSlice, i] - Path.beads[tModPlus, j])^2 )
        end
    end
end

#=
Computes the determinant of the given matrix.
=#
# For recursion reasons, this should probably be changed to
@inbounds function Determinant(Param::Params, matrix::Array{Float64,})
    #@inbounds function Determinant(Param::Params, Path::Paths, tSlice::Int64)
    #Neg1o2tau = -1.0 / (2.0 * Param.tau)
    #tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    if (Param.nPar == 1)
        #return 1.0
        return MathConstants.e
    elseif (Param.nPar > 1)
        return det(matrix)
    else
        println("Something has gone wrong in the determinant function!")
        exit()
    end
end

# For recursion reasons, this should probably be changed to
# @inbounds function Permanant(Param.tau::Float64, Param.nTsl::Int64, Param.nPar::Int64, beads::Array{Float64,}, tSlice::Int64)
@inbounds function Permanant(Param::Params, matrix::Array{Float64,}, tSlice::Int64)
    Neg1o2tau = -1.0 / (2.0 * Param.tau)
    tModPlus = ModTslice(tSlice + 1, Param.nTsl)

    if (Param.nPar == 1)
        return 1.0  #MathConstants.e
    elseif (Param.nPar == 2)
        return ( 
                matrix[tSlice, 1]*matrix[tModPlus, 2] + matrix[tModPlus, 1]*matrix[tSlice, 2]
                #exp(Neg1o2tau * ( (matrix[tSlice, 1] - matrix[tModPlus, 1] )^2 + 
                #                 (matrix[tSlice, 2] - matrix[tModPlus, 2] )^2 ) ) + 
                #    exp(Neg1o2tau * ( (matrix[tSlice, 1] - matrix[tModPlus, 2] )^2 + 
                #                     (matrix[tSlice, 2] - matrix[tModPlus, 1] )^2 ) )
               )
    else
    # TODO: FINISH THE PROPER IMPLEMENTATION OF THIS
        println("This part isn't done yet!")
        exit()
    end
end
#=
Simply returns 1.0, since e^ln(1) = 1, and adding in action calculations
=#
@inline function Boltzmannant(Param::Params, Path::Paths, tSlice::Int64)
    return 1.0  #MathConstants.e
end

#=
Initializes the (computed) determinants vector. Only needs a single index, 
for iterating over time slices. That is, holds the computed determinants or
permanants.
Inputs: permutation function, parameters, particle beads matrix
=#
TODO: REWRITE THIS TO WORK FOR VECTOR Path.determinants
@inbounds function InstantiateManents(Manent::Function, Param::Params, detsMatrix::Array{Float64,3})
    for tSlice = 1:Param.nTsl
        for ptcl = 1:Param.nPar
            Path.determinants[tSlice] = Manent(Param, detsMatrix[:,:,tSlice])
        end
    end
end

#=
Instantiates determinants tensor elements. This tensor is what is fed into
the permutation function to compute the contribution from particle permutations.
Inputs: permutation function, parameters, determinant tensor
=#
@inbounds function InstantiateManentsDets(Manent::Function, Param::Params, dets::Array{Float64, 3})
    for tSlice = 1:Param.nTsl
        for i = 1:Param.nPar
            for j = 1:Param.nPar
                Path.dets[i, j, tSlice] = 
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
                abs(log(Complex((Path.determinants[tSlice]))))
        end
    end
    return action
end

@inbounds function UpdatePotential(Param::Params, Path::Paths, tSlice::Int64, ptcl::Int64)
    tModPlus = ModTSlice(tSlice + 1, Param.nTsl)
    vextT = ExtPotential(Param.lam, Path.beads[tSlice, ptcl])
    vextTPlus = ExtPotential(Param.lam, Path.beads[tModPlus, ptcl])

    #Path.potentials[tSlice,ptcl] = ExtPotential(Param.lam,Path.beads[tSlice,ptcl])
    Path.potentials[tSlice, ptcl] = Param.tau / 2.0 * (vextT + vextTPlus)
end

@inbounds function UpdateManent(Manent::Function, Param::Params, beads::Array{Float64,}, tSlice::Int64, ptcl::Int64)
    tModMinus = ModTslice(tSlice - 1, Param.nTsl)
    
    Path.determinants[tModMinus] = Manent(Param, beads, tModMinus)
    Path.determinants[tSlice] = Manent(Param, beads, tSlice)
end
