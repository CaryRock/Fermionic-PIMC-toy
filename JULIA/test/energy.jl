# This file contains most of the code relevant to the energy aspect of the 
# program, sans "ExtPotential" which remains in the main file (for now).
# Energy estimators are to go here.

# This is the basic method of computing the KE estimator. This shouldn't be 
# used. Instead, the thermodynamic KE esetimator will be used.
@inbounds function KineticEnergy(Param, beads::Matrix{Float64})
    # Computes the KE of the particle(s)
    #dim = 1.0   # Float to aid julia's casting
    tot = 0.0
    norm = 1.0/(4.0 * Param.lam * Param.tau * Param.tau)
    for tSlice = 1:Param.nTsl
        tModPlus = ModTslice(tSlice + 1, Param.nTsl)
        for ptcl = 1:Param.nPar
            delR = beads[tModPlus,ptcl] - beads[tSlice,ptcl]
            
            tot -= norm * delR * delR
        end
    end
    energy = 0.5 * Param.nPar / Param.tau + tot/Param.nTsl
    return energy/Param.nPar
end

@inbounds function PotentialEnergy(Param::Params, potentials::Matrix{Float64})
    # Computes the potential energy of the particle(s)
    pe = 0.0
    
    for tSlice = 1:Param.nTsl
        #tModPlus = ModTslice(tSlice + 1, Param.nTsl)
        for ptcl = 1:Param.nPar
            #pe += potentials[tSlice, ptcl] + potentials[tModPlus,ptcl]
            pe += potentials[tSlice, ptcl]
        end
    end
    return pe / (Param.nPar * Param.nTsl) #(2.0 * Param.nTsl)
end

function Energy(Param::Params, Path::Paths)
    # Returns the KE + PE of the particle(s)

    Path.KE = KineticEnergy(Param, Path.beads)
    Path.PE = PotentialEnergy(Param, Path.potentials)
    
    return (Path.KE + Path.PE)
end

