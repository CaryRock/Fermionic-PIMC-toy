double KineticEnergy(double beads[][], int numTimeSlices, int numParticles, 
        double lam, double tau)
{
    int dim = 1;
    double tot = 0.0;
    double norm = 1/(4 * lam * tau * tau);
    double delR = 0.0;
    
    for(int i = 0; i < numTimeSlices; i++)
    {
        int tModPlus = ModTslice(i + 1, numTimeSlices);
        for(int j = 0; j < numParticles; j++)
        {
            delR = beads[tModPlus][j] - beads[i][j];

            tot -= norm * delR * delR;
        }
    }

    return dim/2 * numParticles/tau + tot/numTimeSlices;
}

double PotentialEnergy(double potentials[][], int numTimeSlices, 
        int numParticles, double lam)
{
    double pe = 0.0;

    for(int i = 0; i < numTimeSlices; i++)
    {
        for(int j = 0; j < numTimeSlices; j++)
        {
            tModPlus = ModTslice(i + 1, numTimeSlices);
            
            pe += potentials[i][j] + potentials[tModPlus][j];
        }
    }

    return pe / (2*numTimeSlices);
}

double Energy(Paths path, double tau)
{
    double KE = KineticEnergy(**beads, numTimeSlices, numParticles, lam, tau);
    double PE = PotentialEnergy(**potentials, numTimeSlices, numParticles, lam);

    Path.KE = KE;
    Path.PE = PE;

    return KE + PE;
}
