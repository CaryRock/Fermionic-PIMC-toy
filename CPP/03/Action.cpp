//double Determinant(double beads[][], int numTimeSlices, int numParticles, int tslice, int iter, double lam)
double Determinant(Paths Path, double tau, int tSlice)
{
    
    return 1    // determinant of the thing - 1 for now
}

void InstantiatePotentials(double beads[][], double potentials[][], 
        int numTimeSlices, int numParticles, double lam, int iter)
{
    int tModPlus = 0;
    double vextT = 0;
    double vexTplus = 0;

    for(int i = 0; i < numTimeSlices; i++)
    {
        for(int j = 0; j < numParticles; j++)
        {
            tModPlus = ModTslice(i + 1, numTimeSlices);
            vextT = ExtPotential(beads[tslice][j], lam);
            vextTplus = ExtPotential(beads[tModPlus][j], lam);

            if (CutOff(beads[tslice][j], beads[tModPlus][j]) == true)
                potentials[i][j] = 0;
            else
                potentials[i][j] = -1/2 * (vextT + vextTplus);
        }
    }
}

double ComputeAction(double beads[][], double potentials[][], 
        double determinants[][], double lam, double tau, int numTimeSlices, 
        int numParticles, int tslice)
{
    int tModPlus = ModTslice(tslice, numTimeSlices);
    double action = 0;
    double vextT = 0;
    double vextTplus = 0;

    for(int i = 0; i < numParticles; i++)
    {
        vextT = potentials[tslice][i];
        vextTplus = potentials[tModPlus][i];

        if ( CutOff(beads[tslice][i], beads[tModPlus][i]) == true)
            action += 0.0;
        else
            action += (vextT + vextTplus) * determinants[tslice][i];
    }

    return action
}

inline void UpdatePotential(double beads[][], double potentials[][], 
        double lam, int tSlice, int ptcl)
{
    potentials[tSlice][ptcl] = ExtPotential(lam, beads[tSlice][ptcl]);
}

void UpdateDeterminant(Paths Path, double tau, int tSlice, int ptcl)
{
    tModMinus = ModTslice(tSlice - 1, numTimeSlices);

    determinants[tModMinus][ptcl] = Determinant(Path, tSlice, tau);
    determinants[tSlice][ptcl] = Determinant(Path, tSlice, tau);
}
