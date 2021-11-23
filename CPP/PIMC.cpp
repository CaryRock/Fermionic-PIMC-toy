#ifndef __RANDOM__
#define __RANDOM__
#endif

void CenterOfMassMove(double beads[][], int numTimeSlices, int numParticles, double tau, double lam)
{
    int i = 0;
//  int j = 0;
    double delta = 0.75;
    double shift = Shift(delta);
    double oldAction = 0.0;
    double newAction = 0.0;

    do
    {
        oldAction += ComputeAction(**beads, lam, tau, numTimeSlices, numParticles, i);
        i++;
    } while (i < numTimeSlices);

    i = 0;

    do
    {
        beads[i][ptcl] += shift;
        i++;
    } while ( i < numTimeSlices);

    i = 0;

    do
    {
        newAction += ComputeAction(**beads, lam, tau, numTimeSlices, numParticles, i);
        i++;
    } while (i < numTimeSlices);
}

void StagingMove(double beads[][], int numTimeSlices, int numParticles, double tau, double lam, int ptcl)
{
    static std::random_device rd2;
    static std::seed_seq seed2{rd2(), rd2(), rd2(), rd2(), rd2(), rd2(), rd2(), rd2()};
    static std::mt19937 mt2(seed2);
    static std::uniform_int_distribution<int> dist2(0,numTimeSlices+1);
    static std::uniform_real_distribution<double> dist3(0,1);
    static std::normal_distribution<double> nDist(0.0,1.0);

    int i = 0;
    int j = 0;

    double oldAction = 0.0;
    double newAction = 0.0;
    double tau1 = 0.0;
    double avex = 0.0;
    double sigma2 = 0.0;

    int edgeExclusion = 2;
    int m = numTimeSlices - 2*edgeExclusion;

    int alpha_start = dist2(mt2);
    int alpha_end = ModTslice((alpha_start + m), numTimeSlices);
    int temp = 0;

    do
    {
        oldAction += ComputeAction(**beads, lam, tau, numTimeSlices, numParticles, i);
        i++;
    } while (i < numTimeSlices);

    i = 1;
    do
    {
        j = ModTslice((alpha_start + a), numTimeSlices);
        tSlicem1 = ModTslice((j - 1), numTimeSlices);
        tau1 = (m - a) * tau;
        
        avex = (tau1 * beads[tSlicem1][ptcl] + tau * beads[alpha_end][ptcl]) /
            (tau + tau1);
        sigma2 = 2.0 * lam / (1.0 / tau + 1.0 / tau1);
        beads[j][ptcl] = avex + sqrt(sigma2) * nDist(mt2);
        newAction += ComputeAction(**beads, lam, tau, numTimeSlices, numParticles, i);
        i++;
    } while (i < m);

    if (dist3(mt2) < exp(-(newAction - oldAction)))
        // in Julia: Path.numAcceptStaging += 1
    else
        i = 1;
        do
        {
            j = ModTslice((alpha_start + a), numTimeSlices);
            beads[j][ptcl] = oldBeads[a-1]; // Cross-reference with Julia and Python
            i++;
        } while (i < m);
}

void PIMC(int numMCsteps, Paths Path, int iter, int numEquilibSteps, 
        int observableSkip, std::string name, double *energyTrace, 
        double *kineticTrace, double *potentialTrace, double x1_ave, 
        double x2_ave)
{

}
