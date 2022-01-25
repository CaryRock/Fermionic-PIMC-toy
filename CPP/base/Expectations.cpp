double X1_expectation(double beads[][], int numTimeSlices, int numParticles)
{
    double x1_av = 0.0;

    for(int i = 0; i < numTimeSlices; i++)
    {
        for(int j = 0; j < numParticles; j++)
        {
            x1_av += beads[i][j];
        }
    }

    x1_av /= numTimeSlices;

    return x1_av;
}

double X2_expectation(double beads[][], int numTimeSlices, int numparticles)
{
    double x2_av = 0.0;

    for(int i = 0; i < numTimeSlices; i++)
    {
        for(int j = 0; j < numParticles; j++)
        {
            x2_av += (beads[i][j] * beads[i][j]);
        }
    }

    x2_av /= numTimeSlices;
    
    return x2_av;
}
