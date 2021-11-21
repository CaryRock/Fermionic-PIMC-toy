#define __SHO_H__

struct Params;
struct Paths;

double ExtPotential(double R, double lam);
uint64_t ModTslice(uint64_t tSlice, uint64_t numSlices);
bool CutOff(double beadA, double beadB, double cutOff);
double Shift(double delta);
// TODO: Implement command line parsing's prototype(s)
//void CopyBeads(double **source, double **dest, int rows, int cols);
//double PotentialAction(double beadVector[], double tau, int numParticles);
//double KineticEnergy(Paths &Path, double tau);
//double PotentialEnergy(Paths &Path);
//double Energy(Paths &Path, double tau);
//double Determinant(double beads[][], int numTimeSlices, int numParticles, int tslice, int iter, double lam);
//void InstantiatePotentials(double beads[][], double potentials[][], int numTimeSlices, int numParticles, double lam, int iter);
//double ComputeAction(double beads[][], double lam, double tau, int numTimeSlices, int numParticles, int tslice);
//void PIMC(Paths &Path, int numSteps, int iter, int numEquilStep, 
//          int observSkip, std::string name);
//void CenterOfMassMove(double beads[][], int numTimeSlices, int numParticles, 
//          double tau, double lam);
//void BinData(double *energyTrace, int binSize, int numTimeSlices, 
//          int numParticles, double temp, int equil, int obs, int numTimeSlices, std::string file_name);
