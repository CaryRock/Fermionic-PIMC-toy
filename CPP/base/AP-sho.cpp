#include "common_headers.hpp"
#include <sstream>  // Do I need this?
#include "AP_header.hpp"

#include <boost/program_options.hpp>

using po = boost::program_options;
using br = boost::random;
using bu = boost::uuids;

using std::string;

struct Params
{
    uint64_t nPar;              // Number of particles
    uint64_t nTsl;              // Number of time slices
    double lam;                 // System parameters - hbar^2/(2mk_B)  = 1/2
    double tau;                 // beta/J (beta/M in Ceperley-ese)
    double x_a;                 // Arbitrary left edge
    double x_b;                 // Arbitrary right edge
    double delta;               // Width of possible shift for a bead
    uint64_t numHistBins;       // Number of bins for histrogram
    uint64_t numBucks;          // Number of MC steps to bin (number of buckets)
    // The following two should be specific to relevant estimators, not global
//    double binWidth;
//    uint64_t numMCbins;
    uint64_t numEquilibSteps;   // Number of steps to skip for equilibriation
    uint64_t observableSkip;    // Number of MC steps to skip between observations
    uint64_t numSamples;        // Sets the # of MC steps in total
    string baseName;            // Base name of the data file(s)
    string uuid;                // UUID for unique identification

    br::mt19937_64 rng;
    br::uniform_real_distribution dist01;

    Params(uint64_t nPar, uint64_t nTsl, double lam, double tau, double x_a, 
            double x_b, double delta, uint64_t numHistBins, uint64_t numBucks, 
            uint64_t numEquilibSteps, uint64_t observableSkip, 
            uint64_t numSamples, string baseName, string uuid);
    
    ~Params();

    double GetRNG01(double delta=0.75);

};

Params::Params(uint64_t nP, uint64_t nT, double lm, double tu, double xa, 
                double xb, double dlt, uint64_t nHstBns, uint64_t nmBcks,
                uint64_t nmEqlbStps, uint64_t obsrvblSkp, uint64_t nmSmpls, 
                string bsNm, string uid)
{
    nPar = nP;
    nTsl = nT;
    lam = lm;
    tau = tu;
    x_a = xa;
    x_b = xb;
    delta = dlt;
    numHistBins = nHstBns;
    numBucks = nmBcks;
    numEquilibSteps = nmEqlbStps;
    numSamples = nmSmpls;
    baseName = bsNm;
    uuid = uid;

    rng = engine;
    dist01 = dist;
}

Params::~Params() {}

inline double Params::GetRNG01(double delta=0.75)
{
    return dist01(rng);
}

struct Paths
{
    double **beads;
    double **determinants;  // Shape is (numTimeSlices, numParticles)
    double **potentials;
    double KE = 0;
    double PE = 0;
    uint64_t numAcceptCOM = 0;
    uint64_t numAcceptStaging = 0;

    Paths(double **beads, double **determinants, double **potentials);
};

// This is wrong - check out the Phys642_Final_Project code for reference
Paths::Paths(double **bds, double **determs, double **potens)
{
    **beads = **bds;
    **determinants = **determs;
    **potentials = **potens;
}

Paths::~Paths() {}

inline double ExtPotential(double R, double lam)
{
    return 1/(4 * lam) * R * R;
}

inline uint64_t ModTslice(uint64_t tSlice, uint64_t numSlices)
{
    return ( ( tSlice - 1 + numSlices ) % numSlices );
}

inline bool CutOff(double beadA, double beadB, double cutoff=pow(10.0,-12))
{
    return (abs(beadA - beadB) < cutoff);
}

double Shift(double delta=0.75)
{
    static std::random_device rd;
    static std::seed_seq seed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    static std::mt19937 mt(seed);
    static std::uniform_real_distribution<double> dist(-1,1);

    return delta*dist(mt);
}

// TODO: Implement command line parsing - use Boost


int WriteHeader(string fileName, string writeOut)
{
    int closed = -1;

    const char *name = fileName.c_str();
    const char *data = writeOut.c_str();

    FILE *outputFile;
    outputFile = fopen(name, "a");
    fprintf(outputFile, "%s", data);
    closed = fclose(outputFile);

    return closed;
}

// TODO: Implement WriteMat()   // I don't think I use this, though
// Or don't. If I need it, it'll get added.

// === End function definitions ===============================================

// === Begin main implementation ==============================================

//#include "./Energy.cpp"
//#include "./Action.cpp"
//#include "./PIMC.cpp"
//#include "./Expectations.cpp"

int main(int argc, char *argv[])
{
    std::stringstream tempArg{ argv[1] };   // Temperature (Kelvin)
    std::stringstream nParArg{ argv[2] };   // numParticles
    std::stringstream nTslArg{ argv[3] };   // numTimeSlices
    std::stringstream nEquArg{ argv[4] };   // numEquilibSteps
    std::stringstream obSkArg{ argv[5] };   // observableSkip
    std::stringstream nSamArg{ argv[6] };   // numSamples

    boost::uuids::uuid u = boost::uuids::random_generator()();
    std::stringstream ss;
    ss << u;    // Maybe I should just use cout instead...
    string uuidName = ss.str();     // UUID for run

    double lam = 1/2;                   // m = hbar^2/k_b
    double temp << tempArg;             // Simulation temperature
    uint64_t nPar << nParArg;           // Number of particles in simulation
    uint64_t nTsl << nTslArg;           // Number of time slices to use
    uint64_t numEquilibSteps << nEquArg;// Number of steps to equilibriate over
    uint64_t observableSkip << obSkArg; // Number of MC steps between observations
    uint64_t numSamples << nSamArg;     // Number of samples to write out by end
    double delta = 1.0;                 // TODO: MAKE THIS READ FROM COMMAND LINE ARGUMENTS
    double x_min = -10.0;               // TODO: MAKE THIS READ FROM COMMAND LINE ARGUMENTS
    double x_max = 10.0;                // TODO: MAKE THIS READ FROM COMMAND LINE ARGUMENTS

    uint64_t numMCsteps = numEquilibSteps + observableSkip * numSamples;
    double tau          = 1/(nTsl * temp);
    int binSize         = 50;   // # of MC steps to bin simulation over TODO: MAKE THIS READ FROM COMMAND LINE ARGUMENTS

    printf("Simulation parameters:\n");
    printf("N\t\t= %d\n", nPar);
    printf("lambda\t\t= %f\n", lam);
    printf("Temperature\t= %f\n", temp);
    printf("numMCsteps\t= %ld\n", numMCsteps);

    return 0;
}
