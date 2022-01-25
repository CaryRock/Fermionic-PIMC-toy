#include <cstdlib>
#include <cstdio>
#include <vector>

#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <chrono>

static double CalcPermanant(std::vector<std::vector<double>> Matrix)
{
    double perm = 0;

    if (Matrix.size() == 1) {return Matrix[0][0];}
    
    else if (Matrix.size() == 2)
    {
        return (Matrix[0][0] * Matrix[1][1] + Matrix[0][1] * Matrix[1][0]);
    }

    else
    {
        for (uint64_t p = 0; p < Matrix[0].size(); p++)
        {
            std::vector<std::vector<double>> tempMatrix;

            for (uint64_t i = 1; i < Matrix.size(); i++)
            {
                std::vector<double> tempRow;

                for (uint64_t j = 0; j < Matrix[i].size(); j++)
                {
                    if (j != p) { tempRow.push_back(Matrix[i][j]); }
                }

                if (tempRow.size() > 0) { tempMatrix.push_back(tempRow); }
            }
        
            perm += Matrix[0][p] * CalcPermanant(tempMatrix);
        }
    }

    return perm;
}


int main()
{
    uint64_t N = 2;

    boost::random::random_device rd;
    boost::random::mt19937_64 rng(rd);
    boost::random::uniform_real_distribution<double> dist(0,1);

    std::vector<std::vector<double>> matrix(N, std::vector<double>(N));

    for (uint64_t i = 0; i < N; i++)
    {
        for (uint64_t j = 0; j < N; j++)
        {
            matrix[i][j] = dist(rng);
        }
    }
    
    auto t1 = std::chrono::high_resolution_clock::now();
    double result = CalcPermanant(matrix);
    auto t2 = std::chrono::high_resolution_clock::now();

    printf("The result is: %f\n", result);
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;
    printf("The execuation time of the permanent was: %f ms\n", ms_double.count());

    return EXIT_SUCCESS;
}

