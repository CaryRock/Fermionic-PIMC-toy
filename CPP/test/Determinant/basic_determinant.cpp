#include <cstdlib>
#include <cstdio>
#include <vector>

#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

static double CalcDeterminant(std::vector<std::vector<double>> Matrix)
{
    double det = 0;

    if (Matrix.size() == 1) {return Matrix[0][0];}
    
    else if (Matrix.size() == 2)
    {
        return (Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0]);
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
        
            det += Matrix[0][p] * pow(-1, p) * CalcDeterminant(tempMatrix);
        }
    }

    return det;
}


int main()
{
    uint64_t N = 10;

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

    double result = CalcDeterminant(matrix);

    printf("The result is: %f\n", result);

    return EXIT_SUCCESS;
}

