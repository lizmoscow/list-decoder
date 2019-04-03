//
// Created by Moskovskaya Elizaveta on 13.01.18.
//

#include <cmath>
#include <vector>
#include <random>
#include <chrono>

#define vec2 std::vector<bool>
#define vec std::vector<double>

void addNoise(double standartDeviation, bool *codeword, double *wordWithNoise, unsigned long n) {
    /*std::default_random_engine generator;
    std::normal_distribution<double> distribution(0, standartDeviation);*/
    unsigned int seed = (unsigned int) std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0.0, standartDeviation);
    double g;
    for (int i = 0; i < n; ++i) {
        g = distribution(generator);
        wordWithNoise[i] = (codeword[i] ? 1 : -1) + g;
        /*uint8_t newVal;
        if (normalDistribution(1 + val, standartDeviation) > normalDistribution(val - 1, standartDeviation)) {
            newVal = 1;
        }
        else {
            newVal = 0;
        }
        if (newVal != wordWithNoise[i]) {
            wordWithNoise[i] = 2;
        }*/
    }
}