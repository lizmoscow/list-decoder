#include <iostream>
#include <cmath>
#include <fstream>
#include "ListDecoder.h"

#define GAUSS
//#define PLOT

void z(unsigned long n, unsigned long k, unsigned long *a, double erasureProb, bool test = false);
void findTrustworthyChannels(unsigned long m, unsigned long k, unsigned long *a,
                             double standartDeviation, bool test = false);
void createWord(bool *word, unsigned long k);
void createNewWord(const bool *word, const unsigned long *index, bool *newWord, unsigned long k, unsigned long m);
void encode(const bool *newWord, bool *encoded, unsigned long m);
void addNoise(double standartDeviation, bool *codeword, double *wordWithNoise, unsigned long m);
double function(unsigned long n, unsigned long m, unsigned long k, unsigned long l, double standartDiviation, double erasureProb,
                unsigned long e, unsigned long t, const unsigned long *a, bool test);
template <class T>
void printVec(T *v, unsigned long s);

int main(int argc, const char *argv[]) {
    try {
        if (argc == 9) {
            //reading arguments from command line
            //if arguments are positive, assigning them to 0
            unsigned long m = (std::stoi(argv[1]) > 0) ? std::stoi(argv[1]) : 0;
            unsigned long k = (std::stoi(argv[2]) > 0) ? std::stoi(argv[2]) : 0;
            unsigned long l = (std::stoi(argv[3]) > 0) ? std::stoi(argv[3]) : 0;
            double signalToNoiseRatio = (std::stof(argv[4]));
            double erasureProb = (std::stof(argv[5]));
            unsigned long e = (std::stoi(argv[6]) > 0) ? std::stoi(argv[6]) : 0;
            unsigned long t = (std::stoi(argv[7]) > 0) ? std::stoi(argv[7]) : 0;
            bool test = std::stoi(argv[8]);


            //checking if the arguments are positive
            //if no, quitting the program
            if (m == 0 || k == 0 || l == 0 || signalToNoiseRatio < 0 ||
                    erasureProb > 1 || erasureProb < 0 || e <= 0 || t <= 0) {
                throw("Error: parameters should be positive!");
            }

            unsigned long n = 1 << m;

            unsigned long *a = new unsigned long[k]; //secure channels numbers



#ifndef PLOT
            double standardDiviation = sqrt(1 / (pow(10, signalToNoiseRatio / 10) * 2 * k / n));
            if (!test) std::cout << "The most trustworthy channels: ";
#ifndef GAUSS
            z(m, k, a, erasureProb, test); //initializing a
#endif
#ifdef GAUSS
            findTrustworthyChannels(m, k, a, standardDiviation, test);
#endif


            double error = function(n, m, k, l, standardDiviation, erasureProb, e, t, a, test);
            std::cout << "Standard Deviation: " << standardDiviation;
            if (!test) std::cout << std::endl << std::string(40, '*');
            std::cout << std::endl << "Error frequency: " << error << std::endl;
#endif PLOT

#ifdef PLOT

#ifndef GAUSS
            if (!test) std::cout << "The most trustworthy channels: ";
            z(m, k, a, erasureProb, test); //initializing a
#endif

#ifdef GAUSS
            if (!test) std::cout << "The most trustworthy channels: ";
            findTrustworthyChannels(m, k, a, sqrt(1 / (pow(10, 2 / 10) * 2 * k / n)), test);
#endif
            std::fstream fout;
            fout.open("plot.csv");
            fout << "ratio,errors\n";
            double error;
            double r = 0;
            double standardDiviation;
            for (unsigned long i = 0; i <= 50; i++) {
                standardDiviation = sqrt(1 / (pow(10, r / 10) * 2 * k / n));
                error = function(n, m, k, l, standardDiviation, erasureProb, e, t, a, test);
                std::cout << "Standard Deviation: " << standardDiviation;
                if (!test) std::cout << std::endl << std::string(40, '*');
                std::cout << std::endl << "Error frequency: " << error << std::endl;
                fout << r << "," << error << "\n";
                r += 0.1;
            }
            fout.close();
#endif PLOT
            delete[] a;
        }
    }
    catch(const char* e) {
        std::cerr << e;
    }

    return 0;
}


template <class T>
void printVec(T *v, unsigned long s) {
    for (unsigned long i = 0; i < s; ++i) {
        std::cout << v[i] << ' ';
    }
    std::cout << '\n';
}


double function(unsigned long n, unsigned long m, unsigned long k, unsigned long l, double standartDiviation, double erasureProb,
                unsigned long e, unsigned long t, const unsigned long *a, bool test) {
    unsigned long count = 0;
    unsigned long errCount = 0;
    bool *word = new bool[k];
    bool *newWord = new bool[n];
    bool *codeword = new bool[n];
    double *codewordWithNoise = new double[n];
    bool *result = new bool[n];
    ListDecoder decoder(m, k, l, standartDiviation, a);
    while (count != t && (errCount != e || test)) {
        createWord(word, k);
        createNewWord(word, a, newWord, k, n);
        for (unsigned long i = 0; i < n; ++i) {
            codeword[i] = newWord[i];
        }
        encode(newWord, codeword, n);

        addNoise(standartDiviation, codeword, codewordWithNoise, n);
        decoder.decode(codewordWithNoise, result);
        if (!test) {
            std::cout << "Random word: ";
            printVec(word, k);
            std::cout << "Passed word: ";
            printVec(newWord, n);
            std::cout << "Encoded word: ";
            printVec(codeword, n);
            std::cout << "Word with added noise: ";
            printVec(codewordWithNoise, n);
            std::cout << "Decoded word: ";
            printVec(result, n);
            std::cout << std::endl;
        }

        for (unsigned long i = 0; i < n; ++i) {
            if (result[i] != codeword[i]) {
                ++errCount;
                break;
            }
        }

        decoder.restart();
        ++count;
    }
    delete[] word;
    delete[] newWord;
    delete[] codeword;
    delete[] codewordWithNoise;
    delete[] result;
    return ((double)errCount) / ((double)count);
}