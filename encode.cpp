//
// Created by Moskovskaya Elizaveta on 13.01.18.
//
#include <vector>
#include <random>

#define vec2 std::vector<bool>

void w(unsigned long i, unsigned long j, vec2& encoded);

void createWord(bool *word, unsigned long k) {
    /*unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(-1, 1);
    for (int i = 0; i < word.size(); ++i) {
        if (distribution(generator) <= 0) {
            word[i] = 0;
        }
        else {
            word[i] = 1;
        }
    }*/
    static std::default_random_engine generator;
    static std::uniform_int_distribution<int> rv(0, 1);
    for(unsigned int i=0; i< k; ++i)
    {
        word[i] = (bool) rv(generator);
    }
}

void createNewWord(const bool *word, const unsigned long *index,
                   bool *newWord, unsigned long k, unsigned long n) {
    int count = 0;
    for (unsigned long i = 0; i < n; ++i) {
        if (i == index[count]) {
            newWord[i] = word[count];
            ++count;
        }
        else {
            newWord[i] = 0;
        }
    }
}

void encode(const bool *newWord, bool *encoded, unsigned long n) {
    unsigned long start;
    unsigned long end;
    for (unsigned long i = 2; i <= n; i *= 2) {
        for (unsigned long j = 0; j < n; j += i) {
            start = j + i - 1;
            long end = j + i / 2 - 1;
            for (long k = start; k > end; --k) {
                encoded[k - i / 2] = encoded[k] ^ encoded[k - i / 2];
            }
        }
    }
}

void w(unsigned long n, unsigned long m, vec2& encoded) {
    unsigned long start = m + n - 1;
    unsigned long end = m + n / 2 - 1;
    for (long i = start; i > end; --i) {
        encoded[i - n / 2] = encoded[i] ^ encoded[i - n / 2];
    }
}