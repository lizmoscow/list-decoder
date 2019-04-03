#include <utility>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#define pairid std::pair<double, unsigned long>

long partition(pairid *array, long p, long r);
void quicksubsort(pairid* array, long p, long r);
void quicksort(pairid* array, long n);

int comp (const void *a, const void *b);

double bigXi(double x);

void z(unsigned long loga, unsigned long k, unsigned long *a, double erasureProb, bool test = false) {
    unsigned long n = (unsigned long)pow(2, loga);
    unsigned long powerof2;
    double *zw2 = new double[n];
    double *zw = new double[n];
    zw[0] = erasureProb;
    for (unsigned long i = 1; i <= loga; ++i) {
        powerof2 = (unsigned long)pow(2, i);
        for (unsigned long j = 0; j < powerof2; ++j) {
            if (j % 2) {
                zw2[j] = pow(zw[(j - 1) / 2], 2);
            }
            else {
                zw2[j] = 2 * zw[j / 2] - pow(zw[j / 2], 2);
            }
        }
        std::swap(zw, zw2);
    }
    pairid *zwArray = new pairid[n];
    for (unsigned long i = 0; i < n; ++i) {
        zwArray[i] = pairid(zw[i], i);
    }
    //std::sort(zwArray, zwArray + n, comp);
    quicksort(zwArray, n);

    for (unsigned long i = 0; i < k; ++i) {
        a[i] = zwArray[i].second;
        if (!test) {
            std::cout << zwArray[i].second << " (" << zwArray[i].first << ")";
            if (i != k - 1) {
                std::cout << ", ";
            }
        }
    }
    std::sort(a, a + k);
    if (!test) {
        std::cout << '\n';
    }
    delete[] zw;
    delete[] zw2;
    delete[] zwArray;
}


int comp (const void *a, const void *b) {
    return ( (*(const pairid*)a).first < (*(const pairid*)b).first );
}


long partition(pairid *array, long p, long r) {
    double x = array[r].first;
    long i = p - 1;
    for (long j = p; j < r; j++)
        if (array[j].first < x)
            std::swap(array[++i], array[j]);
    std::swap(array[i + 1], array[r]);
    return i + 1;
}


void quicksubsort(pairid* array, long p, long r) {
    if (p < r)
    {
        long q = partition(array, p, r);
        quicksubsort(array, p, q - 1);
        quicksubsort(array, q + 1, r);
    }
    return;
}


void quicksort(pairid* array, long n) {
    quicksubsort(array, 0, n - 1);
    return;
}


void findTrustworthyChannels(unsigned long m, unsigned long k, unsigned long *a,
                             double standartDeviation, bool test = false) {
    unsigned long n = (1 << m);
    unsigned long powerof2 = 1;
    double *l2 = new double[n];
    double *l = new double[n];
    l[0] = 2 / standartDeviation / standartDeviation;
    for (unsigned long i = 1; i <= m; ++i) {
        powerof2 <<= 1;
        for (unsigned long j = 0; j < powerof2; ++j) {
            if (j % 2) {
                l2[j] = 2 * l[(j - 1) / 2];
            }
            else {
                l2[j] = bigXi(l[j / 2]);
            }
        }
        std::swap(l, l2);
    }
    pairid *larray = new pairid[n];
    for (unsigned long i = 0; i < n; ++i) {
        larray[i] = pairid(l[i], i);
    }
    quicksort(larray, n);

    for (unsigned long i = 0; i < k; ++i) {
        a[i] = larray[n - i - 1].second;
        if (!test) {
            std::cout << larray[n - i - 1].second << " (" << larray[n - i - 1].first << ")";
            if (i != k - 1) {
                std::cout << ", ";
            }
        }
    }
    std::sort(a, a + k);
    if (!test) {
        std::cout << '\n';
    }
    delete[] l;
    delete[] l2;
    delete[] larray;
}


double bigXi(double x) {
    if (x > 12) {
        return 0.9861 * x - 2.3152;
    }
    if ( x > 3.5 && x <= 12) {
        return x * (9.005 * 0.001 * x + 0.7694) - 0.9507;
    }
    if (x > 1 && x <= 3.5) {
        return x * (0.062883 * x + 0.3678) - 0.1627;
    }
    return x * (0.2202 * x + 0.06448);
}