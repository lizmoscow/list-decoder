#define HIDE

#include <cmath>
#include <list>
#include "ListDecoder.h"
#include "probabilities.h"


template <class T>
void copy(T** one, T** two, unsigned long length, unsigned long width);
double normalDistribution(double x, double dispersion, double expValue);
unsigned long pop(std::stack<unsigned long>& stack);
void bucketSort (double **array, std::list<probs> *brray, unsigned int n);


ListDecoder::ListDecoder(unsigned long m,
                         unsigned long k,
                         unsigned long l,
                         double standartDeviation,
                         bool *frozens,
                         bool frozenVal = false):
        m(m), k(k), l(l), standartDeviation(standartDeviation), frozenVal(frozenVal) {

    inactivePathIndices = std::stack<unsigned long>();
    activePath = new bool[l];
    arrayPointer_p = new double***[m + 1];
    arrayPointer_c = new bool***[m + 1];
    pathIndexToArrayIndex = new unsigned long*[m + 1];
    inactiveArrayIndices = new std::stack<unsigned long>[m + 1];
    arrayReferenceCount = new unsigned long*[m + 1];

    //initialization
    long length;
    for (unsigned int i = 0; i <= m; ++i) {
        arrayPointer_p[i] = new double**[l];
        arrayPointer_c[i] = new bool**[l];
        pathIndexToArrayIndex[i] = new unsigned long[l];
        arrayReferenceCount[i] = new unsigned long[l];
        length = 1 << (m - i);
        for (unsigned int j = 0; j < l; ++j) {
            arrayPointer_p[i][j] = new double*[length];
            arrayPointer_c[i][j] = new bool*[length];

            for (unsigned long t = 0; k < length; ++k) {
                arrayPointer_p[i][j][t] = new double[2];
                arrayPointer_c[i][j][t] = new bool[2];
            }

            arrayReferenceCount[i][j] = 0;
            inactiveArrayIndices[i].push(j);

            activePath[j] = false;
            inactivePathIndices.push(j);
        }
    }

    for (unsigned long i = 0; i < 1 << m; ++i) {
        this -> frozens[i] = frozens[i];
    }
}


ListDecoder::~ListDecoder() {
    long length;
    for (unsigned int i = 0; i <= m; ++i) {
        length = 1 << (m - i);
        for (unsigned int j = 0; j < l; ++j) {
            for (unsigned long t = 0; t < length; ++t) {
                delete[] arrayPointer_p[i][j][t];
                delete[] arrayPointer_c[i][j][t];
            }
            delete[] arrayPointer_p[i][j];
            delete[] arrayPointer_c[i][j];
        }
        delete[] arrayPointer_p[i];
        delete[] arrayPointer_c[i];
        delete[] pathIndexToArrayIndex[i];
        delete[] arrayReferenceCount[i];
    }
    delete[] arrayPointer_p;
    delete[] arrayPointer_c;
    delete[] pathIndexToArrayIndex;
    delete[] arrayReferenceCount;
    delete[] activePath;
    delete[] inactiveArrayIndices;
}


//takes inactive path, activates it and returns its index
unsigned long ListDecoder::assignInitialPath() {
    unsigned long t = pop(inactivePathIndices);
    activePath[t] = true;

    unsigned long s;
    for (unsigned long i = 0; i <= m; ++i) {
        s = pop(inactiveArrayIndices[i]);
        pathIndexToArrayIndex[i][t] = s;
        arrayReferenceCount[i][s] = 1;
    }
    return t;
}


//takes existing path and duplicates it
unsigned long ListDecoder::clonePath(unsigned long t) {
    unsigned long tt = pop(inactivePathIndices);
    activePath[tt] = true;

    unsigned long s;
    for (unsigned long i = 0; i <= m; ++i) {
        s = pathIndexToArrayIndex[i][t];
        pathIndexToArrayIndex[i][tt] = s;
        ++arrayReferenceCount[i][s];
    }
    return tt;
}


void ListDecoder::killPath(unsigned long t) {
    activePath[t] = false;
    inactivePathIndices.push(t);

    unsigned long s;
    for (unsigned long i = 0; i <= m; ++i) {
        s = pathIndexToArrayIndex[i][t];
        --arrayReferenceCount[i][s];
        if (!arrayReferenceCount[i][s]) {
            inactiveArrayIndices[i].push(s);
        }
    }
}


//returns a pointer to a probability array
//in case of multiple paths sharing probability array, method creates a copy of array and returns a pointer to it
double** ListDecoder::getArrayPointer_p_write(unsigned long lambda, unsigned long t) {
    unsigned long ss, s = pathIndexToArrayIndex[lambda][t];
    if (arrayReferenceCount[lambda][s] == 1){
        ss = s;
    }
    else {
        ss = pop(inactiveArrayIndices[lambda]);
        copy<double>(arrayPointer_p[lambda][s], arrayPointer_p[lambda][ss], 1 << (m - lambda), 2);
        --arrayReferenceCount[lambda][s];
        arrayReferenceCount[lambda][ss] = 1;
        pathIndexToArrayIndex[lambda][t] = ss;
    }
    return arrayPointer_p[lambda][ss];
}

//returns a pointer to a bit array
//in case of multiple paths sharing bit array, method creates a copy of array and returns a pointer to it
bool** ListDecoder::getArrayPointer_c_write(unsigned long lambda, unsigned long t) {
    unsigned long ss, s = pathIndexToArrayIndex[lambda][t];
    if (arrayReferenceCount[lambda][s] == 1){
        ss = s;
    }
    else {
        ss = pop(inactiveArrayIndices[lambda]);
        copy<bool>(arrayPointer_c[lambda][s], arrayPointer_c[lambda][ss], 1 << (m - lambda), 2);
        --arrayReferenceCount[lambda][s];
        arrayReferenceCount[lambda][ss] = 1;
        pathIndexToArrayIndex[lambda][t] = ss;
    }
    return arrayPointer_c[lambda][ss];
}

//returns a pointer to a probability array
double** ListDecoder::getArrayPointer_p_read(unsigned long lambda, unsigned long t) {
    return arrayPointer_p[lambda][pathIndexToArrayIndex[lambda][t]];
}

//returns a pointer to a bit array
bool** ListDecoder::getArrayPointer_c_read(unsigned long lambda, unsigned long t) {
    return arrayPointer_c[lambda][pathIndexToArrayIndex[lambda][t]];
}


void ListDecoder::recursivelyCalc_p(unsigned long lambda, unsigned long fi) {
    if (!lambda) {
        return;
    }
    unsigned long xi = fi / 2;
    if (!(fi % 2)) {
        recursivelyCalc_p(lambda - 1, xi);
    }
    double sigma = 0;
    double **pLambda1, **pLambda2;
    bool **cLambda;
    unsigned long border;
    bool u, uu;
    for (unsigned long i = 0; i < l; ++i) {
        if (!activePath[i]) {
            continue;
        }
        pLambda1 = getArrayPointer_p_write(lambda, i);
        pLambda2 = getArrayPointer_p_read(lambda - 1, i);
        cLambda = getArrayPointer_c_read(lambda, i);
        border = 1 << (m - lambda);
        for (unsigned long j = 0; j < border; ++j) {
            if (fi % 2) {
                u = cLambda[j][0];
                pLambda1[i][false] = pLambda2[2 * i][u] * pLambda2[2 * i + 1][false] / 2;
                pLambda1[i][true] = pLambda2[2 * i][!u] * pLambda2[2 * i + 1][true] / 2;
                sigma = std::max(sigma, pLambda1[i][false], pLambda1[i][true]);
            }
            else {
                pLambda1[j][false] = pLambda2[2 * j][false] * pLambda2[2 * j + 1][false] / 2
                                     + pLambda2[2 * j][true] * pLambda2[2 * j + 1][true] / 2;
                pLambda1[j][true] = pLambda2[2 * j][true] * pLambda2[2 * j + 1][false] / 2
                                    + pLambda2[2 * i][false] * pLambda2[2 * i + 1][true] / 2;
                sigma = std::max(sigma, pLambda1[j][false], pLambda1[j][true]);
            }
        }
    }

    //normalisation
    double **pLambda;
    long length = 1 << (m - lambda);
    for (unsigned long i = 0; i < l; ++i) {
        if (!activePath[i]) {
            continue;
        }
        pLambda = getArrayPointer_p_write(lambda, i);
        for (unsigned long j = 0; j < length; ++j) {
            pLambda[j][0] = pLambda[j][0] / sigma;
            pLambda[j][1] = pLambda[j][1] / sigma;
        }
    }
}


void ListDecoder::recursivelyUpdate_c(unsigned long lambda, unsigned long fi) {
    unsigned long xi = fi / 2;
    bool **c1, **c2;
    long length = 1 << (m - lambda);
    for (unsigned long i = 0; i < l; ++i) {
        if (!activePath[i]) {
            continue;
        }
        c1 = getArrayPointer_c_read(lambda, i);
        c2 = getArrayPointer_c_write(lambda - 2, i);
        for (unsigned long j = 0; j < length; ++j) {
            c2[2 * j][xi % 2] = c1[j][0] ^ c1[j][1];
            c2[2 * j + 1][xi % 2] = c1[j][1];
        }
    }

    if (xi % 2) {
        recursivelyUpdate_c(lambda - 1, xi);
    }
}


void ListDecoder::continuePaths_UnfrozenBit(unsigned long fi) {
    double **probForks = new double*[l];
    bool **contForks = new bool*[l];
    for (unsigned long i = 0; i < l; ++i) {
        probForks[i] = new double[2];
        contForks[i] = new bool[2];
        for (unsigned long j = 0; j < 2; ++j) {
            contForks[i][j] = false;
        }
    }
    double **p;
    unsigned long count = 0;
    for (unsigned long i = 0; i < l; ++i) {
        if (activePath[i]) {
            p = getArrayPointer_p_read(m, i);
            probForks[i][0] = p[0][0];
            probForks[i][1] = p[0][1];
            ++count;
        }
        else {
            probForks[i][0] = -1;
            probForks[i][1] = -1;
        }
    }

    unsigned long ro = std::min(2 * count, l);

#ifndef HIDE
    unsigned long **mins = new unsigned long*[ro];
    double *minProbs = new double[ro];
    for (unsigned long i = 0; i < ro; ++i) {
        mins[i] = new unsigned long[2];
        minProbs[i] = 0;
        for(unsigned long j = 0; j < 2; ++j) {
            mins[i][j] = 0;
        }
    }
    double diff;
    long place;
    for(unsigned long i = 0; i < l; ++i) {
        for (unsigned long j = 0; j < 2; ++j) {
            diff = 0;
            place = -1;
            for (unsigned long t = 0; t < ro; ++t) {
                if (probForks[i][j] - minProbs[t] > diff) {
                    diff = probForks[i][j] - minProbs[t];
                    place = t;
                }
            }
            if (place >= 0) {
                minProbs[place] = probForks[i][j];
                mins[place][0] = i;
                mins[place][1] = j;
            }
        }
    }
    for (unsigned long i = 0; i < ro; ++i) {
        contForks[mins[i][0]][mins[i][1]] = true;
    }
#endif

#ifdef HIDE
    std::list<probs>* brray = new std::list<probs>[2 * l];
    bucketSort(probForks, brray, 2 * l);
    count = 0;
    unsigned long index = 0;
    probs temp;
    while (count != ro) {
        if (brray[index].empty() && index != 2 * l - 1) {
            ++index;
        }
        temp = brray[index].front();
        brray[index].pop_front();
        contForks[temp.getI()][temp.getJ()] = true;
        ++count;
    }
#endif

    //getting rid of non-continuing paths
    for (unsigned long i = 0; i < l; ++i) {
        if(!activePath[i]) {
            continue;
        }
        if (!contForks[i][0] && !contForks[i][1]) {
            killPath(i);
        }
    }

    bool **c;
    unsigned long ii;
    for (unsigned long i = 0; i < l; ++i) {
        if (!contForks[i][0] && !contForks[i][1]) {
            continue;
        }
        c = getArrayPointer_c_write(m, i);
        if (contForks[i][0] && contForks[i][1]) {
            c[0][fi % 2] = false;
            ii = clonePath(i);
            c = getArrayPointer_c_write(m, ii);
            c[0][fi % 2] = true;
        }
        else {
            if (contForks[i][0]) {
                c[0][fi % 2] = false;
            }
            else {
                c[0][fi % 2] = true;
            }
        }
    }
}


bool* ListDecoder::decode(bool *y) {
    unsigned long t = assignInitialPath();
    double **p = getArrayPointer_p_write(0, t);
    bool **c;
    for (unsigned long i = 0; i < 1 << m; ++i) {
        p[i][0] = normalDistribution(y[i], standartDeviation, -1.0);
        p[i][1] = normalDistribution(y[i], standartDeviation, 1.0);
    }
    for (unsigned long fi = 0; fi < 1 << m; ++fi) {
        recursivelyCalc_p(m, fi);
        if (frozens[fi]) {
            for (unsigned long j = 0; j < l; ++j) {
                if (!activePath[j]) {
                    continue;
                }
                c = getArrayPointer_c_write(m, j);
                c[0][fi % 2] = frozenVal;
            }
        }
        else {
            continuePaths_UnfrozenBit(fi);
        }
        if (fi % 2) {
            recursivelyUpdate_c(m, fi);
        }
    }

    unsigned long tt = 0, pp = 0;
    for (unsigned long t = 0; t < l; ++t) {
        if (!activePath[t]) {
            continue;
        }
        c = getArrayPointer_c_read(m, t);
        p = getArrayPointer_p_read(m, t);
        if (pp < p[0][c[0][1]]) {
            tt = t;
            pp = p[0][c[0][1]];
        }
    }

    c = getArrayPointer_c_read(0, tt);
    bool *result = new bool[1 << m];
    for (unsigned long i = 0; i < 1 << m; ++i) {
        result[i] = c[i][0];
    }
    return result;
}


double normalDistribution(double x, double dispersion, double expValue) {
    static const double pi = 3.14159265359;
    return 1 / (sqrt(2 * pi) * dispersion) * exp(-pow((x - expValue) / dispersion, 2) / 2);
}


template <class T>
void copy(T** one, T** two, unsigned long length, unsigned long width) {
    for (unsigned long i = 0; i < length; ++i) {
        for (unsigned long j = 0; j < width; ++j) {
            two[i][j] = one[i][j];
        }
    }
}


unsigned long pop(std::stack<unsigned long>& stack) {
    unsigned long k = stack.top();
    stack.pop();
    return k;
}


void bucketSort (double **array, std::list<probs> *brray, unsigned long n) {
    for (unsigned int i = 0; i < n; i++) {
        brray[(unsigned int)(n * array[i][0])].push_back(probs(array[i][0], i, 0));
        brray[(unsigned int)(n * array[i][1])].push_back(probs(array[i][1], i, 1));
    }
    for (unsigned int i = 0; i < n; i++) {
        brray[i].sort();
    }
}