//
// Created by Moskovskaya Elizaveta on 28.02.18.
//


#include <cmath>
#include <list>
#include <algorithm>
#include "ListDecoder.h"
#include "probabilities.h"


template <class T>
void copy(T* one, T* two, unsigned long length);
double normalDistribution(double x, double dispersion, double expValue);
unsigned long myPop(std::stack<unsigned long> &stack);
void bucketSort (double **array, std::list<probabilities> *brray, unsigned long n);


ListDecoder::ListDecoder(unsigned long m,
                         unsigned long k,
                         unsigned long l,
                         double standartDeviation,
                         const unsigned long /*const*/ *unfrozens,
                         bool frozenVal):
        m(m), k(k), l(l), standartDeviation(standartDeviation), frozenVal(frozenVal) {

    unsigned long size = ((1 << (m + 1)) - 1) * l * 2;

    inactivePathIndices = std::stack<unsigned long>();
    activePath = new bool[l];
    arrayPointer_p = new double[size];
    arrayPointer_c = new bool[size];
    pathIndexToArrayIndex = new unsigned long*[m + 1];
    inactiveArrayIndices = new std::stack<unsigned long>[m + 1];
    arrayReferenceCount = new unsigned long[(m + 1) * l];
    layerBeginnings = new unsigned long[m + 1];
    probForks = new probabilities[2 * l];
    contForks = new int[l];

    //initialization
    unsigned long length;
    for (int i = 0; i <= m; ++i) {
        pathIndexToArrayIndex[i] = new unsigned long[l];
        length = 1 << (m - i + 1);
        for (unsigned int j = 0; j < l; ++j) {
            arrayReferenceCount[i * l + j] = 0;
            inactiveArrayIndices[i].push(j * length + ((i - 1 < 0) ? 0 : layerBeginnings[i - 1]));
        }
        layerBeginnings[i] = l * length + ((i - 1 < 0) ? 0 : layerBeginnings[i - 1]);
    }

    for (unsigned long i = 0; i < size; ++i) {
        arrayPointer_c[i] = false;
    }

    for (unsigned long i = 0; i < l; ++i) {
        activePath[i] = false;
        inactivePathIndices.push(i);
        for (unsigned long p = 0; p < (1 << m); ++p) {
            arrayPointer_c[i * l + 2 * p + 1] = 0;
        }
    }

    this -> unfrozens = new unsigned long[k];
    for (unsigned long i = 0; i < k; ++i) {
        this -> unfrozens[i] = unfrozens[i];
    }
}


ListDecoder::~ListDecoder() {
    for (unsigned int i = 0; i <= m; ++i) {
        delete[] pathIndexToArrayIndex[i];
    }
    delete[] unfrozens;
    delete[] layerBeginnings;
    delete[] arrayPointer_p;
    delete[] arrayPointer_c;
    delete[] pathIndexToArrayIndex;
    delete[] arrayReferenceCount;
    delete[] activePath;
    delete[] inactiveArrayIndices;
}


void ListDecoder::restart() {
    unsigned long size = ((1 << (m + 1)) - 1) * l * 2;
    unsigned long length;
    for (int i = 0; i <= m; ++i) {
        length = 1 << (m - i + 1);
        while (!inactiveArrayIndices[i].empty()) {
            inactiveArrayIndices[i].pop();
        }
        for (unsigned int j = 0; j < l; ++j) {
            arrayReferenceCount[i * l + j] = 0;
            inactiveArrayIndices[i].push(j * length + ((i - 1 < 0) ? 0 : layerBeginnings[i - 1]));
        }
        layerBeginnings[i] = l * length + ((i - 1 < 0) ? 0 : layerBeginnings[i - 1]);
    }

    for (unsigned long i = 0; i < size; ++i) {
        arrayPointer_c[i] = false;
    }

    while (!inactivePathIndices.empty()) {
        inactivePathIndices.pop();
    }
    for (unsigned long i = 0; i < l; ++i) {
        activePath[i] = false;
        inactivePathIndices.push(i);
        for (unsigned long p = 0; p < (1 << m); ++p) {
            arrayPointer_c[i * l + 2 * p + 1] = 0;
        }
    }
}


//takes inactive path, activates it and returns its index
unsigned long ListDecoder::assignInitialPath() {
    unsigned long t = myPop(inactivePathIndices);
    activePath[t] = true;

    unsigned long s;
    for (unsigned long i = 0; i <= m; ++i) {
        s = myPop(inactiveArrayIndices[i]);
        pathIndexToArrayIndex[i][t] = s;
        arrayReferenceCount[i * l + getIndex(s, i)] = 1;
    }
    return t;
}


//takes existing path and duplicates it
unsigned long ListDecoder::clonePath(unsigned long t) {
    unsigned long tt = myPop(inactivePathIndices);
    activePath[tt] = true;

    unsigned long s;
    for (unsigned long i = 0; i <= m; ++i) {
        s = pathIndexToArrayIndex[i][t];
        pathIndexToArrayIndex[i][tt] = s;
        ++arrayReferenceCount[i * l + getIndex(s, i)];
    }
    return tt;
}


void ListDecoder::killPath(unsigned long t) {
    activePath[t] = false;
    inactivePathIndices.push(t);

    unsigned long s;
    for (unsigned long i = 0; i <= m; ++i) {
        s = pathIndexToArrayIndex[i][t];
        --arrayReferenceCount[i * l + getIndex(s, i)];
        if (!arrayReferenceCount[i * l + getIndex(s, i)]) {
            inactiveArrayIndices[i].push(s);
        }
    }
}


//returns a pointer to a probability array
//in case of multiple paths sharing probability array, method creates a copy of array and returns a pointer to it
double* ListDecoder::getArrayPointer_p_write(unsigned long lambda, unsigned long t) {
    unsigned long ss, s = pathIndexToArrayIndex[lambda][t];
    if (arrayReferenceCount[lambda * l + getIndex(s, lambda)] == 1){
        ss = s;
    }
    else {
        ss = myPop(inactiveArrayIndices[lambda]);
        copy<double>(arrayPointer_p + s, arrayPointer_p + ss, 1 << (m - lambda + 1));
        --arrayReferenceCount[lambda * l + getIndex(s, lambda)];
        arrayReferenceCount[lambda * l + getIndex(ss, lambda)] = 1;
        pathIndexToArrayIndex[lambda][t] = ss;
    }
    return arrayPointer_p + ss;
}


//returns a pointer to a bit array
//in case of multiple paths sharing bit array, method creates a copy of array and returns a pointer to it
bool* ListDecoder::getArrayPointer_c_write(unsigned long lambda, unsigned long t) {
    unsigned long ss, s = pathIndexToArrayIndex[lambda][t];
    if (arrayReferenceCount[lambda * l + getIndex(s, lambda)] == 1){
        ss = s;
    }
    else {
        ss = myPop(inactiveArrayIndices[lambda]);
        copy<bool>(arrayPointer_c + s, arrayPointer_c + ss, 1 << (m - lambda + 1));
        --arrayReferenceCount[lambda * l + getIndex(s, lambda)];
        arrayReferenceCount[lambda * l + getIndex(ss, lambda)] = 1;
        pathIndexToArrayIndex[lambda][t] = ss;
    }
    return arrayPointer_c + ss;
}


//returns a pointer to a probability array
double* ListDecoder::getArrayPointer_p_read(unsigned long lambda, unsigned long t) {
    return arrayPointer_p + pathIndexToArrayIndex[lambda][t];
}


//returns a pointer to a bit array
bool* ListDecoder::getArrayPointer_c_read(unsigned long lambda, unsigned long t) {
    return arrayPointer_c + pathIndexToArrayIndex[lambda][t];
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
    double *pLambda1, *pLambda2;
    bool *cLambda;
    unsigned long border;
    int u;
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
                u = (!cLambda[2 * j]) ? 0 : 1;
                pLambda1[2 * j] = pLambda2[2 * j + u] * pLambda2[2 * (j + border)];

                pLambda1[2 * j + 1] = pLambda2[2 * j + (u + 1) % 2] * pLambda2[2 * (j + border) + 1];
                sigma = std::max(sigma, std::max(pLambda1[2 * j], pLambda1[2 * j + 1]));
            }
            else {
                pLambda1[2 * j] = pLambda2[2 * j] * pLambda2[2 * (j + border)]
                                 + pLambda2[2 * j + 1] * pLambda2[2 * (j + border) + 1];

                pLambda1[2 * j + 1] = pLambda2[2 * j + 1] * pLambda2[2 * (j + border)]
                                 + pLambda2[2 * j] * pLambda2[2 * (j + border) + 1];
                sigma = std::max(sigma, std::max(pLambda1[2 * j], pLambda1[2 * j + 1]));
            }
        }
    }

    //normalisation
    double *pLambda;
    long length = 1 << (m - lambda);
    for (unsigned long i = 0; i < l; ++i) {
        if (!activePath[i]) {
            continue;
        }
        pLambda = getArrayPointer_p_write(lambda, i);
        for (unsigned long j = 0; j < length; ++j) {
            pLambda[2 * j] /= sigma;
            pLambda[2 * j + 1] /= sigma;
        }
    }
}


void ListDecoder::recursivelyUpdate_c(unsigned long lambda, unsigned long fi) {
    unsigned long xi = fi / 2;
    bool *c1, *c2;
    long length = 1 << (m - lambda);
    for (unsigned long i = 0; i < l; ++i) {
        if (!activePath[i]) {
            continue;
        }

        c1 = getArrayPointer_c_read(lambda, i);

        c2 = getArrayPointer_c_write(lambda - 1, i);
        for (unsigned long j = 0; j < length; ++j) {
            c2[2 * j + xi % 2] = c1[2 * j] ^ c1[2 * j + 1];
            c2[2 * (j + length) + xi % 2] = c1[2 * j + 1];
        }
    }

    if (xi % 2) {
        recursivelyUpdate_c(lambda - 1, xi);
    }
}


void ListDecoder::continuePaths_UnfrozenBit(unsigned long fi) {
    for (unsigned long i = 0; i < l; ++i) {
        contForks[i] = 0;
    }
    double *p;
    unsigned long count = 0;
    for (unsigned long i = 0; i < l; ++i) {
        if (activePath[i]) {

            p = getArrayPointer_p_read(m, i);

            probForks[i * 2].set(p[0], i, 1);
            probForks[i * 2 + 1].set(p[1], i, 2);
            ++count;
        }
        else {
            probForks[i * 2].setNull();
            probForks[i * 2 + 1].setNull();
        }
    }
    unsigned long ro = std::min(2 * count, l);

    std::sort(probForks, probForks + 2 * l);
    unsigned long index = 2 * l - 1;
    probabilities temp;
    while (index != 2 * l - 1 - ro) {
        temp = probForks[index];
        contForks[temp.getI()] |= temp.getJ();
        --index;
    }

    //getting rid of non-continuing paths
    for (unsigned long i = 0; i < l; ++i) {
        if(!activePath[i]) {
            continue;
        }
        if (!(contForks[i] & 1) && !(contForks[i] & 2)) {
            killPath(i);
        }
    }

    bool *c;
    unsigned long ii;
    for (unsigned long i = 0; i < l; ++i) {
        if (!(contForks[i] & 1) && !(contForks[i] & 2)) {
            continue;
        }
        c = getArrayPointer_c_write(m, i);
        if (contForks[i] & 1 && contForks[i] & 2) {
            c[fi % 2] = false;
            ii = clonePath(i);
            c = getArrayPointer_c_write(m, ii);
            c[fi % 2] = true;
        }
        else {
            if (contForks[i] & 1) {
                c[fi % 2] = false;
            }
            else {
                c[fi % 2] = true;
            }
        }
    }
}


unsigned long ListDecoder::getIndex(unsigned long s, unsigned long i) const {
    return (s - (((int)i - 1 < 0) ? 0 : layerBeginnings[i - 1])) / (1 << (m - i + 1));
}


void ListDecoder::decode(double *y, bool *res) {
    unsigned long t = assignInitialPath();
    double *p = getArrayPointer_p_write(0, t);
    bool *c;
    int count = 0;
    for (unsigned long i = 0; i < 1 << m; ++i) {
        p[2 * i] = normalDistribution(y[i], standartDeviation, -1.0);
        p[2 * i + 1] = normalDistribution(y[i], standartDeviation, 1.0);
    }
    for (unsigned long fi = 0; fi < 1 << m; ++fi) {
        recursivelyCalc_p(m, fi);
        if (fi != unfrozens[count] && count < k) {
            for (unsigned long j = 0; j < l; ++j) {
                if (!activePath[j]) {
                    continue;
                }
                c = getArrayPointer_c_write(m, j);
                c[fi % 2] = frozenVal;
            }
        }
        else if (count >= k) {
           throw "An error occured with an index of unfrozens array";
        }
        else {
            continuePaths_UnfrozenBit(fi);
            ++count;
        }
        if (fi % 2) {
            recursivelyUpdate_c(m, fi);
        }
    }

    unsigned long tt = 0;
    double pp = 0;
    for (unsigned long t = 0; t < l; ++t) {
        if (!activePath[t]) {
            continue;
        }

        c = getArrayPointer_c_read(m, t);
        p = getArrayPointer_p_read(m, t);

        if (pp < p[c[1]]) {
            tt = t;
            pp = p[c[1]];
        }
    }

    c = getArrayPointer_c_read(0, tt);

    for (unsigned long i = 0; i < (1 << m); ++i) {
        res[i] = c[2 * i];
    }
}


double normalDistribution(double x, double dispersion, double expValue) {
    static const double pi = 3.14159265359;
    return 1 / (sqrt(2 * pi) * dispersion) * exp(-pow((x - expValue) / dispersion, 2) / 2);
}


template <class T>
void copy(T* one, T* two, unsigned long length) {
    for (unsigned long i = 0; i < length; ++i) {
            two[i] = one[i];
    }
}


unsigned long myPop(std::stack<unsigned long> &stack) {
    unsigned long k = stack.top();
    stack.pop();
    return k;
}