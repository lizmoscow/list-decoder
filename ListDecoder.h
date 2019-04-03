//
// Created by Moskovskaya Elizaveta on 28.02.18.
//

#ifndef LIST_DECODER_LISTDECODER_H
#define LIST_DECODER_LISTDECODER_H


#include <stack>
#include <list>
#include "probabilities.h"

class ListDecoder {

private:
    unsigned long m;
    unsigned long k;
    unsigned long l;
    double standartDeviation;
    unsigned long *unfrozens; //a sorted array of unfrozen bit indexes
    bool frozenVal;

    std::stack<unsigned long> inactivePathIndices;
    bool* activePath;
    double *arrayPointer_p;
    bool *arrayPointer_c;
    unsigned long **pathIndexToArrayIndex;
    std::stack<unsigned long> *inactiveArrayIndices;
    unsigned long *arrayReferenceCount;
    probabilities *probForks;
    int *contForks;

    unsigned long *layerBeginnings;

    unsigned long assignInitialPath();
    unsigned long clonePath(unsigned long t);
    void killPath(unsigned long t);
    double* getArrayPointer_p_read(unsigned long lambda, unsigned long t);
    bool* getArrayPointer_c_read(unsigned long lambda, unsigned long t);
    double* getArrayPointer_p_write(unsigned long lambda, unsigned long t);
    bool* getArrayPointer_c_write(unsigned long lambda, unsigned long t);
    void recursivelyCalc_p(unsigned long lambda, unsigned long fi);
    void recursivelyUpdate_c(unsigned long lambda, unsigned long fi);
    void continuePaths_UnfrozenBit(unsigned long fi);

    unsigned long getIndex(unsigned long s, unsigned long i) const;

public:
    ListDecoder(unsigned long m, unsigned long k, unsigned long l, double standartDeviation,
                const unsigned long /*const*/ *unfrozens, bool frozenVal = false);
    ~ListDecoder();
    void decode(double *y, bool *res);
    void restart();
};


#endif //LIST_DECODER_LISTDECODER_H
