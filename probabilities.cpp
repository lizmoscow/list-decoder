//
// Created by Moskovskaya Elizaveta on 01.03.18.
//

#include "probabilities.h"

probabilities::probabilities(double x, unsigned long y, unsigned long z) :
        prob(x), i(y), j(z) { };

probabilities::probabilities(const probabilities &other) {
    this -> prob = other.prob;
    this -> i = other.i;
    this -> j = other.j;
}

probabilities::probabilities() : prob(0), i(0), j(0) { };

bool probabilities::operator==(const probabilities& other) const {
    return this == &other;
};

bool probabilities::operator!=(const probabilities& other) const {
    return this != &other;
};

bool probabilities::operator<(const probabilities& other) const {
    return this -> prob < other.prob;
};

bool probabilities::operator<=(const probabilities& other) const {
    return this -> prob <= other.prob;
};

bool probabilities::operator>(const probabilities& other) const {
    return this -> prob > other.prob;
};

bool probabilities::operator>=(const probabilities& other) const {
    return this -> prob >= other.prob;
};

probabilities& probabilities::operator=(const probabilities& other) {
    this -> prob = other.prob;
    this -> i = other.i;
    this -> j = other.j;
    return *this;
};

double probabilities::getProb() const {
    return prob;
};

unsigned long probabilities::getI() const {
    return i;
};

unsigned long probabilities::getJ() const {
    return j;
};

void probabilities::setNull() {
    prob = 0;
    i = 0;
    j = 0;
};

void probabilities::set(double prob, unsigned long i, unsigned long j) {
    this -> prob = prob;
    this -> i = i;
    this -> j = j;
}