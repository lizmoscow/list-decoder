//
// Created by Moskovskaya Elizaveta on 01.03.18.
//

#ifndef LIST_DECODER_PROBS_H
#define LIST_DECODER_PROBS_H


class probabilities {
private:
    double prob;
    unsigned long i;
    unsigned long j;
public:
    probabilities(double, unsigned long, unsigned long);
    probabilities(const probabilities&);
    probabilities();
    ~probabilities() = default;
    bool operator==(const probabilities&) const;
    bool operator!=(const probabilities&) const;
    bool operator<(const probabilities&) const;
    bool operator<=(const probabilities&) const;
    bool operator>(const probabilities&) const;
    bool operator>=(const probabilities&) const;
    probabilities& operator=(const probabilities&);

    double getProb() const;
    unsigned long getI() const;
    unsigned long getJ() const;
    void setNull();
    void set(double prob, unsigned long i, unsigned long j);
};


#endif //LIST_DECODER_PROBS_H
