#ifndef TOM_IMPLEMENTATIONS_H
#define TOM_IMPLEMENTATIONS_H

#include "tom.h"

namespace tom {

double Oom::crossEntropyOfKOrderMarkovApproximation(int k, const Sequence &sequence) {
    auto hist = wordsOverAlphabet(nOutputSymbols(), 0, k, k);
    auto symb = wordsOverAlphabet(nOutputSymbols(), 0);
    reset();
    auto P = f(*symb, *hist);
    normalizeCols(P);
    reset();
    double l = -log2(f(sequence.sub(0, k)));
    for (long t = k; t<sequence.length(); ++t)
        l -= log2(P(sequence.u(t) * sequence.nInputSymbols() + sequence.o(t), sequence.sub(t - k, k).lexicographicIndex(true)));
    return l / sequence.length();
}

} //namespace tom

#endif //TOM_IMPLEMENTATIONS_H
