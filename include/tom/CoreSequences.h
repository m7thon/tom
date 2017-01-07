/**
 * @file   CoreSequences.h
 * @author Michael Thon <mthon@jacobs-university.de>
 *
 * @brief  Provides functions to obtain the required core sequences from data.
 */

#ifndef CORESEQUENCES_H
#define CORESEQUENCES_H

#include "tom.h"

namespace tom {

/** Reverse the given `words` in-place. */
void reverseWords(std::shared_ptr<Sequences> words) {
    for (unsigned long i = 0; i < words->size(); ++i)
        (*words)[i] = (*words)[i].reverse();
}

/** Sort the given `words` in-place by length and then lexicographically.*/
void sortWords(std::shared_ptr<Sequences> words) {
    std::sort(words->begin(), words->end());
}

SWIGCODE(%feature ("kwargs") wordsOverAlphabet;)
/** Return in lexicographic order all words of length between `minLength` and `maxLength` over the alphabet with `nOutputSymbols` output symbols if `nInputSymbols` is zero, or otherwise over the alphabet of input-output symbol pairs with `nOutputSymbols` output symbols and `nInputSymbols` input symbols.
 *
 * @param nOutputSymbols the number of output symbols
 * @param nInputSymbols the number of input symbols (default 0)
 * @param minLength the minimum length for returned words (default 1)
 * @param maxLength the maximum length for returned words (default 1)
 * @return an array of words
 */
std::shared_ptr<Sequences> wordsOverAlphabet(int nOutputSymbols, int nInputSymbols = 0, int minLength = 1, int maxLength = 1) {
    auto words = std::make_shared<Sequences>();
    struct WordIterator {
        Sequence word;
        void increment(long position = 0, bool outputSymbol = true) {
            if (position == word.rawSize()) {
                word = Sequence(word.isIO() ? position/2 + 1 : position + 1, word.nOutputSymbols(), word.nInputSymbols());
                return;
            }
            int new_value_at_pos = word.rawAt(-1-position) + 1;
            if (new_value_at_pos == (outputSymbol ? word.nOutputSymbols() : word.nInputSymbols())) {
                word.rawAt(-1-position, 0);
                increment(position+1, not word.isIO() or not outputSymbol);
            } else {
                word.rawAt(-1-position, new_value_at_pos);
            }
        }
        Sequence operator()() { Sequence thisWord = word.copy(); increment(); return thisWord; }
    };
    WordIterator nextWord{ Sequence(minLength, nOutputSymbols, nInputSymbols) };
    int nSymbols = nOutputSymbols * std::max(nInputSymbols, 1);
    for (int length = minLength; length <= maxLength; ++length) {
        long nSequencesForLength = (long)std::pow((double)nSymbols, (double)length);
        for (long count = 0; count < nSequencesForLength; ++count) {
            words->push_back(nextWord());
        }
    }
    return words;
}


SWIGCODE(%feature ("kwargs") wordsFromData;)
/** Find words occurring in a sequence according to given constraints. This can be used to find sets of indicative or characteristic words from training data that are appropriate for the OOM learning algorithms.
 *
 * From a suffix tree representation of the data sequence given by `dataSuffixTree`, this function efficiently computes an array of words with length between `minLength` and `maxLength` that occur at least `minCount` times as subsequences in the data sequence, sorted by occurrence count and then lexicographically. If `uniquePositions` is `true`, then only the first (shortest) word is retained for words occurring at the same set of positions in the data sequence. Finally, only the first at most `maxWords` resulting words are returned.
 *
 * @param dataSuffixTree a suffix tree representation of the data sequence
 * @param minLength the minimum length for returned words (default 0)
 * @param maxLength the maximum length for returned words, or `0` (default) for no limit
 * @param minRelevance the minimum "relevance" in the data sequence for returned words (default 1.0)
 * @param maxWords the maximum number of returned words, or `0` (default) for no limit
 * @param prefixUnique if set to `true`, then only retain the shortest prefix among words occurring at the same set of positions in the data sequence (default `false`). This may be useful when selecting characteristic sequences.
 * @param suffixUnique if set to `true`, then only retain the shortest suffix among words ending at the same set of positions in the data sequence (default `false`).This may be useful when selecting indicative sequences.
 * @param relevance a `PositionRelevance` object to compute the *relevance value* for words, which defaults to their occurrence count in the data (plus a fraction to penalize longer words). This can be customized by passing a `relevance` inherited from `PositionRelevance` with an overwritten `.compute()` method.
 *
 * @return an array of words occurring in the data sequence according to the given constraints, and sorted by occurrence count (descending) and then lexicographically.
 *
 * Notes
 * -----
 * If two words a and b with |a| < |b| always occur at the same positions in the data sequence, then a must be a prefix of b, i.e., b = ac for some non-empty word c. Furthermore, this means that for any word x, the occurrence counts for the words xa and xb will be the same, and moreover these counts will be based on exactly the same occurrences in the data sequence. When collecting such occurrence statistics, it may therefore be useful to use only one (e.g., the shortest) of the words with identical occurrence positions. This is what `uniquePositions` is for.
 *
 * Note that in the context of OOM learning, the above means that `uniquePositions` makes sense when finding **characteristic** words from the training data. To achieve the same (avoid exact duplicate data exploitation) for **indicative** words, one needs to use this function with `uniquePositions` for a **reversed** training sequence, and then reverse each word in the resulting set.
 *
 */
std::shared_ptr<Sequences> wordsFromData(const std::shared_ptr<const stree::STree>& dataSuffixTree,
                                                long minLength = 0, long maxLength = 0, double minRelevance = 1.0, long maxWords = 0,
                                                bool prefixUnique = false, bool suffixUnique = false, const stree::PositionRelevance& relevance = stree::PositionRelevance()) throw(std::invalid_argument) {
    if (prefixUnique and suffixUnique) throw std::invalid_argument("`prefixUnique` and `suffixUnique` cannot both be set");
    stree::PositionRelevance& key = const_cast<stree::PositionRelevance&>(relevance);
    auto sequence = dataSuffixTree->sequence();
    int nSymbols = sequence.nInputSymbols() * std::min(sequence.nOutputSymbols(), 1);
    int IO = sequence.isIO() ? 2 : 1;
    if (maxLength == 0) maxLength = sequence.length();
    auto words = std::make_shared<Sequences>();
    auto comparePositions = [&key] (const stree::Position& a, const stree::Position& b) {
        return key.compute(a) < key.compute(b);
    };
    std::priority_queue<stree::Position, std::vector<stree::Position>, decltype(comparePositions) > positionQueue(comparePositions);

    auto p = stree::Position(dataSuffixTree);
    if (maxLength < minLength or not p.isValid() or key.compute(p) < minRelevance) { return words; }
    positionQueue.push(p);

    while ((not positionQueue.empty()) and (maxWords == 0 or words->size() < maxWords)) {
        p.set(positionQueue.top());
        positionQueue.pop();
        if (p.depth() >= IO * minLength) {
            if (!suffixUnique or (p.depth() == IO * minLength) or (p.suffix().count() > p.count())) {
                words->push_back(p.sequence());
                if (prefixUnique) { p.toExplicit(); }
            }
        } else {
            p.toDepth(std::min((long)p.edge().depth(), IO * (minLength - 1)));
        }
        if (p.depth() < IO * maxLength) {
            if (p.depth() % IO == 0) { p.toChild(); }
            while (p.isValid()) {
                if (IO == 2) {
                    auto grandchild = p.child();
                    while (grandchild.isValid()) {
                        if (key.compute(grandchild) >= minRelevance) positionQueue.push(grandchild);
                        grandchild.toSibling();
                    }
                } else {
                    if (key.compute(p) >= minRelevance) positionQueue.push(p);
                }
                p.toSibling();
            }
        }
    }
    return words;
}

SWIGCODE(%feature ("kwargs") wordsFromModel;)
/** Return all words satisfying the given constraints sorted (descending) by their probability according to the given `oom`.
 *
 * @param oom the `Oom` model from which to compute the word probabilities
 * @param minLength the minimum length for returned words (default 0)
 * @param maxLength the maximum length for returned words, or `0` (default) for no limit
 * @param minProbability the minimum probability for returned words (default 1e-5)
 * @param maxWords the maximum number of returned words, or `0` (default) for no limit
 */
std::shared_ptr<Sequences> wordsFromModel(Oom &oom, long minLength = 0, long maxLength = 0, double minProbability = 1e-5, long maxWords = 0) {
    auto words = std::make_shared<Sequences>();
    auto compare = [] (const std::pair<double, Sequence>& a, const std::pair<double, Sequence>& b) {
        return a.first < b.first or (a.first == b.first and a.second.length() > b.second.length());
    };
    std::priority_queue<std::pair<double, Sequence>, std::vector<std::pair<double, Sequence> >, decltype(compare) > wordQueue (compare);
    wordQueue.push(std::make_pair(1.0, Sequence(0, oom.nOutputSymbols(), oom.nInputSymbols())));
    while ((not wordQueue.empty()) and (maxWords == 0 or words->size() < maxWords)) {
        double p_word = wordQueue.top().first;
        Sequence word = wordQueue.top().second;
        wordQueue.pop();
        if (word.length() >= minLength) { words->push_back(word); }
        if (word.length() < maxLength or maxLength == 0) {
            for (Symbol u = 0; u < (oom.isIO() ? oom.nInputSymbols() : 1); ++u) {
                for (Symbol o = 0; o < oom.nOutputSymbols(); ++o) {
                    Sequence nextWord(word);
                    if (oom.isIO()) { nextWord = nextWord + u; }
                    nextWord = nextWord + o;
                    double p_nextWord = oom.f(nextWord);
                    if (oom.isIO()) { p_nextWord *= std::pow(1.0/oom.nInputSymbols(), nextWord.length()); }
                    if (p_nextWord > minProbability) { wordQueue.push(std::make_pair(p_nextWord, nextWord)); }
                }
            }
        }
    }
    return words;
}

} // namespace tom

#endif // CORESEQUENCES_H
