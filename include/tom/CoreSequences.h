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

/** Reverse all words in the \c words vector in-place. */
void reverseWords(std::shared_ptr<Sequences> words) {
    for (unsigned long i = 0; i < words->size(); ++i)
        (*words)[i] = (*words)[i].reverse();
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
 * @param minCount the minimum number of occurrences in the data sequence for returned words (default 1)
 * @param maxWords the maximum number of returned words, or `0` (default) for no limit
 * @param uniquePositions if set to `true`, then only retain the shortest word for words occurring at the same set of positions in the data sequence (default `false`)
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
                                                long minLength = 0, long maxLength = 0, long minCount = 1, long maxWords = 0,
                                                bool uniquePositions = false) {
    auto sequence = dataSuffixTree->sequence();
    int IO = sequence.isIO() ? 2 : 1;
    if (maxLength == 0) maxLength = sequence.length();
    auto words = std::make_shared<Sequences>();
    auto compareEdgeNodes = [] (const stree::EdgeNode& a, const stree::EdgeNode& b) { return a.count() < b.count(); };
    std::priority_queue<stree::EdgeNode, std::vector<stree::EdgeNode>, decltype(compareEdgeNodes) > nodeQueue (compareEdgeNodes);

    auto node = stree::EdgeNode(dataSuffixTree);
    auto child = stree::EdgeNode(dataSuffixTree); // used later
    if (maxLength < minLength or not node.isValid() or node.count() < minCount) return words;
    nodeQueue.push(node);
    if (minLength == 0) { words->push_back(sequence.rawSub(0, 0)); minLength = 1; }

    while (not nodeQueue.empty() and (maxWords == 0 or words->size() < maxWords)) {
        node = nodeQueue.top();
        nodeQueue.pop();
        stree::nidx_t node_depth = node.depth();
        if (node_depth < IO * maxLength) {
            child = node.child();
            while (child.isValid()) {
                if (child.count() >= minCount) nodeQueue.push(child);
                child.toSibling();
            }
        }
        if (node_depth < IO * minLength) { continue; }
        if (node_depth > IO * maxLength) { node_depth = IO * maxLength; }
        stree::nidx_t node_min_depth = node.parent().depth() + 1;
        if ((IO == 2) and (node_min_depth & 1)) ++node_min_depth;
        if (node_min_depth < IO * minLength) node_min_depth = IO * minLength;
        if (uniquePositions) {
            if (node_min_depth <= node_depth) words->push_back(sequence.rawSub(node.headIndex(), node_min_depth));
        } else {
            stree::nidx_t head_index = node.headIndex();
            while (node_min_depth <= node_depth and (maxWords == 0 or words->size() < maxWords)) {
                words->push_back(sequence.rawSub(head_index, node_min_depth));
                node_min_depth += IO;
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

std::shared_ptr<std::vector<stree::nidx_t> > getIndicativeSequenceNodes(
        const std::shared_ptr<stree::STree> reverseDataSuffixTree, int minIndCount, int maxIndLen) {
    auto indNodes = std::make_shared<std::vector<stree::nidx_t> >();
    for (auto node = stree::DFSIterator(reverseDataSuffixTree); node.isValid(); node.toNext()) {
        if (node.isFirstVisit()) {
            if (node.count() < minIndCount or node.depth() > maxIndLen) {
                node.setUpPass();
            } else {
                indNodes->push_back(node.nidx());
            }
        }
    }
    return indNodes;
}

} // namespace tom

#endif // CORESEQUENCES_H
