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

/** reverse all the sequences in the \c Sequences vector \c seqs. */
void reverseSequences(std::shared_ptr<Sequences> seqs) {
	for (unsigned long i = 0; i < seqs->size(); ++i)
		(*seqs).at(i).reverse();
}

/**
 * Find core sequences to be used in the OOM learning algorithms. Using a suffix tree representation of the (reverse) training sequence given by \c sfxTree, this function returns all sequences of the given length between \c minSeqLen and \c maxSeqLen that occur at least \c minCounts times as subsequences in the (reverse) training sequence. It may happen that two (or more) sequences x and xy are found such that x is always followed by y in the training sequence. In this case, if \c unique is set to \c true, only the shorter sequence x will be returned. The sequences are returned sorted in descending order according to their occurence counts in the training sequence. If \c maxCoreSeqs is set, then only the first \c maxCoreSeqs number of sequences are returned.
 *
 * @param sfxTree a suffix tree representation of the (reverse) training sequence
 * @param minSeqLen the minimum length for core sequences
 * @param maxSeqLen the maximum length for core sequences or -1, if \c maxSeqLen should equal \c minSeqLen
 * @param minCounts the minimum number of occurrence counts in the training sequence for core sequences
 * @param maxCoreSeq the maximum number of core sequences to return, or -1 if unlimited
 * @param unique in case two seqences x and xy qualify such that x is always followed by y in the training sequence, then only the shorter sequence x is returned if \c unique is set to \c true.
 *
 * @return the core sequences
 */
std::shared_ptr<Sequences> coreSequences(const stree::STree* sfxTree,
                                         int minSeqLen = 0,
                                         int maxSeqLen = -1,
                                         int minCounts = 1,
                                         int maxCoreSeq = -1,
                                         bool unique = true) {
	int IO = sfxTree->symbolSize_;
	if (maxSeqLen <= minSeqLen) maxSeqLen = minSeqLen;
	std::shared_ptr<Sequences> coreSeqs(new Sequences());

	stree::STreeEdge node = stree::STreeEdge(sfxTree);
	std::priority_queue<stree::STreeEdge> nodeQueue;
	nodeQueue.push(node);
	if (minSeqLen == 0) coreSeqs->push_back(sfxTree->text_.rawSub(0,0)); // handles root case
	while (!nodeQueue.empty() and (maxCoreSeq == -1 or coreSeqs->size() < maxCoreSeq)) {
		node = nodeQueue.top();
		nodeQueue.pop();
		stree::Idx node_depth = node.depth();
		if (node_depth < IO * maxSeqLen) {
			stree::STreeEdge child = node.getChild();
			while (child.isValid()) {
				if (child.count() >= minCounts) nodeQueue.push(child);
				child.sibling();
			}
		}
		if (node_depth < IO * minSeqLen) { continue; }
		if (node_depth >= IO * maxSeqLen) node_depth = IO * maxSeqLen;
		stree::Idx d = node.parentDepth() + 1; if ((IO == 2) and (d & 1)) d++;
		if (d < IO * minSeqLen) d = IO * minSeqLen;
		if (unique) {
			if (d <= node_depth) coreSeqs->push_back(sfxTree->text_.rawSub(node.headIndex(), d));
		}
		else {
			stree::Idx head_index = node.headIndex();
			while (d <= node_depth and (maxCoreSeq == -1 or coreSeqs->size() < maxCoreSeq)) {
				coreSeqs->push_back(sfxTree->text_.rawSub(head_index, d));
				d += IO;
			}
		}
	}
	return coreSeqs;
}

std::shared_ptr<std::vector<stree::NodeId> > getIndicativeSequenceNodes(stree::STree* reverseSTree, int minIndCount, int maxIndLen) {
  std::shared_ptr<std::vector<stree::NodeId> > indNodes(new std::vector<stree::NodeId>);
  for (stree::DFSIterator node = stree::DFSIterator(reverseSTree); node.isValid(); node.next()) {
		if (node.isFirstVisit()) {
			if (node.count() < minIndCount) { node.setUpPass(); continue; }
			long nodeDepth = node.depth();
			if (nodeDepth > maxIndLen) { node.setUpPass(); continue; }
      indNodes->push_back(node.nidx());
    }
  }
  return indNodes;
}

} // namespace tom

#endif // CORESEQUENCES_H
