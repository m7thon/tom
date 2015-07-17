#include "stree.h"

namespace stree {

NodeId & STree::c(const NodeId node, Char chr) {
  NodeId *child = &(c(node));
  if (key(*child, d(node)) == chr) return *child;
  child = &(l(*child));
  if (!(*child & VALID)) return *child;
  if (key(*child, d(node)) == chr) return *child;
  return rbt::find(l(*child), chr, RBTreeNodeTraits(this, d(node)));
}

const NodeId & STree::c(const NodeId node, Char chr) const {
	return const_cast<STree*>(this)->c(node, chr);
}


void STree::addChild(const NodeId node, NodeId newChild, Char chr) {
  RBTreeNodeTraits rbnt(this, d(node));
  NodeId * x = &(c(node));
  if (rbnt.less(chr, *x)) {
    l(newChild) = l(*x);
    r(newChild) = r(*x);
    NodeId tmp = *x;
    *x = newChild;
    newChild = tmp;
  }
  x = &(l(*x));
  if (*x & VALID) {
    if (rbnt.less(chr, *x)) {
      l(newChild) = l(*x);
      r(newChild) = r(*x);
      NodeId tmp = *x;
      *x = newChild | COLOR;
      newChild = tmp;
    }
    x = &(l(*x));
		if (!(*x & VALID)) { // pass the suffix-link to the newChild
			l(newChild) = *x & ~COLOR;
			r(newChild) = *x & ~COLOR;
		}
    rbt::insert(*x, newChild, chr, rbnt);
  }
  else { // This case is *only* called when adding the second child to the root
    *x = newChild | COLOR;
    l(newChild) = NODE; // suffix link
    r(newChild) = 0; // d;
  }
}


NodeId STree::sib(const NodeId node) const {
  const NodeId * x = &(r(node));
  if (*x & VALID) { // find next sibling in rb-tree
    for (const NodeId * tmp = &(l(*x)); *tmp & VALID; tmp = &(l(*x)))
      x = tmp;
    return *x;
  }
	if (*x & COLOR) { // normal case: follow rb-tree thread to next sibling
    return *x | VALID;
	}
  // special case: first, second, or last sibling
	if (*x & NODE) { // must be last sibling;
		return 0;
	}
	else { // first or second sibling
		x = &(l(node));
		if (*x & COLOR) // node is first sibling
			return *x;
		else { // node is second sibling: return leftmost in rb-tree
			for (const NodeId * tmp = x; *tmp & VALID; tmp = &(l(*x)))
				x = tmp;
			return *x;
		}
	}
}

void STree::createNewLeaf(bool swap) {
	// NOTE: this function invalidates the edgePtr of the currentPos
  if (!currentPos.isExplicit()) {
    const NodeId oldEdge = *currentPos.edgePtr;
    const NodeId newNode = (VALID | NODE | (*currentPos.edgePtr & COLOR) | (Idx)(nodes.size()));
    *currentPos.edgePtr = newNode;
    currentPos.edgePtr = NULL; // ... since it is invalidated anyway by the following:
    if (swap and (at(pos) < at(currentPos.hIndex + currentPos.depth))) {
      nodes.push_back(InternalNode(l(oldEdge), r(oldEdge), (VALID | leaves.size())));
      leaves.push_back(LeafNode((oldEdge | COLOR), (pos - currentPos.depth)));
      l(oldEdge) = 0;
      r(oldEdge) = currentPos.depth; // sets depth of the new internal node
    }
    else {
      nodes.push_back(InternalNode(l(oldEdge), r(oldEdge), oldEdge));
      leaves.push_back(LeafNode(0, currentPos.depth));
      l(oldEdge) = VALID | COLOR | (leaves.size() - 1);
      r(oldEdge) = pos - currentPos.depth; // sets hi of the new internal node
    }
		if (r(newNode) & (VALID | NODE | COLOR)) // the new node is organized in an RBTree structure
				rbt::fixThreading(newNode, RBTreeNodeTraits(this));
    // At this point we can set the suffix link leading TO this internal node:
    if (suffixLinkFrom & VALID) sl(suffixLinkFrom, newNode);
    suffixLinkFrom = newNode;
  }
  else {
    leaves.push_back(LeafNode());
    addChild(currentPos.node, ((leaves.size()-1) | VALID), at(pos));
  }
}


void STree::createTemporaryInternalNodes() {
  assert(suffixLinkFrom == 0);
  Position currentPosOld = currentPos;
  nTemporaryInternalNodes = 0;
  while (!currentPos.isExplicit()) {
    createNewLeaf(false);
    l(nodes.back().c) &= ~VALID; // the sentinel is not counted as a "real" leaf
    nTemporaryInternalNodes++;
    currentPos.followSuffixLink(this);
  }
	if (suffixLinkFrom & VALID) sl(suffixLinkFrom, currentPos.node);
	suffixLinkFrom = 0;
  // Reset the current tree position:
  currentPos = currentPosOld;
  if (nTemporaryInternalNodes > 0)
    currentPos.edgePtr = NULL; // since it will have been invalidated!
}

void STree::removeTemporaryInternalNodes() {
  if (nTemporaryInternalNodes > 0) {
    Position currentPosOld = currentPos;
    currentPos.preCanonize(this);
    for (Idx nToDo = nTemporaryInternalNodes; nToDo > 0; nToDo--) {
      NodeId newEdge = *(currentPos.edgePtr);
      NodeId oldEdge = c(newEdge);
      *(currentPos.edgePtr) = oldEdge;
      l(oldEdge) = l(newEdge);
      r(oldEdge) = r(newEdge);
			if (r(oldEdge) & (VALID | NODE | COLOR)) // the oldEdge is organized in an RBTree structure
				rbt::fixThreading(oldEdge, RBTreeNodeTraits(this));
			currentPos.node = sl(currentPos.node);
			if (symbolSize_ > currentPos.depth) {
				currentPos.hIndex += currentPos.depth;
				currentPos.depth = 0;
			}
			else {
				currentPos.hIndex += symbolSize_;
				currentPos.depth -= symbolSize_;
			}
      currentPos.preCanonize(this);
    }
    leaves.resize(leaves.size()-nTemporaryInternalNodes);
    nodes.resize(nodes.size()-nTemporaryInternalNodes);
    currentPos = currentPosOld;
    nTemporaryInternalNodes = 0;
  }
}

void STree::extendTo(Idx size) {
	assert(size >= size_);
  size_ = size;
	removeTemporaryInternalNodes();
  leaves.reserve(size_); // This may invalidate the current Position:
  currentPos.canonize(this);
  while (pos < size_) {
    if ((currentPos.depth != 0) or (pos % symbolSize_ == 0)) {
      while (!currentPos.followChar(this, at(pos))) {
				createNewLeaf();
				if (currentPos.depth == 0) break;
				currentPos.followSuffixLink(this);
				if ((suffixLinkFrom & VALID) and currentPos.isExplicit()) {
					sl(suffixLinkFrom, currentPos.node);
					suffixLinkFrom = 0;
				}
				if ((currentPos.depth == 0) and (pos % symbolSize_ != 0)) break;
      }
    }
		pos++;
  }
	createTemporaryInternalNodes();
	if (annotated_) annotate();
}

void  STree::initialize(const String& text, Idx size, unsigned int symbolSize, bool annotated) {
	size_ = size; annotated_ = annotated; symbolSize_ = symbolSize;
	text_ = text;
	if (size_ == 0) size_ = text_.rawSize();
  leaves.push_back(LeafNode((NODE), (0)));
  nodes.push_back(InternalNode((NODE), (0), (VALID)));
  suffixLinkFrom = 0;
  pos = 1;
  nTemporaryInternalNodes = 0;
  extendTo(size_);
}

bool STree::Position::followChar(STree* stree, const Char chr) {
  if (isExplicit()) {
    NodeId * nextNode = &(stree->c(node, chr));
    if (!(*nextNode & VALID)) return false;
    edgePtr = nextNode;
    hIndex = stree->hi(*edgePtr);
  }
  else {
    if (chr != (stree->at(hIndex + depth)))
      return false;
  }
  depth++;
  if ((*edgePtr & NODE) and (stree->d(*edgePtr) == depth))
    node = *edgePtr;
  return true;
}


void STree::Position::preCanonize(STree* stree) {
  edgePtr = &node;
  while (stree->d(*edgePtr) < depth) {
    node = *edgePtr;
    edgePtr = &(stree->c(node, stree->at(hIndex + stree->d(node))));
  }
  hIndex = stree->hi(*edgePtr);
}


void STree::Position::canonize(STree* stree) {
  preCanonize(stree);
  if ((*edgePtr & NODE) and (stree->d(*edgePtr) == depth))
    node = *edgePtr;
}


void STree::Position::followSuffixLink(STree* stree) {
  node = stree->sl(node);
  if (stree->symbolSize_ > depth) {
    depth = 0;
    hIndex += depth;
  }
  else {
    hIndex += stree->symbolSize_;
    depth -= stree->symbolSize_;
  }
  canonize(stree);
}

void STree::annotate() {
	nOccurrences.assign(nodes.size(), 0);
	if (currentPos.depth > 0) {
		Position tempIntNode = currentPos;
		tempIntNode.canonize(this);
		for (NodeId tnode = tempIntNode.node; (tnode & INDEX) != 0; tnode = sl(tnode))
			nOccurrences[tnode & INDEX]++;
	}
	for (PostfixIterator it = PostfixIterator(this); it.isValid(); it.next()) {
		if (it.isLeaf())
			nOccurrences[it.getParent().index()]++;
		else if (it.nAncestors() > 0) {
			nOccurrences[it.getParent().index()] += nOccurrences[it.index()];
		}
	}
}

STreeNode STree::getDeepestVirtualLeafBranch() {
  Position deepestVirtualLeafBranch = currentPos;
  deepestVirtualLeafBranch.canonize(this);
  return STreeNode(this, deepestVirtualLeafBranch.node | VALID);
}

} // namespace stree
