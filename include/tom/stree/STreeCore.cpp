#include "stree.h"

namespace stree {

nidx_t & STree::c(const nidx_t node, Symbol chr) {
    nidx_t *child = &(c(node));
    if (edgeKey(*child, d(node)) == chr) return *child;
    child = &(l(*child));
    if (!(*child & VALID)) return *child;
    if (edgeKey(*child, d(node)) == chr) return *child;
    return internal::RBSiblingTree::find(l(*child), chr, internal::RBTreeNodeTraits(this, d(node)));
}

const nidx_t & STree::c(const nidx_t node, Symbol chr) const {
    return const_cast<STree*>(this)->c(node, chr);
}

void STree::addChild(const nidx_t node, nidx_t newChild, Symbol chr) {
    internal::RBTreeNodeTraits rbnt(this, d(node));
    nidx_t * x = &(c(node));
    if (rbnt.less(chr, *x)) {
        l(newChild) = l(*x);
        r(newChild) = r(*x);
        nidx_t tmp = *x;
        *x = newChild;
        newChild = tmp;
    }
    x = &(l(*x));
    if (*x & VALID) {
        if (rbnt.less(chr, *x)) {
            l(newChild) = l(*x);
            r(newChild) = r(*x);
            nidx_t tmp = *x;
            *x = newChild | COLOR;
            newChild = tmp;
        }
        x = &(l(*x));
        if (!(*x & VALID)) { // pass the suffix-link to the newChild
            l(newChild) = *x & ~COLOR;
            r(newChild) = *x & ~COLOR;
        }
        internal::RBSiblingTree::insert(*x, newChild, chr, rbnt);
    } else { // This case is *only* called when adding the second child to the root
        *x = newChild | COLOR;
        l(newChild) = INTERNAL; // suffix link
        r(newChild) = 0; // d;
    }
}

nidx_t STree::sib(const nidx_t node) const {
    const nidx_t * x = &(r(node));
    if (*x & VALID) { // find next sibling in rb-tree
        for (const nidx_t * tmp = &(l(*x)); *tmp & VALID; tmp = &(l(*x)))
            x = tmp;
        return *x;
    }
    if (*x & COLOR) { // normal case: follow rb-tree thread to next sibling
        return *x | VALID;
    }
    // special case: first, second, or last sibling
    if (*x & INTERNAL) { // must be last sibling;
        return 0;
    } else { // first or second sibling
        x = &(l(node));
        if (*x & COLOR) // node is first sibling
        return *x;
        else { // node is second sibling: return leftmost in rb-tree
            for (const nidx_t * tmp = x; *tmp & VALID; tmp = &(l(*x)))
                x = tmp;
            return *x;
        }
    }
}

void STree::createNewLeaf(bool swap) {
    // NOTE: this function invalidates the edgePtr of the currentPos
    if (!currentPos_.isExplicit()) {
        const nidx_t oldEdge = *currentPos_.edgePtr_;
        const nidx_t newNode = (VALID | INTERNAL | (*currentPos_.edgePtr_ & COLOR) | (nidx_t)(nodes_.size()));
        *currentPos_.edgePtr_ = newNode;
        currentPos_.edgePtr_ = NULL; // ... since it is invalidated anyway by the following:
        if (swap and (sequence_.rawAt(pos_) < sequence_.rawAt(currentPos_.hIndex_ + currentPos_.depth_))) {
            nodes_.push_back( { l(oldEdge), r(oldEdge), VALID | (nidx_t)(leaves_.size()) } );
            leaves_.push_back( { oldEdge | COLOR, pos_ - currentPos_.depth_ } );
            l(oldEdge) = 0;
            r(oldEdge) = currentPos_.depth_; // sets depth of the new internal node
        } else {
            nodes_.push_back( { l(oldEdge), r(oldEdge), oldEdge } );
            leaves_.push_back( { 0, currentPos_.depth_ } );
            l(oldEdge) = VALID | COLOR | (leaves_.size() - 1);
            r(oldEdge) = pos_ - currentPos_.depth_; // sets hi of the new internal node
        }
        if (r(newNode) & (VALID | INTERNAL | COLOR)) // the new node is organized in an RBTree structure
        internal::RBSiblingTree::fixThreading(newNode, internal::RBTreeNodeTraits(this));
        // At this point we can set the suffix link leading TO this internal node:
        if (suffixLinkFrom_ & VALID) sl(suffixLinkFrom_, newNode);
        suffixLinkFrom_ = newNode;
    } else {
        leaves_.push_back( { 0, 0 } );
        addChild(currentPos_.node_, ((leaves_.size() - 1) | VALID), sequence_.rawAt(pos_));
    }
}

void STree::createTemporaryInternalNodes() {
    assert(suffixLinkFrom_ == 0);
    internal::Pos currentPosOld = currentPos_;
    nTemporaryInternalNodes_ = 0;
    while (!currentPos_.isExplicit()) {
        createNewLeaf(false);
        l(nodes_.back().c_) &= ~VALID; // the sentinel is not counted as a "real" leaf
        nTemporaryInternalNodes_++;
        currentPos_.followSuffixLink(this);
    }
    if (suffixLinkFrom_ & VALID) sl(suffixLinkFrom_, currentPos_.node_);
    suffixLinkFrom_ = 0;
    // Reset the current tree position:
    currentPos_ = currentPosOld;
    if (nTemporaryInternalNodes_ > 0) currentPos_.edgePtr_ = NULL; // since it will have been invalidated!
}

void STree::removeTemporaryInternalNodes() {
    if (nTemporaryInternalNodes_ > 0) {
        internal::Pos currentPosOld = currentPos_;
        currentPos_.preCanonize(this);
        for (nidx_t nToDo = nTemporaryInternalNodes_; nToDo > 0; nToDo--) {
            nidx_t newEdge = *(currentPos_.edgePtr_);
            nidx_t oldEdge = c(newEdge);
            *(currentPos_.edgePtr_) = oldEdge;
            l(oldEdge) = l(newEdge);
            r(oldEdge) = r(newEdge);
            if (r(oldEdge) & (VALID | INTERNAL | COLOR)) // the oldEdge is organized in an RBTree structure
            internal::RBSiblingTree::fixThreading(oldEdge, internal::RBTreeNodeTraits(this));
            currentPos_.node_ = sl(currentPos_.node_);
            if (symbolSize_ > currentPos_.depth_) {
                currentPos_.hIndex_ += currentPos_.depth_;
                currentPos_.depth_ = 0;
            } else {
                currentPos_.hIndex_ += symbolSize_;
                currentPos_.depth_ -= symbolSize_;
            }
            currentPos_.preCanonize(this);
        }
        leaves_.resize(leaves_.size() - nTemporaryInternalNodes_);
        nodes_.resize(nodes_.size() - nTemporaryInternalNodes_);
        currentPos_ = currentPosOld;
        nTemporaryInternalNodes_ = 0;
    }
}

void STree::extendTo(const Sequence& sequence) CHECK(throw (std::invalid_argument)) {
#ifdef TOM_CHECK
    if (sequence_.rawSize() > sequence.rawSize() or not (sequence.rawSub(0,sequence_.rawSize()) == sequence_)) {
        throw std::invalid_argument("Current sequence is not a prefix of the given sequence!");
    }
#endif
    if (sequence.rawSize() == 0) { return; }
    if (sequence_.rawSize() == 0) {
        leaves_.push_back( { INTERNAL, 0 } );
        nodes_.push_back( { INTERNAL, 0, VALID } );
    }
    sequence_ = sequence;
    removeTemporaryInternalNodes();
    leaves_.reserve(sequence_.rawSize()); // This may invalidate the current Position:
    currentPos_.canonize(this);
    while (pos_ < sequence_.rawSize()) {
        if ((currentPos_.depth_ != 0) or (pos_ % symbolSize_ == 0)) {
            while (!currentPos_.followSymbol(this, sequence_.rawAt(pos_))) {
                createNewLeaf();
                if (currentPos_.depth_ == 0) break;
                currentPos_.followSuffixLink(this);
                if ((suffixLinkFrom_ & VALID) and currentPos_.isExplicit()) {
                    sl(suffixLinkFrom_, currentPos_.node_);
                    suffixLinkFrom_ = 0;
                }
                if ((currentPos_.depth_ == 0) and (pos_ % symbolSize_ != 0)) break;
            }
        }
        pos_++;
    }
    createTemporaryInternalNodes();
    annotate();
}

bool internal::Pos::followSymbol(STree* stree, const Symbol chr) {
    if (isExplicit()) {
        nidx_t * nextNode = &(stree->c(node_, chr));
        if (!(*nextNode & VALID)) return false;
        edgePtr_ = nextNode;
        hIndex_ = stree->hi(*edgePtr_);
    } else {
        if (chr != (stree->sequence_.rawAt(hIndex_ + depth_))) return false;
    }
    depth_++;
    if ((*edgePtr_ & INTERNAL) and (stree->d(*edgePtr_) == depth_)) node_ = *edgePtr_;
    return true;
}

void internal::Pos::preCanonize(STree* stree) {
    edgePtr_ = &node_;
    while (stree->d(*edgePtr_) < depth_) {
        node_ = *edgePtr_;
        edgePtr_ = &(stree->c(node_, stree->sequence_.rawAt(hIndex_ + stree->d(node_))));
    }
    hIndex_ = stree->hi(*edgePtr_);
}

void internal::Pos::canonize(STree* stree) {
    preCanonize(stree);
    if ((*edgePtr_ & INTERNAL) and (stree->d(*edgePtr_) == depth_)) node_ = *edgePtr_;
}

void internal::Pos::followSuffixLink(STree* stree) {
    node_ = stree->sl(node_);
    if (stree->symbolSize_ > depth_) {
        depth_ = 0;
        hIndex_ += depth_;
    } else {
        hIndex_ += stree->symbolSize_;
        depth_ -= stree->symbolSize_;
    }
    canonize(stree);
}

void STree::annotate() {
    nOccurrences_.assign(nodes_.size(), 0);
    if (currentPos_.depth_ > 0) {
        internal::Pos tempIntNode = currentPos_;
        tempIntNode.canonize(this);
        for (nidx_t tnode = tempIntNode.node_; (tnode & INDEX) != 0; tnode = sl(tnode))
            nOccurrences_[tnode & INDEX]++;
    }
    for (auto it = PostfixIterator(std::shared_ptr<STree>(this, nodelete())); it.isValid(); it.toNext()) {
        if (it.isLeaf()) nOccurrences_[it.parent().index()]++;
        else {
            auto parent = it.parent();
            if (parent.isValid()) {
                nOccurrences_[parent.index()] += nOccurrences_[it.index()];
            }
        }
    }
}

nidx_t STree::deepestInternalSuffixNidx() const {
    internal::Pos deepestVirtualLeafBranch = currentPos_;
    deepestVirtualLeafBranch.canonize(const_cast<STree*>(this));
    nidx_t ret = deepestVirtualLeafBranch.node_ | VALID;
    if (not validate(ret)) ret &= ~VALID;
    return ret;
}

} // namespace stree
