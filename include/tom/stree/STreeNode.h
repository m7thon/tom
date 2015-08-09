#ifndef STREE_NODE_H
#define STREE_NODE_H

#include "stree.h"

namespace stree {

//TODO: Add checks to functions where required. These are currently interface functions that can crash tom with improper use!

//TODO: Clean up iterator constructors. Both from stree and from objects should give same result.

//TODO: Complete documentation

//TODO: Write some tests to cover the various iterators.

SWIGCODE(%feature("python:slot", "tp_repr", functype="reprfunc") Node::repr;)
/** This class represents a node in the suffix tree and can be used for extracting information or navigating the suffix tree. It contains a pointer internally to the \c STree that it belongs to. Use `STree.node()` to construct a `Node`.
 *
 * The `to...()` methods are provided to navigate the suffix tree structure:
 *
 * - `to...()` sets this node to its e.g. child, suffix, sibling, etc.
 * - If no such exists, then this node is simply marked as invalid.
 * - For an invalid node, the `to...()` methods have no effect.
 * - Calling `setValid()` after the node has been marked as invalid by a `to...()` method will reset this node to the last valid node during the traversal.
 */
class Node {
    friend class STree;
    friend class EdgeNode;
    friend class Position;
protected:
    const STree* stree_; /**< the `STree` that this node refers to */
    nidx_t nidx_; /**< the internal node reference */

public:
    /** Construct a `Node` for the given `stree` corresponding to the given `nidx`. If no `nidx` is given, it defaults to the root of the suffix tree.
     */
    Node(const STree* stree, nidx_t nidx = ROOT) : stree_(stree), nidx_(nidx) {
        if (!stree_->validate(nidx_)) setValid(false);
    }

	/** Return \c true if valid, otherwise return \c false.
	 */
	bool isValid() const { return stree_ and (nidx_ & VALID); }

	/** Mark this \c Node as \c valid (or \c invalid, if \c valid is \c false).
	 */
	void setValid(bool valid = true) CHECK(throw(std::invalid_argument)) {
	    if (valid) {
	        CHECK(if (!stree_->validate(nidx_)) throw std::invalid_argument("Node does not correspond to a node in the suffix tree. Cannot set valid.");)
	        nidx_ |= VALID;
	    } else {
	        nidx_ &= ~VALID;
	    }
	}

	/** Return `true` if this is an internal node.
	 */
	bool isInternal() const { return (nidx_ & INTERNAL); }

    /** Return `true` if this is a leaf node.
     */
	bool isLeaf() const { return !isInternal(); }

	/** Return `true` if this is the root node.
	 */
	bool isRoot() const { return isValid() and isInternal() and index() == 0; }

	/** The `index` of a valid leaf or a valid internal node is a unique number between 0 and `STree.nLeafNodes()` or between 0 and `STree.nInternalNodes()`, respectively.
	 */
	nidx_t index() const { return nidx_ & INDEX; }

	/** Return the `nidx_t` corresponding to this `Node`.
	 */
	nidx_t nidx() const { return nidx_; }

	/** Return `true` if this `Node` is the same as the given `other`.
	 */
	bool operator ==(const Node& other) const { return ((nidx_ | COLOR) == (other.nidx() | COLOR)); }
    
	/** Return `true` if the `other` node occurs more often in the represented sequence than this node.
	 */
	bool operator <(const Node& other) const CHECK(throw(std::invalid_argument)) {
	    CHECK(if (!isValid() or !other.isValid()) throw std::invalid_argument("Cannot compare invalid nodes.");)
	    return count() < other.count();
	}

	/** Return the "depth" of the node in the suffix tree, which is the size of the represented (sub-)sequence.
	 */
	nidx_t depth() const CHECK(throw(std::invalid_argument)) {
	    CHECK(if (!isValid()) throw std::invalid_argument(".depth() called on invalid node.");)
	    return stree_->d(nidx_);
	}

	/** Return the "headindex" of this node, which is an index in the sequence represented by the suffix tree where the (sub-)sequence represented by this node occurs. I.e., the (sub-)sequence represented by this node is `seq.rawSub(headindex(), depth())`, where `seq` is the sequence represented by the suffix tree.
	 */
	nidx_t headIndex() const CHECK(throw(std::invalid_argument)) {
        CHECK(if (!isValid()) throw std::invalid_argument(".headIndex() called on invalid node.");)
	    return stree_->hi(nidx_);
	}

	/** Return the number of occurrences of the sequence represented by this node in the sequence represented by the suffix tree.
	 */
	nidx_t count() const CHECK(throw(std::invalid_argument)) {
        CHECK(if (!isValid()) throw std::invalid_argument(".count() called on invalid node.");)
	    return stree_->n(nidx_);
	}

	/** Return the first child node of this node. If no such node exists, a `Node` marked as invalid is returned. Note that the children are ordered lexicographically according to their edge labels.
	 */
	Node child() const { if (isValid() and isInternal()) return Node(stree_, stree_->c(nidx_)); else return Node(stree_, nidx_ & ~VALID); }

	/** Set this node to its first child if such a node exists, otherwise mark this node as invalid. Note that the children are ordered lexicographically according to their edge labels.
	 */
	void toChild() { if (isValid() and isInternal()) nidx_ = stree_->c(nidx_); else setValid(false); }

	/** Return the child node along the edge leading away whose label begins with the given `symbol`. If no such node exists, a `Node` marked as invalid is returned.
	 */
	Node child(Symbol symbol) const { if (isValid() and isInternal()) return Node(stree_, stree_->c(nidx_, symbol)); else return Node(stree_, nidx_ & ~VALID); }

	/** Set this node to the child node along the edge leading away whose label begins with the given `symbol`\. If no such node exists, mark this node as invalid instead.
	 */
	void toChild(Symbol symbol) {
		if (isValid() and isInternal()) { nidx_t chld = stree_->c(nidx_, symbol); if (chld & VALID) nidx_ = chld; else setValid(false); }
		else setValid(false);
	}

	/** Return the next sibling of this node. If no such node exists, a `Node` marked as invalid is returned. Note that the siblings are ordered lexicographically according to their edge labels.
	 */
	Node sibling() const { if (isValid()) return Node(stree_, stree_->sib(nidx_)); else return Node(stree_, nidx_ & ~VALID); }

	/** Set this node to its next sibling if a next sibling exists, otherwise mark this node as invalid. Note that the siblings are ordered lexicographically according to their edge labels.
	 */
	void toSibling() { if (isValid()) { nidx_t sbl = stree_->sib(nidx_); if (sbl & VALID) nidx_ = sbl; else setValid(false); } }

	/** Return the node corresponding to the first suffix of the represented sequence. This follows the "suffix link" of the suffix tree. If no such node exists, a `Node` marked as invalid is returned.
	 */
	Node suffix() const {
		if (isValid() and isInternal() and (index() != 0)) return Node(stree_, stree_->sl(nidx_) | VALID);
		else return Node(stree_, nidx_ & ~VALID);
	}

	/** Set this node to the node corresponding to the first suffix of the represented sequence\. If no such node exists, mark this node as invalid instead. This follows the "suffix link" of the suffix tree.
	 */
	void toSuffix() {
		if (isValid() and isInternal() and (index() != 0)) { nidx_ = stree_->sl(nidx_) | VALID; }
		else { setValid(false); }
	}

	/** Return the (sub-)sequence represented by this node. Note that this is `seq.rawSub(headIndex(), depth())`, where `seq` is the sequence represented by the suffix tree.
	 */
	Sequence sequence() const CHECK(throw(std::invalid_argument)) {
        CHECK(if (!isValid()) throw std::invalid_argument(".sequence() called on invalid node.");)
	    return stree_->sequence_.rawSub(headIndex(), depth());
	}

	/** Return a string representation of the data of this node\. This is useful for debugging or understanding the suffix tree structure.
	 */
	std::string dataStr(int width = 5) const {
        std::stringstream s;
		if (not isValid()) {
		    s << "[ " << nidxStr(width) << " ]";
		    return s.str();
		}
		s << "[ " << nidxStr(width)
          << " , hIdx = " << std::setw(width) << headIndex()
          << " , size = "  << std::setw(width) << depth()
          << " , nOcc = "  << std::setw(width) << count();
		if (isInternal()) {
			s << " , chld = "  << child().nidxStr(width)
              << " , sfxL = " << suffix().nidxStr(width);
		}
		s << " ]";
		return s.str();
	}

	/** Return a string representation of the underlying `nidx_t`.
	 */
	std::string nidxStr(int width = 3) const {
		std::stringstream s;
		s << ( (nidx_ & VALID)    ? '+' : '-' )
          << ( (nidx_ & INTERNAL) ? '+' : '-' )
          << ( (nidx_ & COLOR)    ? '+' : '-' )
          << std::setw(width) << std::setfill('0') << (nidx_ & INDEX);
		return s.str();
	}

	/** Return a string representation to display in python.
	 */
	std::string repr() const {
        std::stringstream s;
        s << ( (nidx_ & VALID)    ? "" : "invalid " )
          << ( (nidx_ & INTERNAL) ? "internal " : "leaf ")
          << "node "
          << (nidx_ & INDEX);
        return s.str();
	}
}; // class STreeNode


/** This class represents a node together with the edge leading to it in the suffix tree and can be used for extracting information or navigating the suffix tree. Note that an `EdgeNode` can have a parent marked as invalid, and then contains no edge information. Such an `EdgeNode` is called "degenerate". This is normally only the case for the root, but can also happen when providing incorrect parent information when constructing an `EdgeNode`.
 *
 * The `to...()` methods are provided to navigate the suffix tree structure:
 *
 * - `to...()` sets this node (and edge leading to it) to its e.g. child, suffix, sibling, etc.
 * - If no such exists, then this node is simply marked as invalid.
 * - For an invalid node, the `to...()` methods have no effect.
 * - Calling `setValid()` after the node has been marked as invalid by a `to...()` method will reset this node (and edge) to the last valid node during the traversal.
 */
class EdgeNode : public Node {
    friend class STree;
    friend class Position;
protected:
    nidx_t parent_; /**< The parent node, which defines the edge. If this is invalid, the `EdgeNode` is degenerate and has no edge information. */
public:
    /** Construct an `EdgeNode` for the given `stree` corresponding to the given `nidx` and parent information in `parent`. If no `nidx` is given, it defaults to the root of the suffix tree. This will search for the parent starting form the given `parent` (or root by default). If no parent can be found (e.g. for the root), this `EdgeNode` will be "degenerate" (have an invalid parent).
     */
    EdgeNode(const STree* stree, nidx_t nidx = ROOT, nidx_t parent = ROOT) : Node(stree, nidx), parent_(nidx &~ VALID) { findParent(Node(stree_, parent)); }

    /** Construct an `EdgeNode` from the given `Node`. This will search for the parent starting form the given `parent` (or root by default). If no parent can be found (e.g. for the root), this `EdgeNode` will be "degenerate" (have an invalid parent).
     */
    EdgeNode(const Node& node, nidx_t parent = ROOT) : Node(node), parent_(parent &~ VALID) { findParent(Node(stree_, parent)); }

    /** Construct an `EdgeNode` from the given `Node`. This will search for the parent starting form the given `parent`. If no parent can be found (e.g. for the root), this `EdgeNode` will be "degenerate" (have an invalid parent).
     */
    EdgeNode(const Node& node, const Node& parent) : Node(node), parent_(node.nidx() &~ VALID) { findParent(parent); }

    /** Return the parent `Node` of this `EdgeNode`. If this `EdgeNode` is degenerate, i.e., no parent node exists or is known, a `Node` marked as invalid is returned.
     */
	Node parent() const { return Node(stree_, parent_); }

	/** Return the first child `EdgeNode` of this `EdgeNode`. If no child exists, an `EdgeNode` marked as invalid is returned. Note that the children are ordered lexicographically according to their edge labels.
     */
 	EdgeNode child() const { return EdgeNode(Node::child(), nidx_); }

    /** Set this `EdgeNode` to its first child (with edge leading to it) if such a node exists, otherwise mark this `EdgeNode` as invalid. Note that the children are ordered lexicographically according to their edge labels.
     */
 	void toChild() { nidx_t par = nidx_; Node::toChild(); if (isValid()) parent_ = par; }

    /** Return the child `EdgeNode` leading away whose label begins with the given `symbol`. If no such `EdgeNode` exists, an `EdgeNode` marked as invalid is returned.
     */
 	EdgeNode child(Symbol chr) const { return EdgeNode(Node::child(chr), nidx_); }

    /** Set this `EdgeNode` to the child `EdgeNode` leading away whose label begins with the given `symbol`\. If no such `EdgeNode` exists, mark this `EdgeNode` as invalid instead.
     */
 	void toChild(Symbol chr) { nidx_t par = nidx_; Node::toChild(chr); if (isValid()) parent_ = par; }

    /** Return the next sibling of this `EdgeNode`. If none exists, an `EdgeNode` marked as invalid is returned. Note that the siblings are ordered lexicographically according to their edge labels.
     */
 	EdgeNode sibling() const { return EdgeNode(Node::sibling(), parent_); }

    /** Return the `EdgeNode` corresponding to the first suffix of the represented sequence. This follows the "suffix link" of the suffix tree (and finds the new parent accordingly). If no such node exists, a `Node` marked as invalid is returned.
     */
 	EdgeNode suffix() const { return EdgeNode(Node::suffix(), Node(stree_, parent_).suffix()); }

    /** Set this to the `EdgeNode` corresponding to the first suffix of the represented sequence\. If no such node exists, mark this node as invalid instead. This follows the "suffix link" of the suffix tree (and finds the new parent accordingly).
     */
 	void toSuffix() { Node::toSuffix(); if (isValid()) findParent(Node(stree_, parent_).suffix()); }

 	/** Return the edge label. Note that this `EdgeNode` must not be degenerate.
 	 */
	Sequence edgeLabel() const CHECK(throw(std::invalid_argument)) {
	    Node par = parent();
	    CHECK(if (!isValid()) throw std::invalid_argument(".edgeLabel() called on invalid EdgeNode.");)
	    CHECK(if (!par.isValid()) throw std::invalid_argument(".edgeLabel() called on degenerate EdgeNode (no valid parent information).");)
	    return stree_->sequence_.rawSub(headIndex() + par.depth(), depth() - par.depth());
	}

    /** Return a string representation to display in python.
     */
    std::string repr() const {
        std::stringstream s;
        s << ( (nidx_ & VALID)    ? "" : "invalid " )
          << ( (nidx_ & INTERNAL) ? "internal " : "leaf ")
          << "node "
          << ( (parent().isValid()) ? "with" : "without" ) << " edge information"
          << (nidx_ & INDEX);
        return s.str();
    }

    /** Return a string representation of the data of this `EdgeNode`\. This is useful for debugging or understanding the suffix tree structure.
     */
    std::string dataStr(int width = 5) const {
        std::stringstream s;
        s << Node::dataStr();
        Node par = parent();
        if (par.isValid()) s << "[ parent = " << par.nidxStr() << " ]";
        return s.str();
    }

private:
    /** Find the parent for this `EdgeNode` by searching the sub-nodes of the given `Node` `searchBelow`. If no parent can be found, the parent will be marked invalid.
     */
    void findParent(Node searchBelow) {
        if (not isValid()) { return; }
        nidx_t d = depth();
        nidx_t par_depth;
        while (searchBelow.isValid() and (par_depth = searchBelow.depth()) < d) {
            Node par_child = searchBelow.child(sequence().rawAt(par_depth));
            if (par_child == *this) { parent_ = par_child.nidx(); return; }
            searchBelow = par_child;
        }
        parent_ = nidx_ &~VALID;
        return;
    }
}; // class STreeEdge



class Path : public Node {
public:
    Path(const STree* stree) : Node(stree), path() {}
	void toChild() { nidx_t current = nidx_; Node::toChild(); if (isValid()) path.push_back(current); }
	void toChild(Symbol chr) { nidx_t current = nidx_; Node::toChild(chr); if (isValid()) path.push_back(current); }
	Node parent() const { if (path.empty()) return stree_->node(INTERNAL); return stree_->node(path.back()); }
	Node ancestor(nidx_t generations) const {
		if (generations > path.size()) return stree_->node(INTERNAL);
		if (generations == 0) return *this;
		return stree_->node(path.at(path.size()-generations));
	}
	void toParent() { if (isValid()) {
			if (path.empty()) setValid(false);
			else { nidx_ = path.back(); path.pop_back(); }
		}}
	nidx_t nAncestors() const { return path.size(); }
	nidx_t parentDepth() const { if (path.empty()) return 0; else return stree_->d(path.back()); }

	Sequence edgeLabel() const { return stree_->sequence_.rawSub(headIndex() + parentDepth(), depth() - parentDepth()); }

private:
	using Node::suffix;
	using Node::toSuffix;
protected:
	std::deque<nidx_t> path;
}; // class STreePath



class Position {
private:
    EdgeNode edge_;
    nidx_t depth_;
public:
	Position(const STree* stree) : edge_(stree), depth_(0) {}

	void setRoot() { edge_ = EdgeNode(edge_.stree_); depth_ = 0; }

	bool isValid() const { return edge_.isValid(); }
	void setValid(bool valid = true) { edge_.setValid(valid); }
	bool isExplicit() const { return (depth_ == edge_.depth()); }
    bool isInternal() const { return !isLeaf(); }
	bool isLeaf() const { return (isExplicit() and edge_.isLeaf()); }

	nidx_t count() const { return edge_.count(); }
	nidx_t headIndex() const { return edge_.headIndex(); }
	nidx_t depth() const { return depth_; }
	nidx_t parentDepth() const { return edge_.parent().depth(); }

  EdgeNode& edge() { return edge_; }

	void toSuffix() { if (isValid()) {
			if (edge_.stree_->symbolSize_ > depth_) { setValid(false); return; }
			depth_-= edge_.stree_->symbolSize_;
			nidx_t hi = headIndex() + edge_.stree_->symbolSize_;
			edge_.parent_ = edge_.stree_->sl(edge_.parent_) | VALID;
			edge_.nidx_ = edge_.parent_;
			edge_.setValid();
			while (edge_.isValid() and (edge_.depth() < depth_))
				edge_.toChild(edge_.stree_->sequence_.rawAt(hi + edge_.depth()));
		}}

	void toSymbol(Symbol chr) {
		if (!isValid()) return;
		if (isExplicit()) {
			edge_.toChild(chr);
			if (isValid()) depth_++;
		}
		else
			if (edge_.stree_->sequence_.rawAt(edge_.headIndex() + depth_) == chr) depth_++;
			else edge_.setValid(false);
	}
	void toSequence(const Sequence& str) {
		for (nidx_t pos = 0; ((pos < str.rawSize()) and (edge_.isValid())); ++pos)
			toSymbol(str.rawAt(pos));
	}

	Sequence sequence() const { return edge_.stree_->sequence_.rawSub(headIndex(), depth_); }
	Sequence edgeLabel() const { return edge_.stree_->sequence_.rawSub(headIndex() + parentDepth(), depth_ - parentDepth()); }

}; // class STreePos

} // namespace stree

#endif // STREE_NODE_H
