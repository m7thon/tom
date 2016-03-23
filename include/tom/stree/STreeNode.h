#ifndef STREE_NODE_H
#define STREE_NODE_H

#include "stree.h"

namespace stree {

//TODO: Write some tests to cover the various classes.

SWIGCODE(%feature("python:slot", "tp_repr", functype="reprfunc") Node::repr;)
/** This class represents a node in the suffix tree and can be used for extracting information or navigating the suffix tree. It contains a pointer internally to the \c STree that it belongs to.
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
    std::shared_ptr<const STree> stree_; /**< the `STree` that this node refers to */
    nidx_t nidx_; /**< the internal node reference */

public:
    /** Construct a `Node` for the given `stree` corresponding to the given `nidx`. If no `nidx` is given, it defaults to the root of the suffix tree.
     */
    Node(const std::shared_ptr<const STree>& stree, nidx_t nidx = ROOT) : stree_(stree), nidx_(nidx) {
        if (!stree_->validate(nidx_)) setValid(false);
    }

	/** Reset this `Node` to the root of the suffix tree.
	 */
	void setRoot() {
		nidx_ = ROOT;
        if (!stree_->validate(nidx_)) setValid(false);
	}

    /** Set this `Node` to the given `node`, which must belong to the same suffix tree. This is a faster version of `(*this) = node)`.
     */
    void set(const Node& node) { nidx_ = node.nidx_; }

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

	/** Return `true` if the subsequence represented by this node is a suffix of the underlying sequence. Note that this does not imply that this is a leaf.
	 */
	bool isSuffix() const { return isValid() and (stree_->hi(nidx_) + stree_->d(nidx_)  == stree_->sequence_.rawSize()); }

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

	/** Return the number of occurrences of the sequence represented by this node in the sequence represented by the suffix tree. For an invalid node, zero is returned.
	 */
	nidx_t count() const { return isValid() ? stree_->n(nidx_) : 0; }

	/** Return the first child node of this node. If no such node exists, a `Node` marked as invalid is returned. Note that the children are ordered lexicographically according to their edge labels.
	 */
	Node child() const { Node ret = Node(*this); ret.toChild(); return ret; }

	/** Set this node to its first child if such a node exists, otherwise mark this node as invalid. Note that the children are ordered lexicographically according to their edge labels.
	 */
	void toChild() { if (isValid() and isInternal()) nidx_ = stree_->c(nidx_); else setValid(false); }

	/** Return the child node along the edge leading away whose label begins with the given `symbol`. If no such node exists, a `Node` marked as invalid is returned.
	 */
	Node child(Symbol symbol) const { Node ret = Node(*this); ret.toChild(symbol); return ret; }

	/** Set this node to the child node along the edge leading away whose label begins with the given `symbol`\. If no such node exists, mark this node as invalid instead.
	 */
	void toChild(Symbol symbol) {
		if (isValid() and isInternal()) { nidx_t chld = stree_->c(nidx_, symbol); if (chld & VALID) nidx_ = chld; else setValid(false); }
		else setValid(false);
	}

	/** Return the next sibling of this node. If no such node exists, a `Node` marked as invalid is returned. Note that the siblings are ordered lexicographically according to their edge labels.
	 */
	Node sibling() const { Node ret = Node(*this); ret.toSibling(); return ret; }

	/** Set this node to its next sibling if a next sibling exists, otherwise mark this node as invalid. Note that the siblings are ordered lexicographically according to their edge labels.
	 */
	void toSibling() { if (isValid()) { nidx_t sbl = stree_->sib(nidx_); if (sbl & VALID) nidx_ = sbl; else setValid(false); } }

	/** Return the node corresponding to the first suffix of the represented sequence. This follows the "suffix link" of the suffix tree. If no such node exists, a `Node` marked as invalid is returned.
	 */
	Node suffix() const { Node ret = Node(*this); ret.toSuffix(); return ret; }

	/** Set this node to the node corresponding to the first suffix of the represented sequence\. If no such node exists, mark this node as invalid instead. This follows the "suffix link" of the suffix tree.
	 */
    void toSuffix();

	/** Return the (sub-)sequence represented by this node. Note that this is `seq.rawSub(headIndex(), depth())`, where `seq` is the sequence represented by the suffix tree.
	 */
	Sequence sequence() const CHECK(throw(std::invalid_argument)) {
        CHECK(if (!isValid()) throw std::invalid_argument(".sequence() called on invalid node.");)
	    return stree_->sequence().rawSub(headIndex(), depth());
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



SWIGCODE(%feature("python:slot", "tp_repr", functype="reprfunc") EdgeNode::repr;)
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
    EdgeNode(const std::shared_ptr<const STree>& stree, nidx_t nidx = ROOT, nidx_t parent = ROOT) : Node(stree, nidx), parent_(nidx &~ VALID) { findParent(Node(stree_, parent)); }

    /** Construct an `EdgeNode` from the given `Node`. This will search for the parent starting form the given `parent` (or root by default). If no parent can be found (e.g. for the root), this `EdgeNode` will be "degenerate" (have an invalid parent).
     */
    EdgeNode(const Node& node, nidx_t parent = ROOT) : Node(node), parent_(parent &~ VALID) { findParent(Node(stree_, parent)); }

    /** Construct an `EdgeNode` from the given `Node`. This will search for the parent starting form the given `parent`. If no parent can be found (e.g. for the root), this `EdgeNode` will be "degenerate" (have an invalid parent).
     */
    EdgeNode(const Node& node, const Node& parent) : Node(node), parent_(node.nidx() &~ VALID) { findParent(parent); }

    /** Reset this `EdgeNode` to the root of the suffix tree.
     */
    void setRoot() {
        Node::setRoot();
        parent_ = ROOT &~ VALID;
    }

    /** Set this `EdgeNode` to the given `edgeNode`, which must belong to the same suffix tree. This is a faster version of `(*this) = edgeNode)`.
    */
    void set(const EdgeNode& edgeNode) { Node::set(edgeNode); parent_ = edgeNode.parent_; }

    /** Return the parent `Node` of this `EdgeNode`. If this `EdgeNode` is invalid or degenerate, i.e., no parent node exists or is known, a `Node` marked as invalid is returned.
     */
	Node parent() const {
	    if (isValid()) return Node(stree_, parent_);
	    return Node(stree_, nidx_);
	}

	/** Return the first child `EdgeNode` of this `EdgeNode`. If no child exists, an `EdgeNode` marked as invalid is returned. Note that the children are ordered lexicographically according to their edge labels.
     */
 	EdgeNode child() const { EdgeNode ret = EdgeNode(*this); ret.toChild(); return ret; }

    /** Set this `EdgeNode` to its first child (with edge leading to it) if such a node exists, otherwise mark this `EdgeNode` as invalid. Note that the children are ordered lexicographically according to their edge labels.
     */
 	void toChild() { nidx_t par = nidx_; Node::toChild(); if (isValid()) parent_ = par; }

    /** Return the child `EdgeNode` leading away whose label begins with the given `symbol`. If no such `EdgeNode` exists, an `EdgeNode` marked as invalid is returned.
     */
 	EdgeNode child(Symbol symbol) const { EdgeNode ret = EdgeNode(*this); ret.toChild(symbol); return ret; }

    /** Set this `EdgeNode` to the child `EdgeNode` leading away whose label begins with the given `symbol`\. If no such `EdgeNode` exists, mark this `EdgeNode` as invalid instead.
     */
 	void toChild(Symbol symbol) { nidx_t par = nidx_; Node::toChild(symbol); if (isValid()) parent_ = par; }

    /** Return the next sibling of this `EdgeNode`. If none exists, an `EdgeNode` marked as invalid is returned. Note that the siblings are ordered lexicographically according to their edge labels.
     */
 	EdgeNode sibling() const { EdgeNode ret = EdgeNode(*this); ret.toSibling(); return ret; }

    /** Return the `EdgeNode` corresponding to the first suffix of the represented sequence. This follows the "suffix link" of the suffix tree (and finds the new parent accordingly). If no such node exists, a `Node` marked as invalid is returned.
     */
 	EdgeNode suffix() const { EdgeNode ret = EdgeNode(*this); ret.toSuffix(); return ret; }

    /** Set this to the `EdgeNode` corresponding to the first suffix of the represented sequence\. If no such node exists, mark this node as invalid instead. This follows the "suffix link" of the suffix tree (and finds the new parent accordingly).
     */
 	void toSuffix() { Node::toSuffix(); if (isValid()) findParent(Node(stree_, parent_).suffix()); }

 	/** Return the edge label.
 	 */
	Sequence label() const CHECK(throw(std::invalid_argument)) {
	    CHECK(if (!isValid()) throw std::invalid_argument(".label() called on invalid EdgeNode.");)
        Node par = parent();
	    if (!par.isValid()) return stree_->sequence().rawSub(headIndex() + depth(), 0);
	    return stree_->sequence().rawSub(headIndex() + par.depth(), depth() - par.depth());
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
        Sequence seq = sequence();
        if (not searchBelow.isValid()) searchBelow = Node(stree_, ROOT);
        while (searchBelow.isValid() and (par_depth = searchBelow.depth()) < d) {
            Node par_child = searchBelow.child(seq.rawAt(par_depth));
            if (par_child == *this) { parent_ = par_child.nidx(); return; }
            searchBelow = par_child;
        }
        parent_ = nidx_ &~VALID;
        return;
    }
}; // class STreeEdge


SWIGCODE(%feature("python:slot", "tp_repr", functype="reprfunc") PathNode::repr;)
/** This class represents a node (together with the path of nodes leading to it) in the suffix tree and can be used for extracting information or navigating the suffix tree.
 *
 * The `to...()` methods are provided to navigate the suffix tree structure:
 *
 * - `to...()` sets this node (with path leading to it) to its e.g. child, suffix, sibling, etc.
 * - If no such exists, then this node is simply marked as invalid.
 * - For an invalid node, the `to...()` methods have no effect.
 * - Calling `setValid()` after the node has been marked as invalid by a `to...()` method will reset this node (and path) to the last valid node during the traversal.
 *
 * Note that the plain methods (`child(), sibling(), parent()`, etc.) return `Node`s and not `PathNode`s. If `PathNode`s are required, use instead:
 *
 *     PathNode child = PathNode(thisPathNode);
 *     child.toChild();
 */
class PathNode : public Node {
protected:
    using Node::suffix;
    using Node::toSuffix;
    std::deque<nidx_t> path_;
public:
    /** Construct a `Node` for the given `stree` corresponding to the given `nidx`. If no `nidx` is given, it defaults to the root of the suffix tree.
     */
    PathNode(const std::shared_ptr<const STree>& stree, nidx_t nidx = ROOT) : Node(stree, nidx), path_() { findPath(); }

    /** Construct a `PathNode` corresponding to the given `node`. */
    PathNode(const Node& node) : Node(node), path_() { findPath(); }

    /** Reset this `PathNode` to the root of the suffix tree.
     */
    void setRoot() {
        Node::setRoot();
        path_.clear();
    }

    /** Set this `PathNode` to the given `pathNode`, which must belong to the same suffix tree.
    */
    void set(const PathNode& pathNode) { Node::set(pathNode); path_ = pathNode.path_; }

    /** Extend this `PathNode` to its first child if such a node exists, otherwise mark this `PathNode` as invalid instead. Note that the children are ordered lexicographically according to their edge labels.
     */
    void toChild() { nidx_t current = nidx_; Node::toChild(); if (isValid()) path_.push_back(current); }

    /** Extend this `PathNode` to the child node along the edge leading away whose label begins with the given `symbol`\. If no such node exists, mark this `PathNode` as invalid instead.
     */
	void toChild(Symbol chr) { nidx_t current = nidx_; Node::toChild(chr); if (isValid()) path_.push_back(current); }

	/** Return the parent `Node`. If none exists, return a `Node` marked as invalid.
	 */
	Node parent() const {
	    if (path_.empty()) return Node(stree_, INTERNAL);
	    return Node(stree_, path_.back());
	}

	/** Set this `PathNode` to the path to the parent of the current node, if such exists\. Otherwise, just mark this `PathNode` as invalid.
	 */
    void toParent() {
        if (isValid()) {
            if (path_.empty()) setValid(false);
            else { nidx_ = path_.back(); path_.pop_back(); }
        }
    }

    /** Return the edge label for the edge leading to the current node. If no edge exists, this will be an empty `Sequence`.
     */
	Sequence label() const CHECK(throw(std::invalid_argument)) {
        CHECK(if (!isValid()) throw std::invalid_argument(".label() called on invalid PathNode.");)
        Node par = parent();
        if (!par.isValid()) return stree_->sequence().rawSub(headIndex(), 0);
        return stree_->sequence_.rawSub(headIndex() + par.depth(), depth() - par.depth()); }

    /** Return a string representation to display in python.
     */
    std::string repr() const {
        std::stringstream s;
        s << ( (nidx_ & VALID)    ? "" : "invalid " )
          << ( (nidx_ & INTERNAL) ? "internal " : "leaf ")
          << "node with path information";
        return s.str();
    }

    /** Set this to the `PathNode` corresponding to the first suffix of the represented sequence\. If no suffix exists (i.e., this is the root) mark this `PathNode` as invalid instead. This uses the "suffix link" of the suffix tree, but needs to recompute the path.
     */
    void toSuffix() {
        Node::toSuffix();
        findPath();
    }

private:
    /** Find and set the path leading to this node. */
    void findPath() {
        if (not isValid()) { return; }
        nidx_t d = depth();
        Sequence seq = sequence();
        Node ancestor = Node(stree_, ROOT);
        path_.clear();
        while (ancestor.isValid() and (ancestor.depth()) < d) {
            path_.push_back(ancestor.nidx());
            ancestor.toChild(seq.rawAt(ancestor.depth()));
        }
    }
}; // class STreePath



SWIGCODE(%feature("python:slot", "tp_repr", functype="reprfunc") Position::repr;)
/** This class represents any position in the suffix tree. A `Position` that is a node (leaf or internal) of the suffix tree is considered "explicit", while a `Position` on some edge is called "implicit".
 *
 * The `to...()` methods are provided to navigate the suffix tree structure:
 *
 * - `to...()` sets this `Position` to its suffix, extension by a sequence, etc.
 * - If no such exists, then it is simply marked as invalid.
 * - For an invalid `Position`, the `to...()` methods have no effect.
 * - Calling `setValid()` after the `Position` has been marked as invalid by a `to...()` method will reset the `Position` to the last valid one during the traversal.
 */
class Position {
private:
    EdgeNode edge_; /**< the edge on which this position lies */
    nidx_t depth_; /**< the depth of the position, i.e., the size of the represented subsequence */
public:
    /** Construct the `Position` in the given `stree` corresponding to the given `sequence`. If no corresponding position exists, this `Position` will correspond to the longest possible prefix of the provided `sequence` and will be marked as invalid.
     */
	Position(const std::shared_ptr<const STree>& stree, Sequence sequence = Sequence()) : edge_(stree), depth_(0) { toSequence(sequence); }

	/** Construct a `Position` from the given `node` of type `EdgeNode`.
	 */
	Position(const EdgeNode& node) : edge_(node), depth_(0) {
	    if (isValid()) depth_ = node.depth();
	}

	/** Reset this `Position` to the root of the suffix tree.
	 */
	void setRoot() { edge_.setRoot(); depth_ = 0; }

    /** Set this `Position` to the given `position`, which must belong to the same suffix tree. This is a faster version of `(*this) = position)`.
    */
    void set(const Position& position) { edge_.set(position.edge_); depth_ = position.depth_; }

    /** Return `true` if valid, otherwise return `false`.
     */
	bool isValid() const { return edge_.isValid(); }

    /** Mark this `Position` as `valid` (default: `true`).
     */
	void setValid(bool valid = true) { edge_.setValid(valid); }

    /** Return `true` if this is a node (internal or leaf). Otherwise, this `Position` is "implicit", i.e., it lies on some edge.
     */
	bool isExplicit() const { return (depth_ == edge_.depth()); }

    /** Return `true` if this is an internal node or an implicit position.
     */
	bool isInternal() const { return !isLeaf(); }

    /** Return `true` if this is a leaf node.
     */
	bool isLeaf() const { return (isExplicit() and edge_.isLeaf()); }

    /** Return `true` if this is the root node.
     */
    bool isRoot() const { return edge_.isRoot(); }

	/** Return `true` if the represented subsequence is a suffix of the underlying sequence. Note that this does not imply that this is a leaf.
 	*/
	bool isSuffix() const { return isValid() and (edge_.headIndex() + depth_ == edge_.stree_->sequence_.rawSize()); }

	//<editor-fold desc="Comparison operators">
	/** Return `true` if this `Position` is the same as the given `other`.
     */
    bool operator ==(const Position& other) const { return edge_ == other.edge_ and depth_ == other.depth_; }

    /** Return `true` if this `Position` is the same as the given `other`.
     */
    bool operator ==(const Node& other) const { return edge_ == other and isExplicit(); }

    /** Return `true` if the `other` `Position` occurs more often in the represented sequence than this node.
     */
    bool operator <(const Position& other) const CHECK(throw(std::invalid_argument)) {
        CHECK(if (!isValid() or !other.isValid()) throw std::invalid_argument("Cannot compare invalid nodes.");)
        return count() < other.count();
    }

    /** Return `true` if the `other` `Node` occurs more often in the represented sequence than this node.
     */
    bool operator <(const Node& other) const CHECK(throw(std::invalid_argument)) {
        CHECK(if (!isValid() or !other.isValid()) throw std::invalid_argument("Cannot compare invalid nodes.");)
        return count() < other.count();
    }

    /** Return the "depth" of this position, which is the size of the represented (sub-)sequence.
     */
    nidx_t depth() const CHECK(throw(std::invalid_argument)) {
        CHECK(if (!isValid()) throw std::invalid_argument(".depth() called on invalid position.");)
        return depth_;
    }

    /** Return the "headindex" of this node, which is an index in the sequence represented by the suffix tree where the (sub-)sequence represented by this node occurs. I.e., the (sub-)sequence represented by this position is `seq.rawSub(headindex(), depth())`, where `seq` is the sequence represented by the suffix tree.
     */
    nidx_t headIndex() const CHECK(throw(std::invalid_argument)) { return edge_.headIndex(); }

    /** Return the number of occurrences of the sequence represented by this `Position` in the sequence represented by the suffix tree. For an invalid `Position`, zero is returned.
     */
    nidx_t count() const { return edge_.count(); }

	/** Return the `EdgeNode` that this `Position` lies on. If this `Position` is explicit, then the `EdgeNode` will be the node (with edge leading to it) of the position.
	 */
    EdgeNode edge() { return edge_; }

    /** Return the `Position` corresponding to the first suffix of the represented sequence. This uses the "suffix link" of the suffix tree. If no suffix exists (i.e., this is the root), a `Position` marked as invalid is returned.
     */
    Position suffix() const { Position ret = Position(*this); ret.toSuffix(); return ret; }

    /** Set this to the `Position` corresponding to the first suffix of the represented sequence\. If no suffix exists (i.e., this is the root) mark this `Position` as invalid instead. This uses the "suffix link" of the suffix tree.
     */
	void toSuffix() {
	    if (not isValid()) { return; }
	    Node par = edge_.parent();
	    if (not par.isValid()) { /* root case */ setValid(false); return; }
        par.toSuffix(); if (not par.isValid()) par = Node(edge_.stree_, ROOT);
        Sequence seq = sequence().slice(1, tom::NoIndex);
        depth_ = seq.rawSize();
        edge_ = EdgeNode(par, par);
		while (isValid() and (edge_.depth() < depth_)) {
			edge_.toChild(seq.rawAt(edge_.depth()));
		}
		assert(isValid());
	}

	/** Update this `Position` such that it represents a subsequence extended by the given `symbol`. If no such position exists, this `Position` is unchanged but marked as invalid. For an invalid `Position` this function has no effect.
	 */
	void toSymbol(Symbol symbol) {
		if (!isValid()) return;
		if (isExplicit()) {
			edge_.toChild(symbol);
			if (isValid()) depth_++;
		} else {
			if (edge_.stree_->sequence().rawAt(edge_.headIndex() + depth_) == symbol) ++depth_;
			else edge_.setValid(false);
		}
	}

	/** If this `Position` is not explicit, i.e., it lies on an edge of the suffix tree structure, then set this `Position` to the deeper node end-point of that edge, making this position explicit.
	 */
	void toExplicit() {
		if (!isValid()) return;
		depth_ = edge_.depth();
	}

	/** Return the first child `Position` in the suffix tree structure viewed as a *suffix trie*, i.e., where all positions are seen as nodes and all edges have length one. If no child exists, a `Position` marked as invalid is returned. Note that the children are ordered lexicographically according to their edge symbols.
	 */
	Position child() const { Position ret = Position(*this); ret.toChild(); return ret; }

	/** Set this `Position` to its first child position in the suffix tree structure viewed as a *suffix trie*, i.e., where all positions are seen as nodes and all edges have length one. If no child exists, mark this `Position` as invalid instead. Note that the children are ordered lexicographically according to their edge symbols.
	 */
	void toChild() {
		if (!isValid()) return;
		if (isExplicit()) {
			edge_.toChild();
			if (isValid()) depth_++;
		} else { ++depth_; }
	}

	/** Return the next sibling `Position` in the suffix tree structure viewed as a *suffix trie*, i.e., where all positions are seen as nodes and all edges have length one. If no sibling exists, a `Position` marked as invalid is returned. Note that the siblings are ordered lexicographically according to their edge symbols.
 	*/
	Position sibling() const { Position ret = Position(*this); ret.toSibling(); return ret; }

	/** Set this `Position` to its next sibling position in the suffix tree structure viewed as a *suffix trie*, i.e., where all positions are seen as nodes and all edges have length one. If no sibling exists, mark this `Position` as invalid instead. Note that the siblings are ordered lexicographically according to their edge symbols.
 	*/
	void toSibling() {
		if (!isValid()) return;
		if (depth_ == edge_.parent().depth() + 1) { edge_.toSibling(); }
		else setValid(false);
	}

    /** Update this `Position` such that it represents a subsequence extended by the given `sequence`. If no such position exists, this `Position` is updated symbol-wise according to the given `sequence` as far as possible and then marked as invalid. For an invalid `Position` this function has no effect.
     */
	void toSequence(const Sequence& sequence) {
		for (nidx_t idx = 0; ((idx < sequence.rawSize()) and (edge_.isValid())); ++idx) { toSymbol(sequence.rawAt(idx)); }
	}

    /** Return the (sub-)sequence represented by this `Position`. Note that this is `seq.rawSub(headIndex(), depth())`, where `seq` is the sequence represented by the suffix tree.
     */
	Sequence sequence() const { return edge_.stree_->sequence().rawSub(headIndex(), depth_); }

    /** Return the sub-sequence of the edge label up to this position. For an explicit position this is just the edge label (see `edge.label()`).
     */
    Sequence label() const CHECK(throw(std::invalid_argument)) {
        CHECK(if (!isValid()) throw std::invalid_argument(".label() called on invalid Position.");)
        Node par = edge_.parent();
        if (!par.isValid()) return edge_.stree_->sequence().rawSub(headIndex(), depth_);
        return edge_.stree_->sequence().rawSub(headIndex() + par.depth(), depth_ - par.depth());
    }

    /** Return a string representation to display in python.
     */
    std::string repr() const {
        std::stringstream s;
        s << ( isValid()  ? "" : "invalid " )
          << "suffix tree position";
        return s.str();
    }

}; // class STreePos



/* Implementation */
#ifndef SWIG

inline void Node::toSuffix() {
    if (isValid()) {
        if (isInternal()) {
            if ((index() != 0)) {
                nidx_ = stree_->sl(nidx_) | VALID;
            } else {
                nidx_ = ROOT & ~VALID;
            }
        } else {
            /* Leaf case */
            if (index() < stree_->nLeafNodes() - 1) {
                ++nidx_;
            } else {
                /* "last" leaf. Suffix will be an internal node. */
                Position pos = Position(stree_);
                pos.toSequence(sequence().slice(1, tom::NoIndex));
                nidx_ = pos.edge().nidx();
            }
        }
    } else {
        nidx_ = nidx_ & ~VALID;
    }
}

#endif /* SWIG */

} // namespace stree

#endif // STREE_NODE_H
