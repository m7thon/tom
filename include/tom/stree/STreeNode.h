#ifndef STREE_NODE_H
#define STREE_NODE_H

#include "stree.h"

namespace stree {

/** This class represents a node in the suffix tree\. It contains a pointer internally to the \c STree that it belongs to.
 */
class Node {
public:
	/** construct an uninitialized and invalid \c Node */
	Node() : stree_(NULL), nidx_(INTERNAL) {}
	/** construct a \c Node initialized as the root of the given \c stree
	 */
	Node(const STree* stree) : stree_(stree), nidx_(ROOT) {}
	/** construct a \c Node initialized to the given \c nidx.
	 */
	Node(const STree* stree, nidx_t nidx) : stree_(stree), nidx_(nidx) {}

	/** return \c true if valid, otherwise return \c false.
	 */
	bool isValid() const { return stree_ and (nidx_ & VALID); }

	/** mark this \c Node as \c valid (or invalid, if \c valid is \c false).
	 */
	void setValid(bool valid = true) { if (valid) nidx_ |= VALID; else nidx_ &= ~VALID;}
	bool isInternal() const { return (nidx_ & INTERNAL); }
	bool isLeaf() const { return !isInternal(); }
	bool isRoot() const { return isInternal() and index() == 0; }

	nidx_t index() const { return nidx_ & INDEX; }
	nidx_t nodeIndex() const { return nidx_ & (INDEX | INTERNAL); }
	nidx_t nidx() const { return nidx_; }

	bool operator ==(const Node& other) const { return ((nidx_ | COLOR) == (other.nidx() | COLOR)); }
	bool operator <(const Node& other) const { return count() < other.count(); }

	nidx_t depth() const { return stree_->d(nidx_); }
	nidx_t headIndex() const { return stree_->hi(nidx_); }
	nidx_t count() const { return stree_->n(nidx_); }

	Node child() const { if (isValid() and isInternal()) return Node(stree_, stree_->c(nidx_)); else return Node(); }
	void toChild() { if (isValid() and isInternal()) nidx_ = stree_->c(nidx_); else setValid(false); }
	Node child(Symbol chr) const { if (isValid() and isInternal()) return Node(stree_, stree_->c(nidx_, chr)); else return Node(); }
	void toChild(Symbol chr) {
		if (isValid() and isInternal()) { nidx_t chld = stree_->c(nidx_, chr); if (chld & VALID) nidx_ = chld; else setValid(false); }
		else setValid(false);
	}

	Node sibling() const { if (isValid()) return Node(stree_, stree_->sib(nidx_)); else return Node(); }
	void toSibling() { if (isValid()) { nidx_t sbl = stree_->sib(nidx_); if (sbl & VALID) nidx_ = sbl; else setValid(false); } }

	Node suffixlink() const {
		if (isValid() and isInternal() and (index() != 0)) return Node(stree_, stree_->sl(nidx_) | VALID);
		else return Node();
	}
	void toSuffixlink() {
		if (isValid() and isInternal() and (index() != 0)) { nidx_ = stree_->sl(nidx_) | VALID; }
		else { setValid(false); }
	}

	Sequence asSequence() const { return stree_->sequence_.rawSub(headIndex(), depth()); }

	std::string dataStr(int width = 5) const {
        std::stringstream s;
		if (not isValid()) {
		    s << "[ " << indexStr(width) << " ]";
		    return s.str();
		}
		s << "[ " << indexStr(width)
          << " | hIdx = " << std::setw(width) << headIndex()
          << " | size = "  << std::setw(width) << depth()
          << " | nOcc = "  << std::setw(width) << count();
		if (isInternal()) {
			s << " | chld = "  << child().indexStr(width)
              << " | sfxL = " << suffixlink().indexStr(width);
		}
		s << " ]";
		return s.str();
	}

	std::string indexStr(int width = 3) const {
		std::stringstream s;
		s << ( (nidx_ & VALID)    ? '+' : '-' )
          << ( (nidx_ & INTERNAL) ? '+' : '-' )
          << ( (nidx_ & COLOR)    ? '+' : '-' )
          << std::setw(width) << std::setfill('0') << (nidx_ & INDEX);
		return s.str();
	}

protected:
	const STree* stree_;
	nidx_t nidx_;
}; // class STreeNode



class Edge : public Node {
	friend class Pos;
public:
	Edge() : Node(), parent_(INTERNAL) {}
	Edge(const STree* stree) : Node(stree), parent_(INTERNAL) {}
	Edge(const STree* stree, nidx_t nidx, nidx_t parent) : Node(stree, nidx), parent_(parent) {}
	Edge(const Node& node, nidx_t parent) : Node(node), parent_(parent) {}

	nidx_t parentDepth() const { return stree_->d(this->parent_); }
	Node parent() const { return Node(stree_, parent_); }
 	Edge child() const { return Edge(Node::child(), nidx_); }
	void toChild() { nidx_t par = nidx_; Node::toChild(); if (isValid()) parent_ = par; }
	Edge child(Symbol chr) const { return Edge(Node::child(chr), nidx_); }
	void toChild(Symbol chr) { nidx_t par = nidx_; Node::toChild(chr); if (isValid()) parent_ = par; }
	Edge sibling() const { return Edge(Node::sibling(), parent_); }

	Sequence edgeLabel() const { return stree_->sequence_.rawSub(headIndex() + parentDepth(), depth() - parentDepth()); }

private:
	using Node::suffixlink;
	using Node::toSuffixlink;
protected:
	nidx_t parent_;
}; // class STreeEdge



class Path : public Node {
public:
	Path() : Node(), path() {}
	Path(const STree* stree) : Node(stree), path() {}

	void toChild() { nidx_t current = nidx_; Node::toChild(); if (isValid()) path.push_back(current); }
	void toChild(Symbol chr) { nidx_t current = nidx_; Node::toChild(chr); if (isValid()) path.push_back(current); }
	Node parent() const { if (path.empty()) return Node(stree_, INTERNAL); return Node(stree_, path.back()); }
	Node ancestor(nidx_t generations) const {
		if (generations > path.size()) return Node(stree_, INTERNAL);
		if (generations == 0) return *this;
		return Node(stree_, path.at(path.size()-generations));
	}
	void toParent() { if (isValid()) {
			if (path.empty()) setValid(false);
			else { nidx_ = path.back(); path.pop_back(); }
		}}
	nidx_t nAncestors() const { return path.size(); }
	nidx_t parentDepth() const { if (path.empty()) return 0; else return stree_->d(path.back()); }

	Sequence edgeLabel() const { return stree_->sequence_.rawSub(headIndex() + parentDepth(), depth() - parentDepth()); }

private:
	using Node::suffixlink;
	using Node::toSuffixlink;
protected:
	std::deque<nidx_t> path;
}; // class STreePath



class Pos {
public:
	Pos() : edge_(), depth_(0) {}
	Pos(const STree* stree) : edge_(stree), depth_(0) {}

	void setRoot() { edge_ = Edge(edge_.stree_); depth_ = 0; }

	bool isValid() const { return edge_.isValid(); }
	void setValid(bool valid = true) { edge_.setValid(valid); }
	bool isExplicit() const { return (depth_ == edge_.depth()); }
    bool isInternal() const { return !isLeaf(); }
	bool isLeaf() const { return (isExplicit() and edge_.isLeaf()); }

	nidx_t count() const { return edge_.count(); }
	nidx_t headIndex() const { return edge_.headIndex(); }
	nidx_t depth() const { return depth_; }
	nidx_t parentDepth() const { return edge_.parentDepth(); }

  Edge& edge() { return edge_; }

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

	Sequence asSequence() const { return edge_.stree_->sequence_.rawSub(headIndex(), depth_); }
	Sequence edgeLabel() const { return edge_.stree_->sequence_.rawSub(headIndex() + parentDepth(), depth_ - parentDepth()); }

private:
	Edge edge_;
	nidx_t depth_;

}; // class STreePos

} // namespace stree

#endif // STREE_NODE_H
