/**
 * @file   STreeNode.h
 * @brief  This file provides several objects for navigating within the suffix tree structure.
 */

#ifndef STREE_NODE_H
#define STREE_NODE_H

namespace stree {

/**
 * This class represents a node in the suffix tree\. It contains a pointer internally to the \a STree that it belongs to.
 */
class STreeNode {
public:
  /** construct an uninitialized and invalid \a STreeNode */
	STreeNode() : stree_(NULL), nidx_(NODE) {}
	/** construct a \a STreeNode initialized as the root of the given \a stree */
	STreeNode(const STree* stree) : stree_(stree), nidx_(ROOT) {}
	/** construct a \a STreeNode initialized to the given \a node \a Nidx. */
	STreeNode(const STree* stree, Nidx node) : stree_(stree), nidx_(node) {}

	/** return \c true if valid, otherwise return \c false. */
	bool isValid() const { return stree_ and (nidx_ & VALID); }
	void setValid(bool valid = true) { if (valid) nidx_ |= VALID; else nidx_ &= ~VALID;}
	bool isNode() const { return (nidx_ & NODE); }
	bool isLeaf() const { return !isNode(); }
  bool isRoot() const { return isNode() and index() == 0; }

	Idx index() const { return nidx_ & INDEX; }
	Idx nodeIndex() const { return nidx_ & (INDEX | NODE); }
	Nidx nidx() const { return nidx_; }

	bool operator ==(const STreeNode& other) const { return ((nidx_ | COLOR) == (other.nidx() | COLOR)); }

	bool operator <(const STreeNode& other) const { return count() < other.count(); }

	Idx depth() const { return stree_->d(nidx_); }
	Idx headIndex() const { return stree_->hi(nidx_); }
	Idx count() const { return stree_->n(nidx_); }

	STreeNode getChild() const { if (isValid() and isNode()) return STreeNode(stree_, stree_->c(nidx_)); else return STreeNode(); }
	void child() { if (isValid() and isNode()) nidx_ = stree_->c(nidx_); else setValid(false); }
	STreeNode getChild(Char chr) const { if (isValid() and isNode()) return STreeNode(stree_, stree_->c(nidx_, chr)); else return STreeNode(); }
	void child(Char chr) {
		if (isValid() and isNode()) { Nidx chld = stree_->c(nidx_, chr); if (chld & VALID) nidx_ = chld; else setValid(false); }
		else setValid(false);
	}

	STreeNode getSibling() const { if (isValid()) return STreeNode(stree_, stree_->sib(nidx_)); else return STreeNode(); }
	void sibling() { if (isValid()) { Nidx sbl = stree_->sib(nidx_); if (sbl & VALID) nidx_ = sbl; else setValid(false); } }

	STreeNode getSuffixLink() const {
		if (isValid() and isNode() and (index() != 0)) return STreeNode(stree_, stree_->sl(nidx_) | VALID);
		else return STreeNode();
	}
	void suffixLink() {
		if (isValid() and isNode() and (index() != 0)) { nidx_ = stree_->sl(nidx_) | VALID; }
		else { setValid(false); }
	}

	String string() const { return stree_->text_.rawSub(headIndex(), depth()); }
	String label(Idx parentDepth) const {
		assert(parentDepth <= depth());
		return stree_->text_.rawSub(headIndex() + parentDepth, depth() - parentDepth);
	}

	std::string dataStr(int width = 5) const {
		assert(isValid());
		std::stringstream s;
		s << "[ " << indexStr(width)
			<< " | hi = " << std::setw(4) << headIndex()
			<< " | d = " << std::setw(4) << depth()
			<< " | n = " << std::setw(4) << count();
		if (isNode()) {
			s << " | c = " << getChild().indexStr(width)
				<< " | sl = " << getSuffixLink().indexStr(width);
		}
		s << " ]";
		return s.str();
	}

	std::string indexStr(int width = 5) const {
		std::stringstream s;
		s << (isValid() ? (isNode() ? "N" : "L") : (isNode() ? "n" : "l"))
			<< ((nidx_ & COLOR) ? std::setfill(':') : std::setfill('.'))
			<< std::setw(width-2) << (nidx_ & INDEX);
		return s.str();
	}

protected:
	const STree* stree_;
	Nidx nidx_;
}; // class STreeNode



class STreeEdge : public STreeNode {
	friend class STreePos;
public:
	STreeEdge() : STreeNode(), parent_(NODE) {}
	STreeEdge(const STree* stree) : STreeNode(stree), parent_(NODE) {}
	STreeEdge(const STree* stree, Nidx nidx, Nidx parent) : STreeNode(stree, nidx), parent_(parent) {}
	STreeEdge(const STreeNode& node, Nidx parent) : STreeNode(node), parent_(parent) {}

	Idx parentDepth() const { return stree_->d(this->parent_); }
	STreeNode getParent() const { return STreeNode(stree_, parent_); }
 	STreeEdge getChild() const { return STreeEdge(STreeNode::getChild(), nidx_); }
	void child() { Nidx par = nidx_; STreeNode::child(); if (isValid()) parent_ = par; }
	STreeEdge getChild(Char chr) const { return STreeEdge(STreeNode::getChild(chr), nidx_); }
	void child(Char chr) { Nidx par = nidx_; STreeNode::child(chr); if (isValid()) parent_ = par; }
	STreeEdge getSibling() const { return STreeEdge(STreeNode::getSibling(), parent_); }

	String label() const { return STreeNode::label(parentDepth()); }

private:
	using STreeNode::getSuffixLink;
	using STreeNode::suffixLink;
protected:
	Nidx parent_;
}; // class STreeEdge



class STreePath : public STreeNode {
public:
	STreePath() : STreeNode(), path() {}
	STreePath(const STree* stree) : STreeNode(stree), path() {}

	void child() { Nidx current = nidx_; STreeNode::child(); if (isValid()) path.push_back(current); }
	void child(Char chr) { Nidx current = nidx_; STreeNode::child(chr); if (isValid()) path.push_back(current); }
	STreeNode getParent() const { if (path.empty()) return STreeNode(stree_, NODE); return STreeNode(stree_, path.back()); }
	STreeNode getAncestor(Idx generations) const {
		if (generations > path.size()) return STreeNode(stree_, NODE);
		if (generations == 0) return *this;
		return STreeNode(stree_, path.at(path.size()-generations));
	}
	void parent() { if (isValid()) {
			if (path.empty()) setValid(false);
			else { nidx_ = path.back(); path.pop_back(); }
		}}
	Idx nAncestors() const { return path.size(); }
	Idx parentDepth() const { if (path.empty()) return 0; else return stree_->d(path.back()); }

	String label() const { return STreeNode::label(parentDepth()); }

private:
	using STreeNode::getSuffixLink;
	using STreeNode::suffixLink;
protected:
	std::deque<Nidx> path;
}; // class STreePath



class STreePos {
public:
	STreePos() : edge_(), depth_(0) {}
	STreePos(const STree* stree) : edge_(stree), depth_(0) {}

	void setRoot() { edge_ = STreeEdge(edge_.stree_); depth_ = 0; }

	bool isValid() const { return edge_.isValid(); }
	void setValid(bool valid = true) { edge_.setValid(valid); }
	bool isExplicit() const { return (depth_ == edge_.depth()); }
	bool isLeaf() const { return (isExplicit() and edge_.isLeaf()); }

	Idx count() const { return edge_.count(); }
	Idx headIndex() const { return edge_.headIndex(); }
	Idx depth() const { return depth_; }
	Idx parentDepth() const { return edge_.parentDepth(); }

  STreeEdge& edge() { return edge_; }

	void suffixLink() { if (isValid()) {
			if (edge_.stree_->symbolSize_ > depth_) { setValid(false); return; }
			depth_-= edge_.stree_->symbolSize_;
			Idx hi = headIndex() + edge_.stree_->symbolSize_;
			edge_.parent_ = edge_.stree_->sl(edge_.parent_) | VALID;
			edge_.nidx_ = edge_.parent_;
			edge_.setValid();
			while (edge_.isValid() and (edge_.depth() < depth_))
				edge_.child(edge_.stree_->at(hi + edge_.depth()));
		}}

	void addChar(Char chr) {
		if (!isValid()) return;
		if (isExplicit()) {
			edge_.child(chr);
			if (isValid()) depth_++;
		}
		else
			if (edge_.stree_->at(edge_.headIndex() + depth_) == chr) depth_++;
			else edge_.setValid(false);
	}
	void addString(const String& str) {
		for (Idx pos = 0; ((pos < str.rawSize()) and (edge_.isValid())); ++pos)
			addChar(str.rawAt(pos));
	}

	String string() const { return edge_.stree_->text_.rawSub(headIndex(), depth_); }
	String label() const { return edge_.stree_->text_.rawSub(headIndex() + parentDepth(), depth_ - parentDepth()); }

private:
	STreeEdge edge_;
	Idx depth_;

}; // class STreePos

} // namespace stree

#endif // STREE_NODE_H
