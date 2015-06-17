/*
 * @file   STreeIterators.h
 * @author Michael Thon
 */

#ifndef _STREE_ITERATORS_H_
#define _STREE_ITERATORS_H_

#include "STreeNode.h"

namespace stree {

class PrefixIterator : public STreePath {
public:
	PrefixIterator(const STree* stree) : STreePath(stree) {}
	void next() { assert(isValid());
		child();
		if (isValid()) return;
		setValid(); sibling();
		if (isValid()) return;
		while(!isValid()) {
			setValid(); parent();
			if (!isValid()) return;
			sibling();
		}
		return;
	}
};

class PostfixIterator : public STreePath {
public:
	PostfixIterator(const STree* stree) : STreePath(stree) {
		while(isValid()) child();
		setValid();
	}
	void next() { assert(isValid());
		sibling();
		if (!isValid()) { setValid(); parent(); return; }
		while(isValid()) child();
		setValid();
		return;
	}
};

class DFSIterator : public STreePath {
public:
	DFSIterator(const STree* stree) : STreePath(stree), firstVisit(true) {}
	void next() {
		assert(isValid());
		if (firstVisit) {
			child();
			if (isValid()) return;
			setValid(); sibling();
			if (isValid()) return;
			firstVisit = false;
			setValid(); parent();
			return;
		}
		else {
			sibling();
			if (isValid()) { firstVisit = true; return; }
			setValid(); parent();
			return;
		}
	}
	bool isFirstVisit() { return firstVisit; }
	void setUpPass() { firstVisit = false; }
private:
	bool firstVisit;
};

} // namespace stree


#endif /* _STREE_ITERATORS_H_ */
