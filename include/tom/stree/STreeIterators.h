/*
 * @file   STreeIterators.h
 * @author Michael Thon
 */

#ifndef STREE_ITERATORS_H
#define STREE_ITERATORS_H

#include "stree.h"

namespace stree {

class PrefixIterator : public Path {
public:
	PrefixIterator(const STree* stree) : Path(stree) {}
	void next() { assert(isValid());
		toChild();
		if (isValid()) return;
		setValid(); toSibling();
		if (isValid()) return;
		while(!isValid()) {
			setValid(); toParent();
			if (!isValid()) return;
			toSibling();
		}
		return;
	}
};

class PostfixIterator : public Path {
public:
	PostfixIterator(const STree* stree) : Path(stree) {
		while(isValid()) toChild();
		setValid();
	}
	void next() { assert(isValid());
		toSibling();
		if (!isValid()) { setValid(); toParent(); return; }
		while(isValid()) toChild();
		setValid();
		return;
	}
};

class DFSIterator : public Path {
public:
	DFSIterator(const STree* stree) : Path(stree), firstVisit(true) {}
	void next() {
		assert(isValid());
		if (firstVisit) {
			toChild();
			if (isValid()) return;
			setValid(); toSibling();
			if (isValid()) return;
			firstVisit = false;
			setValid(); toParent();
			return;
		}
		else {
			toSibling();
			if (isValid()) { firstVisit = true; return; }
			setValid(); toParent();
			return;
		}
	}
	bool isFirstVisit() { return firstVisit; }
	void setUpPass() { firstVisit = false; }
private:
	bool firstVisit;
};

} // namespace stree


#endif // STREE_ITERATORS_H
