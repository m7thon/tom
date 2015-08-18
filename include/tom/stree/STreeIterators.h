#ifndef STREE_ITERATORS_H
#define STREE_ITERATORS_H

#include "stree.h"

namespace stree {


/** An iterator to traverse the nodes (leaf and internal) of the suffix tree in *prefix* order, which means traversing the tree depth first and visiting the nodes on the downward pass. A `PrefixIterator` is a `PathNode` with a `toNext()` method to move to the next node in prefix order. The end of iteration is signaled by marking the iterator as invalid, which can be checked by the inherited method `isValid()`.

Note that this iterator is wrapped as a native python iterator, i.e., the following is possible:

    for node in PrefixIterator(stree):
        print(node.nidxStr())
 */
class PrefixIterator : public PathNode {
public:
    /** Create a `PrefixIterator` for the given `stree`.
     */
	PrefixIterator(const std::shared_ptr<const STree>& stree) : PathNode(stree) {}

	/** Set this `PrefixIterator` to the next node in prefix order. If none exists, this `PrefixIterator` will be marked as invalid. Calling `toNext()` on an invalid `PrefixIterator` has no effect.
	 */
	void toNext() {
	    if (not isValid()) return;
	    toChild();
		if (isValid()) return;
		setValid(); toSibling();
		if (isValid()) return;
		while(!isValid()) {
			setValid(); toParent();
			if (!isValid()) return;
			toSibling();
		}
	}
};

#ifdef SWIG
%extend PrefixIterator {
    %feature("python:slot", "tp_iter", functype="getiterfunc") __iter__;
    stree::PrefixIterator __iter__() { return stree::PrefixIterator(*$self); }
    %feature("python:slot", "tp_iternext", functype="iternextfunc") __next__;
    Node __next__() throw(swig::stop_iteration) {
        if (not $self->isValid()) { throw swig::stop_iteration(); }
        stree::Node currentItem(*$self);
        $self->toNext();
        return currentItem;
    }
};
#endif // SWIG



/** An iterator to traverse the nodes (leaf and internal) of the suffix tree in *postfix* order, which means traversing the tree depth first and visiting the nodes on the upward pass. A `PostfixIterator` is a `PathNode` with a `toNext()` method to move to the next node in postfix order. The end of iteration is signaled by marking the iterator as invalid, which can be checked by the inherited method `isValid()`.

Note that this iterator is wrapped as a native python iterator, i.e., the following is possible:

    for node in PostfixIterator(stree):
        print(node.nidxStr())
 */
class PostfixIterator : public PathNode {
public:
    /** Create a `PostfixIterator` for the given `stree`.
     */
	PostfixIterator(const std::shared_ptr<const STree>& stree) : PathNode(stree) {
		while(isValid()) toChild();
		setValid();
	}

    /** Set this `PostfixIterator` to the next node in postfix order. If none exists, this `PostfixIterator` will be marked as invalid. Calling `toNext()` on an invalid `PostfixIterator` has no effect.
     */
	void toNext() {
	    if (not isValid()) return;
		toSibling();
		if (!isValid()) { setValid(); toParent(); return; }
		while(isValid()) toChild();
		setValid();
		return;
	}
};

#ifdef SWIG
%extend PostfixIterator {
    %feature("python:slot", "tp_iter", functype="getiterfunc") __iter__;
    stree::PostfixIterator __iter__() { return stree::PostfixIterator(*$self); }
    %feature("python:slot", "tp_iternext", functype="iternextfunc") __next__;
    Node __next__() throw(swig::stop_iteration) {
        if (not $self->isValid()) { throw swig::stop_iteration(); }
        stree::Node currentItem(*$self);
        $self->toNext();
        return currentItem;
    }
};
#endif // SWIG



/** An iterator to traverse the nodes (leaf and internal) of the suffix tree depth first. This means that every node (leaf and internal) is visited exactly twice (once on the downward pass, and once on the upward pass). The method `isFirstVisit()` can be used to check if the current node is being visited for the first or second time. Furthermore, the method `setUpPass()` can be used to stop traversing deeper into the current tree branch. A `DFSIterator` is a `PathNode` with a `toNext()` method to move to the next node of the depth first traversal. The end of iteration is signaled by marking the iterator as invalid, which can be checked by the inherited method `isValid()`.

Note that this iterator is not wrapped as a native python iterator, i.e., the following is *not* possible:

    for node in DFSIterator(stree):
        print(node.nidxStr())
 */
class DFSIterator : public PathNode {
public:
    /** Create a `DFSIterator` for the given `stree`.
     */
	DFSIterator(const std::shared_ptr<const STree>& stree) : PathNode(stree), firstVisit(true) {}

    /** Set this `DFSIterator` to the next node in the depth first traversal. If the current node has been visited for the first time and `setUpPass()` was not called, the next node will be the next in prefix order (or the same leaf node, if the current node is a leaf node), otherwise the next in postfix order. Note that every node (including leaf nodes) is visited exactly twice, i.e., once on the downwards pass and once on the upwards pass. If no next node exists, this `DFSIterator` will be marked as invalid. Calling `toNext()` on an invalid `DFSIterator` has no effect.
     */
	void toNext() {
        if (not isValid()) return;
		if (firstVisit) {
		    if (isLeaf()) { firstVisit = false; return; }
			toChild();
			if (isValid()) return;
			setValid(); toSibling();
			if (isValid()) return;
			firstVisit = false;
			setValid(); toParent();
		} else {
			toSibling();
			if (isValid()) { firstVisit = true; return; }
			setValid(); toParent();
		}
		return;
	}

	/** Return `true` if the current node is being visited for the first time (i.e., on the downward pass).
	 */
	bool isFirstVisit() { return firstVisit; }

	/** Calling this method on the first visit of a node causes the traversal of nodes below the current node to be skipped, i.e., iteration continues as if this node has been visited for the second time (on the upward pass). Calling this on the second visit has no effect.
	 */
	void setUpPass() { firstVisit = false; }
private:
	bool firstVisit;
};

} // namespace stree


#endif // STREE_ITERATORS_H
