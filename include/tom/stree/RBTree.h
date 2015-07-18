#ifndef RBTREE_H
#define RBTREE_H

#include "stree.h"

namespace stree {

/**
 * An example of which functions an \c RBTreeNodeTraits object must provide to be used with \c RBTree. It is assumed that each node has a left and right pointer to child nodes, forming a binary tree. These pointers may be any user defined data structure used to index the nodes. For this reason, all needed functions to manipulate the tree structure must me provided in an \c RBTreeNodeTraits object. Note that nodes are referred to only by \c RBNodePtr objects and never directly.
 */
template< typename NodeType, typename KeyType >
class RBTreeNodeTraitsTemplate {
	/** the type of pointer to a node: These may be simple pointers, but may also be more sophisticated structures used to index nodes. */
  typedef NodeType* RBNodePtr;
	/** the key type used for comparison operations. */
  typedef KeyType RBKey;
	/**
	 * \c true if the given \c key is less than the key of the node corresponding to \c n
	 */
  bool less(RBKey& k, RBNodePtr& n);
	/**
	 * \c true if the given \c key is equal to the key of the node corresponding to \c n
	 */
  bool equals(RBKey& k, RBNodePtr& n);
	/**
	 * \c true if \c n is a \c NULL pointer
	 */
  bool isNull(RBNodePtr& n);
  /**
	 * the left \c RBNodePtr of the node corresponding to \c n
	 */
  RBNodePtr& left(RBNodePtr& n);
  /**
	 * the right \c RBNodePtr of the node corresponding to \c n
	 */
  RBNodePtr& right(RBNodePtr& n);
	/**
	 * set the \c RBNodePtr \c n to the given \c RBNodePtr \c nNew. Note that if color, nullness and threading properties are stored in the \c RBNodePtr structure, then these must be preserved by the set operation. This will generaly just amount to n = nNew.
	 */
  void set(RBNodePtr& n, RBNodePtr& nNew);
	/**
	 * mark the \c RBNodePtr \c n as \c NULL and indicate a threading if possible.
	 */
  void setThread(RBNodePtr& n);
	/**
	 * set the color of the node corresponding to \c n to the given color \c c.
	 */
  void setColor(RBNodePtr& n, bool c);
	/**
	 * return the color of the node corresponding to \c n.
	 */
	bool getColor(RBNodePtr& n) const;
};


/**
 * An implementation of left-leaning red-black trees following "Left-Leaning Red-Black Trees" by Robert Sedgewick from 2008. The functions can be used with generic kinds of nodes that need to be specified by providing an \c RBTreeNodeTraits object when calling the functions.
 *
 * Note that the red-black tree will be threaded if a \c NULL \c RBNodePtr can still address a node. This allows iterating over the stored values in their order merely by following the left or right pointers (even if \c NULL). The left-and rightmost \c RBNodePtr will be inherited from the left and right \c RBNodePtr of the original root, i.e., the first node \c RBNodePtr used in the first insertion operation when constructing the red-black tree, while every threaded right or left \c RBNodePtr will be set to the according \c RBNodePtr and then marked as a thread by calling the function \c setThread. It is up to the user to implement a suitable \c RBNodePtr structure to be able to distinguish the left-/rightmost \c RBNodePtr from a threaded \c RBNodePtr (e.g., if the color of nodes is stored in the \c RBNodePtr, even if \c NULL, then a threaded \c RBNodePtr may be \c NULL and colored red, while the left-/rightmost \c RBNodePtr may be \c NULL and black).
 */
template< typename RBTreeNodeTraits >
class RBTree {
  typedef RBTreeNodeTraits RBNT;
  typedef typename RBNT::RBKey RBKey;
  typedef typename RBNT::RBNodePtr RBNodePtr;
public:
  /**
   * search for a node matching a given \c key in the binary search tree below a given node \c h.
   *
   * @param h the node below which (and including) to search for the key
   * @param key the key value to be searched for
   * @param rbnt an \c RBTreeNodeTraits object
   *
   * @return the \c RBNodePtr& of the parent node to the node found, or some \c NULL \c RBNodePtr, if no matching node found.
   */
  static RBNodePtr&
  find(RBNodePtr& h, const RBKey& key, const RBNT& rbnt) {
    RBNodePtr* hPtr = &h;
    while (!rbnt.isNull(*hPtr)) {
      if (rbnt.less(key, *hPtr)) hPtr = &(rbnt.left(*hPtr));
			else if (rbnt.equals(key, *hPtr)) return *hPtr;
      else hPtr = &(rbnt.right(*hPtr));
    }
    return *hPtr;
  }

  /**
   * insert a node into the red-black tree according to a given key. Please see the general remarks about threading.
   *
   * @param h the \c RBNodePtr to the root of the red-black tree (this will point to the new root after the insertion operation)
   * @param n the node to be inserted
   * @param key the key value of the new node
   * @param rbnt an \c RBTreeNodeTraits object
   */
  static void
  insert(RBNodePtr& h, RBNodePtr n, const RBKey& key, const RBNT& rbnt) {
		if (rbnt.isNull(h)) { // initial construction of the red-black tree
      // you may want to insure that the leftmost and rightmost (null)-pointers in the tree can be recognized.
			rbnt.set(h, n);
			rbnt.setColor(h, false);
			return;
		}

    RBNodePtr* hPtr = &h;
    std::stack<RBNodePtr*> stack;

    // The pass down the tree...
    while (true) {
      stack.push(hPtr);
      if (!rbnt.isNull(rbnt.left(*hPtr)) && !rbnt.isNull(rbnt.right(*hPtr)) &&
					rbnt.getColor(rbnt.left(*hPtr)) && rbnt.getColor(rbnt.right(*hPtr)) )
				split(*hPtr, rbnt);
      if (rbnt.less(key, *hPtr))
				if (rbnt.isNull(rbnt.left(*hPtr))) { addNodeLeft(*hPtr, n, rbnt); break; }
				else hPtr = &(rbnt.left(*hPtr));
      else
				if (rbnt.isNull(rbnt.right(*hPtr))) { addNodeRight(*hPtr, n, rbnt); break; }
				else hPtr = &(rbnt.right(*hPtr));
    }

    // The pass up the tree...
    while (!stack.empty()) {
      hPtr = stack.top();
      stack.pop();
      if ( !rbnt.isNull(rbnt.right(*hPtr)) && rbnt.getColor(rbnt.right(*hPtr)) &&
	   (rbnt.isNull(rbnt.left(*hPtr)) || !rbnt.getColor(rbnt.left(*hPtr))) )
	rotateLeft(*hPtr, rbnt);
      if ( !rbnt.isNull(rbnt.left(*hPtr)) && rbnt.getColor(rbnt.left(*hPtr)) &&
	   !rbnt.isNull(rbnt.left(rbnt.left(*hPtr))) && rbnt.getColor(rbnt.left(rbnt.left(*hPtr))) )
	rotateRight(*hPtr, rbnt);
    }

    rbnt.setColor(h, false);
  }


  /**
   * fixes the threads leading to the node \c n in the red-black tree; this function needs to be called after replacing a node in the red-black tree structure.
   *
   * @param n the node whose threading needs to be fixed
   * @param rbnt an \c RBTreeNodeTraits object
   */
	static void
	fixThreading(const RBNodePtr& n, const RBNT& rbnt) {
		RBNodePtr* xPtr = &(rbnt.left(n));
		if (!rbnt.isNull(*xPtr)) {
			while (!rbnt.isNull(*xPtr)) xPtr = &(rbnt.right(*xPtr));
			rbnt.set(*xPtr, n);
			rbnt.setThread(*xPtr);
		}
		xPtr = &(rbnt.right(n));
		if (!rbnt.isNull(*xPtr)) {
			while (!rbnt.isNull(*xPtr)) xPtr = &(rbnt.left(*xPtr));
			rbnt.set(*xPtr, n);
			rbnt.setThread(*xPtr);
		}
	}

private:
	static void
	addNodeLeft(RBNodePtr& h, RBNodePtr& n, const RBNT& rbnt) {
		rbnt.set(rbnt.left(n), rbnt.left(h)); // inherit left-threading
		rbnt.set(rbnt.right(n), h); // setup right-threading
		rbnt.setThread(rbnt.right(n));
		rbnt.setColor(n, true); // new nodes are initially red
		rbnt.set(rbnt.left(h), n); // attach new node
	}
	static void
	addNodeRight(RBNodePtr& h, RBNodePtr& n, const RBNT& rbnt) {
		rbnt.set(rbnt.right(n), rbnt.right(h)); // inherit right-threading
		rbnt.set(rbnt.left(n), h); // setup left-threading
		rbnt.setThread(rbnt.left(n));
		rbnt.setColor(n, true); // new nodes are initially red
		rbnt.set(rbnt.right(h), n); // attach new node
	}
  static void
  rotateLeft(RBNodePtr& h, const RBNT& rbnt) {
    RBNodePtr x;
    rbnt.set(x, rbnt.right(h));
		if (rbnt.isNull(rbnt.left(x))) { // we need to take care of the threading
			rbnt.setThread(rbnt.right(h));
		}
		else { // normal case
			rbnt.set(rbnt.right(h), rbnt.left(x));
		}
    rbnt.setColor(x, rbnt.getColor(h));
    rbnt.setColor(h, true);
    rbnt.set(rbnt.left(x), h);
    rbnt.set(h, x);
  }

  static void
  rotateRight(RBNodePtr& h, const RBNT& rbnt) {
    RBNodePtr x;
    rbnt.set(x, rbnt.left(h));
		if (rbnt.isNull(rbnt.right(x))) { // we need to take care of the threading
			rbnt.setThread(rbnt.left(h));
		}
		else { // normal case
			rbnt.set(rbnt.left(h), rbnt.right(x));
		}
    rbnt.setColor(x, rbnt.getColor(h));
    rbnt.setColor(h, true);
    rbnt.set(rbnt.right(x), h);
    rbnt.set(h, x);
  }

  static void
  split(RBNodePtr& h, const RBNT& rbnt) {
    rbnt.setColor(h, true);
    rbnt.setColor(rbnt.left(h), false);
    rbnt.setColor(rbnt.right(h), false);
  }

}; // class RBTree

} // namespace stree

#endif // RBTREE_H
