// -*- C++ -*-

/**
 * @file   STreeCore.h
 * @brief  An implementation of suffix trees for generic string types.
 */

#ifndef STREE_CORE_H
#define STREE_CORE_H

namespace stree {

/**
 * An implementation of suffix trees for generic string types. For string types other than \a std::string, one must define \a stree::STREE_STRING_TYPE which must conform to a standard string interface. (Note that template arguments are avoided on purpose to insure better compatibility to a python interface). The children of an internal node are stored in a self-balancing binary tree, so that this implementation is suitable especially for strings over large alphabets. The space requirement of the index structure (including suffix links) is at most 5 * 4 bytes / input symbol, and the maximum string length is 2**29-1.
 *
 * In the following comes a brief description of the internal structure of this suffix tree implementation:
 * - The suffix tree has two types of nodes (internal and leaf nodes), which are stored in vectors \a nodes and \a leaves respectively.
 * - Nodes are addressed by a 32-bit \a Nidx ("node index"). The first three bits of a \a Nidx indicate whether the \a Nidx is valid (1) or nil (0), addresses an (internal) node (1) or a leaf (0), and is "colored" (1) or not (0), respectively. In the case of a valid \a Nidx, the remaining 29 bits form the index value of the addressed node in the \a nodes or \a leaves vector. Otherwise, these can have some other significance.
 * - Every node (internal or leaf) has a left and right \a Nidx. Internal nodes have an additional child \a Nidx, which addresses the first child in the suffix tree structure.
 * - The children of any node (which are siblings) are basically organized in a self-balancing binary tree structure formed by the left and right \a Nidx "pointers". For this, a left-leaning red-black tree is used (which uses the color flag of the \a Nidx). However, the first two children of any node have special roles:
 * - The right \a Nidx of the first child encodes the headindex of the parent node, while the right \a Nidx of the second child encodes the depth of the parent node. Both right \a Nidx addresses of the first two children have their first three bits set to zero.
 * - The left \a Nidx of the first child addresses the second child and has its color flag set. The left \a Nidx of the second child addresses the root of the red-black tree in which all further siblings are organized, and has its color flag unset (to be able to distinguish first and second children). However, if there are only two siblings, then the left \a Nidx of the second child will be marked as invalid, but will address the suffix-link of the parent. Still, the color flag will be unset.
 * - Suffix-links of any (internal) node are stored in the right \a Nidx of the rightmost (in the binary tree) child, or alternatively in the left \a Nidx of the second child if there are only two children. This \a Nidx will always be marked as invalid and uncolored, but as addressing a node.
 * - The left \a Nidx of the leftmost sibling in the binary tree will always be marked as invalid and uncolored, but has no further meaning.
 * - The binary tree of siblings is threaded. This means that invalid left and right \a Nidx entries address previous and next siblings respectively. This is true for all siblings in the red-black tree, except for the left-and rightmost. To distinguish these, all invalid \a Nidx that indicate a thread are marked as colored (as opposed to the left-and rightmost, which are always marked as uncolored).
 * - Siblings are stored in lexicographic order with respect to their edge labels. It is ensured that the first and second siblings are always the lexicographically smallest.
 *
 */
class STree {
	friend class Position;
public:
	/** Create an uninitialized \a STree object. */
	STree() {}

	/** Create a suffix tree for a given \a text.
	 *
	 * \param text The text that should be represented by the suffix tree
	 * \param size The size of the \a text that should be represented, or 0 if the entire text should be represented
	 * \param symbolSize If set to 1, every suffix of the \a text is represented; if set to 2, only every second suffix is represented (so in a sense one symbol consists of two characters), and so on
	 * \param annotated \c true it the nodes should be annotated by the occurrence counts of the corresponding substrings in the \a text
	 */
  STree(const String& text, Idx size = 0, unsigned int symbolSize = 1, bool annotated = true) {
		initialize(text, size, symbolSize, annotated);
	};

	/** (Re-)Initialize the \a STree object.
	 *
	 * \param text The text that should be represented by the suffix tree
	 * \param size The size of the \a text that should be represented, or 0 if the entire text should be represented
	 * \param symbolSize If set to 1, every suffix of the \a text is represented; if set to 2, only every second suffix is represented (so in a sense one symbol consists of two characters), and so on
	 * \param annotated \c true it the nodes should be annotated by the occurrence counts of the corresponding substrings in the \a text
	 */
  void initialize(const String& text, Idx size = 0, unsigned int symbolSize = 1, bool annotated = true);

	/** Extend the suffix tree representation of the underlying \a text to the given \a size\. Note that the given \a size must be smaller or equal to the size of the underlying \a text\. If the string does not end with a unique terminal symbol, temporary internal nodes (that have only one child!) will be inserted into the edges for positions corresponding to suffixes of the input string. */
  void extendTo(Idx size);

  /** @name High-level suffix tree interface functions */ //@{

  /** Return the number of leaf nodes in the suffix tree. */
	Idx nLeaves() const { return leaves.size(); }

	/** Return the number of internal nodes in the suffix tree. */
	Idx nInternalNodes() const { return nodes.size(); }

  /** Return the number of nodes (internal and leaves) in the suffix tree. */
	Idx nNodes() const { return leaves.size() + nodes.size(); }
  //@}

  /** return the deepest node in the SuffixTree having a virtual / temporary leaf attached\. This occurs when the input text does not terminate with a unique symbol, and corresponds to the "active position" in the suffix tree construction\. The other virtual leaf branches can be found by traversing the suffix links until reaching the root node. */
  STreeNode getDeepestVirtualLeafBranch();

  /* Return the number of leaves below this \a node, which is the same as the number of occurrences of the substring corresponding to this \a node in the \a text. */
	Idx n(const Nidx node) const { if (!((annotated_) and (node & VALID))) return 0;
		if (node & NODE) return nOccurrences[node & INDEX]; else return 1;
	}

	/** Return the depth of the \a node, i.e., the size of the substring represented by the \a node. */
  Idx d(const Nidx node) const { return (node & NODE) ? r(l(c(node))) : size_ - ((node & INDEX) * symbolSize_); }
	/** Return the head index of the \a node, i.e., a position in the underlying \a text where the substring represented by the \a node can be found. */
  Idx hi(const Nidx node) const { return (node & NODE) ? r(c(node)) : (node & INDEX) * symbolSize_; }
	/** Return the left \a Nidx\& of the given \a node. */
  const Nidx& l(const Nidx node) const { return (node & NODE) ? nodes[node & INDEX].l : leaves[node & INDEX].l; }
	/** Return the right \a Nidx\& of the given \a node. */
  const Nidx& r(const Nidx node) const { return (node & NODE) ? nodes[node & INDEX].r : leaves[node & INDEX].r; }
	/** Return the child \a Nidx\& of the given \a node\. Note that \a node must specify an internal node and not a leaf. */
  const Nidx& c(const Nidx node) const { return nodes[node & INDEX].c; }
	/** Return the child \a Nidx\& of the given \a node corresponding to the given \a chr (which is the first character of the edge leading away from the \a node)\. If no corresponding child is found, the (null) \a Nidx\& will be returned that corresponds to the place where the according child node would need to be inserted. */
	const Nidx& c(const Nidx node, Char chr) const;
	/** Return the suffix link \a Nidx\& of the given \a node. */
	const Nidx& sl(const Nidx node) const { const Nidx* x = &(l(l(c(node)))); while (*x & VALID) x = &(r(*x)); return *x; }
	/** Return the next sibling (according to lexicographic ordering of edge labels) of the given \a node, or a null \a Nidx if no further sibling exists. */
  Nidx sib(const Nidx node) const;

	/** Return the character at index \a pos in the string represented by this suffix tree. */
	Char at(Idx pos) const { return text_.rawAt(pos); }

  String text_;            //< a copy of the underlying text for which the suffix tree is built. This must not be changed during the lifetime of the suffix tree.
  Idx size_;                         //< the size of the represented text
  unsigned int symbolSize_;          //< 1 if every suffix of the text is represented; 2 if only every second suffix is represented (so in a sense one symbol consists of two characters), and so on
	bool annotated_;                   //< true it this suffix tree is annotated (i.e., the number of occurrences of substrings are computed and stored in \a nOccurrences)

private:
	/** Return the \b non-const left \a Nidx\& of the given \a node. */
  Nidx& l(const Nidx node) { return (node & NODE) ? nodes[node & INDEX].l : leaves[node & INDEX].l; }
	/** Return the \b non-const right \a Nidx\& of the given \a node. */
  Nidx& r(const Nidx node) { return (node & NODE) ? nodes[node & INDEX].r : leaves[node & INDEX].r; }
	/** Return the \b non-const child \a Nidx\& of the given \a node\. Note that \a node must specify an internal node and not a leaf. */
  Nidx& c(const Nidx node) { return nodes[node & INDEX].c; }
	/** Return the \b non-const child \a Nidx\& of the given \a node corresponding to the given \a chr (which is the first character of the edge leading away from the \a node)\. If no corresponding child is found, the (null) \a Nidx\& will be returned that corresponds to the place where the according child node would need to be inserted. */
	Nidx& c(const Nidx node, Char chr);
  /** Set the depth of the given \a node to the given \a depth. */
  void d(const Nidx node, const Idx depth) { r(l(c(node))) = depth; }
  /** Set the head index of the given \a node to the given \a headIndex. */
  void hi(const Nidx node, const Idx headIndex) { r(c(node)) = headIndex; }
  /** Set the suffix link of the given \a node to the given \a suffixLink. */
  void sl(const Nidx node, Nidx suffixLink) { Nidx* x = &(l(l(c(node)))); while (*x & VALID) x = &(r(*x)); *x = suffixLink & ~(VALID | COLOR); }

  /** Add the node \a newChild as a new child according to the given \a chr (the first character of the edge leading from \a node to \a newChild) to the given \a node. */
  void addChild(const Nidx node, Nidx newChild, Char chr);

	/** Return the first character of the edge label leading from its parent to the \a node, where \a parentDepth specifies the depth of the parent node (this needs to be given since parent information is not stored in the suffix tree). */
	Char key(const Nidx node, const Idx parentDepth) const { return at(hi(node) + parentDepth); }

  /** Create a new leaf at the current position according to \a currentPos in the suffix tree\. If the \a currentPos is not an internal node, then a new internal node is created also\. The \a swap parameter only plays a role in the last case, and means, when set to \c false, that the new leaf will always be the second child of the new internal node, regardless of the correct ordering (this is used for \a temporary \a internal \a nodes)\. Finally, the suffix-link leading to this (new) node is set\. Note that this function invalidates the \a edgePtr of the \a currentPos ! */
  void createNewLeaf(bool swap = true);

  /** Create internal nodes as if a terminal symbol was appended to the \a text at the current position\. These are needed to be able to count the number of occurrences of a substring by counting the number of leaves below the corresponding node. */
  void createTemporaryInternalNodes();
	/** Remove the temporary internal nodes created by \ref createTemporaryInternalNodes()\. After calling this function, the suffix tree can be extended further. */
  void removeTemporaryInternalNodes();

  /** Annotate the suffix tree with leaf counts, i.e., for each node, count the number of leafs (= number of occurrences of the corresponding substring) and store this in \a nOccurrences. */
	void annotate();

private:
	/** A \a LeafNode consists of a left and right \a Nidx. */
  class LeafNode {
	public:
		/** Create a \a LeafNode with zero (hence null) left and right \a Nidx. */
    LeafNode() : l(0), r(0) {}
		/** Create a \a LeafNode with the given left (\a l_) and right (\a r_) \a Nidx. */
    LeafNode(Nidx l_, Nidx r_) : l(l_), r(r_) {}
    Nidx l; //< the left \a Nidx.
		Nidx r; //< the right \a Nidx.
  };

	/** An \a InternalNode consists of a left, right and child \a Nidx. */
  class InternalNode
  { public:
		/** Create a \a LeafNode with zero (hence null) left, right and child \a Nidx. */
    InternalNode() : l(0), r(0), c(0) {}
		/** Create a \a LeafNode with the given left (\a l_), right (\a r_) and child (\a c_) \a Nidx. */
    InternalNode(Nidx l_, Nidx r_, Nidx c_) : l(l_), r(r_), c(c_) {}
    Nidx l; //< the left \a Nidx.
		Nidx r; //< the right \a Nidx.
    Nidx c; //< the child \a Nidx.
  };

	/** An object specifying the node traits for the underlying red-black tree implementation. */
  class RBTreeNodeTraits
  {
  public:
    typedef Nidx RBNodePtr;
    typedef Char RBKey;
    RBTreeNodeTraits(STree* bst_, Idx parentDepth_ = 0) : bst(bst_), parentDepth(parentDepth_) {}
    bool less(RBKey k, RBNodePtr n) const { return (k < bst->key(n, parentDepth)); }
    bool equals(RBKey k, RBNodePtr n) const { return (k == bst->key(n, parentDepth)); }
    bool isNull(RBNodePtr n) const { return !(n & VALID); }
    void setThread(RBNodePtr& n) const { n &= ~VALID; n |= COLOR;} // threaded nodes are colored red
    RBNodePtr& left(RBNodePtr n) const { return bst->l(n); }
    RBNodePtr& right(RBNodePtr n) const { return bst->r(n); }
    void set(RBNodePtr& n, RBNodePtr nNew) const { n = nNew; }
    bool getColor(RBNodePtr n) const { return (n & COLOR); }
    void setColor(RBNodePtr& n, bool c) const { if (c) n |= COLOR; else n &= ~COLOR; }
  private:
    STree* bst;
    Idx parentDepth;
  }; // class RBTreeNodeTraits

  typedef RBTree<RBTreeNodeTraits> rbt;

  /** A \a Position refers to the location in the suffix tree that corresponds to some substring of the represented \a String\. This position is unique, but can either be an expicit node (internal or leaf) or an implicit internal node\. This class is only for internal use in the suffix tree construction\. Please use the class \a STreePos instead! */
	class Position {
	public:
		/** Create a \a Position corresponding to the root of the suffix tree */
		Position() : node(ROOT), hIndex(0), depth(0) { edgePtr = &node; }

		Nidx node;     //< The last / deepest node on the path from the root to the represented position in the suffix tree
		Nidx* edgePtr; //< If the represented position lies on an edge (i.e., it is an implicit internal node), then this is  a pointer from the parent \a node to the next node, which defines this edge. Otherwise, this is a pointer from the parent to the current \a node. Note that this pointer is invalidated by any change to the data structures underlying the suffix tree (the \a nodes and \a leaves arrays).
		Idx hIndex;    //< The index in the underlying \a String coresponding to the represented substring
		Idx depth;     //< The length of the represented substring = depth in the (uncompressed) suffix trie

		/** If the current STreePosition corresponds to the substring \a str, then attempt to move to \a str concatenated with \a chr\. Return true if this is successful, i.e., if \a str + \a chr is also a substring of the represented \a text. */
		bool followChar(STree* stree, const Char chr);
		inline void preCanonize(STree* stree);
		inline void canonize(STree* stree);
		void followSuffixLink(STree* stree);
		/** Return \c true if the \a STreePosition corresponds to an explicit node. */
		bool isExplicit() const {return (*edgePtr & (NODE | INDEX)) == (node & (NODE | INDEX));}
	}; // class Position

  std::vector<InternalNode> nodes;   //< the vector of internal nodes
  std::vector<LeafNode> leaves;      //< the vector of leaf nodes
  Idx nTemporaryInternalNodes;       //< the number of temporary internal nodes created by \ref createTemporaryInternalNodes().
	std::vector<Idx> nOccurrences;     //< for each internal node, the number of occurrences of the corresponding substring in the represented text.

	Position currentPos;               //< the current position in the suffix tree construction
  Idx pos;                           //< the position in the \a text corresponding to the current step in the construction process
  Nidx suffixLinkFrom;

public:
	std::string serialize() {
		std::stringstream os;
		os << "Nodes:" << std::endl;
		for (int i = 0; i < nodes.size(); ++i)
			os << nodes[i].l << " " << nodes[i].r << " " << nodes[i].c << " " << nOccurrences[i] << std::endl;
		os << "Leaves:" << std::endl;
		for (int i = 0; i < leaves.size(); ++i)
			os << leaves[i].l << " " << leaves[i].r << std::endl;
		os << "nTemporaryInternalNodes: " << nTemporaryInternalNodes << std::endl;
		return os.str();
	}
}; // class STree

} // namespace stree

#endif // STREE_CORE_H
