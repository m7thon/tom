#ifndef STREE_CORE_H
#define STREE_CORE_H

#include "stree.h"

namespace stree {
namespace internal {
    
/** A \c LeafNode consists of a left and right \c nidx_t. */
struct LeafNode {
    nidx_t l_; /**< the left \c nidx_t. */
    nidx_t r_; /**< the right \c nidx_t. */
};

/** An \c InternalNode consists of a left, right and child \c nidx_t. */
struct InternalNode {
    nidx_t l_; //< the left \c nidx_t.
    nidx_t r_; //< the right \c nidx_t.
    nidx_t c_; //< the child \c nidx_t.
};

SWIGCODE(%ignore Pos;)
/** A \c Pos refers to the location in the suffix tree that corresponds to some substring of the represented \c Sequence\. This position is unique, but can either be an explicit node (internal or leaf) or an implicit internal node\. This class is only for internal use in the suffix tree construction\. Please use the class \c stree::Position instead! */
class Pos {
public:
    /** Create a \c Pos corresponding to the root of the suffix tree */
    Pos() : node_(ROOT), hIndex_(0), depth_(0) { edgePtr_ = &node_; }
    
    nidx_t node_;     //< The last / deepest node on the path from the root to the represented position in the suffix tree
    nidx_t * edgePtr_; //< If the represented position lies on an edge (i.e., it is an implicit internal node), then this is  a pointer from the parent \c node to the next node, which defines this edge. Otherwise, this is a pointer from the parent to the current \c node. Note that this pointer is invalidated by any change to the data structures underlying the suffix tree (the \c nodes and \c leaves arrays).
    nidx_t hIndex_;    //< The index in the underlying `Sequence` corresponding to the represented substring
    nidx_t depth_;     //< The length of the represented substring = depth in the (uncompressed) suffix trie
    
    /** If the current `Pos` corresponds to the subsequence \c seq, then attempt to move to \c seq concatenated with \c chr\. Return true if this is successful, i.e., if \c seq + \c chr is also a substring of the represented \c sequence. */
    bool followSymbol(STree* stree, const Symbol chr);
    inline void preCanonize(STree* stree);
    inline void canonize(STree* stree);
    void followSuffixLink(STree* stree);
    /** Return \c true if the \c Pos corresponds to an explicit node. */
    bool isExplicit() const {return (*edgePtr_ & (INTERNAL | INDEX)) == (node_ & (INTERNAL | INDEX));}
    }; // class Pos

SWIGCODE(%ignore RBTreeNodeTraits;)
class RBTreeNodeTraits;
typedef RBTree<RBTreeNodeTraits> RBSiblingTree;

} //namespace internal

/**
 * An implementation of suffix trees for the `Sequence` type, which can represent a plain or io-sequence.
 *
 * Features
 * --------
 *
 * - Sequences are not required to end in a unique termination symbol:
 *   - Every suffix is guaranteed to end in a leaf *or* internal node. The suffixes ending in internal nodes are given by `internalSuffixes()`.
 *   - `.count()` always means the number of occurrences as a subsequence in the represented sequence. This is *not* the same as the leaf count.
 * - `extendTo(seq)` can be used to efficiently "extend" a suffix tree representation for a sequence `s` to a suffix tree for the sequence `seq` if `s` is a prefix of `seq`.
 * - Io-sequences are basically treated as raw sequences, but only every second suffix (that is io-aligned) is represented in the suffix tree. Therefore:
 *   - Occurrence counts are supported also for io-sequences of the form $u_0o_0, ..., u_k$ (not back aligned).
 *   - The suffix of $u_0o_0...u_ko_k$ is considered to be $u_1o_1...u_ko_k$.
 * - Sequences of size up to $2^29 - 1$ are supported.
 * - The space requirement is at most 24 bytes / symbol, i.e., less than 12 GB for the maximum supported sequence size.
 *
 * Details
 * -------
 * In the following comes a brief description of the internal structure of this suffix tree implementation:
 *
 * - The suffix tree has two types of nodes (internal and leaf nodes), which are stored in vectors \c nodes and \c leaves respectively.
 * - Nodes are addressed by a 32-bit \c nidx_t ("node index"). The first three bits of a \c nidx_t indicate whether the \c nidx_t is valid (1) or nil (0), addresses an (internal) node (1) or a leaf (0), and is "colored" (1) or not (0), respectively. In the case of a valid \c nidx_t, the remaining 29 bits form the index value of the addressed node in the \c nodes or \c leaves vector. Otherwise, these can have some other significance.
 * - Every node (internal or leaf) has a left and right \c nidx_t. Internal nodes have an additional child \c nidx_t, which addresses the first child in the suffix tree structure.
 * - The children (which are siblings) of any internal node are basically organized in a self-balancing binary tree structure formed by the left and right \c nidx_t "pointers". For this, a left-leaning red-black tree is used (which uses the color flag of the \c nidx_t). However, the first two children of any node have special roles:
 *   - The right \c nidx_t of the first child encodes the headindex of the parent node, while the right \c nidx_t of the second child encodes the depth of the parent node. Both right \c nidx_t addresses of the first two children have their first three bits set to zero.
 *   - The left \c nidx_t of the first child addresses the second child and has its color flag set. The left \c nidx_t of the second child addresses the root of the red-black tree in which all further siblings are organized, and has its color flag unset (to be able to distinguish first and second children). However, if there are only two siblings, then the left \c nidx_t of the second child will be marked as invalid, but will address the suffix-link of the parent. Still, the color flag will be unset.
 * - Suffix-links of any (internal) node are stored in the right \c nidx_t of the rightmost (in the binary tree) child, or alternatively in the left \c nidx_t of the second child if there are only two children. This \c nidx_t will always be marked as invalid and uncolored, but as addressing a node.
 * - The left \c nidx_t of the leftmost sibling in the binary tree will always be marked as invalid and uncolored, but has no further meaning.
 * - The binary tree of siblings is threaded. This means that invalid left and right \c nidx_t entries address previous and next siblings respectively. This is true for all siblings in the red-black tree, except for the left-and rightmost. To distinguish these, all invalid \c nidx_t that indicate a thread are marked as colored (as opposed to the left-and rightmost, which are always marked as uncolored).
 * - Siblings are stored in lexicographic order with respect to their edge labels. It is ensured that the first and second siblings are always the lexicographically smallest.
 *
 */
class STree {
    friend class internal::Pos;
    friend class internal::RBTreeNodeTraits;
    friend class Node;
    friend class EdgeNode;
    friend class PathNode;
    friend class Position;
public:
    /** Create a suffix tree for the given `sequence`. Note that the sequence must have size at least one.
     */
    STree(const Sequence& sequence) {
        sequence_ = sequence.rawSub(0,0);
        symbolSize_ = sequence_.isIO() ? 2 : 1;
        extendTo(sequence);
    }
    
    /** Extend the current suffix tree representation to a representation for a given `sequence`, which requires that the current suffix tree represents a prefix of the given `sequence`.
     */
    void extendTo(const Sequence& sequence, bool checkExtendability = true) CHECK(throw (std::invalid_argument));

    /** Return the represented sequence.
     */
    const Sequence sequence() const { return sequence_; }

    /** Return the number of leaf nodes in the suffix tree. */
    nidx_t nLeafNodes() const { return leaves_.size() - nTemporaryInternalNodes_; }
    
    /** Return the number of internal nodes in the suffix tree. */
    nidx_t nInternalNodes() const { return nodes_.size(); }
    
    /** Return the number of nodes (internal and leaves) in the suffix tree. */
    nidx_t nNodes() const { return leaves_.size() + nodes_.size(); }
    
    /** Return the node index (`nidx_t`) of the deepest internal node in the suffix tree representing a suffix of the represented sequence. This can and should be converted to a `Node` by calling `Node(stree, stree.deepestInternalSuffix())`. [This is required for technical reasons of memory management safety].
     *
     * If the input sequence terminates with a unique symbol, then this will always be the root of the suffix tree, corresponding to the empty suffix. If the input sequence does not terminate with a unique symbol, this corresponds to the "active position" in the suffix tree construction, which will be an internal node. The other internal nodes representing a suffix can be found by traversing the suffix link (calling `toSuffix()` on the converted returned `Node`) until reaching the root node (to exclude the empty suffix) or until `toSuffix()` results in an invalid `Node` (to include the empty suffix). The remaining (longer) suffixes correspond to the leaf nodes.
     */
    nidx_t deepestInternalSuffixNidx() const;

private:
    Sequence sequence_; /**< a copy of the underlying sequence for which the suffix tree is built. This must not be changed during the lifetime of the suffix tree. */
    unsigned int symbolSize_; /**< 1 if every suffix of the sequence is represented; 2 if only every second suffix is represented (so in a sense one symbol consists of two characters), and so o.n */
    
    std::vector<internal::InternalNode> nodes_;   //< the vector of internal nodes
    std::vector<internal::LeafNode> leaves_;      //< the vector of leaf nodes
    nidx_t nTemporaryInternalNodes_ = 0;       //< the number of temporary internal nodes created by \ref createTemporaryInternalNodes().
    std::vector<nidx_t> nOccurrences_;     //< for each internal node, the number of occurrences of the corresponding substring in the represented sequence.
    internal::Pos currentPos_;               //< the current position in the suffix tree construction
    nidx_t pos_ = 1;                           //< the position in the \c sequence corresponding to the current step in the construction process
    nidx_t suffixLinkFrom_ = 0;
    
    /** Return the number of occurrences of the substring corresponding to this \c node in the represented \c sequence.
     */
    nidx_t n(const nidx_t node) const { if (!(node & VALID)) return 0;
        if (node & INTERNAL) return nOccurrences_[node & INDEX]; else return 1;
    }

    /** Return `true` if the given `nidx` can refer to a valid node in the suffix tree. Note that the `valid` and `color` flags are ignored.
     */
    bool validate(nidx_t nidx) const {
        if (nidx & INTERNAL) {
            if ((nidx & INDEX) < nInternalNodes()) { return true; }
        } else {
            if ((nidx & INDEX) < nLeafNodes()) { return true; }
        }
        return false;
    }

#ifndef FOLDABLE_CODE_GROUP_HACK /** @name Node data manipulation */
    /**@{*/ //MARK: Node data manipulation

    /** Return the depth of the \c node, i.e., the size of the substring represented by the \c node. */
    nidx_t d(const nidx_t node) const { return (node & INTERNAL) ? r(l(c(node))) : sequence_.rawSize() - ((node & INDEX) * symbolSize_); }

    /** Set the depth of the given \c node to the given \c depth. */
    void d(const nidx_t node, const nidx_t depth) { r(l(c(node))) = depth; }
    
    /** Return the head index of the \c node, i.e., a position in the underlying \c sequence where the substring represented by the \c node can be found. */
    nidx_t hi(const nidx_t node) const { return (node & INTERNAL) ? r(c(node)) : (node & INDEX) * symbolSize_; }

    /** Set the head index of the given \c node to the given \c headIndex. */
    void hi(const nidx_t node, const nidx_t headIndex) { r(c(node)) = headIndex; }

    /** Return the left \c nidx_t\& of the given \c node. */
    const nidx_t & l(const nidx_t node) const { return (node & INTERNAL) ? nodes_[node & INDEX].l_ : leaves_[node & INDEX].l_; }

    /** Return the \b non-const left \c nidx_t\& of the given \c node. */
    nidx_t & l(const nidx_t node) { return (node & INTERNAL) ? nodes_[node & INDEX].l_ : leaves_[node & INDEX].l_; }

    /** Return the right \c nidx_t\& of the given \c node. */
    const nidx_t & r(const nidx_t node) const { return (node & INTERNAL) ? nodes_[node & INDEX].r_ : leaves_[node & INDEX].r_; }
    
    /** Return the \b non-const right \c nidx_t\& of the given \c node. */
    nidx_t & r(const nidx_t node) { return (node & INTERNAL) ? nodes_[node & INDEX].r_ : leaves_[node & INDEX].r_; }

    /** Return the child \c nidx_t\& of the given \c node\. Note that \c node must specify an internal node and not a leaf. */
    const nidx_t & c(const nidx_t node) const { return nodes_[node & INDEX].c_; }

    /** Return the \b non-const child \c nidx_t\& of the given \c node\. Note that \c node must specify an internal node and not a leaf. */
    nidx_t & c(const nidx_t node) { return nodes_[node & INDEX].c_; }

    /** Return the child \c nidx_t\& of the given \c node corresponding to the given \c chr (which is the first character of the edge leading away from the \c node)\. If no corresponding child is found, the (null) \c nidx_t\& will be returned that corresponds to the place where the according child node would need to be inserted. */
    const nidx_t & c(const nidx_t node, Symbol chr) const;
    
    /** Return the \b non-const child \c nidx_t\& of the given \c node corresponding to the given \c chr (which is the first character of the edge leading away from the \c node)\. If no corresponding child is found, the (null) \c nidx_t\& will be returned that corresponds to the place where the according child node would need to be inserted. */
    nidx_t & c(const nidx_t node, Symbol chr);

    /** Return the next sibling (according to lexicographic ordering of edge labels) of the given \c node, or a null \c nidx_t if no further sibling exists. */
    nidx_t sib(const nidx_t node) const;

    /** Return the suffix link \c nidx_t\& of the given internal `node`. */
    const nidx_t & sl(const nidx_t node) const { const nidx_t * x = &(l(l(c(node)))); while (*x & VALID) x = &(r(*x)); return *x; }
    
    /** Set the suffix link of the given internal \c node to the given \c suffixLink. */
    void sl(const nidx_t node, nidx_t suffixLink) { nidx_t * x = &(l(l(c(node)))); while (*x & VALID) x = &(r(*x)); *x = suffixLink & ~(VALID | COLOR); }

    //@}
#endif

    /** Add the node \c newChild as a new child according to the given \c chr (the first character of the edge leading from \c node to \c newChild) to the given \c node. */
    void addChild(const nidx_t node, nidx_t newChild, Symbol chr);
    
    /** Return the first character of the edge label leading from its parent to the \c node, where \c parentDepth specifies the depth of the parent node (this needs to be given since parent information is not stored in the suffix tree). */
    Symbol edgeKey(const nidx_t node, const nidx_t parentDepth) const { return sequence_.rawAt(hi(node) + parentDepth); }
    
    /** Create a new leaf at the current position according to \c currentPos in the suffix tree\. If the \c currentPos is not an internal node, then a new internal node is created also\. The \c swap parameter only plays a role in the last case, and means, when set to \c false, that the new leaf will always be the second child of the new internal node, regardless of the correct ordering (this is used for \c temporary \c internal \c nodes)\. Finally, the suffix-link leading to this (new) node is set\. Note that this function invalidates the \c edgePtr of the \c currentPos ! */
    void createNewLeaf(bool swap = true);
    
    /** Create internal nodes as if a terminal symbol was appended to the \c sequence at the current position\. These are needed to be able to count the number of occurrences of a substring by counting the number of leaves below the corresponding node. */
    void createTemporaryInternalNodes();

    /** Remove the temporary internal nodes created by `createTemporaryInternalNodes()`\. After calling this function, the suffix tree can be extended further. */
    void removeTemporaryInternalNodes();
    
    /** Annotate the suffix tree with leaf counts, i.e., for each node, count the number of leafs (= number of occurrences of the corresponding substring) and store this in \c nOccurrences. */
    void annotate();
    
}; // class STree



namespace internal {
    /** An object specifying the node traits for the underlying red-black tree implementation. */
    class RBTreeNodeTraits
    {
    public:
        typedef nidx_t RBNodePtr;
        typedef Symbol RBKey;
        RBTreeNodeTraits(STree* bst, nidx_t parentDepth = 0) : bst_(bst), parentDepth_(parentDepth) {}
        bool less(RBKey k, RBNodePtr n) const { return (k < bst_->edgeKey(n, parentDepth_)); }
        bool equals(RBKey k, RBNodePtr n) const { return (k == bst_->edgeKey(n, parentDepth_)); }
        bool isNull(RBNodePtr n) const { return !(n & VALID); }
        void setThread(RBNodePtr& n) const { n &= ~VALID; n |= COLOR;} // threaded nodes are colored red
        RBNodePtr& left(RBNodePtr n) const { return bst_->l(n); }
        RBNodePtr& right(RBNodePtr n) const { return bst_->r(n); }
        void set(RBNodePtr& n, RBNodePtr nNew) const { n = nNew; }
        bool getColor(RBNodePtr n) const { return (n & COLOR); }
        void setColor(RBNodePtr& n, bool c) const { if (c) n |= COLOR; else n &= ~COLOR; }
    private:
        STree* bst_;
        nidx_t parentDepth_;
    }; // class RBTreeNodeTraits
    
} //namespace internal

} // namespace stree

#endif // STREE_CORE_H
