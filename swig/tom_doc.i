
// File: index.xml

// File: class_linear_membership_function.xml
%feature("docstring") LinearMembershipFunction "

This class implements membership functions that are either piecewise
linear of simply a single dirac delta function.

C++ includes: MembershipFunctions.h ";

/*  IO-functions  */

%feature("docstring")  LinearMembershipFunction::from_string "

read the parameters from the given string.

The format must correspond to what the output functions produce. ";

%feature("docstring")  LinearMembershipFunction::to_string "

output the parameters as a string. ";

%feature("docstring")  LinearMembershipFunction::show "

write the continuous OOM parameters to std::cout ";

%feature("docstring")
LinearMembershipFunction::LinearMembershipFunction "

constructs an uninitialized (!) LinearMembershipFunction object ";

%feature("docstring")
LinearMembershipFunction::LinearMembershipFunction "

constructs a dirac delta membership function at x ";

%feature("docstring")  LinearMembershipFunction::setPoints "

setup a LinearMembershipFunction from the given points. ";

%feature("docstring")  LinearMembershipFunction::init "

initialize the LinearMembershipFunction from given nP and points. ";

%feature("docstring")  LinearMembershipFunction::value "

evaluates this function at x ";

%feature("docstring")  LinearMembershipFunction::mean "

the mean of this membership function viewed as a distribution ";

%feature("docstring")  LinearMembershipFunction::sample "

sample from this membership function viewed as a distribution ";


// File: classstd_1_1allocator.xml
%feature("docstring") std::allocator "

STL class. ";


// File: classstd_1_1auto__ptr.xml
%feature("docstring") std::auto_ptr "

STL class. ";


// File: classstd_1_1bad__alloc.xml
%feature("docstring") std::bad_alloc "

STL class. ";


// File: classstd_1_1bad__cast.xml
%feature("docstring") std::bad_cast "

STL class. ";


// File: classstd_1_1bad__exception.xml
%feature("docstring") std::bad_exception "

STL class. ";


// File: classstd_1_1bad__typeid.xml
%feature("docstring") std::bad_typeid "

STL class. ";


// File: classstd_1_1basic__fstream.xml
%feature("docstring") std::basic_fstream "

STL class. ";


// File: classstd_1_1basic__ifstream.xml
%feature("docstring") std::basic_ifstream "

STL class. ";


// File: classstd_1_1basic__ios.xml
%feature("docstring") std::basic_ios "

STL class. ";


// File: classstd_1_1basic__iostream.xml
%feature("docstring") std::basic_iostream "

STL class. ";


// File: classstd_1_1basic__istream.xml
%feature("docstring") std::basic_istream "

STL class. ";


// File: classstd_1_1basic__istringstream.xml
%feature("docstring") std::basic_istringstream "

STL class. ";


// File: classstd_1_1basic__ofstream.xml
%feature("docstring") std::basic_ofstream "

STL class. ";


// File: classstd_1_1basic__ostream.xml
%feature("docstring") std::basic_ostream "

STL class. ";


// File: classstd_1_1basic__ostringstream.xml
%feature("docstring") std::basic_ostringstream "

STL class. ";


// File: classstd_1_1basic__string.xml
%feature("docstring") std::basic_string "

STL class. ";


// File: classstd_1_1basic__string_1_1const__iterator.xml
%feature("docstring") std::basic_string::const_iterator "

STL iterator class. ";


// File: classstd_1_1basic__string_1_1const__reverse__iterator.xml
%feature("docstring") std::basic_string::const_reverse_iterator "

STL iterator class. ";


// File: classstd_1_1basic__string_1_1iterator.xml
%feature("docstring") std::basic_string::iterator "

STL iterator class. ";


// File: classstd_1_1basic__string_1_1reverse__iterator.xml
%feature("docstring") std::basic_string::reverse_iterator "

STL iterator class. ";


// File: classstd_1_1basic__stringstream.xml
%feature("docstring") std::basic_stringstream "

STL class. ";


// File: classstd_1_1bitset.xml
%feature("docstring") std::bitset "

STL class. ";


// File: classstd_1_1complex.xml
%feature("docstring") std::complex "

STL class. ";


// File: classstd_1_1deque.xml
%feature("docstring") std::deque "

STL class. ";


// File: classstd_1_1deque_1_1const__iterator.xml
%feature("docstring") std::deque::const_iterator "

STL iterator class. ";


// File: classstd_1_1deque_1_1const__reverse__iterator.xml
%feature("docstring") std::deque::const_reverse_iterator "

STL iterator class. ";


// File: classstd_1_1deque_1_1iterator.xml
%feature("docstring") std::deque::iterator "

STL iterator class. ";


// File: classstd_1_1deque_1_1reverse__iterator.xml
%feature("docstring") std::deque::reverse_iterator "

STL iterator class. ";


// File: classstd_1_1domain__error.xml
%feature("docstring") std::domain_error "

STL class. ";


// File: classstd_1_1exception.xml
%feature("docstring") std::exception "

STL class. ";


// File: classstd_1_1fstream.xml
%feature("docstring") std::fstream "

STL class. ";


// File: classstd_1_1ifstream.xml
%feature("docstring") std::ifstream "

STL class. ";


// File: classstd_1_1invalid__argument.xml
%feature("docstring") std::invalid_argument "

STL class. ";


// File: classstd_1_1ios.xml
%feature("docstring") std::ios "

STL class. ";


// File: classstd_1_1ios__base.xml
%feature("docstring") std::ios_base "

STL class. ";


// File: classstd_1_1ios__base_1_1failure.xml
%feature("docstring") std::ios_base::failure "

STL class. ";


// File: classstd_1_1istream.xml
%feature("docstring") std::istream "

STL class. ";


// File: classstd_1_1istringstream.xml
%feature("docstring") std::istringstream "

STL class. ";


// File: classstd_1_1length__error.xml
%feature("docstring") std::length_error "

STL class. ";


// File: classstd_1_1list.xml
%feature("docstring") std::list "

STL class. ";


// File: classstd_1_1list_1_1const__iterator.xml
%feature("docstring") std::list::const_iterator "

STL iterator class. ";


// File: classstd_1_1list_1_1const__reverse__iterator.xml
%feature("docstring") std::list::const_reverse_iterator "

STL iterator class. ";


// File: classstd_1_1list_1_1iterator.xml
%feature("docstring") std::list::iterator "

STL iterator class. ";


// File: classstd_1_1list_1_1reverse__iterator.xml
%feature("docstring") std::list::reverse_iterator "

STL iterator class. ";


// File: classstd_1_1logic__error.xml
%feature("docstring") std::logic_error "

STL class. ";


// File: classstd_1_1map.xml
%feature("docstring") std::map "

STL class. ";


// File: classstd_1_1map_1_1const__iterator.xml
%feature("docstring") std::map::const_iterator "

STL iterator class. ";


// File: classstd_1_1map_1_1const__reverse__iterator.xml
%feature("docstring") std::map::const_reverse_iterator "

STL iterator class. ";


// File: classstd_1_1map_1_1iterator.xml
%feature("docstring") std::map::iterator "

STL iterator class. ";


// File: classstd_1_1map_1_1reverse__iterator.xml
%feature("docstring") std::map::reverse_iterator "

STL iterator class. ";


// File: classstd_1_1multimap.xml
%feature("docstring") std::multimap "

STL class. ";


// File: classstd_1_1multimap_1_1const__iterator.xml
%feature("docstring") std::multimap::const_iterator "

STL iterator class. ";


// File: classstd_1_1multimap_1_1const__reverse__iterator.xml
%feature("docstring") std::multimap::const_reverse_iterator "

STL iterator class. ";


// File: classstd_1_1multimap_1_1iterator.xml
%feature("docstring") std::multimap::iterator "

STL iterator class. ";


// File: classstd_1_1multimap_1_1reverse__iterator.xml
%feature("docstring") std::multimap::reverse_iterator "

STL iterator class. ";


// File: classstd_1_1multiset.xml
%feature("docstring") std::multiset "

STL class. ";


// File: classstd_1_1multiset_1_1const__iterator.xml
%feature("docstring") std::multiset::const_iterator "

STL iterator class. ";


// File: classstd_1_1multiset_1_1const__reverse__iterator.xml
%feature("docstring") std::multiset::const_reverse_iterator "

STL iterator class. ";


// File: classstd_1_1multiset_1_1iterator.xml
%feature("docstring") std::multiset::iterator "

STL iterator class. ";


// File: classstd_1_1multiset_1_1reverse__iterator.xml
%feature("docstring") std::multiset::reverse_iterator "

STL iterator class. ";


// File: classstd_1_1ofstream.xml
%feature("docstring") std::ofstream "

STL class. ";


// File: classstd_1_1ostream.xml
%feature("docstring") std::ostream "

STL class. ";


// File: classstd_1_1ostringstream.xml
%feature("docstring") std::ostringstream "

STL class. ";


// File: classstd_1_1out__of__range.xml
%feature("docstring") std::out_of_range "

STL class. ";


// File: classstd_1_1overflow__error.xml
%feature("docstring") std::overflow_error "

STL class. ";


// File: classstd_1_1priority__queue.xml
%feature("docstring") std::priority_queue "

STL class. ";


// File: classstd_1_1queue.xml
%feature("docstring") std::queue "

STL class. ";


// File: classstd_1_1range__error.xml
%feature("docstring") std::range_error "

STL class. ";


// File: classstd_1_1runtime__error.xml
%feature("docstring") std::runtime_error "

STL class. ";


// File: classstd_1_1set.xml
%feature("docstring") std::set "

STL class. ";


// File: classstd_1_1set_1_1const__iterator.xml
%feature("docstring") std::set::const_iterator "

STL iterator class. ";


// File: classstd_1_1set_1_1const__reverse__iterator.xml
%feature("docstring") std::set::const_reverse_iterator "

STL iterator class. ";


// File: classstd_1_1set_1_1iterator.xml
%feature("docstring") std::set::iterator "

STL iterator class. ";


// File: classstd_1_1set_1_1reverse__iterator.xml
%feature("docstring") std::set::reverse_iterator "

STL iterator class. ";


// File: classstd_1_1smart__ptr.xml
%feature("docstring") std::smart_ptr "

STL class. ";


// File: classstd_1_1stack.xml
%feature("docstring") std::stack "

STL class. ";


// File: classstd_1_1string.xml
%feature("docstring") std::string "

STL class. ";


// File: classstd_1_1string_1_1const__iterator.xml
%feature("docstring") std::string::const_iterator "

STL iterator class. ";


// File: classstd_1_1string_1_1const__reverse__iterator.xml
%feature("docstring") std::string::const_reverse_iterator "

STL iterator class. ";


// File: classstd_1_1string_1_1iterator.xml
%feature("docstring") std::string::iterator "

STL iterator class. ";


// File: classstd_1_1string_1_1reverse__iterator.xml
%feature("docstring") std::string::reverse_iterator "

STL iterator class. ";


// File: classstd_1_1stringstream.xml
%feature("docstring") std::stringstream "

STL class. ";


// File: classstd_1_1underflow__error.xml
%feature("docstring") std::underflow_error "

STL class. ";


// File: classstd_1_1unique__ptr.xml
%feature("docstring") std::unique_ptr "

STL class. ";


// File: classstd_1_1valarray.xml
%feature("docstring") std::valarray "

STL class. ";


// File: classstd_1_1vector.xml
%feature("docstring") std::vector "

STL class. ";


// File: classstd_1_1vector_1_1const__iterator.xml
%feature("docstring") std::vector::const_iterator "

STL iterator class. ";


// File: classstd_1_1vector_1_1const__reverse__iterator.xml
%feature("docstring") std::vector::const_reverse_iterator "

STL iterator class. ";


// File: classstd_1_1vector_1_1iterator.xml
%feature("docstring") std::vector::iterator "

STL iterator class. ";


// File: classstd_1_1vector_1_1reverse__iterator.xml
%feature("docstring") std::vector::reverse_iterator "

STL iterator class. ";


// File: classstd_1_1weak__ptr.xml
%feature("docstring") std::weak_ptr "

STL class. ";


// File: classstd_1_1wfstream.xml
%feature("docstring") std::wfstream "

STL class. ";


// File: classstd_1_1wifstream.xml
%feature("docstring") std::wifstream "

STL class. ";


// File: classstd_1_1wios.xml
%feature("docstring") std::wios "

STL class. ";


// File: classstd_1_1wistream.xml
%feature("docstring") std::wistream "

STL class. ";


// File: classstd_1_1wistringstream.xml
%feature("docstring") std::wistringstream "

STL class. ";


// File: classstd_1_1wofstream.xml
%feature("docstring") std::wofstream "

STL class. ";


// File: classstd_1_1wostream.xml
%feature("docstring") std::wostream "

STL class. ";


// File: classstd_1_1wostringstream.xml
%feature("docstring") std::wostringstream "

STL class. ";


// File: classstd_1_1wstring.xml
%feature("docstring") std::wstring "

STL class. ";


// File: classstd_1_1wstring_1_1const__iterator.xml
%feature("docstring") std::wstring::const_iterator "

STL iterator class. ";


// File: classstd_1_1wstring_1_1const__reverse__iterator.xml
%feature("docstring") std::wstring::const_reverse_iterator "

STL iterator class. ";


// File: classstd_1_1wstring_1_1iterator.xml
%feature("docstring") std::wstring::iterator "

STL iterator class. ";


// File: classstd_1_1wstring_1_1reverse__iterator.xml
%feature("docstring") std::wstring::reverse_iterator "

STL iterator class. ";


// File: classstd_1_1wstringstream.xml
%feature("docstring") std::wstringstream "

STL class. ";


// File: classstree_1_1_d_f_s_iterator.xml
%feature("docstring") stree::DFSIterator "C++ includes:
STreeIterators.h ";

%feature("docstring")  stree::DFSIterator::DFSIterator "";

%feature("docstring")  stree::DFSIterator::next "";

%feature("docstring")  stree::DFSIterator::isFirstVisit "";

%feature("docstring")  stree::DFSIterator::setUpPass "";

%feature("docstring")  stree::DFSIterator::child "";

%feature("docstring")  stree::DFSIterator::child "";

%feature("docstring")  stree::DFSIterator::getParent "";

%feature("docstring")  stree::DFSIterator::getAncestor "";

%feature("docstring")  stree::DFSIterator::parent "";

%feature("docstring")  stree::DFSIterator::nAncestors "";

%feature("docstring")  stree::DFSIterator::parentDepth "";

%feature("docstring")  stree::DFSIterator::label "";

%feature("docstring")  stree::DFSIterator::label "";

%feature("docstring")  stree::DFSIterator::isValid "";

%feature("docstring")  stree::DFSIterator::setInvalid "";

%feature("docstring")  stree::DFSIterator::setValid "";

%feature("docstring")  stree::DFSIterator::isNode "";

%feature("docstring")  stree::DFSIterator::isLeaf "";

%feature("docstring")  stree::DFSIterator::depth "";

%feature("docstring")  stree::DFSIterator::headIndex "";

%feature("docstring")  stree::DFSIterator::index "";

%feature("docstring")  stree::DFSIterator::nidx "";

%feature("docstring")  stree::DFSIterator::indexStr "";

%feature("docstring")  stree::DFSIterator::count "";

%feature("docstring")  stree::DFSIterator::dataStr "";

%feature("docstring")  stree::DFSIterator::getChild "";

%feature("docstring")  stree::DFSIterator::getChild "";

%feature("docstring")  stree::DFSIterator::getSibling "";

%feature("docstring")  stree::DFSIterator::sibling "";

%feature("docstring")  stree::DFSIterator::getSuffixLink "";

%feature("docstring")  stree::DFSIterator::suffixLink "";

%feature("docstring")  stree::DFSIterator::string "";


// File: classstree_1_1internal_1_1_s_tree_position.xml
%feature("docstring") stree::internal::STreePosition "

A  STreePosition refers to the location in the suffix tree that
corresponds to some substring of the represented String.

This position is unique, but can either be an expicit node (internal
or leaf) or an implicit internal node. This class is only for internal
use in the suffix tree construction. Please use the class  STreePos
instead!

C++ includes: STreePosition.h ";

%feature("docstring")  stree::internal::STreePosition::STreePosition "

Create a  STreePosition corresponding to the root of the suffix tree.
";

%feature("docstring")  stree::internal::STreePosition::followChar "

If the current STreePosition corresponds to the substring str, then
attempt to move to str concatenated with chr.

Return true if this is successful, i.e., if str + chr is also a
substring of the represented text. ";

%feature("docstring")  stree::internal::STreePosition::preCanonize "";

%feature("docstring")  stree::internal::STreePosition::canonize "";

%feature("docstring")
stree::internal::STreePosition::followSuffixLink "";

%feature("docstring")  stree::internal::STreePosition::isExplicit "

Return true if the  STreePosition corresponds to an explicit node. ";


// File: classstree_1_1_postfix_iterator.xml
%feature("docstring") stree::PostfixIterator "C++ includes:
STreeIterators.h ";

%feature("docstring")  stree::PostfixIterator::PostfixIterator "";

%feature("docstring")  stree::PostfixIterator::next "";

%feature("docstring")  stree::PostfixIterator::child "";

%feature("docstring")  stree::PostfixIterator::child "";

%feature("docstring")  stree::PostfixIterator::getParent "";

%feature("docstring")  stree::PostfixIterator::getAncestor "";

%feature("docstring")  stree::PostfixIterator::parent "";

%feature("docstring")  stree::PostfixIterator::nAncestors "";

%feature("docstring")  stree::PostfixIterator::parentDepth "";

%feature("docstring")  stree::PostfixIterator::label "";

%feature("docstring")  stree::PostfixIterator::label "";

%feature("docstring")  stree::PostfixIterator::isValid "";

%feature("docstring")  stree::PostfixIterator::setInvalid "";

%feature("docstring")  stree::PostfixIterator::setValid "";

%feature("docstring")  stree::PostfixIterator::isNode "";

%feature("docstring")  stree::PostfixIterator::isLeaf "";

%feature("docstring")  stree::PostfixIterator::depth "";

%feature("docstring")  stree::PostfixIterator::headIndex "";

%feature("docstring")  stree::PostfixIterator::index "";

%feature("docstring")  stree::PostfixIterator::nidx "";

%feature("docstring")  stree::PostfixIterator::indexStr "";

%feature("docstring")  stree::PostfixIterator::count "";

%feature("docstring")  stree::PostfixIterator::dataStr "";

%feature("docstring")  stree::PostfixIterator::getChild "";

%feature("docstring")  stree::PostfixIterator::getChild "";

%feature("docstring")  stree::PostfixIterator::getSibling "";

%feature("docstring")  stree::PostfixIterator::sibling "";

%feature("docstring")  stree::PostfixIterator::getSuffixLink "";

%feature("docstring")  stree::PostfixIterator::suffixLink "";

%feature("docstring")  stree::PostfixIterator::string "";


// File: classstree_1_1_prefix_iterator.xml
%feature("docstring") stree::PrefixIterator "C++ includes:
STreeIterators.h ";

%feature("docstring")  stree::PrefixIterator::PrefixIterator "";

%feature("docstring")  stree::PrefixIterator::next "";

%feature("docstring")  stree::PrefixIterator::child "";

%feature("docstring")  stree::PrefixIterator::child "";

%feature("docstring")  stree::PrefixIterator::getParent "";

%feature("docstring")  stree::PrefixIterator::getAncestor "";

%feature("docstring")  stree::PrefixIterator::parent "";

%feature("docstring")  stree::PrefixIterator::nAncestors "";

%feature("docstring")  stree::PrefixIterator::parentDepth "";

%feature("docstring")  stree::PrefixIterator::label "";

%feature("docstring")  stree::PrefixIterator::label "";

%feature("docstring")  stree::PrefixIterator::isValid "";

%feature("docstring")  stree::PrefixIterator::setInvalid "";

%feature("docstring")  stree::PrefixIterator::setValid "";

%feature("docstring")  stree::PrefixIterator::isNode "";

%feature("docstring")  stree::PrefixIterator::isLeaf "";

%feature("docstring")  stree::PrefixIterator::depth "";

%feature("docstring")  stree::PrefixIterator::headIndex "";

%feature("docstring")  stree::PrefixIterator::index "";

%feature("docstring")  stree::PrefixIterator::nidx "";

%feature("docstring")  stree::PrefixIterator::indexStr "";

%feature("docstring")  stree::PrefixIterator::count "";

%feature("docstring")  stree::PrefixIterator::dataStr "";

%feature("docstring")  stree::PrefixIterator::getChild "";

%feature("docstring")  stree::PrefixIterator::getChild "";

%feature("docstring")  stree::PrefixIterator::getSibling "";

%feature("docstring")  stree::PrefixIterator::sibling "";

%feature("docstring")  stree::PrefixIterator::getSuffixLink "";

%feature("docstring")  stree::PrefixIterator::suffixLink "";

%feature("docstring")  stree::PrefixIterator::string "";


// File: classstree_1_1_r_b_tree.xml
%feature("docstring") stree::RBTree "

An implementation of red-black tree algorithms.

The functions can be used with generic kinds of nodes that need to be
specified by providing an RBTreeNodeTraits object when calling the
functions.

Note that the red-black tree will be threaded if a NULL RBNodePtr can
still address a node. This allows iterating over the stored values in
their order merely by following the left or right pointers (even if
NULL). The left-and rightmost RBNodePtr will be inherited from the
left and right RBNodePtr of the original root, i.e., the first node
RBNodePtr used in the first insertion operation when constructing the
red-black tree, while every threaded right or left RBNodePtr will be
set to the according RBNodePtr and then marked as a thread by calling
the function setThread. It is up to the user to implement a suitable
RBNodePtr structure to be able to distinguish the left-/rightmost
RBNodePtr from a threaded RBNodePtr (e.g., if the color of nodes is
stored in the RBNodePtr, even if NULL, then a threaded RBNodePtr may
be NULL and colored red, while the left-/rightmost RBNodePtr may be
NULL and black).

C++ includes: RBTree.h ";


// File: classstree_1_1_r_b_tree_node_traits_template.xml
%feature("docstring") stree::RBTreeNodeTraitsTemplate "

An example of which functions an RBTreeNodeTraits object must provide
to be used with  RBTree.

It is assumed that each node has a left and right pointer to child
nodes, forming a binary tree. These pointers may be any user defined
data structure used to index the nodes. For this reason, all needed
functions to manipulate the tree structure must me provided in an
RBTreeNodeTraits object. Note that nodes are referred to only by
RBNodePtr objects and never directly.

C++ includes: RBTree.h ";


// File: classstree_1_1_s_tree.xml
%feature("docstring") stree::STree "

An implementation of suffix trees for generic string types.

For string types other than std::string, one must define
STREE_STRING_TRAITS and define an appropriate  STreeStringTraits
object placed in the namespace stree. (Note that template arguments
are avoided on purpose to insure better compatibility to a python
interface). The children of an internal node are stored in a self-
balancing binary tree, so that this implementation is suitable
especially for strings over large alphabets. The space requirement of
the index structure (including suffix links) is at most 5 * 4 bytes /
input symbol, and the maximum string length is 2**29-1.

In the following comes a brief description of the internal structure
of this suffix tree implementation: The suffix tree has two types of
nodes (internal and leaf nodes), which are stored in vectors nodes and
leaves respectively.

Nodes are addressed by a 32-bit Nidx (\"node index\"). The first three
bits of a Nidx indicate whether the Nidx is valid (1) or nil (0),
addresses an (internal) node (1) or a leaf (0), and is \"colored\" (1)
or not (0), respectively. In the case of a valid Nidx, the remaining
29 bits form the index value of the addressed node in the nodes or
leaves vector. Otherwise, these can have some other significance.

Every node (internal or leaf) has a left and right Nidx. Internal
nodes have an additional child Nidx, which addresses the first child
in the suffix tree structure.

The children of any node (which are siblings) are basically organized
in a self-balancing binary tree structure formed by the left and right
Nidx \"pointers\". For this, a left-leaning red-black tree is used
(which uses the color flag of the Nidx). However, the first two
children of any node have special roles:

The right Nidx of the first child encodes the headindex of the parent
node, while the right Nidx of the second child encodes the depth of
the parent node. Both right Nidx addresses of the first two children
have their first three bits set to zero.

The left Nidx of the first child addresses the second child and has
its color flag set. The left Nidx of the second child addresses the
root of the red-black tree in which all further siblings are
organized, and has its color flag unset (to be able to distinguish
first and second children). However, if there are only two siblings,
then the left Nidx of the second child will be marked as invalid, but
will address the suffix-link of the parent. Still, the color flag will
be unset.

Suffix-links of any (internal) node are stored in the right Nidx of
the rightmost (in the binary tree) child, or alternatively in the left
Nidx of the second child if there are only two children. This Nidx
will always be marked as invalid and uncolored, but as addressing a
node.

The left Nidx of the leftmost sibling in the binary tree will always
be marked as invalid and uncolored, but has no further meaning.

The binary tree of siblings is threaded. This means that invalid left
and right Nidx entries address previous and next siblings
respectively. This is true for all siblings in the red-black tree,
except for the left-and rightmost. To distinguish these, all invalid
Nidx that indicate a thread are marked as colored (as opposed to the
left- and rightmost, which are always marked as uncolored).

Siblings are stored in lexicographic order with respect to their edge
labels. It is ensured that the first and second siblings are always
the lexicographically smallest.

C++ includes: STreeCore.h ";

%feature("docstring")  stree::STree::STree "";

%feature("docstring")  stree::STree::extendTo "

Extend the suffix tree representation of the underlying text to the
given size.

Note that the given size must be smaller or equal to the size of the
underlying text. ";

%feature("docstring")  stree::STree::~STree "";

%feature("docstring")  stree::STree::nLeaves "

Return the number of leaf nodes in the suffix tree. ";

%feature("docstring")  stree::STree::nInternalNodes "

Return the number of internal nodes in the suffix tree. ";

%feature("docstring")  stree::STree::nNodes "

Return the number of nodes (internal and leaves) in the suffix tree.
";

%feature("docstring")  stree::STree::n "";

%feature("docstring")  stree::STree::d "

Return the depth of the node, i.e., the size of the substring
represented by the node. ";

%feature("docstring")  stree::STree::hi "

Return the head index of the node, i.e., a position in the underlying
text where the substring represented by the node can be found. ";

%feature("docstring")  stree::STree::l "

Return the left Nidx& of the given node. ";

%feature("docstring")  stree::STree::r "

Return the right Nidx& of the given node. ";

%feature("docstring")  stree::STree::c "

Return the child Nidx& of the given node.

Note that node must specify an internal node and not a leaf. ";

%feature("docstring")  stree::STree::c "

Return the child Nidx& of the given node corresponding to the given
chr (which is the first character of the edge leading away from the
node).

If no corresponding child is found, the (null) Nidx& will be returned
that corresponds to the place where the according child node would
need to be inserted. ";

%feature("docstring")  stree::STree::sl "

Return the suffix link Nidx& of the given node. ";

%feature("docstring")  stree::STree::sib "

Return the next sibling (according to lexicographic ordering of edge
labels) of the given node, or a null Nidx if no further sibling
exists. ";

%feature("docstring")  stree::STree::subString "

return a substring of the represented String text of length len
starting at position pos. ";

%feature("docstring")  stree::STree::at "

Return the character at index pos in the string represented by this
suffix tree. ";


// File: classstree_1_1_s_tree_1_1_internal_node.xml


// File: classstree_1_1_s_tree_1_1_leaf_node.xml


// File: classstree_1_1_s_tree_1_1_r_b_tree_node_traits.xml


// File: classstree_1_1_s_tree_edge.xml
%feature("docstring") stree::STreeEdge "C++ includes: STreeNode.h ";

%feature("docstring")  stree::STreeEdge::STreeEdge "";

%feature("docstring")  stree::STreeEdge::STreeEdge "";

%feature("docstring")  stree::STreeEdge::STreeEdge "";

%feature("docstring")  stree::STreeEdge::STreeEdge "";

%feature("docstring")  stree::STreeEdge::parentDepth "";

%feature("docstring")  stree::STreeEdge::getParent "";

%feature("docstring")  stree::STreeEdge::getChild "";

%feature("docstring")  stree::STreeEdge::child "";

%feature("docstring")  stree::STreeEdge::getChild "";

%feature("docstring")  stree::STreeEdge::child "";

%feature("docstring")  stree::STreeEdge::getSibling "";

%feature("docstring")  stree::STreeEdge::label "";

%feature("docstring")  stree::STreeEdge::isValid "";

%feature("docstring")  stree::STreeEdge::setInvalid "";

%feature("docstring")  stree::STreeEdge::setValid "";

%feature("docstring")  stree::STreeEdge::isNode "";

%feature("docstring")  stree::STreeEdge::isLeaf "";

%feature("docstring")  stree::STreeEdge::depth "";

%feature("docstring")  stree::STreeEdge::headIndex "";

%feature("docstring")  stree::STreeEdge::index "";

%feature("docstring")  stree::STreeEdge::nidx "";

%feature("docstring")  stree::STreeEdge::indexStr "";

%feature("docstring")  stree::STreeEdge::count "";

%feature("docstring")  stree::STreeEdge::dataStr "";

%feature("docstring")  stree::STreeEdge::getChild "";

%feature("docstring")  stree::STreeEdge::getChild "";

%feature("docstring")  stree::STreeEdge::getSibling "";

%feature("docstring")  stree::STreeEdge::sibling "";

%feature("docstring")  stree::STreeEdge::getSuffixLink "";

%feature("docstring")  stree::STreeEdge::suffixLink "";

%feature("docstring")  stree::STreeEdge::string "";

%feature("docstring")  stree::STreeEdge::label "";


// File: classstree_1_1_s_tree_node.xml
%feature("docstring") stree::STreeNode "

This class represents a node in the suffix tree.

C++ includes: STreeNode.h ";

%feature("docstring")  stree::STreeNode::STreeNode "";

%feature("docstring")  stree::STreeNode::STreeNode "";

%feature("docstring")  stree::STreeNode::STreeNode "";

%feature("docstring")  stree::STreeNode::isValid "";

%feature("docstring")  stree::STreeNode::setInvalid "";

%feature("docstring")  stree::STreeNode::setValid "";

%feature("docstring")  stree::STreeNode::isNode "";

%feature("docstring")  stree::STreeNode::isLeaf "";

%feature("docstring")  stree::STreeNode::depth "";

%feature("docstring")  stree::STreeNode::headIndex "";

%feature("docstring")  stree::STreeNode::index "";

%feature("docstring")  stree::STreeNode::nidx "";

%feature("docstring")  stree::STreeNode::indexStr "";

%feature("docstring")  stree::STreeNode::count "";

%feature("docstring")  stree::STreeNode::dataStr "";

%feature("docstring")  stree::STreeNode::getChild "";

%feature("docstring")  stree::STreeNode::child "";

%feature("docstring")  stree::STreeNode::getChild "";

%feature("docstring")  stree::STreeNode::child "";

%feature("docstring")  stree::STreeNode::getSibling "";

%feature("docstring")  stree::STreeNode::sibling "";

%feature("docstring")  stree::STreeNode::getSuffixLink "";

%feature("docstring")  stree::STreeNode::suffixLink "";

%feature("docstring")  stree::STreeNode::string "";

%feature("docstring")  stree::STreeNode::label "";


// File: classstree_1_1_s_tree_path.xml
%feature("docstring") stree::STreePath "C++ includes: STreeNode.h ";

%feature("docstring")  stree::STreePath::STreePath "";

%feature("docstring")  stree::STreePath::STreePath "";

%feature("docstring")  stree::STreePath::child "";

%feature("docstring")  stree::STreePath::child "";

%feature("docstring")  stree::STreePath::getParent "";

%feature("docstring")  stree::STreePath::getAncestor "";

%feature("docstring")  stree::STreePath::parent "";

%feature("docstring")  stree::STreePath::nAncestors "";

%feature("docstring")  stree::STreePath::parentDepth "";

%feature("docstring")  stree::STreePath::label "";

%feature("docstring")  stree::STreePath::isValid "";

%feature("docstring")  stree::STreePath::setInvalid "";

%feature("docstring")  stree::STreePath::setValid "";

%feature("docstring")  stree::STreePath::isNode "";

%feature("docstring")  stree::STreePath::isLeaf "";

%feature("docstring")  stree::STreePath::depth "";

%feature("docstring")  stree::STreePath::headIndex "";

%feature("docstring")  stree::STreePath::index "";

%feature("docstring")  stree::STreePath::nidx "";

%feature("docstring")  stree::STreePath::indexStr "";

%feature("docstring")  stree::STreePath::count "";

%feature("docstring")  stree::STreePath::dataStr "";

%feature("docstring")  stree::STreePath::getChild "";

%feature("docstring")  stree::STreePath::getChild "";

%feature("docstring")  stree::STreePath::getSibling "";

%feature("docstring")  stree::STreePath::sibling "";

%feature("docstring")  stree::STreePath::getSuffixLink "";

%feature("docstring")  stree::STreePath::suffixLink "";

%feature("docstring")  stree::STreePath::string "";

%feature("docstring")  stree::STreePath::label "";


// File: classstree_1_1_s_tree_pos.xml
%feature("docstring") stree::STreePos "C++ includes: STreeNode.h ";

%feature("docstring")  stree::STreePos::STreePos "";

%feature("docstring")  stree::STreePos::STreePos "";

%feature("docstring")  stree::STreePos::setRoot "";

%feature("docstring")  stree::STreePos::isValid "";

%feature("docstring")  stree::STreePos::setValid "";

%feature("docstring")  stree::STreePos::setInvalid "";

%feature("docstring")  stree::STreePos::isExplicit "";

%feature("docstring")  stree::STreePos::isLeaf "";

%feature("docstring")  stree::STreePos::count "";

%feature("docstring")  stree::STreePos::headIndex "";

%feature("docstring")  stree::STreePos::depth "";

%feature("docstring")  stree::STreePos::parentDepth "";

%feature("docstring")  stree::STreePos::suffixLink "";

%feature("docstring")  stree::STreePos::addChar "";

%feature("docstring")  stree::STreePos::addString "";

%feature("docstring")  stree::STreePos::string "";

%feature("docstring")  stree::STreePos::label "";


// File: classstree_1_1_s_tree_string_traits.xml
%feature("docstring") stree::STreeStringTraits "

This class specifies the required string traits for std::string.

To use another string type, one must define STREE_STRING_TRAITS and
include an appropriate  STreeStringTraits object in the namespace
stree.

C++ includes: stree.h ";


// File: structtom_1_1_c_i_p__by___s_v_d___e_m.xml
%feature("docstring") tom::CIP_by_SVD_EM "C++ includes:
ContOomTrain.h ";


// File: classtom_1_1_cont_oom.xml
%feature("docstring") tom::ContOom "

This class provides the basic functionality of continuous observable
operator models, which are blended from standard discrete OOM.

The main aspects are extracting information from a given continuous
OOM (sequence generation, prediction, fingerprints, etc.), and

transformation functions (stabilization, equivalence transforms,
dimension reduction, etc.) To estimate (\"learn\") an OOM from data,
use the ContOomTrain class.

C++ includes: ContOom.h ";

%feature("docstring")  tom::ContOom::ContOom "

creates an uninitialized (!) ContOom. ";

%feature("docstring")  tom::ContOom::ContOom "

read a continuous OOM from a file and initialize it.

The file format must correspond to what the output functions produce.

Parameters:
-----------

filename:  the file to read the OOM parameters from ";

%feature("docstring")  tom::ContOom::sample "

sample the next output from the distribution modeled by this
continuous OOM (conditioned on the current in the case of an input-
output model).

Parameters:
-----------

u:  the current input in the case of an input-output model

the sampled output ";

%feature("docstring")  tom::ContOom::predict "

predict the next output.

This gives the expectation of the next output according to the
distribution modeled by this continuous OOM (conditioned on the
current in the case of an input output model).

Parameters:
-----------

u:  the current input in the case of an input-output model

the expectation of the next output ";

%feature("docstring")  tom::ContOom::generate "

generate a sample sequence of a given length according to the
continuous OOM. ";

%feature("docstring")  tom::ContOom::loglik "

robustly compute the average log-likeliehood of this model for the
given sequence seq.

Parameters:
-----------

seq:  ";

%feature("docstring")  tom::ContOom::initializeOomType "";

%feature("docstring")  tom::ContOom::readOomType "";

%feature("docstring")  tom::ContOom::writeOomType "";

%feature("docstring")  tom::ContOom::rawStateUpdateOomType "";


// File: classtom_1_1_estimator.xml
%feature("docstring") tom::Estimator "

This class computes estimates $\\\\hat{f}(\\\\bar{x})$ and
corresponding variance estimates for sequences $\\\\bar{x}$ based on a
suffix tree representation of a sample sequence.

The following parameters control how the estimates and their variances
are computed:  iidPi_ is only relevant for input-output sequences, and
means that the input policy should be assumed to be \"blind\" (iid
inputs). The estimates will be computed accordingly. The input Symbol
probabilities are estimated from the sample sequence and stored in
uProb_, which can be overwritten if the probabilities are known.

useBayes_ determines whether Bayesian estimates are returned (note
that Bayesian estimates are always used for the variance estimation
for numerical reasons).

priorStrength_ and proportionalPrior_ are used to control the Bayesian
priors used for the Bayesian estimates (which are always used for the
variance estimation).

varEst_ determines the method used to compute the variance. If this is
set to zero, an exact formula is used, otherwise a rough approximation
is computed. For the case of default estimates for input-output
sequences the following options for rough variance approximations
exist: 1: v() = $\\\\hat{f}^2$, where $\\\\hat{f}$ is a Bayesian
estimate

2: v() = $\\\\hat{f}$

3: v() = $\\\\#\\\\bar{x}$, where $\\\\#\\\\bar{x}$ are the number of
occurrences of $\\\\bar{x}$ in the sample sequence plus some pseudo-
counts

4: v() = $\\\\#\\\\bar{x}/\\\\hat{\\\\pi}$, where $\\\\hat{\\\\pi}$ is
a Bayesian estimate of the input policy

5: v() = $\\\\hat{f}/\\\\hat{\\\\pi}$

6: v() = $\\\\#\\\\bar{x}/\\\\hat{\\\\pi}^2$

7: v() = $\\\\#\\\\bar{x}/\\\\hat{\\\\pi}^2$, where more pseudo-counts
are added to $\\\\#\\\\bar{x}$

C++ includes: Estimators.h ";

%feature("docstring")  tom::Estimator::Estimator "

create an  Estimator from a given sfxTree -- a suffix tree
representation of a sample sequence ";

%feature("docstring")  tom::Estimator::reset "

reset the  Estimator to the empty sequence. ";

%feature("docstring")  tom::Estimator::f "

return the current estimate ";

%feature("docstring")  tom::Estimator::f "

return the matrix of estimates
$\\\\hat{F}^{I,J}=[\\\\hat{f}(\\\\bar{x}_j\\\\bar{x}_i)]_{i,j}$ for
the given set chaSeqs of characteristic sequences $\\\\bar{x}_j$ and
(optionally) set indSeqs of indicative sequences $\\\\bar{x}_j$. ";

%feature("docstring")  tom::Estimator::f "

return the matrix of estimates
$\\\\hat{F}_z^{I,J}=[\\\\hat{f}(\\\\bar{x}_jz\\\\bar{x}_i)]_{i,j}$ for
the given set chaSeqs of characteristic sequences $\\\\bar{x}_j$, the
set indSeqs of indicative sequences $\\\\bar{x}_j$ and the input-
output pair z = ( u, o). ";

%feature("docstring")  tom::Estimator::v "

return a variance estimate for the current estiamte ";

%feature("docstring")  tom::Estimator::v "

return the matrix of variance estimates
$\\\\hat{V}^{I,J}[\\\\hat{\\\\rm{Var}}(\\\\hat{f}(\\\\bar{x}_j\\\\bar{x}_i))]_{i,j}$
corresponding to the estimates $\\\\hat{f}(\\\\bar{x}_j\\\\bar{x}_i)$
for the given set chaSeqs of characteristic sequences $\\\\bar{x}_j$
and (optionally) set indSeqs of indicative sequences $\\\\bar{x}_j$.
";

%feature("docstring")  tom::Estimator::v "

return the matrix of variance estimates
$\\\\hat{V}^{I,J}[\\\\hat{\\\\rm{Var}}(\\\\hat{f}(\\\\bar{x}_jz\\\\bar{x}_i))]_{i,j}$
corresponding to the estimates $\\\\hat{f}(\\\\bar{x}_jz\\\\bar{x}_i)$
for the given set chaSeqs of characteristic sequences $\\\\bar{x}_j$,
the set indSeqs of indicative sequences $\\\\bar{x}_j$, and the input-
output pair z = ( u, o). ";

%feature("docstring")  tom::Estimator::c "

return the current count number ";


// File: classtom_1_1_learner.xml
%feature("docstring") tom::Learner "C++ includes: OomTrain.h ";


// File: classtom_1_1_oom.xml
%feature("docstring") tom::Oom "

This class provides the basic functionality of observable operator
models.

The main aspects are extracting information from a given OOM (sequence
generation, prediction, fingerprints, etc), and

transformation functions (stabilization, equivalence transforms,
dimension reduction, etc.) To estimate (\"learn\") an OOM from data,
use the OomTrain class.

C++ includes: Oom.h ";

/*  Constructors and initialization  */

%feature("docstring")  tom::Oom::Oom "

Construct an uninitialized (!) Oom. ";

%feature("docstring")  tom::Oom::Oom "

Construct a simple OOM of dimension dim that models a stochastic
process of independently uniformly distributed random variables with
nO number of possible observations (and nU number of possible inputs,
that are ignored).

Parameters:
-----------

nO:  the size of the output alphabet

nU:  the size of the input alphabet (0 for an output-only OOM)

dim:  the dimension of the OOM to construct ";

%feature("docstring")  tom::Oom::Oom "

read an OOM from file. The file format must correspond to what the
output functions produce. ";

%feature("docstring")  tom::Oom::Oom "

Create an OOM with the given ( sig, tauArray, w0) and initialize it.

The tauArray must be a matrix of of size ( nO * dim x nU * dim), where
tau(o,u) is the dim x dim submatrix at position ( nO * dim), nU *
dim). ";

%feature("docstring")  tom::Oom::Oom "

Create an OOM from a given HMM.

The HMM parameters must be specified as follows:

Parameters:
-----------

T:  the state transition probabilities: $T_{i,j} = P(s_{t+1}=j |
s_{t}=i)$

E:  the emission probabilities: $E_{i,j} = P(o_t=j | s_t=i)$

w:  the initial state probabilities: $w_i = P(s_0 = i)$. ";

%feature("docstring")  tom::Oom::setSize "

setup the internal structure for an OOM of the desired size without
performing any initialization. Typically, the parameters sig, tau(o,u)
and w0 will be assigned next, and then  initialize() must be called.
";

%feature("docstring")  tom::Oom::initialize "

initialize the OOM. This assumes that all essential parameters (e.g.,
nU_, nO_, dim_, sig_, tau_, w0_ ) have been set. ";

/*  Accessors  */

%feature("docstring")  tom::Oom::nU "

return the size of the input alphabet ";

%feature("docstring")  tom::Oom::nO "

return the size of the output alphabet ";

%feature("docstring")  tom::Oom::dim "

return the model dimension ";

%feature("docstring")  tom::Oom::sig "

return a reference to the evaluation functional vector $\\\\sigma$ ";

%feature("docstring")  tom::Oom::sig "

set the evaluation functional vector $\\\\sigma$ to s and reset() ";

%feature("docstring")  tom::Oom::tau "

return a reference to the observable operator corresponding to
observation o and input u. ";

%feature("docstring")  tom::Oom::w0 "

return a reference to the initial state ";

%feature("docstring")  tom::Oom::w0 "

set the initial state to w and reset() ";

%feature("docstring")  tom::Oom::wt "

return a reference to the current state ";

%feature("docstring")  tom::Oom::wt "

set the current state to w ";

%feature("docstring")  tom::Oom::prediction "

return a reference to the current prediction vector of the next output
symbol probabilities, i.e., $P(\\\\cdot|u_t, \\\\omega_t)$. ";

%feature("docstring")  tom::Oom::prediction "

return the probability of the next output symbol o, i.e., $P(o|u_t,
\\\\omega_t)$. ";

%feature("docstring")  tom::Oom::maxSetback "

return the maximum number of steps to \"replay\" during a setBack()
operation ";

%feature("docstring")  tom::Oom::maxSetback "

set the maximum number of steps to \"replay\" during a setBack()
operation to value ";

/*  Basic OOM functionality  */

%feature("docstring")  tom::Oom::reset "

reset the OOM to its initial state and reset the error counters. ";

%feature("docstring")  tom::Oom::resetCounters "

reset the error counters ";

%feature("docstring")  tom::Oom::updateState "

update the OOM state according to the input-output pair ( u, o), i.e.,
set $\\\\omega_t = \\\\tau_{u,o} \\\\omega_t$. ";

%feature("docstring")  tom::Oom::evaluate "

return $\\\\sigma\\\\omega_t$, the evaluation of the current state
$\\\\omega_t$ by the evaluation functional $\\\\sigma$. ";

%feature("docstring")  tom::Oom::normalizeState "

attempt to normalize the current state, i.e., set $\\\\omega_t =
\\\\omega_t / \\\\sigma\\\\omega_t$, and return true if successful,
i.e., $\\\\sigma\\\\omega_t$> setbackProbMargin_ + zeroMargin_. ";

%feature("docstring")  tom::Oom::updatePrediction "

update the prediction vector of the next output symbol probabilities
according to the current state $\\\\omega_t$ and input u, i.e.,
compute $P(\\\\cdot|u_t, \\\\omega_t)$. ";

%feature("docstring")  tom::Oom::fixPrediction "

attempt to fix the prediction vector of the next output symbol
probabilities $P(\\\\cdot|u_t, \\\\omega_t)$ such that all
probabilities are at least minPredictionProb_ and the probabilities
sum to one. Return a measure of the required change to the prediction
vector: 1.5 * nO() * squared norm of the difference. ";

%feature("docstring")  tom::Oom::setBack "

attempt to perform a state setback operation for at most maxSetback_
time-steps. Note that calling this method repeatedly will attempt a
setback for a shorter history each time. Return true if a setback
could be performed. ";

%feature("docstring")  tom::Oom::condition "

perform the required stabilization and update the prediction vector of
next output probabilities conditioned on the current state and input
u. This method should be called at each time step before updating the
state or using the prediction.

This method does the following: normalize the current state

update the prediction vector according to u

fix the prediction vector

perform setback operations where required

update the nSetback and nFixPrediction counters. ";

%feature("docstring")  tom::Oom::generate "

generate a sample sequence of a given length according to the OOM and
an input policy (by default an iud-process is used as input) starting
from the current state $\\\\omega_t$. ";

%feature("docstring")  tom::Oom::f "

return the function value (probability) for the input-output pair ( u,
o) given the current state, i.e., $P(o|u, \\\\omega_t)$, and update
the state. ";

%feature("docstring")  tom::Oom::f "

return the function value (probability) for the given sequence seq
given the current state, i.e., $f(seq | \\\\omega_t)$, and update the
state. ";

%feature("docstring")  tom::Oom::log_f "

return the log function value of the given sequence seq given the
current state, i.e.

$\\\\log f(seq | \\\\omega_t)$, and update the state. To deal
gracefully with time-steps n in the sequence seq where $ P(o_n | u_n,
\\\\omega_{t+n}) \\\\approx 0$, all such probabilities are treated as
impossibleProbMargin_ in the computation of the function value. Every
time this happens, the counter nImpossible_ is incremented. Note that
this problem can be avoided by increasing minPredictionProb_. ";

%feature("docstring")  tom::Oom::averageOneStepPredictionError "

calculate the average one-step squared prediction error along the
given sample sequence seq according to a correct  Oom model gen. Both
this  Oom and the given gen are evaluated from their current states,
and their states are updated accordingly. ";

/*  internal functionality  */

%feature("docstring")  tom::Oom::validate "";

/*  IO-functions  */

%feature("docstring")  tom::Oom::from_string "

read the parameters from the given string.

The format must correspond to what the output functions produce. ";

%feature("docstring")  tom::Oom::to_string "

output the parameters as a string. ";

%feature("docstring")  tom::Oom::show "

write the OOM parameters to std::cout ";

%feature("docstring")  tom::Oom::load "

load the OOM parameters from a file and initialize it.

The format must correspond to what the output functions produce.

Parameters:
-----------

filename:  the file to load the OOM parameters from ";

%feature("docstring")  tom::Oom::save "

write the OOM parameters to a file.

Parameters:
-----------

filename:  the file to write the OOM parameters to ";


// File: classtom_1_1_random.xml
%feature("docstring") tom::Random "C++ includes: Random.h ";


// File: classtom_1_1_sequence.xml
%feature("docstring") tom::Sequence "

This class represents a sequence to be used with the OOM algorithms.
It may represent a sequence, subsequence view or io-sequence, and
stores information about the output and possibly input alphabet.

C++ includes: Sequence.h ";

/*  Constructors  */

%feature("docstring")  tom::Sequence::Sequence "

Construct a  Sequence with output alphabet size nO and input alphabet
size nU from a given data vector. The data vector is copied, and the
sequence is viewed as an input-output sequence if nU != 0. ";

%feature("docstring")  tom::Sequence::Sequence "

Construct a zero  Sequence with output alphabet size nO and input
alphabet size nU of a given lengthIO. The size of the  Sequence will
be 2 * lengthIO if it is an input-output sequence, i.e., if nU != 0.
";

%feature("docstring")  tom::Sequence::Sequence "

Load a  Sequence from file filename. The format of the file must
correspond to what the output functions produce. ";

%feature("docstring")  tom::Sequence::Sequence "

Make a copy of the given  Sequence seq. This will be a subsequence if
the subsequence if seq is a subsequence or the given parameter sub is
true. ";

%feature("docstring")  tom::Sequence::~Sequence "

Destructor. Note that you are responsible to insure that no
subsequence views are using memory of this  Sequence. ";

/*  Accessors  */

%feature("docstring")  tom::Sequence::nU "

Return the size of the input alphabet. If this is zero, then the
Sequence is an ordinary sequence, else an io-sequence. ";

%feature("docstring")  tom::Sequence::nO "

Return the size of the output alphabet. ";

%feature("docstring")  tom::Sequence::size "

Return the size of the represented (sub)-sequence. Note that for an
io-sequence $ u_0o_0\\\\ldots u_{N-1}o_{N-1}$ this is 2*N. ";

%feature("docstring")  tom::Sequence::isSub "

Return true if this is a subsequence view, i.e., it uses memory from
another  Sequence object. ";

%feature("docstring")  tom::Sequence::isReversed "

Return true if this is a reversed  Sequence. ";

%feature("docstring")  tom::Sequence::pos "

Return the index of the first symbol of the represented (sub)-sequence
in the underlying data vector. ";

%feature("docstring")  tom::Sequence::at "

Return the (input or output) symbol at position n. ";

%feature("docstring")  tom::Sequence::at "

Access the (input or output) symbol at position n. ";

/*  Functionality  */

%feature("docstring")  tom::Sequence::reverse "

Reverse the Sequnece. Note that does not actually change the
underlying data vector, and therefore does not affect subsequence
views. ";

%feature("docstring")  tom::Sequence::sub "

Return a subsequence view to this  Sequence. ";

%feature("docstring")  tom::Sequence::sub "

Return a subsequence view to the subsequence starting at the given
position pos and of the given size. ";

%feature("docstring")  tom::Sequence::append "

try to append the Symbol or input-output Symbol pair ( u, o) to the
Sequence. For subsequences this will only be successful if the
subsequence is not reversed and the requested extension is the
extension in the underlying sequence. ";

/*  Special functions for input-outpus sequences  */

%feature("docstring")  tom::Sequence::isStrictlyIO "

Return true if this is an input-output sequence, i.e., if the input
alphabet size nU is zero. ";

%feature("docstring")  tom::Sequence::isAlignedIO "

Return true if this sequence is aligned. This is always true for
normal (non-io) sequences. For io-sequences this means that the
sequence begins with an input-output pair, i.e., this  Sequence is
aligned if it is not reversed and begins with an input symbol, or is
reversed and begins with an output symbol. ";

%feature("docstring")  tom::Sequence::isValidIO "

Return true if this sequence is a valid (io)-sequence. This is always
true for normal (non-io) sequences. For io-sequences this means that
the sequence is aligned and has an even number of symbols. ";

%feature("docstring")  tom::Sequence::subIO "

Return a subsequence view to the subsequence starting at the given
index posIO and of the given length, i.e., if this  Sequence is $
u_0o_0\\\\ldots u_{N-1}o_{N-1}$, then return a subsequence view to $
u_{posIO}o_{posIO}\\\\ldots u_{posIO+length-1}o_{posIO+length-1}$. For
normal (non-io) sequences this is the same as sub. ";

%feature("docstring")  tom::Sequence::lengthIO "

Return the length of this sequence, which is its size for normal (non-
io) sequences and half its size for io-sequences. Note that this only
makes sense for valid (io)-sequences. ";

%feature("docstring")  tom::Sequence::o "

Return the n-th output symbol. Note that this is only correct for
aligned (io)-sequences. ";

%feature("docstring")  tom::Sequence::o "

Access the n-th output symbo. Note that this is only correct for
aligned (io)-sequences. ";

%feature("docstring")  tom::Sequence::u "

Return the n-th input symbol. Note that this is only correct for
aligned and strictly io-sequences. ";

%feature("docstring")  tom::Sequence::u "

Access the n-th input symbol. Note that this is only correct for
aligned and strictly io-sequences. ";

/*  IO-functions  */

%feature("docstring")  tom::Sequence::from_string "

read the  Sequence from the given string.

The format must correspond to what the output functions produce. ";

%feature("docstring")  tom::Sequence::to_string "

output the  Sequence as a string. ";

%feature("docstring")  tom::Sequence::show "

write the  Sequence to std::cout ";

%feature("docstring")  tom::Sequence::load "

load the  Sequence from a file and initialize it.

The format must correspond to what the output functions produce.

Parameters:
-----------

filename:  the file to load the  Sequence from ";

%feature("docstring")  tom::Sequence::save "

write the  Sequence to a file.

Parameters:
-----------

filename:  the file to write the  Sequence to ";


// File: namespace_eigen.xml


// File: namespacestd.xml


// File: namespacestree.xml
%feature("docstring")  stree::internal::showTree "

Display the suffix-tree structure; show full detail if verbose is
true. ";


// File: namespacestree_1_1internal.xml


// File: namespacetom.xml
%feature("docstring")  tom::estimate_E "";

%feature("docstring")  tom::compute_T "";

%feature("docstring")  tom::compute_CIP_by_SVD_EM "";

%feature("docstring")  tom::estimateContOom "";

%feature("docstring")  tom::compute_ContOom "";

%feature("docstring")  tom::reverseSequences "

reverse all the sequences in the Sequences vector seqs. ";

%feature("docstring")  tom::coreSequences "

find core sequences to be used in the OOM learning algorithms.

This function can be used in several ways. First, one can give a
suffix tree representation of the training data, and a core sequence
length as minSeqLen, leaving the other parameters at their default;
this will return all sequences of the given length that occurr as
subsequences in the training sequence. Secondly, one can specify
minSeqLen, maxSeqLen and minOccurrences; in this case

Parameters:
-----------

sfxTree:  a suffix tree representation of the training sequence

minSeqLen:  the minimum length for core sequences

maxSeqLen:  the maximum length for core sequences or -1, if maxSeqLen
should equal minSeqLen

minCounts:  the minimum number of occurrence counts in the training
sequence for core sequences

the core sequences ";

%feature("docstring")  tom::ipow "

return $a^b$ for non-negative integer bases and exponents. Note that
$0^0 := 1$. ";

%feature("docstring")  tom::kron "

return the Kronecker-product $A\\\\otimes B$ of the matrices A and B.
";

%feature("docstring")  tom::k_kron_times_vec_inplace "";

%feature("docstring")  tom::pinvQR "

return the pseudoinverse $M^\\\\dagger$ of the given matrix M, which
must have full column or row rank. This is computed via QR
decomposition and is not numerically stable if M is ill-conditioned.
";

%feature("docstring")  tom::pinv "

return the pseudoinverse $M^\\\\dagger$ of the given matrix M,
computed via SVD.

All singular values below tolerance are treated as zero. If tolerance
= -1, the default tolerance is used. ";

%feature("docstring")  tom::mldivideOD "

return $X$ that minimizes ${\\\\|AX - B\\\\|}_{\\\\rm{F}}$, i.e.,
which solves the (overdetermined!) system of linear equations $AX=B$
in the least-squares sense. This function assumes that A is an m-by- n
matrix with $ m \\\\ge n$ that has full column rank and is not too
ill-conditioned. This is computed from the QR decomposition of A with
column pivoting. ";

%feature("docstring")  tom::mrdivideOD "

return $X$ that minimizes ${\\\\|XA - B\\\\|}_{\\\\rm{F}}$, i.e.,
which solves the (overdetermined!) system of linear equations $XA=B$
in the least-squares sense. This function assumes that A is an m-by- n
matrix with $ n \\\\ge m$ that has full row rank and is not too ill-
conditioned. This is computed from the QR decomposition of the
transpose of A with column pivoting. ";

%feature("docstring")  tom::compute_CQ_by_SVD "";

%feature("docstring")  tom::compute_CQ_by_EC "";

%feature("docstring")  tom::learnOom "";


// File: _cont_oom_8h.xml


// File: _cont_oom_train_8h.xml


// File: _core_sequences_8h.xml


// File: _estimators_8h.xml


// File: _linear_algebra_8h.xml


// File: _macros_8h.xml


// File: _membership_functions_8h.xml


// File: _oom_8h.xml


// File: _oom_train_8h.xml


// File: _random_8h.xml


// File: _r_b_tree_8h.xml


// File: _sequence_8h.xml


// File: stree_8h.xml


// File: _s_tree_core_8h.xml


// File: _s_tree_iterators_8h.xml


// File: _s_tree_node_8h.xml


// File: _s_tree_position_8h.xml


// File: _s_tree_tools_8h.xml


// File: _s_v_d___a_c_c_e_l_e_r_a_t_e_8h.xml


// File: tom_8h.xml


// File: dir_68267d1309a1af8e8297ef4c3efbcdba.xml


// File: dir_3ad80c9ea3d4b88a91a0297702954b3c.xml


// File: indexpage.xml

