
// File: index.xml

// File: classstree_1_1_d_f_s_iterator.xml


%feature("docstring") stree::DFSIterator "
``DFSIterator(stree)``  

An iterator to traverse the nodes (leaf and internal) of the suffix tree depth
first.  

This means that every node (leaf and internal) is visited exactly twice (once on
the downward pass, and once on the upward pass). The method ``isFirstVisit()``
can be used to check if the current node is being visited for the first or
second time. Furthermore, the method ``setUpPass()`` can be used to stop
traversing deeper into the current tree branch. A ``DFSIterator`` is a
``PathNode`` with a ``toNext()`` method to move to the next node of the depth
first traversal. The end of iteration is signaled by marking the iterator as
invalid, which can be checked by the inherited method ``isValid()``.  

Note that this iterator is not wrapped as a native python iterator, i.e., the
following is *not* possible:  

    for node in DFSIterator(stree):
        print(node.nidxStr())  

Constructors
------------
* ``DFSIterator(stree)``  
    
    Create a ``DFSIterator`` for the given ``stree``.  

C++ includes: STreeIterators.h
";

%feature("docstring") stree::DFSIterator::sequence "
``sequence() -> Sequence``  

Return the (sub-)sequence represented by this node.  

Note that this is ``seq.rawSub(headIndex(), depth())``, where ``seq`` is the
sequence represented by the suffix tree.  
";

%feature("docstring") stree::DFSIterator::label "
``label() -> Sequence``  

Return the edge label for the edge leading to the current node.  

If no edge exists, this will be an empty ``Sequence``.  
";

%feature("docstring") stree::DFSIterator::toSuffix "
``toSuffix()``  

Set this to the ``PathNode`` corresponding to the first suffix of the
represented sequence.  

If no suffix exists (i.e., this is the root) mark this ``PathNode`` as invalid
instead. This uses the \"suffix link\" of the suffix tree, but needs to
recompute the path.  
";

%feature("docstring") stree::DFSIterator::toNext "
``toNext()``  

Set this ``DFSIterator`` to the next node in the depth first traversal.  

If the current node has been visited for the first time and ``setUpPass()`` was
not called, the next node will be the next in prefix order (or the same leaf
node, if the current node is a leaf node), otherwise the next in postfix order.
Note that every node (including leaf nodes) is visited exactly twice, i.e., once
on the downwards pass and once on the upwards pass. If no next node exists, this
``DFSIterator`` will be marked as invalid. Calling ``toNext()`` on an invalid
``DFSIterator`` has no effect.  
";

%feature("docstring") stree::DFSIterator::isSuffix "
``isSuffix() -> bool``  

Return ``true`` if the subsequence represented by this node is a suffix of the
underlying sequence.  

Note that this does not imply that this is a leaf.  
";

%feature("docstring") stree::DFSIterator::setValid "
``setValid(valid=true)``  

Mark this ``Node`` as ``valid`` (or ``invalid``, if ``valid`` is ``false``).  
";

%feature("docstring") stree::DFSIterator::setRoot "
``setRoot()``  

Reset this ``PathNode`` to the root of the suffix tree.  
";

%feature("docstring") stree::DFSIterator::isLeaf "
``isLeaf() -> bool``  

Return ``true`` if this is a leaf node.  
";

%feature("docstring") stree::DFSIterator::toParent "
``toParent()``  

Set this ``PathNode`` to the path to the parent of the current node, if such
exists.  

Otherwise, just mark this ``PathNode`` as invalid.  
";

%feature("docstring") stree::DFSIterator::isValid "
``isValid() -> bool``  

Return ``true`` if valid, otherwise return ``false``.  
";

%feature("docstring") stree::DFSIterator::index "
``index() -> nidx_t``  

The ``index`` of a valid leaf or a valid internal node is a unique number
between 0 and ``STree.nLeafNodes()`` or between 0 and
``STree.nInternalNodes()``, respectively.  
";

%feature("docstring") stree::DFSIterator::isRoot "
``isRoot() -> bool``  

Return ``true`` if this is the root node.  
";

%feature("docstring") stree::DFSIterator::nidxStr "
``nidxStr(width=3) -> std::string``  

Return a string representation of the underlying ``nidx_t``.  
";

%feature("docstring") stree::DFSIterator::setUpPass "
``setUpPass()``  

Calling this method on the first visit of a node causes the traversal of nodes
below the current node to be skipped, i.e., iteration continues as if this node
has been visited for the second time (on the upward pass).  

Calling this on the second visit has no effect.  
";

%feature("docstring") stree::DFSIterator::headIndex "
``headIndex() -> nidx_t``  

Return the \"headindex\" of this node, which is an index in the sequence
represented by the suffix tree where the (sub-)sequence represented by this node
occurs.  

I.e., the (sub-)sequence represented by this node is ``seq.rawSub(headindex(),
depth())``, where ``seq`` is the sequence represented by the suffix tree.  
";

%feature("docstring") stree::DFSIterator::isFirstVisit "
``isFirstVisit() -> bool``  

Return ``true`` if the current node is being visited for the first time (i.e.,
on the downward pass).  
";

%feature("docstring") stree::DFSIterator::toSibling "
``toSibling()``  

Set this node to its next sibling if a next sibling exists, otherwise mark this
node as invalid.  

Note that the siblings are ordered lexicographically according to their edge
labels.  
";

%feature("docstring") stree::DFSIterator::sibling "
``sibling() -> Node``  

Return the next sibling of this node.  

If no such node exists, a ``Node`` marked as invalid is returned. Note that the
siblings are ordered lexicographically according to their edge labels.  
";

%feature("docstring") stree::DFSIterator::child "
``child() -> Node``  
``child(symbol) -> Node``  

Overloaded function
-------------------
* ``child() -> Node``  
    
    Return the first child node of this node.  

    If no such node exists, a ``Node`` marked as invalid is returned. Note that
    the children are ordered lexicographically according to their edge labels.  

* ``child(symbol) -> Node``  
    
    Return the child node along the edge leading away whose label begins with
    the given ``symbol``.  

    If no such node exists, a ``Node`` marked as invalid is returned.  
";

%feature("docstring") stree::DFSIterator::count "
``count() -> nidx_t``  

Return the number of occurrences of the sequence represented by this node in the
sequence represented by the suffix tree.  

For an invalid node, zero is returned.  
";

%feature("docstring") stree::DFSIterator::DFSIterator "
``DFSIterator(stree)``  

Create a ``DFSIterator`` for the given ``stree``.  
";

%feature("docstring") stree::DFSIterator::depth "
``depth() -> nidx_t``  

Return the \"depth\" of the node in the suffix tree, which is the size of the
represented (sub-)sequence.  
";

%feature("docstring") stree::DFSIterator::set "
``set(pathNode)``  
``set(node)``  

Overloaded function
-------------------
* ``set(pathNode)``  
    
    Set this ``PathNode`` to the given ``pathNode``, which must belong to the
    same suffix tree.  

* ``set(node)``  
    
    Set this ``Node`` to the given ``node``, which must belong to the same
    suffix tree.  

    This is a faster version of ``(*this) = node)``.  
";

%feature("docstring") stree::DFSIterator::repr "
``repr() -> std::string``  

Return a string representation to display in python.  
";

%feature("docstring") stree::DFSIterator::dataStr "
``dataStr(width=5) -> std::string``  

Return a string representation of the data of this node.  

This is useful for debugging or understanding the suffix tree structure.  
";

%feature("docstring") stree::DFSIterator::suffix "
``suffix() -> Node``  

Return the node corresponding to the first suffix of the represented sequence.  

This follows the \"suffix link\" of the suffix tree. If no such node exists, a
``Node`` marked as invalid is returned.  
";

%feature("docstring") stree::DFSIterator::nidx "
``nidx() -> nidx_t``  

Return the ``nidx_t`` corresponding to this ``Node``.  
";

%feature("docstring") stree::DFSIterator::parent "
``parent() -> Node``  

Return the parent ``Node``.  

If none exists, return a ``Node`` marked as invalid.  
";

%feature("docstring") stree::DFSIterator::isInternal "
``isInternal() -> bool``  

Return ``true`` if this is an internal node.  
";

%feature("docstring") stree::DFSIterator::toChild "
``toChild()``  
``toChild(chr)``  

Overloaded function
-------------------
* ``toChild()``  
    
    Extend this ``PathNode`` to its first child if such a node exists, otherwise
    mark this ``PathNode`` as invalid instead.  

    Note that the children are ordered lexicographically according to their edge
    labels.  

* ``toChild(chr)``  
    
    Extend this ``PathNode`` to the child node along the edge leading away whose
    label begins with the given ``symbol``.  

    If no such node exists, mark this ``PathNode`` as invalid instead.  
";

// File: classstree_1_1_edge_node.xml


%feature("docstring") stree::EdgeNode "
``EdgeNode(stree, nidx=ROOT, parent=ROOT)``  
``EdgeNode(node, parent=ROOT)``  
``EdgeNode(node, parent)``  

This class represents a node together with the edge leading to it in the suffix
tree and can be used for extracting information or navigating the suffix tree.  

Note that an ``EdgeNode`` can have a parent marked as invalid, and then contains
no edge information. Such an ``EdgeNode`` is called \"degenerate\". This is
normally only the case for the root, but can also happen when providing
incorrect parent information when constructing an ``EdgeNode``.  

The ``to...()`` methods are provided to navigate the suffix tree structure:  

*   ``to...()`` sets this node (and edge leading to it) to its e.g. child,
    suffix, sibling, etc.  
*   If no such exists, then this node is simply marked as invalid.  
*   For an invalid node, the ``to...()`` methods have no effect.  
*   Calling ``setValid()`` after the node has been marked as invalid by a
    ``to...()`` method will reset this node (and edge) to the last valid node
    during the traversal.  

Constructors
------------
* ``EdgeNode(stree, nidx=ROOT, parent=ROOT)``  
    
    Construct an ``EdgeNode`` for the given ``stree`` corresponding to the given
    ``nidx`` and parent information in ``parent``.  

    If no ``nidx`` is given, it defaults to the root of the suffix tree. This
    will search for the parent starting form the given ``parent`` (or root by
    default). If no parent can be found (e.g. for the root), this ``EdgeNode``
    will be \"degenerate\" (have an invalid parent).  

* ``EdgeNode(node, parent=ROOT)``  
    
    Construct an ``EdgeNode`` from the given ``Node``.  

    This will search for the parent starting form the given ``parent`` (or root
    by default). If no parent can be found (e.g. for the root), this
    ``EdgeNode`` will be \"degenerate\" (have an invalid parent).  

* ``EdgeNode(node, parent)``  
    
    Construct an ``EdgeNode`` from the given ``Node``.  

    This will search for the parent starting form the given ``parent``. If no
    parent can be found (e.g. for the root), this ``EdgeNode`` will be
    \"degenerate\" (have an invalid parent).  

C++ includes: STreeNode.h
";

%feature("docstring") stree::EdgeNode::toSibling "
``toSibling()``  

Set this node to its next sibling if a next sibling exists, otherwise mark this
node as invalid.  

Note that the siblings are ordered lexicographically according to their edge
labels.  
";

%feature("docstring") stree::EdgeNode::toChild "
``toChild()``  
``toChild(symbol)``  

Overloaded function
-------------------
* ``toChild()``  
    
    Set this ``EdgeNode`` to its first child (with edge leading to it) if such a
    node exists, otherwise mark this ``EdgeNode`` as invalid.  

    Note that the children are ordered lexicographically according to their edge
    labels.  

* ``toChild(symbol)``  
    
    Set this ``EdgeNode`` to the child ``EdgeNode`` leading away whose label
    begins with the given ``symbol``.  

    If no such ``EdgeNode`` exists, mark this ``EdgeNode`` as invalid instead.  
";

%feature("docstring") stree::EdgeNode::depth "
``depth() -> nidx_t``  

Return the \"depth\" of the node in the suffix tree, which is the size of the
represented (sub-)sequence.  
";

%feature("docstring") stree::EdgeNode::sequence "
``sequence() -> Sequence``  

Return the (sub-)sequence represented by this node.  

Note that this is ``seq.rawSub(headIndex(), depth())``, where ``seq`` is the
sequence represented by the suffix tree.  
";

%feature("docstring") stree::EdgeNode::suffix "
``suffix() -> EdgeNode``  

Return the ``EdgeNode`` corresponding to the first suffix of the represented
sequence.  

This follows the \"suffix link\" of the suffix tree (and finds the new parent
accordingly). If no such node exists, a ``Node`` marked as invalid is returned.  
";

%feature("docstring") stree::EdgeNode::headIndex "
``headIndex() -> nidx_t``  

Return the \"headindex\" of this node, which is an index in the sequence
represented by the suffix tree where the (sub-)sequence represented by this node
occurs.  

I.e., the (sub-)sequence represented by this node is ``seq.rawSub(headindex(),
depth())``, where ``seq`` is the sequence represented by the suffix tree.  
";

%feature("docstring") stree::EdgeNode::isInternal "
``isInternal() -> bool``  

Return ``true`` if this is an internal node.  
";

%feature("docstring") stree::EdgeNode::setValid "
``setValid(valid=true)``  

Mark this ``Node`` as ``valid`` (or ``invalid``, if ``valid`` is ``false``).  
";

%feature("docstring") stree::EdgeNode::index "
``index() -> nidx_t``  

The ``index`` of a valid leaf or a valid internal node is a unique number
between 0 and ``STree.nLeafNodes()`` or between 0 and
``STree.nInternalNodes()``, respectively.  
";

%feature("docstring") stree::EdgeNode::nidx "
``nidx() -> nidx_t``  

Return the ``nidx_t`` corresponding to this ``Node``.  
";

%feature("docstring") stree::EdgeNode::nidxStr "
``nidxStr(width=3) -> std::string``  

Return a string representation of the underlying ``nidx_t``.  
";

%feature("docstring") stree::EdgeNode::child "
``child() -> EdgeNode``  
``child(symbol) -> EdgeNode``  

Overloaded function
-------------------
* ``child() -> EdgeNode``  
    
    Return the first child ``EdgeNode`` of this ``EdgeNode``.  

    If no child exists, an ``EdgeNode`` marked as invalid is returned. Note that
    the children are ordered lexicographically according to their edge labels.  

* ``child(symbol) -> EdgeNode``  
    
    Return the child ``EdgeNode`` leading away whose label begins with the given
    ``symbol``.  

    If no such ``EdgeNode`` exists, an ``EdgeNode`` marked as invalid is
    returned.  
";

%feature("docstring") stree::EdgeNode::isRoot "
``isRoot() -> bool``  

Return ``true`` if this is the root node.  
";

%feature("docstring") stree::EdgeNode::dataStr "
``dataStr(width=5) -> std::string``  

Return a string representation of the data of this ``EdgeNode``.  

This is useful for debugging or understanding the suffix tree structure.  
";

%feature("docstring") stree::EdgeNode::repr "
``repr() -> std::string``  

Return a string representation to display in python.  
";

%feature("docstring") stree::EdgeNode::set "
``set(edgeNode)``  
``set(node)``  

Overloaded function
-------------------
* ``set(edgeNode)``  
    
    Set this ``EdgeNode`` to the given ``edgeNode``, which must belong to the
    same suffix tree.  

    This is a faster version of ``(*this) = edgeNode)``.  

* ``set(node)``  
    
    Set this ``Node`` to the given ``node``, which must belong to the same
    suffix tree.  

    This is a faster version of ``(*this) = node)``.  
";

%feature("docstring") stree::EdgeNode::count "
``count() -> nidx_t``  

Return the number of occurrences of the sequence represented by this node in the
sequence represented by the suffix tree.  

For an invalid node, zero is returned.  
";

%feature("docstring") stree::EdgeNode::EdgeNode "
``EdgeNode(stree, nidx=ROOT, parent=ROOT)``  
``EdgeNode(node, parent=ROOT)``  
``EdgeNode(node, parent)``  

Overloaded function
-------------------
* ``EdgeNode(stree, nidx=ROOT, parent=ROOT)``  
    
    Construct an ``EdgeNode`` for the given ``stree`` corresponding to the given
    ``nidx`` and parent information in ``parent``.  

    If no ``nidx`` is given, it defaults to the root of the suffix tree. This
    will search for the parent starting form the given ``parent`` (or root by
    default). If no parent can be found (e.g. for the root), this ``EdgeNode``
    will be \"degenerate\" (have an invalid parent).  

* ``EdgeNode(node, parent=ROOT)``  
    
    Construct an ``EdgeNode`` from the given ``Node``.  

    This will search for the parent starting form the given ``parent`` (or root
    by default). If no parent can be found (e.g. for the root), this
    ``EdgeNode`` will be \"degenerate\" (have an invalid parent).  

* ``EdgeNode(node, parent)``  
    
    Construct an ``EdgeNode`` from the given ``Node``.  

    This will search for the parent starting form the given ``parent``. If no
    parent can be found (e.g. for the root), this ``EdgeNode`` will be
    \"degenerate\" (have an invalid parent).  
";

%feature("docstring") stree::EdgeNode::setRoot "
``setRoot()``  

Reset this ``EdgeNode`` to the root of the suffix tree.  
";

%feature("docstring") stree::EdgeNode::toSuffix "
``toSuffix()``  

Set this to the ``EdgeNode`` corresponding to the first suffix of the
represented sequence.  

If no such node exists, mark this node as invalid instead. This follows the
\"suffix link\" of the suffix tree (and finds the new parent accordingly).  
";

%feature("docstring") stree::EdgeNode::isSuffix "
``isSuffix() -> bool``  

Return ``true`` if the subsequence represented by this node is a suffix of the
underlying sequence.  

Note that this does not imply that this is a leaf.  
";

%feature("docstring") stree::EdgeNode::sibling "
``sibling() -> EdgeNode``  

Return the next sibling of this ``EdgeNode``.  

If none exists, an ``EdgeNode`` marked as invalid is returned. Note that the
siblings are ordered lexicographically according to their edge labels.  
";

%feature("docstring") stree::EdgeNode::isValid "
``isValid() -> bool``  

Return ``true`` if valid, otherwise return ``false``.  
";

%feature("docstring") stree::EdgeNode::label "
``label() -> Sequence``  

Return the edge label.  
";

%feature("docstring") stree::EdgeNode::parent "
``parent() -> Node``  

Return the parent ``Node`` of this ``EdgeNode``.  

If this ``EdgeNode`` is invalid or degenerate, i.e., no parent node exists or is
known, a ``Node`` marked as invalid is returned.  
";

%feature("docstring") stree::EdgeNode::isLeaf "
``isLeaf() -> bool``  

Return ``true`` if this is a leaf node.  
";

// File: classtom_1_1_estimator.xml


%feature("docstring") tom::Estimator "
``Estimator(stree)``  

This class computes estimates for :math:`f( x )` and corresponding variance
estimates for sequences :math:`x` based on a suffix tree representation of a
sample sequence.  

Constructors
------------
* ``Estimator(stree)``  
    
    Create an ``Estimator`` for a sample sequence data given by a suffix tree
    representation ``stree``.  

C++ includes: Estimator.h
";

%feature("docstring") tom::Estimator::fv "
``fv(z) -> tuple< double, double >``  
``fv(o, u) -> tuple< double, double >``  
``fv(sequence) -> tuple< double, double >``  
``fv(Y, X) -> tuple< MatrixXd, MatrixXd >``  
``fv(Y, z, X) -> tuple< MatrixXd, MatrixXd >``  
``fv(Y, o, u, X) -> tuple< MatrixXd, MatrixXd >``  
``fv(Y, s, X) -> tuple< MatrixXd, MatrixXd >``  

Overloaded function
-------------------
* ``fv(z) -> tuple< double, double >``  
    
    Return in a tuple (``f``, ``v``) an estimate ``f`` of f( ``z`` ) for the
    given output symbol ``z`` together with the corresponding variance estimate
    ``v``.  

* ``fv(o, u) -> tuple< double, double >``  
    
    Return in a tuple (``f``, ``v``) an estimate ``f`` of f( z ) for the given
    input-output symbol pair z = (``u``, ``o``) together with the corresponding
    variance estimate ``v``.  

    In the case of an output-only system, the input ``u`` is simply ignored.  

* ``fv(sequence) -> tuple< double, double >``  
    
    Return in a tuple (``f``, ``v``) an estimate of f( ``sequence`` ) together
    with the corresponding variance estimate ``v``.  

* ``fv(Y, X) -> tuple< MatrixXd, MatrixXd >``  
    
    Return in a tuple (``F``, ``V``) the matrix ``F`` of estimates for :math:`[
    f( x y ) ]_{y \\in X, x \\in X}` with rows indexed by the given set ``Y`` of
    characteristic sequences and columns indexed by the given set ``X`` of
    indicative sequences, together with the corresponding matrix ``V`` of
    element-wise variance estimates for the estimates returned in ``F``.  

* ``fv(Y, z, X) -> tuple< MatrixXd, MatrixXd >``  
    
    Return in a tuple (``F``, ``V``) the matrix ``F`` of estimates for :math:`[
    f( x z y ) ]_{y \\in X, x \\in X}` with rows indexed by the given set ``Y``
    of characteristic sequences and columns indexed by the given set ``X`` of
    indicative sequences for a given output symbol ``z``, together with the
    corresponding matrix ``V`` of element-wise variance estimates for the
    estimates returned in ``F``.  

* ``fv(Y, o, u, X) -> tuple< MatrixXd, MatrixXd >``  
    
    Return in a tuple (``F``, ``V``) the matrix ``F`` of estimates for :math:`[
    f( x z y ) ]_{y \\in X, x \\in X}` with rows indexed by the given set ``Y``
    of characteristic sequences and columns indexed by the given set ``X`` of
    indicative sequences for a given input-output symbol pair z = (``u``,
    ``o``), together with the corresponding matrix ``V`` of element-wise
    variance estimates for the estimates returned in ``F``.  

    In the case of an output-only system, the input ``u`` is simply ignored.  

* ``fv(Y, s, X) -> tuple< MatrixXd, MatrixXd >``  
    
    Return in a tuple (``F``, ``V``) the matrix ``F`` of estimates for :math:`[
    f( x s y ) ]_{y \\in X, x \\in X}` with rows indexed by the given set ``Y``
    of characteristic sequences and columns indexed by the given set ``X`` of
    indicative sequences for a given ``Sequence`` ``s``, together with the
    corresponding matrix ``V`` of element-wise variance estimates for the
    estimates returned in ``F``.  
";

%feature("docstring") tom::Estimator::sequence "
``sequence() -> Sequence``  

Return the data sequence.  
";

%feature("docstring") tom::Estimator::regularization "
``regularization(vPC=-1, vMin=-1, preset=\"\") -> tuple``  

Set (optional) and then return the regularization parameters for the
``Estimator``.  

This is done as follows:  

*   First, if a ``preset`` (\"none\" / \"default\") is specified, all
    regularization parameters are set accordingly.  
*   Next, any non-default argument causes the corresponding regularization
    parameter to be set to the given value, while any argument left at its
    default value (-1) has no effect.  
*   Finally, the current regularization parameters are returned in a tuple in
    the same order as they appear as function arguments. This allows writing
    python code such as:  

        old_params = estimator.regularization()
        ...
        estimator.regularization(*old_params)  

Parameters
----------
* ``vPC`` :  
    number of pseudo-counts to add for the variance computation  
* ``vMin`` :  
    determines a lower bound for the returned variance, which will simply be
    ``vMin`` if this is < 0.25, or computed as :math:`vMin / (\\pi N)^2`, where
    ``N`` is the length of the sequence data.  
* ``preset`` :  
    a preset that may be specified as explained above { \"default\", \"none\" }  
";

%feature("docstring") tom::Estimator::v "
``v(z) -> double``  
``v(o, u) -> double``  
``v(sequence) -> double``  
``v(Y, X) -> MatrixXd``  
``v(Y, z, X) -> MatrixXd``  
``v(Y, o, u, X) -> MatrixXd``  
``v(Y, s, X) -> MatrixXd``  

Overloaded function
-------------------
* ``v(z) -> double``  
    
    Return a variance estimate for the estimate of f( ``z`` ) for the given
    output symbol ``z``.  

* ``v(o, u) -> double``  
    
    Return a variance estimate for the estimate of f( z ) for the given input-
    output symbol pair z = (``u``, ``o``).  

    In the case of an output-only system, the input ``u`` is simply ignored.  

* ``v(sequence) -> double``  
    
    Return a variance estimate for the estimate of f( ``sequence`` ).  

* ``v(Y, X) -> MatrixXd``  
    
    Return the matrix of element-wise variance estimates corresponding to the
    estimates for :math:`[ f( x y ) ]_{y \\in X, x \\in X}` with rows indexed by
    the given set ``Y`` of characteristic sequences and columns indexed by the
    given set ``X`` of indicative sequences.  

* ``v(Y, z, X) -> MatrixXd``  
    
    Return the matrix of element-wise variance estimates corresponding to the
    estimates for :math:`[ f( x z y ) ]_{y \\in X, x \\in X}` with rows indexed
    by the given set ``Y`` of characteristic sequences and columns indexed by
    the given set ``X`` of indicative sequences for a given output symbol ``z``.  

* ``v(Y, o, u, X) -> MatrixXd``  
    
    Return the matrix of element-wise variance estimates corresponding to the
    estimates for :math:`[ f( x z y ) ]_{y \\in X, x \\in X}` with rows indexed
    by the given set ``Y`` of characteristic sequences and columns indexed by
    the given set ``X`` of indicative sequences for a given input-output symbol
    pair z = (``u``, ``o``).  

    In the case of an output-only system, the input ``u`` is simply ignored.  

* ``v(Y, s, X) -> MatrixXd``  
    
    Return the matrix of element-wise variance estimates corresponding to the
    estimates for :math:`[ f( x s y ) ]_{y \\in X, x \\in X}` with rows indexed
    by the given set ``Y`` of characteristic sequences and columns indexed by
    the given set ``X`` of indicative sequences for a given ``Sequence`` ``s``.  
";

%feature("docstring") tom::Estimator::nInputSymbols "
``nInputSymbols() -> int``  

Return the size of the input alphabet.  
";

%feature("docstring") tom::Estimator::nOutputSymbols "
``nOutputSymbols() -> int``  

Return the size of the output alphabet.  
";

%feature("docstring") tom::Estimator::Estimator "
``Estimator(stree)``  

Create an ``Estimator`` for a sample sequence data given by a suffix tree
representation ``stree``.  
";

%feature("docstring") tom::Estimator::f "
``f(z) -> double``  
``f(o, u) -> double``  
``f(sequence) -> double``  
``f(Y, X) -> MatrixXd``  
``f(Y, z, X) -> MatrixXd``  
``f(Y, o, u, X) -> MatrixXd``  
``f(Y, s, X) -> MatrixXd``  

Overloaded function
-------------------
* ``f(z) -> double``  
    
    Return an estimate of f( ``z`` ) for the given output symbol ``z``.  

* ``f(o, u) -> double``  
    
    Return an estimate of f( z ) for the given input-output symbol pair z =
    (``u``, ``o``).  

    In the case of an output-only system, the input ``u`` is simply ignored.  

* ``f(sequence) -> double``  
    
    Return an estimate of f( ``sequence`` ).  

* ``f(Y, X) -> MatrixXd``  
    
    Return the matrix of estimates for :math:`[ f( x y ) ]_{y \\in X, x \\in X}`
    with rows indexed by the given set ``Y`` of characteristic sequences and
    columns indexed by the given set ``X`` of indicative sequences.  

* ``f(Y, z, X) -> MatrixXd``  
    
    Return the matrix of estimates for :math:`[ f( x z y ) ]_{y \\in X, x \\in
    X}` with rows indexed by the given set ``Y`` of characteristic sequences and
    columns indexed by the given set ``X`` of indicative sequences for a given
    output symbol ``z``.  

* ``f(Y, o, u, X) -> MatrixXd``  
    
    Return the matrix of estimates for :math:`[ f( x z y ) ]_{y \\in X, x \\in
    X}` with rows indexed by the given set ``Y`` of characteristic sequences and
    columns indexed by the given set ``X`` of indicative sequences for a given
    input-output symbol pair z = (``u``, ``o``).  

    In the case of an output-only system, the input ``u`` is simply ignored.  

* ``f(Y, s, X) -> MatrixXd``  
    
    Return the matrix of estimates for :math:`[ f( x s y ) ]_{y \\in X, x \\in
    X}` with rows indexed by the given set ``Y`` of characteristic sequences and
    columns indexed by the given set ``X`` of indicative sequences for a given
    ``Sequence`` ``s``.  
";

// File: classtom_1_1_hmm.xml


%feature("docstring") tom::Hmm "
``Hmm(nStates, nObservations, nInputs=0, exponent=1, rnd=Random())``  

This class provides provides a rudimentry structure for HMMs and POMDPs.  

It purpose is to  

*   create randomly initialized HMMs or POMDPs  
*   learn HMMs / POMDPs from data using EM (Baum-Welch) Further operations are
    available after conversion into an ``Oom``.  

Constructors
------------
* ``Hmm(nStates, nObservations, nInputs=0, exponent=1, rnd=Random())``  

C++ includes: Hmm.h
";

/*
 Accessors 
*/

%feature("docstring") tom::Hmm::toJSON "
``toJSON() -> std::string``  
";

%feature("docstring") tom::Hmm::fromJSON "
``fromJSON(string)``  
";

%feature("docstring") tom::Hmm::pi "
``pi() -> const VectorXd &``  
``pi(_pi)``  

Overloaded function
-------------------
* ``pi() -> const VectorXd &``  

* ``pi(_pi)``  
";

%feature("docstring") tom::Hmm::normalize "
``normalize() -> bool``  
";

%feature("docstring") tom::Hmm::T "
``T(a=0) -> const MatrixXd &``  
``T(_T)``  
``T(a, _Ta)``  

Overloaded function
-------------------
* ``T(a=0) -> const MatrixXd &``  

* ``T(_T)``  

* ``T(a, _Ta)``  
";

%feature("docstring") tom::Hmm::E "
``E(o, a=0) -> const VectorXd &``  
``E(o, a, _Eoa)``  
``E(o, _Eo)``  

Overloaded function
-------------------
* ``E(o, a=0) -> const VectorXd &``  

* ``E(o, a, _Eoa)``  

* ``E(o, _Eo)``  
";

%feature("docstring") tom::Hmm::nInputs "
``nInputs() -> int``  
";

%feature("docstring") tom::Hmm::init "
``init()``  
";

%feature("docstring") tom::Hmm::Hmm "
``Hmm(nStates, nObservations, nInputs=0, exponent=1, rnd=Random())``  
";

%feature("docstring") tom::Hmm::setSize "
``setSize(nStates, nObservations, nInputs, zeroParameters=false)``  
";

%feature("docstring") tom::Hmm::cereal::access "
``cereal::access() -> friend class``  
";

%feature("docstring") tom::Hmm::trainEM "
``trainEM(trainSequence, stopCondition=StopCondition(100, 1e-7, 0)) -> double``  
";

%feature("docstring") tom::Hmm::repr "
``repr() -> std::string``  

return a representation to display in interactive python.  
";

%feature("docstring") tom::Hmm::nStates "
``nStates() -> int``  
";

%feature("docstring") tom::Hmm::randomize "
``randomize(exponent=1, rnd=Random())``  
";

%feature("docstring") tom::Hmm::nObservations "
``nObservations() -> int``  
";

// File: structstree_1_1internal_1_1_internal_node.xml


%feature("docstring") stree::internal::InternalNode "

An ``InternalNode`` consists of a left, right and child ``nidx_t``.  

Attributes
----------
* ``l_`` : ``nidx_t``  

* ``r_`` : ``nidx_t``  

* ``c_`` : ``nidx_t``  

C++ includes: STreeCore.h
";

// File: structstree_1_1internal_1_1_leaf_node.xml


%feature("docstring") stree::internal::LeafNode "

A ``LeafNode`` consists of a left and right ``nidx_t``.  

Attributes
----------
* ``l_`` : ``nidx_t``  
    the left ``nidx_t``.  

* ``r_`` : ``nidx_t``  
    the right ``nidx_t``.  

C++ includes: STreeCore.h
";

// File: classtom_1_1_learner.xml


%feature("docstring") tom::Learner "
``Learner()``  

Constructors
------------
* ``Learner()``  

Attributes
----------
* ``nU_`` : ``int``  
    the size of the input alphabet  

* ``nO_`` : ``int``  
    the size of the output alphabet  

* ``dim_`` : ``int``  
    the target dimension of the *Oom* model to be learnt  

* ``seqLen_`` : ``unsigned long``  
    the length of the sample sequence to use for learning  

* ``weightingScheme_`` : ``int``  

* ``algorithm_`` : ``int``  
    which algorithm to use to estimate the *Oom* model.  

    Currently, the only choices are:  

    *   SVD_ALGO: Spectral learning
        -   EC_ALGO: Same as spectral learning, but the SVD is computed via EM  

* ``cache_`` : ``int``  
    determines which values in the computation to cache  

* ``ec_err_threshold_`` : ``double``  

* ``ec_max_iterations_`` : ``int``  

* ``ec_iterations_`` : ``int``  

* ``ec_err_`` : ``double``  

* ``trainSequence_`` : ``Sequence``  
    the sample sequence from which to learn the *Oom* model.  

* ``suffixTree_`` : ``stree::STree``  
    a suffix tree representation of the *trainSequence_* or a subsequence  

* ``estimator_`` : ``Estimator``  
    the *Estimator* to use to obtain estimates from the sample sequence.  

* ``characteristicSequences_`` : ``SHARED_PTR< Sequences >``  
    the set of characteristic sequences  

* ``indicativeSequences_`` : ``SHARED_PTR< Sequences >``  
    the set of indicative sequences  

C++ includes: Learner.h
";

/*
 Accessors 
*/

/*
 Basic OOM functionality 
*/

%feature("docstring") tom::Learner::weightedSpectral "
``weightedSpectral() -> Oom *``  
";

%feature("docstring") tom::Learner::indicativeSequences "
``indicativeSequences() -> SHARED_PTR< Sequences >``  
``indicativeSequences(indicativeSequences_new)``  

Overloaded function
-------------------
* ``indicativeSequences() -> SHARED_PTR< Sequences >``  

* ``indicativeSequences(indicativeSequences_new)``  
";

%feature("docstring") tom::Learner::F "
``F() -> Eigen::MatrixXd &``  
``F(F_new)``  

Overloaded function
-------------------
* ``F() -> Eigen::MatrixXd &``  

* ``F(F_new)``  
";

%feature("docstring") tom::Learner::W "
``W() -> Eigen::MatrixXd &``  
``W(W_new)``  

Overloaded function
-------------------
* ``W() -> Eigen::MatrixXd &``  

* ``W(W_new)``  
";

%feature("docstring") tom::Learner::WI "
``WI() -> Eigen::MatrixXd &``  
``WI(WI_new)``  

Overloaded function
-------------------
* ``WI() -> Eigen::MatrixXd &``  

* ``WI(WI_new)``  
";

%feature("docstring") tom::Learner::C "
``C() -> Eigen::MatrixXd &``  
``C(C_new)``  

Overloaded function
-------------------
* ``C() -> Eigen::MatrixXd &``  

* ``C(C_new)``  
";

%feature("docstring") tom::Learner::init "
``init()``  
";

%feature("docstring") tom::Learner::characteristicSequences "
``characteristicSequences() -> SHARED_PTR< Sequences >``  
``characteristicSequences(characteristicSequences_new)``  

Overloaded function
-------------------
* ``characteristicSequences() -> SHARED_PTR< Sequences >``  

* ``characteristicSequences(characteristicSequences_new)``  
";

%feature("docstring") tom::Learner::Wz "
``Wz(o, u=0) -> Eigen::MatrixXd &``  
``Wz(o, Wz_new)``  
``Wz(o, u, Wz_new)``  

Overloaded function
-------------------
* ``Wz(o, u=0) -> Eigen::MatrixXd &``  

* ``Wz(o, Wz_new)``  

* ``Wz(o, u, Wz_new)``  
";

%feature("docstring") tom::Learner::computeCQ "
``computeCQ()``  
";

%feature("docstring") tom::Learner::oom "
``oom() -> Oom *``  
";

%feature("docstring") tom::Learner::Learner "
``Learner()``  
";

%feature("docstring") tom::Learner::~Learner "
``~Learner()``  
";

%feature("docstring") tom::Learner::FJ "
``FJ() -> Eigen::MatrixXd &``  
``FJ(FJ_new)``  

Overloaded function
-------------------
* ``FJ() -> Eigen::MatrixXd &``  

* ``FJ(FJ_new)``  
";

%feature("docstring") tom::Learner::FI "
``FI() -> Eigen::MatrixXd &``  
``FI(FI_new)``  

Overloaded function
-------------------
* ``FI() -> Eigen::MatrixXd &``  

* ``FI(FI_new)``  
";

%feature("docstring") tom::Learner::Fz "
``Fz(o, u=0) -> Eigen::MatrixXd &``  
``Fz(o, Fz_new)``  
``Fz(o, u, Fz_new)``  

Overloaded function
-------------------
* ``Fz(o, u=0) -> Eigen::MatrixXd &``  

* ``Fz(o, Fz_new)``  

* ``Fz(o, u, Fz_new)``  
";

%feature("docstring") tom::Learner::WJ "
``WJ() -> Eigen::MatrixXd &``  
``WJ(WJ_new)``  

Overloaded function
-------------------
* ``WJ() -> Eigen::MatrixXd &``  

* ``WJ(WJ_new)``  
";

%feature("docstring") tom::Learner::Q "
``Q() -> Eigen::MatrixXd &``  
``Q(Q_new)``  

Overloaded function
-------------------
* ``Q() -> Eigen::MatrixXd &``  

* ``Q(Q_new)``  
";

%feature("docstring") tom::Learner::clearCQ "
``clearCQ()``  
";

// File: classstree_1_1_node.xml


%feature("docstring") stree::Node "
``Node(stree, nidx=ROOT)``  

This class represents a node in the suffix tree and can be used for extracting
information or navigating the suffix tree.  

It contains a pointer internally to the ``STree`` that it belongs to.  

The ``to...()`` methods are provided to navigate the suffix tree structure:  

*   ``to...()`` sets this node to its e.g. child, suffix, sibling, etc.  
*   If no such exists, then this node is simply marked as invalid.  
*   For an invalid node, the ``to...()`` methods have no effect.  
*   Calling ``setValid()`` after the node has been marked as invalid by a
    ``to...()`` method will reset this node to the last valid node during the
    traversal.  

Constructors
------------
* ``Node(stree, nidx=ROOT)``  
    
    Construct a ``Node`` for the given ``stree`` corresponding to the given
    ``nidx``.  

    If no ``nidx`` is given, it defaults to the root of the suffix tree.  

C++ includes: STreeNode.h
";

%feature("docstring") stree::Node::depth "
``depth() -> nidx_t``  

Return the \"depth\" of the node in the suffix tree, which is the size of the
represented (sub-)sequence.  
";

%feature("docstring") stree::Node::count "
``count() -> nidx_t``  

Return the number of occurrences of the sequence represented by this node in the
sequence represented by the suffix tree.  

For an invalid node, zero is returned.  
";

%feature("docstring") stree::Node::isInternal "
``isInternal() -> bool``  

Return ``true`` if this is an internal node.  
";

%feature("docstring") stree::Node::sequence "
``sequence() -> Sequence``  

Return the (sub-)sequence represented by this node.  

Note that this is ``seq.rawSub(headIndex(), depth())``, where ``seq`` is the
sequence represented by the suffix tree.  
";

%feature("docstring") stree::Node::headIndex "
``headIndex() -> nidx_t``  

Return the \"headindex\" of this node, which is an index in the sequence
represented by the suffix tree where the (sub-)sequence represented by this node
occurs.  

I.e., the (sub-)sequence represented by this node is ``seq.rawSub(headindex(),
depth())``, where ``seq`` is the sequence represented by the suffix tree.  
";

%feature("docstring") stree::Node::setRoot "
``setRoot()``  

Reset this ``Node`` to the root of the suffix tree.  
";

%feature("docstring") stree::Node::toSuffix "
``toSuffix()``  

Set this node to the node corresponding to the first suffix of the represented
sequence.  

If no such node exists, mark this node as invalid instead. This follows the
\"suffix link\" of the suffix tree.  
";

%feature("docstring") stree::Node::isLeaf "
``isLeaf() -> bool``  

Return ``true`` if this is a leaf node.  
";

%feature("docstring") stree::Node::toSibling "
``toSibling()``  

Set this node to its next sibling if a next sibling exists, otherwise mark this
node as invalid.  

Note that the siblings are ordered lexicographically according to their edge
labels.  
";

%feature("docstring") stree::Node::repr "
``repr() -> std::string``  

Return a string representation to display in python.  
";

%feature("docstring") stree::Node::suffix "
``suffix() -> Node``  

Return the node corresponding to the first suffix of the represented sequence.  

This follows the \"suffix link\" of the suffix tree. If no such node exists, a
``Node`` marked as invalid is returned.  
";

%feature("docstring") stree::Node::dataStr "
``dataStr(width=5) -> std::string``  

Return a string representation of the data of this node.  

This is useful for debugging or understanding the suffix tree structure.  
";

%feature("docstring") stree::Node::isValid "
``isValid() -> bool``  

Return ``true`` if valid, otherwise return ``false``.  
";

%feature("docstring") stree::Node::isSuffix "
``isSuffix() -> bool``  

Return ``true`` if the subsequence represented by this node is a suffix of the
underlying sequence.  

Note that this does not imply that this is a leaf.  
";

%feature("docstring") stree::Node::toChild "
``toChild()``  
``toChild(symbol)``  

Overloaded function
-------------------
* ``toChild()``  
    
    Set this node to its first child if such a node exists, otherwise mark this
    node as invalid.  

    Note that the children are ordered lexicographically according to their edge
    labels.  

* ``toChild(symbol)``  
    
    Set this node to the child node along the edge leading away whose label
    begins with the given ``symbol``.  

    If no such node exists, mark this node as invalid instead.  
";

%feature("docstring") stree::Node::setValid "
``setValid(valid=true)``  

Mark this ``Node`` as ``valid`` (or ``invalid``, if ``valid`` is ``false``).  
";

%feature("docstring") stree::Node::index "
``index() -> nidx_t``  

The ``index`` of a valid leaf or a valid internal node is a unique number
between 0 and ``STree.nLeafNodes()`` or between 0 and
``STree.nInternalNodes()``, respectively.  
";

%feature("docstring") stree::Node::nidxStr "
``nidxStr(width=3) -> std::string``  

Return a string representation of the underlying ``nidx_t``.  
";

%feature("docstring") stree::Node::set "
``set(node)``  

Set this ``Node`` to the given ``node``, which must belong to the same suffix
tree.  

This is a faster version of ``(*this) = node)``.  
";

%feature("docstring") stree::Node::isRoot "
``isRoot() -> bool``  

Return ``true`` if this is the root node.  
";

%feature("docstring") stree::Node::nidx "
``nidx() -> nidx_t``  

Return the ``nidx_t`` corresponding to this ``Node``.  
";

%feature("docstring") stree::Node::Node "
``Node(stree, nidx=ROOT)``  

Construct a ``Node`` for the given ``stree`` corresponding to the given
``nidx``.  

If no ``nidx`` is given, it defaults to the root of the suffix tree.  
";

%feature("docstring") stree::Node::sibling "
``sibling() -> Node``  

Return the next sibling of this node.  

If no such node exists, a ``Node`` marked as invalid is returned. Note that the
siblings are ordered lexicographically according to their edge labels.  
";

%feature("docstring") stree::Node::child "
``child() -> Node``  
``child(symbol) -> Node``  

Overloaded function
-------------------
* ``child() -> Node``  
    
    Return the first child node of this node.  

    If no such node exists, a ``Node`` marked as invalid is returned. Note that
    the children are ordered lexicographically according to their edge labels.  

* ``child(symbol) -> Node``  
    
    Return the child node along the edge leading away whose label begins with
    the given ``symbol``.  

    If no such node exists, a ``Node`` marked as invalid is returned.  
";

// File: structstree_1_1nodelete.xml


%feature("docstring") stree::nodelete "

Hack needed to convert a normal pointer to a ``std::shared_ptr`` safely.  

C++ includes: stree.h
";

// File: classtom_1_1_oom.xml


%feature("docstring") tom::Oom "
``Oom()``  
``Oom(dimension, nOutputSymbols, nInputSymbols=0, randomExponent=0, zero_threshold=0, randomSource=Random())``  
``Oom(json_representation)``  
``Oom(hmm)``  

This class provides the basic functionality of observable operator models.  

The main aspects are  

*   extracting information from a given OOM (sequence generation, prediction,
    fingerprints, etc), and  
*   transformation functions (stabilization, equivalence transforms, dimension
    reduction, etc.) To estimate (\"learn\") an OOM from data, use the OomTrain
    class.  

Constructors
------------
* ``Oom()``  
    
    Construct an uninitialized (!) Oom.  

* ``Oom(dimension, nOutputSymbols, nInputSymbols=0, randomExponent=0, zero_threshold=0, randomSource=Random())``  
    
    Construct a simple (random) OOM of dimension ``dimension`` that models a
    stochastic process with ``nOutputSymbols`` number of possible observations
    and ``nInputSymbols`` number of possible inputs.  

    The initial state will be set to the stationary state, assuming iid inputs
    in the case of an input-output OOM.  

    **Parameters:**  

    * ``dimension`` :  
        the dimension of the OOM  
    * ``nOutputSymbols`` :  
        the size of the output alphabet  
    * ``nInputSymbols`` :  
        the size of the input alphabet, or 0 (default) for an output-only
        ``Oom``  
    * ``randomExponent`` :  
        exponent of value distribution of tau operator matrix entries. When set
        to 1, the entries of matrices are sampled from a uniform distribution.
        Higher values lead to sample distributions that are increasingly skewed
        to the right, while a value of zero (default) will lead to an ``Oom``
        that produces iid outputs.  
    * ``zero_threshold`` :  
        set all parameters less than ``zero_threshold`` to zero and renormalize.  
    * ``randomSource`` :  
        the ``Random`` object to use as the source of randomness  

* ``Oom(json_representation)``  
    
    Construct an ``Oom`` corresponding to the given string
    ``json_representation``.  

    The format must correspond to what ``toJSON()`` produces.  

* ``Oom(hmm)``  
    
    Construct an ``Oom`` equivalent to the ``Hmm`` given by ``hmm``.  

Attributes
----------
* ``fixPredictionMargin_`` : ``double``  
    the largest tolerated error of the prediction vector before the
    normalization is considered a fixing event  

* ``nSetback_`` : ``int``  
    the number of times that the OOM state needed to be fixed by a setback
    operation (introduces a considerable error)  

* ``nFixPrediction_`` : ``int``  
    the number of times that the predicted probabilities were adjusted by more
    than ``fixPredictionMargin_``.  

* ``nImpossible_`` : ``int``  
    the number of times that a symbol was encountered that has probability
    smaller than ``impossibleProbMargin_`` according to the OOM  

C++ includes: Oom.h
";

/*
 Constructors and initialization 
*/

/*
 Main Oom parameters 
*/

/*
 Stabilization 
*/

/*
 Basic OOM functionality 
*/

/*
 Transformation functions 
*/

/*
 IO-functions 
*/

/*
 Internalals\\. Use only if you know what you are doing! 
*/

%feature("docstring") tom::Oom::setBack "
``setBack() -> bool``  

Attempt to perform a state setback operation for at most ``maxSetback_`` time-
steps. Note that calling this method repeatedly will attempt a setback for a
shorter history each time. Return ``true`` if a setback could be performed.  
";

%feature("docstring") tom::Oom::minPrediction "
``minPrediction() -> double``  
``minPrediction(new_value)``  

Overloaded function
-------------------
* ``minPrediction() -> double``  
    
    Return the minimum probability for any observation symbol.  

    The ``prediction()`` is normalized at each time-step such that every output
    symbol has at least this probability.  

* ``minPrediction(new_value)``  
    
    Set the minimum probability for any observation symbol to the given
    ``new_value``.  

    The ``prediction()`` is normalized at each time-step such that every output
    symbol has at least this probability.  
";

%feature("docstring") tom::Oom::stabilization "
``stabilization(minPrediction=-1, normalizationTolerance=-1, maxSetback=-1, impossibilityThreshold=-1, preset=\"\") -> tuple``  

Set (optional) and then return the stabilization parameters for the ``Oom``.  

This is done as follows:  

*   First, if a ``preset`` (\"none\" / \"default\") is specified, all
    stabilization parameters are set accordingly.  
*   Next, any non-default argument causes the corresponding stabilization
    parameter to be set to the given value, while any argument left at its
    default value (-1) has no effect.  
*   Finally, the current stabilization parameters are returned in a tuple in the
    same order as they appear as function arguments. This allows writing python
    code such as:  

        old_params = oom.stabilization()
        ...
        oom.stabilization(*old_params)  

Parameters
----------
* ``minPrediction`` :  
    The minimum probability for any observation symbol. The ``prediction()`` is
    normalized at each time-step such that every output symbol has at least this
    probability.  
* ``normalizationTolerance`` :  
    The maximum allowed tolerance between the unnormalized and normalized
    prediction vector before invoking a ``setBack()``. The tolerance is computed
    as 1.5 * ``nOutputSymbols()`` * squared norm of the difference. See also
    ``normalizePrediction()``.  
* ``maxSetback`` :  
    The maximum number of steps to \"replay\" during a ``setBack()`` operation.  
* ``impossibilityThreshold`` :  
    The smallest probability value considered as non-zero. If an observation
    symbol is encountered that this model predicts to have a probability lower
    than ``impossibilityThreshold()``, this is considered as an impossible event
    and a ``setBack()`` must be performed, as the state ``wt()`` can no longer
    be normalized.  
* ``preset`` :  
    A preset for the stabilization parameters to apply before setting the other
    parameters. This can be either \"none\" to disable stabilization, or
    \"default\" to use the default stabilization settings.  
";

%feature("docstring") tom::Oom::f "
``f(z, reset=true) -> double``  
``f(o, u, reset=true) -> double``  
``f(sequence, reset=true) -> double``  
``f(Y, X, reset=true) -> MatrixXd``  
``f(Y, z, X, reset=true) -> MatrixXd``  
``f(Y, o, u, X, reset=true) -> MatrixXd``  
``f(Y, s, X, reset=true) -> MatrixXd``  

Overloaded function
-------------------
* ``f(z, reset=true) -> double``  
    
    If ``reset`` is ``true`` (default), perform a state ``reset()`` first.  

    Then return the value for the given output symbol ``z`` of the prediction
    function for the state ``wt()``, i.e., the \"probability\" P( ``z`` |
    ``wt()`` ), and update the state.  

* ``f(o, u, reset=true) -> double``  
    
    If ``reset`` is ``true`` (default), perform a state ``reset()`` first.  

    Then return the value for the given input-output pair (``u``,``o``) of the
    prediction function for the state ``wt()``, i.e., the \"probability\" P(
    ``o`` | ``u``, ``wt()`` ), and update the state. In the case of an output-
    only ``Oom``, the input ``u`` is simply ignored.  

* ``f(sequence, reset=true) -> double``  
    
    If ``reset`` is ``true`` (default), perform a state ``reset()`` first.  

    Then return the value for the given ``sequence`` of the prediction function
    for the state ``wt()``, i.e., the \"probability\" P( ``sequence`` | ``wt()``
    ), and update the state.  

* ``f(Y, X, reset=true) -> MatrixXd``  
    
    If ``reset`` is ``true`` (default), perform a state ``reset()`` first.  

    Then return the matrix of prediction function values :math:`[ f(x y) ]_{y
    \\in Y, x \\in X}` with rows indexed by the given set ``Y`` of
    characteristic sequences and columns indexed by the given set ``X`` of
    indicative sequences.  

* ``f(Y, z, X, reset=true) -> MatrixXd``  
    
    If ``reset`` is ``true`` (default), perform a state ``reset()`` first.  

    Then return the matrix of prediction function values :math:`[ f(x z y) ]_{y
    \\in Y, x \\in X}` with rows indexed by the given set ``Y`` of
    characteristic sequences and columns indexed by the given set ``X`` of
    indicative sequences for a given output symbol ``z``.  

* ``f(Y, o, u, X, reset=true) -> MatrixXd``  
    
    If ``reset`` is ``true`` (default), perform a state ``reset()`` first.  

    Then return the matrix of prediction function values :math:`[ f(x z y) ]_{y
    \\in Y, x \\in X}` with rows indexed by the given set ``Y`` of
    characteristic sequences and columns indexed by the given set ``X`` of
    indicative sequences for a given input-output symbol pair z = (``u``,
    ``o``). In the case of an output-only ``Oom``, the input ``u`` is simply
    ignored.  

* ``f(Y, s, X, reset=true) -> MatrixXd``  
    
    If ``reset`` is ``true`` (default), perform a state ``reset()`` first.  

    Then return the matrix of prediction function values :math:`[ f(x s y) ]_{y
    \\in Y, x \\in X}` with rows indexed by the given set ``Y`` of
    characteristic sequences and columns indexed by the given set ``X`` of
    indicative sequences for a given ``Sequence`` ``s``.  
";

%feature("docstring") tom::Oom::averageOneStepPredictionError "
``averageOneStepPredictionError(sequence, trueModel) -> double``  

Return the average one-step squared prediction error computed along the given
sample ``sequence`` according to a correct ``Oom`` ``trueModel``, after first
performing a state ``reset()`` operation on both this ``Oom`` and the
``trueModel``.  

Their states are updated accordingly.  
";

%feature("docstring") tom::Oom::sig "
``sig() -> const RowVectorXd &``  
``sig(new_value)``  

Overloaded function
-------------------
* ``sig() -> const RowVectorXd &``  
    
    Return the evaluation functional row vector :math:`\\sigma`.  

* ``sig(new_value)``  
    
    Set the evaluation functional vector :math:`\\sigma` to the given row vector
    ``new_value``.  

    Note that you *must* call ``initialize()`` after re-setting the ``Oom``
    parameters ``sig`` or ``tau(o,u)``.  
";

%feature("docstring") tom::Oom::nInputSymbols "
``nInputSymbols() -> int``  

Return the size of the input alphabet.  

Use ``setSize()`` to modify.  
";

%feature("docstring") tom::Oom::tau "
``tau(o, u=0) -> const MatrixXd &``  
``tau(z) -> const MatrixXd &``  
``tau(z, new_value)``  
``tau(o, u, new_value)``  
``tau(o, new_value)``  

Overloaded function
-------------------
* ``tau(o, u=0) -> const MatrixXd &``  
    
    Return the observable operator corresponding to observation ``o`` and input
    ``u``.  

    The parameter ``u`` defaults to 0 for the case of no inputs.  

* ``tau(z) -> const MatrixXd &``  
    
    Return the observable operator corresponding to the symbol ``z`` given as a
    ``Sequence`` of length one.  

* ``tau(z, new_value)``  
    
    Set the observable operator corresponding to the symbol ``z`` given as a
    ``Sequence`` of length one to the given matrix ``new_value``.  

    Note that you *must* call ``initialize()`` after re-setting the ``Oom``
    parameters ``sig`` or ``tau(o,u)``.  

* ``tau(o, u, new_value)``  
    
    Set the observable operator corresponding to observation ``o`` and input
    ``u`` to the given matrix ``new_value``.  

    Note that you *must* call ``initialize()`` after re-setting the ``Oom``
    parameters ``sig`` or ``tau(o,u)``.  

* ``tau(o, new_value)``  
    
    Set the observable operator corresponding to the observation ``o`` to the
    given matrix ``new_value``  
";

%feature("docstring") tom::Oom::repr "
``repr() -> std::string``  

return a representation to display in interactive python.  
";

%feature("docstring") tom::Oom::reset "
``reset()``  

Reset the ``Oom`` to its initial state and ``resetStabilizationStatistics()``.  
";

%feature("docstring") tom::Oom::w0 "
``w0() -> const VectorXd &``  
``w0(new_value)``  

Overloaded function
-------------------
* ``w0() -> const VectorXd &``  
    
    Return the initial state vector :math:`\\omega_0`.  

* ``w0(new_value)``  
    
    Set the initial state vector :math:`\\omega_0` to the given vector
    ``new_value`` and perform a ``reset()``.  
";

%feature("docstring") tom::Oom::entropy "
``entropy(sample_length, randomSource=Random(), policy=Policy()) -> double``  

Return entropy of this ``Oom`` estimated on a sample sequence of the given
``sample_length``.  
";

%feature("docstring") tom::Oom::cereal::access "
``cereal::access() -> friend class``  
";

%feature("docstring") tom::Oom::maxSetback "
``maxSetback() -> int``  
``maxSetback(new_value)``  

Overloaded function
-------------------
* ``maxSetback() -> int``  
    
    Return the maximum number of steps to \"replay\" during a ``setBack()``
    operation.  

* ``maxSetback(new_value)``  
    
    Set the maximum number of steps to \"replay\" during a ``setBack()``
    operation to ``new_value``.  
";

%feature("docstring") tom::Oom::normalizationTolerance "
``normalizationTolerance() -> double``  
``normalizationTolerance(new_value)``  

Overloaded function
-------------------
* ``normalizationTolerance() -> double``  
    
    Return the maximum allowed tolerance between the unnormalized and normalized
    prediction vector before invoking a ``setBack()``.  

    The tolerance is computed as 1.5 * ``nOutputSymbols()`` * squared norm of
    the difference. See also ``normalizePrediction()``.  

* ``normalizationTolerance(new_value)``  
    
    Set the maximum allowed tolerance between the unnormalized and normalized
    prediction vector before invoking a ``setBack()`` to the given
    ``new_value``.  

    The tolerance is computed as 1.5 * ``nOutputSymbols()`` * squared norm of
    the difference. See also ``normalizePrediction()``.  
";

%feature("docstring") tom::Oom::setSize "
``setSize(dimension, nOutputSymbols, nInputSymbols=0)``  

Set the internal structure for an OOM of the desired size without performing any
initialization. Typically, the parameters ``sig``, ``tau(o,u)`` and ``w0`` will
be assigned next, and then ``initialize()`` must be called.  

Parameters
----------
* ``dimension`` :  
    the dimension of the OOM  
* ``nOutputSymbols`` :  
    the size of the output alphabet  
* ``nInputSymbols`` :  
    the size of the input alphabet, or 0 (default) for an output-only ``Oom``  
";

%feature("docstring") tom::Oom::nOutputSymbols "
``nOutputSymbols() -> int``  

Return the size of the output alphabet.  

Use ``setSize()`` to modify  
";

%feature("docstring") tom::Oom::normalizePrediction "
``normalizePrediction() -> double``  

Attempt to fix the prediction vector of the next output symbol probabilities
:math:`P(\\cdot|u_t, \\omega_t)` such that all probabilities are at least
``minPrediction_`` and the probabilities sum to one. Return a measure of the
required change to the prediction vector: 1.5 * nO() * squared norm of the
difference.  
";

%feature("docstring") tom::Oom::harvestStates "
``harvestStates(sequence, reset=true) -> MatrixXf``  

Return the ``float`` matrix of the ``sequence.length()`` states (in its columns)
occurring during the computation of ``f(sequence, reset)``.  
";

%feature("docstring") tom::Oom::l2l "
``l2l(sequence) -> double``  

Return the log2-likelihood of the ``Oom`` for the given ``sequence``.  

That is, return -``log2_f(sequence, true)`` / ``sequence.length()``.  
";

%feature("docstring") tom::Oom::crossEntropyOfKOrderMarkovApproximation "
``crossEntropyOfKOrderMarkovApproximation(k, sequence) -> double``  

Return the cross-entropy of the best ``k``-order Markov model approximation
estimated on the given sample ``sequence``.  
";

%feature("docstring") tom::Oom::initialize "
``initialize()``  

Initialize the OOM. This assumes that all essential parameters (i.e.,
``dimension``, ``nOutputSymbols``, ``nInputSymbols``, ``sig``, ``tau(o,u)`` and
``w0``) have been set.  
";

%feature("docstring") tom::Oom::sample "
``sample(length, randomSource=Random(), policy=Policy(), exponent=1) -> Sequence``  

Sample, using the given ``randomSource``, a sequence of given ``length`` from
the ``Oom``, in the case of inputs together with the given input ``policy`` (by
default iud inputs), starting from the **current** state ``wt()``.  

For each time-step an observation is sampled from the ``prediction()`` vector
raised element-wise to the power ``exponent`` (default 1). Higher ``exponent``
values add a bias towards the most likely sequences, while lower values bias
towards uniformly distributed sequences.  
";

%feature("docstring") tom::Oom::reverse "
``reverse(normalize=true) -> std::shared_ptr< Oom >``  

Return the \"reverse\" of this ``Oom``.  
";

%feature("docstring") tom::Oom::stationaryState "
``stationaryState(policy=Policy(), maxIterations=10000) -> Eigen::VectorXd``  

Return the stationary state (in the case of an input-output ``Oom`` according to
the given ``policy``), computed by the power method with at most
``maxIterations`` number of iterations.  
";

%feature("docstring") tom::Oom::transform "
``transform(sig, w0=VectorXd::Zero(0))``  

Transform this ``Oom`` to an equivalent ``Oom`` that has given ``sig`` and
``w0`` as parameters for ``sig()`` and ``w0()``.  

This will only yield an (equivalent) ``Oom`` if ``sig`` * ``w0`` = 1.  
";

%feature("docstring") tom::Oom::Oom "
``Oom()``  
``Oom(dimension, nOutputSymbols, nInputSymbols=0, randomExponent=0, zero_threshold=0, randomSource=Random())``  
``Oom(json_representation)``  
``Oom(hmm)``  

Overloaded function
-------------------
* ``Oom()``  
    
    Construct an uninitialized (!) Oom.  

* ``Oom(dimension, nOutputSymbols, nInputSymbols=0, randomExponent=0, zero_threshold=0, randomSource=Random())``  
    
    Construct a simple (random) OOM of dimension ``dimension`` that models a
    stochastic process with ``nOutputSymbols`` number of possible observations
    and ``nInputSymbols`` number of possible inputs.  

    The initial state will be set to the stationary state, assuming iid inputs
    in the case of an input-output OOM.  

    **Parameters:**  

    * ``dimension`` :  
        the dimension of the OOM  
    * ``nOutputSymbols`` :  
        the size of the output alphabet  
    * ``nInputSymbols`` :  
        the size of the input alphabet, or 0 (default) for an output-only
        ``Oom``  
    * ``randomExponent`` :  
        exponent of value distribution of tau operator matrix entries. When set
        to 1, the entries of matrices are sampled from a uniform distribution.
        Higher values lead to sample distributions that are increasingly skewed
        to the right, while a value of zero (default) will lead to an ``Oom``
        that produces iid outputs.  
    * ``zero_threshold`` :  
        set all parameters less than ``zero_threshold`` to zero and renormalize.  
    * ``randomSource`` :  
        the ``Random`` object to use as the source of randomness  

* ``Oom(json_representation)``  
    
    Construct an ``Oom`` corresponding to the given string
    ``json_representation``.  

    The format must correspond to what ``toJSON()`` produces.  

* ``Oom(hmm)``  
    
    Construct an ``Oom`` equivalent to the ``Hmm`` given by ``hmm``.  
";

%feature("docstring") tom::Oom::prediction "
``prediction() -> const VectorXd &``  

Return the current prediction vector of the next output symbol probabilities.  

In the case of an input-output ``Oom`` these probabilities depend on the current
input symbol u_t, so ``condition(u_t)`` *must* have been called first.  
";

%feature("docstring") tom::Oom::log2_f "
``log2_f(sequence, reset=true) -> double``  

If ``reset`` is ``true`` (default), perform a state ``reset()`` first.  

Then return the log_2 of the prediction function for the given ``sequence``
given the state ``wt()``, i.e., log_2 f( ``sequence`` | ``wt()`` ), and update
the state.  

To deal gracefully with observations that have a prediction value below
``impossibilityThreshold()`` at some time step but occur nevertheless in the
``sequence``, the prediction for this occurrence is treated as having a
probability of ``impossibilityThreshold()``. Every time this happens, the
counter ``nImpossible_`` is incremented. Note that this problem can be avoided
by increasing ``minPrediction()`` above ``impossibilityThreshold()``.  
";

%feature("docstring") tom::Oom::impossibilityThreshold "
``impossibilityThreshold() -> double``  
``impossibilityThreshold(new_value)``  

Overloaded function
-------------------
* ``impossibilityThreshold() -> double``  
    
    Return the smallest probability value considered as non-zero.  

    If an observation symbol is encountered that this model predicts to have a
    probability lower than ``impossibilityThreshold()``, this is considered as
    an impossible event and a ``setBack()`` must be performed, as the state
    ``wt()`` can no longer be normalized.  

* ``impossibilityThreshold(new_value)``  
    
    Set the smallest probability value considered as non-zero to the given
    ``new_value``.  

    If an observation symbol is encountered that this model predicts to have a
    probability lower than ``impossibilityThreshold()``, this is considered as
    an impossible event and a ``setBack()`` must be performed, as the state
    ``wt()`` can no longer be normalized.  
";

%feature("docstring") tom::Oom::isIO "
``isIO() -> bool``  

Return ``true`` if this is an input-output sequence, i.e., if the input alphabet
size ``nInputSymbols()`` is non-zero.  
";

%feature("docstring") tom::Oom::update "
``update(o, u=0)``  

Update the ``Oom`` state according to the input-output pair (``u``,``o``),
normalize, and in the case of an output-only ``Oom``, additionally call
``condition()``.  

That is, first set ``wt`` to ``tau(o,u)`` * ``wt()``, and then attempt to
normalize to ``wt()`` / ( ``sig()`` * ``wt()`` ) if ``sig()`` * ``wt()`` is
greater than the ``impossibilityThreshold()``, else perform ``setback()``
operations.  
";

%feature("docstring") tom::Oom::wt "
``wt() -> const VectorXd &``  
``wt(new_value, history=Sequence())``  

Overloaded function
-------------------
* ``wt() -> const VectorXd &``  
    
    Return the (normalized) current state vector :math:`\\omega_t`.  

* ``wt(new_value, history=Sequence())``  
    
    Set the current state to the (normalized) given vector ``new_value``, and
    (optionally) specify a given input-output ``history`` (relevant for
    stabilization).  

    In the case of an output-only ``Oom`` this automatically calls
    ``condition()``.  
";

%feature("docstring") tom::Oom::conjugate "
``conjugate(rho, rhoInv)``  

Conjugate this ``Oom`` by the given matrices ``rho`` and ``rhoInv``.  

That is, set:  

*   ``w0()`` = ``rho`` * ``w0()``  
*   ``tau(o,u)`` = ``rho`` * ``tau(o,u)`` * ``rhoInv``  
*   ``sig()`` = ``sig()`` * ``rhoInv``.  
";

%feature("docstring") tom::Oom::history "
``history() -> Sequence``  

Return the most recent input-output history that is relevant for stabilization
purposes.  
";

%feature("docstring") tom::Oom::condition "
``condition(u=0)``  

Compute and normalize the prediction vector of the next output symbol
probabilities according to the current state ``wt()`` and input ``u``, i.e.,
compute P( · | ``u``, ``wt()`` ).  

Note that this function calls ``setback()`` if the ``normalizationTolerance()``
is exceeded. Therefore, calling e.g., ``condition(0)`` after ``condition(1)``
may not give the same result as calling just ``condition(0)``.  
";

%feature("docstring") tom::Oom::dimension "
``dimension() -> int``  

Return the model dimension.  
";

// File: classstree_1_1_path_node.xml


%feature("docstring") stree::PathNode "
``PathNode(stree, nidx=ROOT)``  
``PathNode(node)``  

This class represents a node (together with the path of nodes leading to it) in
the suffix tree and can be used for extracting information or navigating the
suffix tree.  

The ``to...()`` methods are provided to navigate the suffix tree structure:  

*   ``to...()`` sets this node (with path leading to it) to its e.g. child,
    suffix, sibling, etc.  
*   If no such exists, then this node is simply marked as invalid.  
*   For an invalid node, the ``to...()`` methods have no effect.  
*   Calling ``setValid()`` after the node has been marked as invalid by a
    ``to...()`` method will reset this node (and path) to the last valid node
    during the traversal.  

Note that the plain methods (``child(), sibling(), parent()``, etc.) return
``Node``s and not ``PathNode``s. If ``PathNode``s are required, use instead:  

    PathNode child = PathNode(thisPathNode);
    child.toChild();  

Constructors
------------
* ``PathNode(stree, nidx=ROOT)``  
    
    Construct a ``Node`` for the given ``stree`` corresponding to the given
    ``nidx``.  

    If no ``nidx`` is given, it defaults to the root of the suffix tree.  

* ``PathNode(node)``  
    
    Construct a ``PathNode`` corresponding to the given ``node``.  

C++ includes: STreeNode.h
";

%feature("docstring") stree::PathNode::isRoot "
``isRoot() -> bool``  

Return ``true`` if this is the root node.  
";

%feature("docstring") stree::PathNode::set "
``set(pathNode)``  
``set(node)``  

Overloaded function
-------------------
* ``set(pathNode)``  
    
    Set this ``PathNode`` to the given ``pathNode``, which must belong to the
    same suffix tree.  

* ``set(node)``  
    
    Set this ``Node`` to the given ``node``, which must belong to the same
    suffix tree.  

    This is a faster version of ``(*this) = node)``.  
";

%feature("docstring") stree::PathNode::setValid "
``setValid(valid=true)``  

Mark this ``Node`` as ``valid`` (or ``invalid``, if ``valid`` is ``false``).  
";

%feature("docstring") stree::PathNode::sibling "
``sibling() -> Node``  

Return the next sibling of this node.  

If no such node exists, a ``Node`` marked as invalid is returned. Note that the
siblings are ordered lexicographically according to their edge labels.  
";

%feature("docstring") stree::PathNode::toSibling "
``toSibling()``  

Set this node to its next sibling if a next sibling exists, otherwise mark this
node as invalid.  

Note that the siblings are ordered lexicographically according to their edge
labels.  
";

%feature("docstring") stree::PathNode::dataStr "
``dataStr(width=5) -> std::string``  

Return a string representation of the data of this node.  

This is useful for debugging or understanding the suffix tree structure.  
";

%feature("docstring") stree::PathNode::isSuffix "
``isSuffix() -> bool``  

Return ``true`` if the subsequence represented by this node is a suffix of the
underlying sequence.  

Note that this does not imply that this is a leaf.  
";

%feature("docstring") stree::PathNode::nidxStr "
``nidxStr(width=3) -> std::string``  

Return a string representation of the underlying ``nidx_t``.  
";

%feature("docstring") stree::PathNode::PathNode "
``PathNode(stree, nidx=ROOT)``  
``PathNode(node)``  

Overloaded function
-------------------
* ``PathNode(stree, nidx=ROOT)``  
    
    Construct a ``Node`` for the given ``stree`` corresponding to the given
    ``nidx``.  

    If no ``nidx`` is given, it defaults to the root of the suffix tree.  

* ``PathNode(node)``  
    
    Construct a ``PathNode`` corresponding to the given ``node``.  
";

%feature("docstring") stree::PathNode::child "
``child() -> Node``  
``child(symbol) -> Node``  

Overloaded function
-------------------
* ``child() -> Node``  
    
    Return the first child node of this node.  

    If no such node exists, a ``Node`` marked as invalid is returned. Note that
    the children are ordered lexicographically according to their edge labels.  

* ``child(symbol) -> Node``  
    
    Return the child node along the edge leading away whose label begins with
    the given ``symbol``.  

    If no such node exists, a ``Node`` marked as invalid is returned.  
";

%feature("docstring") stree::PathNode::suffix "
``suffix() -> Node``  

Return the node corresponding to the first suffix of the represented sequence.  

This follows the \"suffix link\" of the suffix tree. If no such node exists, a
``Node`` marked as invalid is returned.  
";

%feature("docstring") stree::PathNode::headIndex "
``headIndex() -> nidx_t``  

Return the \"headindex\" of this node, which is an index in the sequence
represented by the suffix tree where the (sub-)sequence represented by this node
occurs.  

I.e., the (sub-)sequence represented by this node is ``seq.rawSub(headindex(),
depth())``, where ``seq`` is the sequence represented by the suffix tree.  
";

%feature("docstring") stree::PathNode::nidx "
``nidx() -> nidx_t``  

Return the ``nidx_t`` corresponding to this ``Node``.  
";

%feature("docstring") stree::PathNode::setRoot "
``setRoot()``  

Reset this ``PathNode`` to the root of the suffix tree.  
";

%feature("docstring") stree::PathNode::toSuffix "
``toSuffix()``  

Set this to the ``PathNode`` corresponding to the first suffix of the
represented sequence.  

If no suffix exists (i.e., this is the root) mark this ``PathNode`` as invalid
instead. This uses the \"suffix link\" of the suffix tree, but needs to
recompute the path.  
";

%feature("docstring") stree::PathNode::label "
``label() -> Sequence``  

Return the edge label for the edge leading to the current node.  

If no edge exists, this will be an empty ``Sequence``.  
";

%feature("docstring") stree::PathNode::repr "
``repr() -> std::string``  

Return a string representation to display in python.  
";

%feature("docstring") stree::PathNode::count "
``count() -> nidx_t``  

Return the number of occurrences of the sequence represented by this node in the
sequence represented by the suffix tree.  

For an invalid node, zero is returned.  
";

%feature("docstring") stree::PathNode::sequence "
``sequence() -> Sequence``  

Return the (sub-)sequence represented by this node.  

Note that this is ``seq.rawSub(headIndex(), depth())``, where ``seq`` is the
sequence represented by the suffix tree.  
";

%feature("docstring") stree::PathNode::parent "
``parent() -> Node``  

Return the parent ``Node``.  

If none exists, return a ``Node`` marked as invalid.  
";

%feature("docstring") stree::PathNode::isLeaf "
``isLeaf() -> bool``  

Return ``true`` if this is a leaf node.  
";

%feature("docstring") stree::PathNode::isInternal "
``isInternal() -> bool``  

Return ``true`` if this is an internal node.  
";

%feature("docstring") stree::PathNode::isValid "
``isValid() -> bool``  

Return ``true`` if valid, otherwise return ``false``.  
";

%feature("docstring") stree::PathNode::index "
``index() -> nidx_t``  

The ``index`` of a valid leaf or a valid internal node is a unique number
between 0 and ``STree.nLeafNodes()`` or between 0 and
``STree.nInternalNodes()``, respectively.  
";

%feature("docstring") stree::PathNode::depth "
``depth() -> nidx_t``  

Return the \"depth\" of the node in the suffix tree, which is the size of the
represented (sub-)sequence.  
";

%feature("docstring") stree::PathNode::toChild "
``toChild()``  
``toChild(chr)``  

Overloaded function
-------------------
* ``toChild()``  
    
    Extend this ``PathNode`` to its first child if such a node exists, otherwise
    mark this ``PathNode`` as invalid instead.  

    Note that the children are ordered lexicographically according to their edge
    labels.  

* ``toChild(chr)``  
    
    Extend this ``PathNode`` to the child node along the edge leading away whose
    label begins with the given ``symbol``.  

    If no such node exists, mark this ``PathNode`` as invalid instead.  
";

%feature("docstring") stree::PathNode::toParent "
``toParent()``  

Set this ``PathNode`` to the path to the parent of the current node, if such
exists.  

Otherwise, just mark this ``PathNode`` as invalid.  
";

// File: classtom_1_1_policy_1_1_plane.xml

// File: classtom_1_1_policy.xml


%feature("docstring") tom::Policy "
``Policy(nU=0, exploration=1)``  

Constructors
------------
* ``Policy(nU=0, exploration=1)``  

Attributes
----------
* ``nU_`` : ``unsigned int``  

* ``exploration_`` : ``double``  

C++ includes: Policy.h
";

%feature("docstring") tom::Policy::p "
``p(w) -> Eigen::VectorXd``  
";

%feature("docstring") tom::Policy::cereal::access "
``cereal::access() -> friend class``  
";

%feature("docstring") tom::Policy::Policy "
``Policy(nU=0, exploration=1)``  
";

%feature("docstring") tom::Policy::addPlane "
``addPlane(u, indices, vals)``  
";

%feature("docstring") tom::Policy::u "
``u(w, r) -> int``  
";

// File: classstree_1_1internal_1_1_pos.xml


%feature("docstring") stree::internal::Pos "
``Pos()``  

A ``Pos`` refers to the location in the suffix tree that corresponds to some
substring of the represented ``Sequence``.  

This position is unique, but can either be an explicit node (internal or leaf)
or an implicit internal node. This class is only for internal use in the suffix
tree construction. Please use the class ``stree::Position`` instead!  

Constructors
------------
* ``Pos()``  
    
    Create a ``Pos`` corresponding to the root of the suffix tree.  

Attributes
----------
* ``node_`` : ``nidx_t``  

* ``edgePtr_`` : ``nidx_t *``  

* ``hIndex_`` : ``nidx_t``  

* ``depth_`` : ``nidx_t``  

C++ includes: STreeCore.h
";

%feature("docstring") stree::internal::Pos::canonize "
``canonize(stree)``  
";

%feature("docstring") stree::internal::Pos::followSuffixLink "
``followSuffixLink(stree)``  
";

%feature("docstring") stree::internal::Pos::isExplicit "
``isExplicit() -> bool``  

Return ``true`` if the ``Pos`` corresponds to an explicit node.  
";

%feature("docstring") stree::internal::Pos::preCanonize "
``preCanonize(stree)``  
";

%feature("docstring") stree::internal::Pos::followSymbol "
``followSymbol(stree, chr) -> bool``  

If the current ``Pos`` corresponds to the subsequence ``seq``, then attempt to
move to ``seq`` concatenated with ``chr``.  

Return true if this is successful, i.e., if ``seq`` + ``chr`` is also a
substring of the represented ``sequence``.  
";

%feature("docstring") stree::internal::Pos::Pos "
``Pos()``  

Create a ``Pos`` corresponding to the root of the suffix tree.  
";

// File: classstree_1_1_position.xml


%feature("docstring") stree::Position "
``Position(stree, sequence=Sequence())``  
``Position(node)``  

This class represents any position in the suffix tree.  

A ``Position`` that is a node (leaf or internal) of the suffix tree is
considered \"explicit\", while a ``Position`` on some edge is called
\"implicit\".  

The ``to...()`` methods are provided to navigate the suffix tree structure:  

*   ``to...()`` sets this ``Position`` to its suffix, extension by a sequence,
    etc.  
*   If no such exists, then it is simply marked as invalid.  
*   For an invalid ``Position``, the ``to...()`` methods have no effect.  
*   Calling ``setValid()`` after the ``Position`` has been marked as invalid by
    a ``to...()`` method will reset the ``Position`` to the last valid one
    during the traversal.  

Constructors
------------
* ``Position(stree, sequence=Sequence())``  
    
    Construct the ``Position`` in the given ``stree`` corresponding to the given
    ``sequence``.  

    If no corresponding position exists, this ``Position`` will correspond to
    the longest possible prefix of the provided ``sequence`` and will be marked
    as invalid.  

* ``Position(node)``  
    
    Construct a ``Position`` from the given ``node`` of type ``EdgeNode``.  

C++ includes: STreeNode.h
";

%feature("docstring") stree::Position::isSuffix "
``isSuffix() -> bool``  

Return ``true`` if the represented subsequence is a suffix of the underlying
sequence.  

Note that this does not imply that this is a leaf.  
";

%feature("docstring") stree::Position::sibling "
``sibling() -> Position``  

Return the next sibling ``Position`` in the suffix tree structure viewed as a
*suffix trie*, i.e., where all positions are seen as nodes and all edges have
length one.  

If no sibling exists, a ``Position`` marked as invalid is returned. Note that
the siblings are ordered lexicographically according to their edge symbols.  
";

%feature("docstring") stree::Position::repr "
``repr() -> std::string``  

Return a string representation to display in python.  
";

%feature("docstring") stree::Position::suffix "
``suffix() -> Position``  

Return the ``Position`` corresponding to the first suffix of the represented
sequence.  

This uses the \"suffix link\" of the suffix tree. If no suffix exists (i.e.,
this is the root), a ``Position`` marked as invalid is returned.  
";

%feature("docstring") stree::Position::child "
``child() -> Position``  

Return the first child ``Position`` in the suffix tree structure viewed as a
*suffix trie*, i.e., where all positions are seen as nodes and all edges have
length one.  

If no child exists, a ``Position`` marked as invalid is returned. Note that the
children are ordered lexicographically according to their edge symbols.  
";

%feature("docstring") stree::Position::set "
``set(position)``  

Set this ``Position`` to the given ``position``, which must belong to the same
suffix tree.  

This is a faster version of ``(*this) = position)``.  
";

%feature("docstring") stree::Position::isInternal "
``isInternal() -> bool``  

Return ``true`` if this is an internal node or an implicit position.  
";

%feature("docstring") stree::Position::isLeaf "
``isLeaf() -> bool``  

Return ``true`` if this is a leaf node.  
";

%feature("docstring") stree::Position::edge "
``edge() -> EdgeNode``  

Return the ``EdgeNode`` that this ``Position`` lies on.  

If this ``Position`` is explicit, then the ``EdgeNode`` will be the node (with
edge leading to it) of the position.  
";

%feature("docstring") stree::Position::label "
``label() -> Sequence``  

Return the sub-sequence of the edge label up to this position.  

For an explicit position this is just the edge label (see ``edge.label()``).  
";

%feature("docstring") stree::Position::toSuffix "
``toSuffix()``  

Set this to the ``Position`` corresponding to the first suffix of the
represented sequence.  

If no suffix exists (i.e., this is the root) mark this ``Position`` as invalid
instead. This uses the \"suffix link\" of the suffix tree.  
";

%feature("docstring") stree::Position::headIndex "
``headIndex() -> nidx_t``  

Return the \"headindex\" of this node, which is an index in the sequence
represented by the suffix tree where the (sub-)sequence represented by this node
occurs.  

I.e., the (sub-)sequence represented by this position is
``seq.rawSub(headindex(), depth())``, where ``seq`` is the sequence represented
by the suffix tree.  
";

%feature("docstring") stree::Position::toChild "
``toChild()``  

Set this ``Position`` to its first child position in the suffix tree structure
viewed as a *suffix trie*, i.e., where all positions are seen as nodes and all
edges have length one.  

If no child exists, mark this ``Position`` as invalid instead. Note that the
children are ordered lexicographically according to their edge symbols.  
";

%feature("docstring") stree::Position::isValid "
``isValid() -> bool``  

Return ``true`` if valid, otherwise return ``false``.  
";

%feature("docstring") stree::Position::Position "
``Position(stree, sequence=Sequence())``  
``Position(node)``  

Overloaded function
-------------------
* ``Position(stree, sequence=Sequence())``  
    
    Construct the ``Position`` in the given ``stree`` corresponding to the given
    ``sequence``.  

    If no corresponding position exists, this ``Position`` will correspond to
    the longest possible prefix of the provided ``sequence`` and will be marked
    as invalid.  

* ``Position(node)``  
    
    Construct a ``Position`` from the given ``node`` of type ``EdgeNode``.  
";

%feature("docstring") stree::Position::setValid "
``setValid(valid=true)``  

Mark this ``Position`` as ``valid`` (default: ``true``).  
";

%feature("docstring") stree::Position::isExplicit "
``isExplicit() -> bool``  

Return ``true`` if this is a node (internal or leaf).  

Otherwise, this ``Position`` is \"implicit\", i.e., it lies on some edge.  
";

%feature("docstring") stree::Position::toSymbol "
``toSymbol(symbol)``  

Update this ``Position`` such that it represents a subsequence extended by the
given ``symbol``.  

If no such position exists, this ``Position`` is unchanged but marked as
invalid. For an invalid ``Position`` this function has no effect.  
";

%feature("docstring") stree::Position::toSibling "
``toSibling()``  

Set this ``Position`` to its next sibling position in the suffix tree structure
viewed as a *suffix trie*, i.e., where all positions are seen as nodes and all
edges have length one.  

If no sibling exists, mark this ``Position`` as invalid instead. Note that the
siblings are ordered lexicographically according to their edge symbols.  
";

%feature("docstring") stree::Position::toSequence "
``toSequence(sequence)``  

Update this ``Position`` such that it represents a subsequence extended by the
given ``sequence``.  

If no such position exists, this ``Position`` is updated symbol-wise according
to the given ``sequence`` as far as possible and then marked as invalid. For an
invalid ``Position`` this function has no effect.  
";

%feature("docstring") stree::Position::setRoot "
``setRoot()``  

Reset this ``Position`` to the root of the suffix tree.  
";

%feature("docstring") stree::Position::isRoot "
``isRoot() -> bool``  

Return ``true`` if this is the root node.  
";

%feature("docstring") stree::Position::depth "
``depth() -> nidx_t``  

Return the \"depth\" of this position, which is the size of the represented
(sub-)sequence.  
";

%feature("docstring") stree::Position::sequence "
``sequence() -> Sequence``  

Return the (sub-)sequence represented by this ``Position``.  

Note that this is ``seq.rawSub(headIndex(), depth())``, where ``seq`` is the
sequence represented by the suffix tree.  
";

%feature("docstring") stree::Position::toDepth "
``toDepth(depth)``  
";

%feature("docstring") stree::Position::toExplicit "
``toExplicit()``  

If this ``Position`` is not explicit, i.e., it lies on an edge of the suffix
tree structure, then set this ``Position`` to the deeper node end-point of that
edge, making this position explicit.  
";

%feature("docstring") stree::Position::count "
``count() -> nidx_t``  

Return the number of occurrences of the sequence represented by this
``Position`` in the sequence represented by the suffix tree.  

For an invalid ``Position``, zero is returned.  
";

// File: classstree_1_1_position_relevance.xml


%feature("docstring") stree::PositionRelevance "

This class computes a \"relevance\" value for suffix tree ``Positions`` in a way
that can be customized from Python by inheriting and overwriting the ``compute``
method.  

C++ includes: STreeNode.h
";

%feature("docstring") stree::PositionRelevance::compute "
``compute(position) -> double``  

Return the relevance value for the given ``position``  
";

%feature("docstring") stree::PositionRelevance::~PositionRelevance "
``~PositionRelevance()``  
";

// File: classstree_1_1_postfix_iterator.xml


%feature("docstring") stree::PostfixIterator "
``PostfixIterator(stree)``  

An iterator to traverse the nodes (leaf and internal) of the suffix tree in
*postfix* order, which means traversing the tree depth first and visiting the
nodes on the upward pass.  

A ``PostfixIterator`` is a ``PathNode`` with a ``toNext()`` method to move to
the next node in postfix order. The end of iteration is signaled by marking the
iterator as invalid, which can be checked by the inherited method ``isValid()``.  

Note that this iterator is wrapped as a native python iterator, i.e., the
following is possible:  

    for node in PostfixIterator(stree):
        print(node.nidxStr())  

Constructors
------------
* ``PostfixIterator(stree)``  
    
    Create a ``PostfixIterator`` for the given ``stree``.  

C++ includes: STreeIterators.h
";

%feature("docstring") stree::PostfixIterator::isInternal "
``isInternal() -> bool``  

Return ``true`` if this is an internal node.  
";

%feature("docstring") stree::PostfixIterator::setRoot "
``setRoot()``  

Reset this ``PathNode`` to the root of the suffix tree.  
";

%feature("docstring") stree::PostfixIterator::label "
``label() -> Sequence``  

Return the edge label for the edge leading to the current node.  

If no edge exists, this will be an empty ``Sequence``.  
";

%feature("docstring") stree::PostfixIterator::nidx "
``nidx() -> nidx_t``  

Return the ``nidx_t`` corresponding to this ``Node``.  
";

%feature("docstring") stree::PostfixIterator::suffix "
``suffix() -> Node``  

Return the node corresponding to the first suffix of the represented sequence.  

This follows the \"suffix link\" of the suffix tree. If no such node exists, a
``Node`` marked as invalid is returned.  
";

%feature("docstring") stree::PostfixIterator::headIndex "
``headIndex() -> nidx_t``  

Return the \"headindex\" of this node, which is an index in the sequence
represented by the suffix tree where the (sub-)sequence represented by this node
occurs.  

I.e., the (sub-)sequence represented by this node is ``seq.rawSub(headindex(),
depth())``, where ``seq`` is the sequence represented by the suffix tree.  
";

%feature("docstring") stree::PostfixIterator::dataStr "
``dataStr(width=5) -> std::string``  

Return a string representation of the data of this node.  

This is useful for debugging or understanding the suffix tree structure.  
";

%feature("docstring") stree::PostfixIterator::isSuffix "
``isSuffix() -> bool``  

Return ``true`` if the subsequence represented by this node is a suffix of the
underlying sequence.  

Note that this does not imply that this is a leaf.  
";

%feature("docstring") stree::PostfixIterator::child "
``child() -> Node``  
``child(symbol) -> Node``  

Overloaded function
-------------------
* ``child() -> Node``  
    
    Return the first child node of this node.  

    If no such node exists, a ``Node`` marked as invalid is returned. Note that
    the children are ordered lexicographically according to their edge labels.  

* ``child(symbol) -> Node``  
    
    Return the child node along the edge leading away whose label begins with
    the given ``symbol``.  

    If no such node exists, a ``Node`` marked as invalid is returned.  
";

%feature("docstring") stree::PostfixIterator::nidxStr "
``nidxStr(width=3) -> std::string``  

Return a string representation of the underlying ``nidx_t``.  
";

%feature("docstring") stree::PostfixIterator::toNext "
``toNext()``  

Set this ``PostfixIterator`` to the next node in postfix order.  

If none exists, this ``PostfixIterator`` will be marked as invalid. Calling
``toNext()`` on an invalid ``PostfixIterator`` has no effect.  
";

%feature("docstring") stree::PostfixIterator::count "
``count() -> nidx_t``  

Return the number of occurrences of the sequence represented by this node in the
sequence represented by the suffix tree.  

For an invalid node, zero is returned.  
";

%feature("docstring") stree::PostfixIterator::depth "
``depth() -> nidx_t``  

Return the \"depth\" of the node in the suffix tree, which is the size of the
represented (sub-)sequence.  
";

%feature("docstring") stree::PostfixIterator::setValid "
``setValid(valid=true)``  

Mark this ``Node`` as ``valid`` (or ``invalid``, if ``valid`` is ``false``).  
";

%feature("docstring") stree::PostfixIterator::isValid "
``isValid() -> bool``  

Return ``true`` if valid, otherwise return ``false``.  
";

%feature("docstring") stree::PostfixIterator::set "
``set(pathNode)``  
``set(node)``  

Overloaded function
-------------------
* ``set(pathNode)``  
    
    Set this ``PathNode`` to the given ``pathNode``, which must belong to the
    same suffix tree.  

* ``set(node)``  
    
    Set this ``Node`` to the given ``node``, which must belong to the same
    suffix tree.  

    This is a faster version of ``(*this) = node)``.  
";

%feature("docstring") stree::PostfixIterator::toChild "
``toChild()``  
``toChild(chr)``  

Overloaded function
-------------------
* ``toChild()``  
    
    Extend this ``PathNode`` to its first child if such a node exists, otherwise
    mark this ``PathNode`` as invalid instead.  

    Note that the children are ordered lexicographically according to their edge
    labels.  

* ``toChild(chr)``  
    
    Extend this ``PathNode`` to the child node along the edge leading away whose
    label begins with the given ``symbol``.  

    If no such node exists, mark this ``PathNode`` as invalid instead.  
";

%feature("docstring") stree::PostfixIterator::PostfixIterator "
``PostfixIterator(stree)``  

Create a ``PostfixIterator`` for the given ``stree``.  
";

%feature("docstring") stree::PostfixIterator::isLeaf "
``isLeaf() -> bool``  

Return ``true`` if this is a leaf node.  
";

%feature("docstring") stree::PostfixIterator::repr "
``repr() -> std::string``  

Return a string representation to display in python.  
";

%feature("docstring") stree::PostfixIterator::sequence "
``sequence() -> Sequence``  

Return the (sub-)sequence represented by this node.  

Note that this is ``seq.rawSub(headIndex(), depth())``, where ``seq`` is the
sequence represented by the suffix tree.  
";

%feature("docstring") stree::PostfixIterator::toParent "
``toParent()``  

Set this ``PathNode`` to the path to the parent of the current node, if such
exists.  

Otherwise, just mark this ``PathNode`` as invalid.  
";

%feature("docstring") stree::PostfixIterator::isRoot "
``isRoot() -> bool``  

Return ``true`` if this is the root node.  
";

%feature("docstring") stree::PostfixIterator::sibling "
``sibling() -> Node``  

Return the next sibling of this node.  

If no such node exists, a ``Node`` marked as invalid is returned. Note that the
siblings are ordered lexicographically according to their edge labels.  
";

%feature("docstring") stree::PostfixIterator::toSuffix "
``toSuffix()``  

Set this to the ``PathNode`` corresponding to the first suffix of the
represented sequence.  

If no suffix exists (i.e., this is the root) mark this ``PathNode`` as invalid
instead. This uses the \"suffix link\" of the suffix tree, but needs to
recompute the path.  
";

%feature("docstring") stree::PostfixIterator::parent "
``parent() -> Node``  

Return the parent ``Node``.  

If none exists, return a ``Node`` marked as invalid.  
";

%feature("docstring") stree::PostfixIterator::toSibling "
``toSibling()``  

Set this node to its next sibling if a next sibling exists, otherwise mark this
node as invalid.  

Note that the siblings are ordered lexicographically according to their edge
labels.  
";

%feature("docstring") stree::PostfixIterator::index "
``index() -> nidx_t``  

The ``index`` of a valid leaf or a valid internal node is a unique number
between 0 and ``STree.nLeafNodes()`` or between 0 and
``STree.nInternalNodes()``, respectively.  
";

// File: classstree_1_1_prefix_iterator.xml


%feature("docstring") stree::PrefixIterator "
``PrefixIterator(stree)``  

An iterator to traverse the nodes (leaf and internal) of the suffix tree in
*prefix* order, which means traversing the tree depth first and visiting the
nodes on the downward pass.  

A ``PrefixIterator`` is a ``PathNode`` with a ``toNext()`` method to move to the
next node in prefix order. The end of iteration is signaled by marking the
iterator as invalid, which can be checked by the inherited method ``isValid()``.  

Note that this iterator is wrapped as a native python iterator, i.e., the
following is possible:  

    for node in PrefixIterator(stree):
        print(node.nidxStr())  

Constructors
------------
* ``PrefixIterator(stree)``  
    
    Create a ``PrefixIterator`` for the given ``stree``.  

C++ includes: STreeIterators.h
";

%feature("docstring") stree::PrefixIterator::headIndex "
``headIndex() -> nidx_t``  

Return the \"headindex\" of this node, which is an index in the sequence
represented by the suffix tree where the (sub-)sequence represented by this node
occurs.  

I.e., the (sub-)sequence represented by this node is ``seq.rawSub(headindex(),
depth())``, where ``seq`` is the sequence represented by the suffix tree.  
";

%feature("docstring") stree::PrefixIterator::label "
``label() -> Sequence``  

Return the edge label for the edge leading to the current node.  

If no edge exists, this will be an empty ``Sequence``.  
";

%feature("docstring") stree::PrefixIterator::parent "
``parent() -> Node``  

Return the parent ``Node``.  

If none exists, return a ``Node`` marked as invalid.  
";

%feature("docstring") stree::PrefixIterator::isSuffix "
``isSuffix() -> bool``  

Return ``true`` if the subsequence represented by this node is a suffix of the
underlying sequence.  

Note that this does not imply that this is a leaf.  
";

%feature("docstring") stree::PrefixIterator::toSibling "
``toSibling()``  

Set this node to its next sibling if a next sibling exists, otherwise mark this
node as invalid.  

Note that the siblings are ordered lexicographically according to their edge
labels.  
";

%feature("docstring") stree::PrefixIterator::depth "
``depth() -> nidx_t``  

Return the \"depth\" of the node in the suffix tree, which is the size of the
represented (sub-)sequence.  
";

%feature("docstring") stree::PrefixIterator::isInternal "
``isInternal() -> bool``  

Return ``true`` if this is an internal node.  
";

%feature("docstring") stree::PrefixIterator::setRoot "
``setRoot()``  

Reset this ``PathNode`` to the root of the suffix tree.  
";

%feature("docstring") stree::PrefixIterator::toParent "
``toParent()``  

Set this ``PathNode`` to the path to the parent of the current node, if such
exists.  

Otherwise, just mark this ``PathNode`` as invalid.  
";

%feature("docstring") stree::PrefixIterator::set "
``set(pathNode)``  
``set(node)``  

Overloaded function
-------------------
* ``set(pathNode)``  
    
    Set this ``PathNode`` to the given ``pathNode``, which must belong to the
    same suffix tree.  

* ``set(node)``  
    
    Set this ``Node`` to the given ``node``, which must belong to the same
    suffix tree.  

    This is a faster version of ``(*this) = node)``.  
";

%feature("docstring") stree::PrefixIterator::index "
``index() -> nidx_t``  

The ``index`` of a valid leaf or a valid internal node is a unique number
between 0 and ``STree.nLeafNodes()`` or between 0 and
``STree.nInternalNodes()``, respectively.  
";

%feature("docstring") stree::PrefixIterator::toSuffix "
``toSuffix()``  

Set this to the ``PathNode`` corresponding to the first suffix of the
represented sequence.  

If no suffix exists (i.e., this is the root) mark this ``PathNode`` as invalid
instead. This uses the \"suffix link\" of the suffix tree, but needs to
recompute the path.  
";

%feature("docstring") stree::PrefixIterator::nidxStr "
``nidxStr(width=3) -> std::string``  

Return a string representation of the underlying ``nidx_t``.  
";

%feature("docstring") stree::PrefixIterator::sibling "
``sibling() -> Node``  

Return the next sibling of this node.  

If no such node exists, a ``Node`` marked as invalid is returned. Note that the
siblings are ordered lexicographically according to their edge labels.  
";

%feature("docstring") stree::PrefixIterator::toNext "
``toNext()``  

Set this ``PrefixIterator`` to the next node in prefix order.  

If none exists, this ``PrefixIterator`` will be marked as invalid. Calling
``toNext()`` on an invalid ``PrefixIterator`` has no effect.  
";

%feature("docstring") stree::PrefixIterator::nidx "
``nidx() -> nidx_t``  

Return the ``nidx_t`` corresponding to this ``Node``.  
";

%feature("docstring") stree::PrefixIterator::count "
``count() -> nidx_t``  

Return the number of occurrences of the sequence represented by this node in the
sequence represented by the suffix tree.  

For an invalid node, zero is returned.  
";

%feature("docstring") stree::PrefixIterator::isLeaf "
``isLeaf() -> bool``  

Return ``true`` if this is a leaf node.  
";

%feature("docstring") stree::PrefixIterator::sequence "
``sequence() -> Sequence``  

Return the (sub-)sequence represented by this node.  

Note that this is ``seq.rawSub(headIndex(), depth())``, where ``seq`` is the
sequence represented by the suffix tree.  
";

%feature("docstring") stree::PrefixIterator::child "
``child() -> Node``  
``child(symbol) -> Node``  

Overloaded function
-------------------
* ``child() -> Node``  
    
    Return the first child node of this node.  

    If no such node exists, a ``Node`` marked as invalid is returned. Note that
    the children are ordered lexicographically according to their edge labels.  

* ``child(symbol) -> Node``  
    
    Return the child node along the edge leading away whose label begins with
    the given ``symbol``.  

    If no such node exists, a ``Node`` marked as invalid is returned.  
";

%feature("docstring") stree::PrefixIterator::dataStr "
``dataStr(width=5) -> std::string``  

Return a string representation of the data of this node.  

This is useful for debugging or understanding the suffix tree structure.  
";

%feature("docstring") stree::PrefixIterator::isValid "
``isValid() -> bool``  

Return ``true`` if valid, otherwise return ``false``.  
";

%feature("docstring") stree::PrefixIterator::isRoot "
``isRoot() -> bool``  

Return ``true`` if this is the root node.  
";

%feature("docstring") stree::PrefixIterator::PrefixIterator "
``PrefixIterator(stree)``  

Create a ``PrefixIterator`` for the given ``stree``.  
";

%feature("docstring") stree::PrefixIterator::toChild "
``toChild()``  
``toChild(chr)``  

Overloaded function
-------------------
* ``toChild()``  
    
    Extend this ``PathNode`` to its first child if such a node exists, otherwise
    mark this ``PathNode`` as invalid instead.  

    Note that the children are ordered lexicographically according to their edge
    labels.  

* ``toChild(chr)``  
    
    Extend this ``PathNode`` to the child node along the edge leading away whose
    label begins with the given ``symbol``.  

    If no such node exists, mark this ``PathNode`` as invalid instead.  
";

%feature("docstring") stree::PrefixIterator::setValid "
``setValid(valid=true)``  

Mark this ``Node`` as ``valid`` (or ``invalid``, if ``valid`` is ``false``).  
";

%feature("docstring") stree::PrefixIterator::repr "
``repr() -> std::string``  

Return a string representation to display in python.  
";

%feature("docstring") stree::PrefixIterator::suffix "
``suffix() -> Node``  

Return the node corresponding to the first suffix of the represented sequence.  

This follows the \"suffix link\" of the suffix tree. If no such node exists, a
``Node`` marked as invalid is returned.  
";

// File: classtom_1_1_random.xml


%feature("docstring") tom::Random "
``Random()``  
``Random(seed)``  

This class provides basic random number generation functionality.  

Constructors
------------
* ``Random()``  
    
    Create a ``Random`` object initialized with a random seed.  

* ``Random(seed)``  
    
    Create a ``Random`` object initialized with the given ``seed``.  

C++ includes: Random.h
";

%feature("docstring") tom::Random::integer "
``integer(n) -> unsigned int``  

Return a non-negative integer sampled uniformly from the set {0, ..., ``n``-1}.  
";

%feature("docstring") tom::Random::sample "
``sample(probArray) -> unsigned int``  

Return a non-negative integer from the set {0, ..., n-1}, where n is the size of
the given array ``probArray``, distributed according to the discrete
distribution described by the ``probArray``.  

Note that the probabilities in ``probArray`` are assumed to sum to 1.  
";

%feature("docstring") tom::Random::random "
``random() -> double``  
``random(m, n) -> Eigen::MatrixXd``  

Overloaded function
-------------------
* ``random() -> double``  
    
    Return a double value sampled uniformly from in the range [0, 1).  

* ``random(m, n) -> Eigen::MatrixXd``  
    
    Return a matrix of size ``m`` x ``n`` with uniformly random entries from [0,
    1).  
";

%feature("docstring") tom::Random::seed "
``seed() -> unsigned int``  
``seed(seedValue)``  

Overloaded function
-------------------
* ``seed() -> unsigned int``  
    
    Randomly seed the random number generator and return the seed.  

* ``seed(seedValue)``  
    
    Seed the random number generator with the given ``seedValue``.  
";

%feature("docstring") tom::Random::Random "
``Random()``  
``Random(seed)``  

Overloaded function
-------------------
* ``Random()``  
    
    Create a ``Random`` object initialized with a random seed.  

* ``Random(seed)``  
    
    Create a ``Random`` object initialized with the given ``seed``.  
";

// File: classstree_1_1_r_b_tree.xml


%feature("docstring") stree::RBTree "

An implementation of left-leaning red-black trees following \"Left-Leaning Red-
Black Trees\" by Robert Sedgewick from 2008.  

The functions can be used with generic kinds of nodes that need to be specified
by providing an ``RBTreeNodeTraits`` object when calling the functions.  

Note that the red-black tree will be threaded if a ``NULL`` ``RBNodePtr`` can
still address a node. This allows iterating over the stored values in their
order merely by following the left or right pointers (even if ``NULL``). The
left-and rightmost ``RBNodePtr`` will be inherited from the left and right
``RBNodePtr`` of the original root, i.e., the first node ``RBNodePtr`` used in
the first insertion operation when constructing the red-black tree, while every
threaded right or left ``RBNodePtr`` will be set to the according ``RBNodePtr``
and then marked as a thread by calling the function ``setThread``. It is up to
the user to implement a suitable ``RBNodePtr`` structure to be able to
distinguish the left-/rightmost ``RBNodePtr`` from a threaded ``RBNodePtr``
(e.g., if the color of nodes is stored in the ``RBNodePtr``, even if ``NULL``,
then a threaded ``RBNodePtr`` may be ``NULL`` and colored red, while the
left-/rightmost ``RBNodePtr`` may be ``NULL`` and black).  

C++ includes: RBTree.h
";

%feature("docstring") stree::RBTree::insert "
``insert(h, n, key, rbnt)``  

insert a node into the red-black tree according to a given key.  

Please see the general remarks about threading.  

Parameters
----------
* ``h`` :  
    the ``RBNodePtr`` to the root of the red-black tree (this will point to the
    new root after the insertion operation)  
* ``n`` :  
    the node to be inserted  
* ``key`` :  
    the key value of the new node  
* ``rbnt`` :  
    an ``RBTreeNodeTraits`` object  
";

%feature("docstring") stree::RBTree::fixThreading "
``fixThreading(n, rbnt)``  

fixes the threads leading to the node ``n`` in the red-black tree; this function
needs to be called after replacing a node in the red-black tree structure.  

Parameters
----------
* ``n`` :  
    the node whose threading needs to be fixed  
* ``rbnt`` :  
    an ``RBTreeNodeTraits`` object  
";

%feature("docstring") stree::RBTree::find "
``find(h, key, rbnt) -> RBNodePtr &``  

search for a node matching a given ``key`` in the binary search tree below a
given node ``h``.  

Parameters
----------
* ``h`` :  
    the node below which (and including) to search for the key  
* ``key`` :  
    the key value to be searched for  
* ``rbnt`` :  
    an ``RBTreeNodeTraits`` object  

Returns
-------
the ``RBNodePtr&`` of the parent node to the node found, or some ``NULL``
``RBNodePtr``, if no matching node found.  
";

// File: classstree_1_1internal_1_1_r_b_tree_node_traits.xml


%feature("docstring") stree::internal::RBTreeNodeTraits "
``RBTreeNodeTraits(bst, parentDepth=0)``  

An object specifying the node traits for the underlying red-black tree
implementation.  

Constructors
------------
* ``RBTreeNodeTraits(bst, parentDepth=0)``  

C++ includes: STreeCore.h
";

%feature("docstring") stree::internal::RBTreeNodeTraits::left "
``left(n) -> RBNodePtr &``  
";

%feature("docstring") stree::internal::RBTreeNodeTraits::less "
``less(k, n) -> bool``  
";

%feature("docstring") stree::internal::RBTreeNodeTraits::getColor "
``getColor(n) -> bool``  
";

%feature("docstring") stree::internal::RBTreeNodeTraits::isNull "
``isNull(n) -> bool``  
";

%feature("docstring") stree::internal::RBTreeNodeTraits::setThread "
``setThread(n)``  
";

%feature("docstring") stree::internal::RBTreeNodeTraits::setColor "
``setColor(n, c)``  
";

%feature("docstring") stree::internal::RBTreeNodeTraits::right "
``right(n) -> RBNodePtr &``  
";

%feature("docstring") stree::internal::RBTreeNodeTraits::RBTreeNodeTraits "
``RBTreeNodeTraits(bst, parentDepth=0)``  
";

%feature("docstring") stree::internal::RBTreeNodeTraits::set "
``set(n, nNew)``  
";

%feature("docstring") stree::internal::RBTreeNodeTraits::equals "
``equals(k, n) -> bool``  
";

// File: classstree_1_1_r_b_tree_node_traits_template.xml


%feature("docstring") stree::RBTreeNodeTraitsTemplate "

An example of which functions an ``RBTreeNodeTraits`` object must provide to be
used with ``RBTree``.  

It is assumed that each node has a left and right pointer to child nodes,
forming a binary tree. These pointers may be any user defined data structure
used to index the nodes. For this reason, all needed functions to manipulate the
tree structure must me provided in an ``RBTreeNodeTraits`` object. Note that
nodes are referred to only by ``RBNodePtr`` objects and never directly.  

C++ includes: RBTree.h
";

// File: classtom_1_1_sequence.xml


%feature("docstring") tom::Sequence "
``Sequence(symbol_list, nO, nU=0)``  
``Sequence(length=0, nO=0, nU=0)``  
``Sequence(json_representation)``  

This object represents a sequence, subsequence view or io-sequence and stores
the size ``nO()`` of the output and ``nU()`` of the input alphabet.  

If the size of the input alphabet is zero, this is just an ordinary sequence
:math:`o_0...o_{N-1}` of symbols :math:`o_t` with zero-based indexing. An io-
sequence is represented as a simple sequence :math:`u_0o_0...u_{N-1}o_{N-1}` of
input symbols :math:`u_t` and output symbols :math:`o_t`. For io-sequences, we
distinguish *size* and *length*: *size* is always the number of symbols, while
the *length* is the number of io symbol pairs, which is just the *size* for
ordinary sequences, and 2 * *size* for (aligned) io-sequences.  

There are three ways to interact with this sequence:  

1.  The ``raw...()`` methods. These just access the sequence as a raw symbol
    sequence, treating each input or output symbol as a separate symbol. I.e.,
    for the io sequence :math:`u_0o_0...u_{N-1}o_{N-1}`, the ``rawSize()`` is
    :math:`2N`, and ``rawAt(i)`` is :math:`u_{i/2}` if :math:`i` is even or
    :math:`o_{i/2}` if :math:`i` is odd. Etc.  
2.  The methods not prefixed by \"raw\" treat each io-symbol *pair* as one
    symbol: ``at(i)`` is the symbol-pair :math:`(u_i, o_i)` (actually, the
    subsequence at index :math:`i` of length 1: ``sub(i,1)``), the methods
    ``u(i)`` and ``o(i)`` return the input symbol :math:`u_i` or respectively
    output symbol :math:`o_i`, and the ``length()`` is the *length* of the
    sequence -- the number N of io symbol pairs. For ordinary sequences these
    methods are equivalent to the ``raw...()`` ones.  
3.  Access using python ``[]``-syntax (including slicing) and python iteration
    just treats all sequences as raw symbol sequences (as in 1.).  

This object always represents a view to underlying sequence data, i.e., copies,
slices and subsequences always point to the same underlying data. To obtain a
real deep copy, the ``copy()`` member function is provided.  

Constructors
------------
* ``Sequence(symbol_list, nO, nU=0)``  
    
    Construct a ``Sequence`` of given output alphabet size ``nO`` and input
    alphabet size ``nU`` from the given ``symbol_list``.  

    This may be a list ``[u_0, ..., u_{N-1}]`` (or ``[u_0, o_0, ..., u_{N-1},
    o_{N-1}]`` for an io-sequence), or a ``std::vector<int>``. The contents of
    ``symbol_list`` is copied.  

* ``Sequence(length=0, nO=0, nU=0)``  
    
    Construct a ``Sequence`` of given ``length``, output alphabet size ``nO``
    and input alphabet size ``nU`` initialized with zeros.  

    In the case of an io-sequence (if ``nU != 0``), the ``rawSize()`` will be 2
    * ``length``.  

* ``Sequence(json_representation)``  
    
    Construct a ``Sequence`` corresponding to the given string
    ``json_representation``.  

    The format must correspond to what ``toJSON()`` produces.  

C++ includes: Sequence.h
";

/*
 Constructors 
*/

/*
 Accessors and Properties 
*/

/*
 Interface as raw symbol sequence 
*/

/*
 Standard interface 
*/

/*
 Functionality 
*/

/*
 IO-functions 
*/

%feature("docstring") tom::Sequence::isAligned "
``isAligned() -> bool``  

Return true if this sequence is io-aligned at its beginning, i.e., if either:  

*   it is not reversed and front aligned  
*   it is reversed and back aligned  

This basically says that the first symbol of this sequence is what it is
supposed to be: an input symbol if it is not reversed, or an output symbol
otherwise.  
";

%feature("docstring") tom::Sequence::nOutputSymbols "
``nOutputSymbols() -> Symbol``  

Return the size of the output alphabet.  
";

%feature("docstring") tom::Sequence::slice "
``slice(begin, end=NoIndex, forwards=true) -> Sequence``  

Return a subsequence from the ``begin`` index up to (and not including) the
``end`` index, or if ``forwards`` is set to ``false``, a reverse sequence from
the ``begin`` position up to (and not including) the ``end`` position, where
each index covers one io-pair.  

Negative indexing is supported, and ``begin`` and ``end`` may be set to
``NoIndex``, and then extend to the beginning or end of the sequence, depending
on ``reverse``.  

Note that the requested slice must define a sub-sequence of this sequence or a
sequence of size zero.  
";

%feature("docstring") tom::Sequence::rawSlice "
``rawSlice(begin=NoIndex, end=NoIndex, forwards=true) -> Sequence``  

Return a subsequence from the ``begin`` index up to (and not including) the
``end`` index, or if ``forwards`` is set to ``false``, a reverse sequence from
the ``begin`` position up to (and not including) the ``end`` position.  

Negative indexing is supported, and ``begin`` and ``end`` may be set to
``NoIndex``, and then extend to the beginning or end of the sequence, depending
on ``reverse``.  

Note that the requested slice must define a sub-sequence of this sequence or a
sequence of size zero.  
";

%feature("docstring") tom::Sequence::u "
``u(idx) -> Symbol``  
``u(idx, u)``  

Overloaded function
-------------------
* ``u(idx) -> Symbol``  
    
    Return the input symbol at index ``idx``, where each index covers one io-
    pair, or return zero if this is a plain (non-io) sequence.  

    Negative indexing is supported.  

* ``u(idx, u)``  
    
    Set the input symbol at index ``idx`` to ``u``, where each index covers one
    io-pair, or do nothing if this is a plain (non-io) sequence.  

    Negative indexing is supported.  
";

%feature("docstring") tom::Sequence::isIO "
``isIO() -> bool``  

Return ``true`` if this is an input-output sequence, i.e., if the input alphabet
size ``nInputSymbols()`` is non-zero.  
";

%feature("docstring") tom::Sequence::rawAt "
``rawAt(idx) -> Symbol``  
``rawAt(idx, x)``  

Overloaded function
-------------------
* ``rawAt(idx) -> Symbol``  
    
    Return the symbol at index ``idx``, treating io-sequences as raw sequences.  

* ``rawAt(idx, x)``  
    
    Set the symbol at index ``idx`` to ``x``, treating io-sequences as raw
    sequences.  

    Negative indexing is supported.  
";

%feature("docstring") tom::Sequence::lexicographicIndex "
``lexicographicIndex(withinSequencesOfSameLength=false) -> unsigned long``  
";

%feature("docstring") tom::Sequence::rawSize "
``rawSize() -> long``  

Return the *size* of the represented sequence, i.e., the raw symbol count,
counting each input and output as one symbol.  

Generally, use ``length()`` instead.  
";

%feature("docstring") tom::Sequence::at "
``at(idx) -> Sequence``  

Return the io symbol pair at index ``idx``, where each index covers one io-pair.  

This returns ``sub(idx, 1)``, so even for plain sequences, the return value is
not a symbol. For plain sequences, generally use ``rawAt(idx)`` or ``o(idx)``
instead.  

Negative indexing is supported.  
";

%feature("docstring") tom::Sequence::count "
``count(seq) -> unsigned int``  

Count the number of occurrences of the given ``Sequence`` ``seq`` as a sub-
sequence of this ``Sequence``.  
";

%feature("docstring") tom::Sequence::toJSON "
``toJSON() -> std::string``  
";

%feature("docstring") tom::Sequence::length "
``length() -> long``  

Return the *length* of this sequence.  

For plain sequences this is the same as ``rawSize()``. For (aligned) io-
sequences this is the number of io-symbol pairs, which is half the *size*.  

For unaligned io-sequences, the *length* means the number of covered io-
sequence-pair indices. Example:  
For the io-sequence :math:`o_0 u_1o_1 ... u_{N-2}o_{N-2} u_{N-1}`, which is
neither front nor back aligned, the *length* is :math:`N`, since :math:`N` io-
pair indices are covered, but the *size* -- the number of raw symbols -- is only
:math:`2N - 2`.  
";

%feature("docstring") tom::Sequence::rawSub "
``rawSub(idx, size) -> Sequence``  

Return a subsequence starting at the given position index ``idx`` and of the
given ``size``, treating io-sequences as raw sequences.  

Negative indexing is supported, and if ``size`` is negative, a reverse sequence
starting at the ``idx`` is returned.  

The given ``idx`` must be a valid position index (i.e., not out of bounds),
unless the requested ``size`` is zero, in which case always a ``Sequence`` of
size zero is returned.  
";

%feature("docstring") tom::Sequence::nInputSymbols "
``nInputSymbols() -> Symbol``  

Return the size of the input alphabet.  

If this is zero, then this is an ordinary sequence, else an io-sequence.  
";

%feature("docstring") tom::Sequence::isReversed "
``isReversed() -> bool``  

Return ``true`` if this is a reversed sequence.  

This is relevant for io-sequences, since the order of input and outputs is then
also reversed.  
";

%feature("docstring") tom::Sequence::hasPrefix "
``hasPrefix(sequence, withSameAlphabet=false) -> bool``  

Return ``true`` if the given ``sequence`` is a prefix of this sequence.  

If ``withSameAlphabet`` is set to ``true`` (default ``false``), then a prefix
must have the same alphabet.  
";

%feature("docstring") tom::Sequence::fromJSON "
``fromJSON(string)``  
";

%feature("docstring") tom::Sequence::o "
``o(idx) -> Symbol``  
``o(idx, o)``  

Overloaded function
-------------------
* ``o(idx) -> Symbol``  
    
    Return the output symbol at indx ``idx``, where each index covers one io-
    pair.  

    Negative indexing is supported.  

* ``o(idx, o)``  
    
    Set the output symbol at index ``idx`` to ``o``, where each index covers one
    io-pair.  

    Negative indexing is supported.  
";

%feature("docstring") tom::Sequence::incr_as_python_iterator_only "
``incr_as_python_iterator_only()``  

Increment the first position, which allows using a ``Sequence`` as an iterator.  
";

%feature("docstring") tom::Sequence::isFrontAligned "
``isFrontAligned() -> bool``  

Return ``true`` if this sequence is io-aligned with respect to its beginning in
the underlying data, i.e., if either:  

*   this is a plain (non-io) sequence  
*   this io-sequence is not reversed and begins with an input symbol  
*   this io-sequence is reversed and ends with an input symbol  
";

%feature("docstring") tom::Sequence::reverse "
``reverse() -> Sequence``  

Return the reverse view of this ``Sequence``.  
";

%feature("docstring") tom::Sequence::isBackAligned "
``isBackAligned() -> bool``  

Return ``true`` if this sequence is io-aligned with respect to its end in the
underlying data, i.e., if either:  

*   this is a plain (non-io) sequence  
*   this io-sequence is not reversed and ends with an output symbol  
*   this io-sequence is reversed and begins with an output symbol  
";

%feature("docstring") tom::Sequence::Sequence "
``Sequence(symbol_list, nO, nU=0)``  
``Sequence(length=0, nO=0, nU=0)``  
``Sequence(json_representation)``  

Overloaded function
-------------------
* ``Sequence(symbol_list, nO, nU=0)``  
    
    Construct a ``Sequence`` of given output alphabet size ``nO`` and input
    alphabet size ``nU`` from the given ``symbol_list``.  

    This may be a list ``[u_0, ..., u_{N-1}]`` (or ``[u_0, o_0, ..., u_{N-1},
    o_{N-1}]`` for an io-sequence), or a ``std::vector<int>``. The contents of
    ``symbol_list`` is copied.  

* ``Sequence(length=0, nO=0, nU=0)``  
    
    Construct a ``Sequence`` of given ``length``, output alphabet size ``nO``
    and input alphabet size ``nU`` initialized with zeros.  

    In the case of an io-sequence (if ``nU != 0``), the ``rawSize()`` will be 2
    * ``length``.  

* ``Sequence(json_representation)``  
    
    Construct a ``Sequence`` corresponding to the given string
    ``json_representation``.  

    The format must correspond to what ``toJSON()`` produces.  
";

%feature("docstring") tom::Sequence::repr "
``repr() -> std::string``  

Return a string representation to display in python.  
";

%feature("docstring") tom::Sequence::copy "
``copy() -> Sequence``  

Return a deep copy of this ``Sequence``, i.e., the copy will use its own memory.  
";

%feature("docstring") tom::Sequence::sub "
``sub(idx, length) -> Sequence``  
``sub(length) -> Sequence``  

Overloaded function
-------------------
* ``sub(idx, length) -> Sequence``  
    
    Return a subsequence starting at the given position index ``idx`` and of the
    given ``length``, where each index covers one io-pair.  

    Negative indexing is supported, and if ``length`` is negative, a reverse
    sequence starting at the ``idx`` is returned.  

    The given ``idx`` must be a valid position index (i.e., not out of bounds),
    unless the requested ``length`` is zero, in which case always a ``Sequence``
    of length zero is returned.  

* ``sub(length) -> Sequence``  
    
    Return a subsequence of the given ``length``, starting at the beginning of
    this sequence if the given ``length`` is positive, else a reverse substring
    starting at the end of this sequence.  
";

// File: classtom_1_1_sequence_data.xml


%feature("docstring") tom::SequenceData "
``SequenceData(nO=0, nU=0)``  
``SequenceData(seq, nO=0, nU=0)``  
``SequenceData(size, nO=0, nU=0)``  

Constructors
------------
* ``SequenceData(nO=0, nU=0)``  

* ``SequenceData(seq, nO=0, nU=0)``  

* ``SequenceData(size, nO=0, nU=0)``  

Attributes
----------
* ``nO_`` : ``Symbol``  
    The size of the output alphabet.  

* ``nU_`` : ``Symbol``  
    The size of the input alphabet, or 0 if there are no inputs.  

* ``seq_`` : ``std::vector< Symbol >``  
    the underlying sequence data  

C++ includes: Sequence.h
";

%feature("docstring") tom::SequenceData::SequenceData "
``SequenceData(nO=0, nU=0)``  
``SequenceData(seq, nO=0, nU=0)``  
``SequenceData(size, nO=0, nU=0)``  

Overloaded function
-------------------
* ``SequenceData(nO=0, nU=0)``  

* ``SequenceData(seq, nO=0, nU=0)``  

* ``SequenceData(size, nO=0, nU=0)``  
";

%feature("docstring") tom::SequenceData::cereal::access "
``cereal::access() -> friend class``  
";

// File: classtom_1_1_estimator_1_1_state.xml

// File: classtom_1_1_stop_condition.xml


%feature("docstring") tom::StopCondition "
``StopCondition(maxIterations=100, relativeImprovementThreshold=1e-7, absoluteImprovementThreshold=1e-12)``  

This class allows to control the stopping condition for iterative algorithms,
and at the same time provides a callback mechanism between iterations.  

Constructors
------------
* ``StopCondition(maxIterations=100, relativeImprovementThreshold=1e-7, absoluteImprovementThreshold=1e-12)``  
    
    Construct a ``StopCondition`` while setting the key parameters.  

Attributes
----------
* ``maxIterations_`` : ``int``  
    the maximum number of iterations to perform  

* ``relativeImprovementThreshold_`` : ``double``  
    stop the iteration if the relative improvement falls below this threshold  

* ``absoluteImprovementThreshold_`` : ``double``  
    stop the iteration if the absolute improvement falls below this threshold  

* ``iteration_`` : ``int``  
    the current iteration  

* ``lastValue_`` : ``double``  
    the most recent \"value\" achieved by the iteration  

C++ includes: StopCondition.h
";

%feature("docstring") tom::StopCondition::stop "
``stop(currentValue) -> bool``  

Return ``true`` if the iteration should be terminated.  

This method gets called before starting a new iteration. The termination
condition can be customized by inheriting from ``StopCondition`` and overwriting
this method (even from Python).  

By default, the iteration is terminated if either of the following conditions is
met:  

*   the ``iteration_`` count reaches ``maxIterations_``  
*   the relative improvement |(``currentValue`` - ``lastValue_``) /
    ``lastValue_``| is less than the ``relativeImprovementThreshold_``  
*   the absolute improvement |(``currentValue`` - ``lastValue_``)| is less than
    the ``absoluteImprovementThreshold_``  
";

%feature("docstring") tom::StopCondition::reset "
``reset()``  

Reset the ``iteration_`` count and the ``lastValue_``.  
";

%feature("docstring") tom::StopCondition::StopCondition "
``StopCondition(maxIterations=100, relativeImprovementThreshold=1e-7, absoluteImprovementThreshold=1e-12)``  

Construct a ``StopCondition`` while setting the key parameters.  
";

%feature("docstring") tom::StopCondition::~StopCondition "
``~StopCondition()``  
";

%feature("docstring") tom::StopCondition::callback "
``callback()``  

This method gets called before every iteration and also right before
termination.  

It is a no-op by default, but can be used as a callback between iterations by
inheriting from ``StopCondition`` and overwriting this method (even from
Python).  
";

// File: classstree_1_1_s_tree.xml


%feature("docstring") stree::STree "
``STree(sequence)``  

An implementation of suffix trees for the ``Sequence`` type, which can represent
a plain or io-sequence.  

Features
--------  

*   Sequences are not required to end in a unique termination symbol:
    -   Every suffix is guaranteed to end in a leaf *or* internal node. The
        suffixes ending in internal nodes are given by ``internalSuffixes()``.  
    -   ``.count()`` always means the number of occurrences as a subsequence in
        the represented sequence. This is *not* the same as the leaf count.  
*   ``extendTo(seq)`` can be used to efficiently \"extend\" a suffix tree
    representation for a sequence ``s`` to a suffix tree for the sequence
    ``seq`` if ``s`` is a prefix of ``seq``.  
*   Io-sequences are basically treated as raw sequences, but only every second
    suffix (that is io-aligned) is represented in the suffix tree. Therefore:
    -   Occurrence counts are supported also for io-sequences of the form
        $u_0o_0, ..., u_k$ (not back aligned).  
    -   The suffix of $u_0o_0...u_ko_k$ is considered to be $u_1o_1...u_ko_k$.  
*   Sequences of size up to $2^29 - 1$ are supported.  
*   The space requirement is at most 24 bytes / symbol, i.e., less than 12 GB
    for the maximum supported sequence size.  

Details
-------  
In the following comes a brief description of the internal structure of this
suffix tree implementation:  

*   The suffix tree has two types of nodes (internal and leaf nodes), which are
    stored in vectors ``nodes`` and ``leaves`` respectively.  
*   Nodes are addressed by a 32-bit ``nidx_t`` (\"node index\"). The first three
    bits of a ``nidx_t`` indicate whether the ``nidx_t`` is valid (1) or nil
    (0), addresses an (internal) node (1) or a leaf (0), and is \"colored\" (1)
    or not (0), respectively. In the case of a valid ``nidx_t``, the remaining
    29 bits form the index value of the addressed node in the ``nodes`` or
    ``leaves`` vector. Otherwise, these can have some other significance.  
*   Every node (internal or leaf) has a left and right ``nidx_t``. Internal
    nodes have an additional child ``nidx_t``, which addresses the first child
    in the suffix tree structure.  
*   The children (which are siblings) of any internal node are basically
    organized in a self-balancing binary tree structure formed by the left and
    right ``nidx_t`` \"pointers\". For this, a left-leaning red-black tree is
    used (which uses the color flag of the ``nidx_t``). However, the first two
    children of any node have special roles:
    -   The right ``nidx_t`` of the first child encodes the headindex of the
        parent node, while the right ``nidx_t`` of the second child encodes the
        depth of the parent node. Both right ``nidx_t`` addresses of the first
        two children have their first three bits set to zero.  
    -   The left ``nidx_t`` of the first child addresses the second child and
        has its color flag set. The left ``nidx_t`` of the second child
        addresses the root of the red-black tree in which all further siblings
        are organized, and has its color flag unset (to be able to distinguish
        first and second children). However, if there are only two siblings,
        then the left ``nidx_t`` of the second child will be marked as invalid,
        but will address the suffix-link of the parent. Still, the color flag
        will be unset.  
*   Suffix-links of any (internal) node are stored in the right ``nidx_t`` of
    the rightmost (in the binary tree) child, or alternatively in the left
    ``nidx_t`` of the second child if there are only two children. This
    ``nidx_t`` will always be marked as invalid and uncolored, but as addressing
    a node.  
*   The left ``nidx_t`` of the leftmost sibling in the binary tree will always
    be marked as invalid and uncolored, but has no further meaning.  
*   The binary tree of siblings is threaded. This means that invalid left and
    right ``nidx_t`` entries address previous and next siblings respectively.
    This is true for all siblings in the red-black tree, except for the left-and
    rightmost. To distinguish these, all invalid ``nidx_t`` that indicate a
    thread are marked as colored (as opposed to the left-and rightmost, which
    are always marked as uncolored).  
*   Siblings are stored in lexicographic order with respect to their edge
    labels. It is ensured that the first and second siblings are always the
    lexicographically smallest.  

Constructors
------------
* ``STree(sequence)``  
    
    Create a suffix tree for the given ``sequence``.  

    Note that the sequence must have size at least one.  

C++ includes: STreeCore.h
";

/*
 Node data manipulation 
*/

%feature("docstring") stree::STree::nLeafNodes "
``nLeafNodes() -> nidx_t``  

Return the number of leaf nodes in the suffix tree.  
";

%feature("docstring") stree::STree::extendTo "
``extendTo(sequence, checkExtendability=true)``  

Extend the current suffix tree representation to a representation for a given
``sequence``, which requires that the current suffix tree represents a prefix of
the given ``sequence``.  
";

%feature("docstring") stree::STree::nNodes "
``nNodes() -> nidx_t``  

Return the number of nodes (internal and leaves) in the suffix tree.  
";

%feature("docstring") stree::STree::nInternalNodes "
``nInternalNodes() -> nidx_t``  

Return the number of internal nodes in the suffix tree.  
";

%feature("docstring") stree::STree::internal::Pos "
``internal::Pos() -> friend class``  
";

%feature("docstring") stree::STree::sequence "
``sequence() -> const Sequence``  

Return the represented sequence.  
";

%feature("docstring") stree::STree::internal::RBTreeNodeTraits "
``internal::RBTreeNodeTraits() -> friend class``  
";

%feature("docstring") stree::STree::STree "
``STree(sequence)``  

Create a suffix tree for the given ``sequence``.  

Note that the sequence must have size at least one.  
";

%feature("docstring") stree::STree::deepestInternalSuffixNidx "
``deepestInternalSuffixNidx() -> nidx_t``  

Return the node index (``nidx_t``) of the deepest internal node in the suffix
tree representing a suffix of the represented sequence.  

This can and should be converted to a ``Node`` by calling ``Node(stree,
stree.deepestInternalSuffix())``. [This is required for technical reasons of
memory management safety].  

If the input sequence terminates with a unique symbol, then this will always be
the root of the suffix tree, corresponding to the empty suffix. If the input
sequence does not terminate with a unique symbol, this corresponds to the
\"active position\" in the suffix tree construction, which will be an internal
node. The other internal nodes representing a suffix can be found by traversing
the suffix link (calling ``toSuffix()`` on the converted returned ``Node``)
until reaching the root node (to exclude the empty suffix) or until
``toSuffix()`` results in an invalid ``Node`` (to include the empty suffix). The
remaining (longer) suffixes correspond to the leaf nodes.  
";

// File: namespacecereal.xml

%feature("docstring") cereal::load "
``load(ar, m)``  
``load(ar, m)``  

Overloaded function
-------------------
* ``load(ar, m)``  

* ``load(ar, m)``  
";

%feature("docstring") cereal::save "
``save(ar, m)``  
``save(ar, m)``  

Overloaded function
-------------------
* ``save(ar, m)``  

* ``save(ar, m)``  
";

// File: namespace_eigen.xml

// File: namespacestd.xml

// File: namespacestree.xml

// File: namespacestree_1_1internal.xml

// File: namespacetom.xml

%feature("docstring") tom::normalize "
``normalize(matrix) -> double``  

Devide the given ``matrix`` by its element-sum, i.e., normalize the matrix to
have an element-sum of one, and return the element-sum.  
";

%feature("docstring") tom::kron "
``kron(A, B) -> MatrixXd``  

Return the Kronecker-product :math:`A\\otimes B` of the matrices ``A`` and
``B``.  
";

%feature("docstring") tom::solveLS "
``solveLS(A, M, W=MatrixXd(), transposed=false, method=\"LDLT\") -> MatrixXd``  

Return the generalized weighted least-squares (GLS) solution to the
overdetermined problem ``A`` * ``X`` = ``M`` (or to the problem ``X`` * ``A`` =
``M`` if ``transposed``) with the given weights ``W`` using a ``method`` from
{\"Cholesky\", \"LDLT\" (default), \"QR\", \"SVD\", \"JacobiSVD\"}.  

This computes ``X`` that minimizes the weighted norm |``A`` * ``X`` - ``M``|_W
(or |``X`` * ``A`` - ``M``|_W if ``transposed``), utilizing the structure of W
that depends on the size of the supplied weights ``W`` in the following way,
assuming ``M`` is of size m x n:  

*   If ``W`` is of size zero (default), then no weights are used and the
    ordinary least squares (OLS) solution minimizing the Frobenius norm is
    returned.  
*   If ``W`` is of size m x n, then element-wise weights are assumed, i.e., W =
    D(``W``), resulting in the weighted least squares (WLS) solution.  
*   If ``W`` is of size m x nm, then a block-diagonal structure for W is
    assumed, i.e., W = D(``W``_1, ..., ``W``_m), where ``W``_j is the j-th (m x
    m)-block of ``W`` corresponding to the weights for [``M``]_j, which must be
    symmetric an positive definite.  
*   If ``W`` is of size m x 1, then these are regarded as row-weights, which
    only make sense if not ``transposed``.  
*   If ``W`` is of size 1 x n, then these are regarded as column-weights, which
    only make sense if ``transposed``.  
*   If ``W`` is of size m+n x 1 then ``W``[:m] are regarded as row-weights and
    ``W``[m:] as column-weights. If ``transposed`` the row-weights, else the
    column-weights, have no effect.  

Note that solving GLS with full weight matrix is expensive, and therefore only
solving with block-diagonal structured weight matrix is supported, and then only
the methods \"Cholesky\" and \"LDLT\" which can make use of extra
simplifications.  

The computation is done by reducing the problem to a set of OLS problems that
are then solved according to the given ``method`` as detailed below. Note that
the weights in ``W`` must be strictly greater than zero.  

The \"Cholesky\" method solves the normal equations using a Cholesky
decomposition. This is the fastest method, but loses most precision and requires
the problem to be overdetermined and ``A`` to have full rank.  

The \"LDLT\" method is essentially the same as \"Cholesky\", but uses a more
robust Cholesky decomposition with pivoting that also avoids taking a square
root. This method is recommended over \"Cholesky\" by Eigen3.  

The \"QR\" method uses a QR decomposition. This is slower than \"Cholesky\", but
gives more precision. The marix ``A`` should have full rank.  

The \"SVD\" uses an SVD decomposition. This is the slowest, but gives best
precision. Also, the matrix ``A`` does not need to have full rank in which case
the least-squares solution with the smallest norm is returned.  

The \"JacobiSVD\" method is similar to the \"SVD\" method, but uses a different
(slower, but potentially more accurate) svd algorithm.  
";

%feature("docstring") tom::normalizeCols "
``normalizeCols(matrix) -> bool``  

Devide each column of the given ``matrix`` by its sum, i.e., normalize the
columns to have column-sum one.  

Return ``true`` if successful, or ``false`` if a column could not be normalized
due to a zero column-sum.  
";

%feature("docstring") tom::solveRowColWLS "
``solveRowColWLS(A, M, W, transposed=false, method=\"LDLT\") -> MatrixXd``  

Return the row or column weighted least-squares solution to the problem ``A`` *
``X`` = ``M`` with row-weights given in the column vector ``W`` (or if
``transposed`` to ``X`` * ``A`` = ``M`` with column-weights given in the row
vector ``W``) using a ``method`` from {\"Cholesky\", \"LDLT\" (default), \"QR\",
\"SVD\", \"JacobiSVD\"}.  

This computes ``X`` that minimizes |D(sqrt_W) * (``A`` * ``X`` - ``M``)|_F (or
|(``X`` * ``A`` - ``M``) * D(sqrt_W)|_F if ``transposed``), where ``sqrt_W`` is
the element-wise square-root of ``W``, i.e., ``W`` = ``sqrt_W`` .* ``sqrt_W``,
and ``.*`` denotes the element-wise product. The computation is done by reducing
the problem to an OLS problem that is then solved according to the given
``method`` as detailed below (see also ``solveOLS()``). Note that the weights in
``W`` must be strictly greater than zero.  

Note that column weights have no effect in the default case, and row weight have
no effect if ``transposed``, and are therefore ommitted.  

The \"Cholesky\" method solves the normal equations using a Cholesky
decomposition. This is the fastest method, but loses most precision and requires
the problem to be overdetermined and ``A`` to have full rank.  

The \"LDLT\" method is essentially the same as \"Cholesky\", but uses a more
robust Cholesky decomposition with pivoting that also avoids taking a square
root. This method is recommended over \"Cholesky\" by Eigen3.  

The \"QR\" method uses a QR decomposition. This is slower than \"Cholesky\", but
gives more precision. The marix ``A`` should have full rank.  

The \"SVD\" uses an SVD decomposition. This is the slowest, but gives best
precision. Also, the matrix ``A`` does not need to have full rank, and in the
case of an underdetermined problem, the least-squares solution with the smallest
norm is returned.  

The \"JacobiSVD\" method is similar to the \"SVD\" method, but uses a different
(slower, but potentially more accurate) svd algorithm.  
";

%feature("docstring") tom::pinv "
``pinv(M, method=\"SVD\") -> MatrixXd``  

Return the pseudo-inverse of the given matrix ``M`` computed according to the
given ``method`` from {\"Cholesky\", \"QR\", \"SVD\" (default), \"JacobiSVD\"}.  

If ``method`` is \"SVD\" or \"JacobiSVD\", the classical pseudo-inverse is
computed from the svd of ``M``.  

If ``method`` is \"QR\", the pseudo-inverse is computed from the QR-
factorization of ``M``, which requires ``M`` to have full rank.  

If ``method`` is \"Cholesky\" or \"LDLT\", the pseudo-inverse is computed as
:math:`(M^\\top M)^{-1} M^\\top` or :math:`M^\\top (M M^\\top)^{-1}` depending
on the size of ``M``, which requires ``M`` to have full rank.  
";

%feature("docstring") tom::reverseWords "
``reverseWords(words)``  

Reverse the given ``words`` in-place.  
";

%feature("docstring") tom::solveWLS "
``solveWLS(A, M, W, transposed=false, method=\"LDLT\") -> MatrixXd``  

Return the (element-wise) D(``W``)-weighted least-squares (WLS) solution to the
problem ``A`` * ``X`` = ``M`` (or to ``X`` * ``A`` = ``M`` if ``transposed``)
using a ``method`` from {\"Cholesky\", \"LDLT\" (default), \"QR\", \"SVD\",
\"JacobiSVD\"}.  

This computes ``X`` that minimizes |``A`` * ``X`` - ``M``|_D(``W``) (or |``X`` *
``A`` - ``M``|_D(``W``) if ``transposed``). The computation is done by reducing
the problem to a set of OLS problems that are then solved according to the given
``method`` as detailed below (see also ``solveOLS()``). Note that the weights in
``W`` must be strictly greater than zero.  

The \"Cholesky\" method solves the normal equations using a Cholesky
decomposition. This is the fastest method, but loses most precision and requires
the problem to be overdetermined and ``A`` to have full rank.  

The \"LDLT\" method is essentially the same as \"Cholesky\", but uses a more
robust Cholesky decomposition with pivoting that also avoids taking a square
root. This method is recommended over \"Cholesky\" by Eigen3.  

The \"QR\" method uses a QR decomposition. This is slower than \"Cholesky\", but
gives more precision. The marix ``A`` should have full rank.  

The \"SVD\" uses an SVD decomposition. This is the slowest, but gives best
precision. Also, the matrix ``A`` does not need to have full rank, and in the
case of an underdetermined problem, the least-squares solution with the smallest
norm is returned.  

The \"JacobiSVD\" method is similar to the \"SVD\" method, but uses a different
(slower, but potentially more accurate) svd algorithm.  
";

%feature("docstring") tom::hmmToOom "
``hmmToOom(T, E, w, transition_first=false) -> SHARED_PTR< Oom >``  
";

%feature("docstring") tom::colwiseMean "
``colwiseMean(matrix, p=1.0) -> RowVectorXd``  

Return the column-wise generalized mean with exponent ``p`` (default 1) of the
given ``matrix``.  

For ``p`` = 1, 0, -1 this is the arithmetic, geometric and harmonic mean,
respectively.  

Note that for values of ``p`` other than {1, 2k} this requires all matrix
entries to be positive.  
";

%feature("docstring") tom::computeWLRA "
``computeWLRA(M, W, B_init, stopCondition=StopCondition(50, 1e-5, 1e-12), method=\"Cholesky\") -> tuple< MatrixXd, MatrixXd >``  

Return in a tuple [B,A] (an approximation to) the best weighted rank-d
approximation to ``M`` with element-wise weights ``W`` such that |``B`` * ``A``
- ``M``|_D(vect(W)) is minimized.  

This is computed iteratively starting from an initial approximation given by
``B_init`` using \"alternating projections\", which are in turn solved via the
given ``method``. The termination of the iteration is controlled by the given
``stopCondition``. See ``tom.util.StopCondition``.  
";

%feature("docstring") tom::sortWords "
``sortWords(words)``  

Sort the given ``words`` in-place by length and then lexicographically.  
";

%feature("docstring") tom::sharpenEfficiency "
``sharpenEfficiency(oom, rStree, indNodes) -> std::shared_ptr< Oom >``  
";

%feature("docstring") tom::improveWLRA "
``improveWLRA(B, A, M, W, stopCondition=StopCondition(50, 1e-5, 1e-12), method=\"Cholesky\")``  

Compute in the arguments ``B`` and ``A`` (an approximation to) the best weighted
rank-d approximation to ``M`` with element-wise weights ``W`` such that |``B`` *
``A`` - ``M``|_D(W) is minimized.  

This is computed iteratively starting from an initial approximation given by
``B`` * ``A`` using \"alternating projections\" solved via the given ``method``.
The termination of the iteration is controlled by the given ``stopCondition``.
See ``tom.util.StopCondition``.  
";

%feature("docstring") tom::wordsOverAlphabet "
``wordsOverAlphabet(nOutputSymbols, nInputSymbols=0, minLength=1, maxLength=1) -> std::shared_ptr< Sequences >``  

Return in lexicographic order all words of length between ``minLength`` and
``maxLength`` over the alphabet with ``nOutputSymbols`` output symbols if
``nInputSymbols`` is zero, or otherwise over the alphabet of input-output symbol
pairs with ``nOutputSymbols`` output symbols and ``nInputSymbols`` input
symbols.  

Parameters
----------
* ``nOutputSymbols`` :  
    the number of output symbols  
* ``nInputSymbols`` :  
    the number of input symbols (default 0)  
* ``minLength`` :  
    the minimum length for returned words (default 1)  
* ``maxLength`` :  
    the maximum length for returned words (default 1)  

Returns
-------
an array of words  
";

%feature("docstring") tom::wordsFromData "
``wordsFromData(dataSuffixTree, minLength=0, maxLength=0, minRelevance=1.0, maxWords=0, prefixUnique=false, suffixUnique=false, relevance=stree::PositionRelevance()) -> std::shared_ptr< Sequences >``  

Find words occurring in a sequence according to given constraints.  

This can be used to find sets of indicative or characteristic words from
training data that are appropriate for the OOM learning algorithms.  

From a suffix tree representation of the data sequence given by
``dataSuffixTree``, this function efficiently computes an array of words with
length between ``minLength`` and ``maxLength`` that occur at least ``minCount``
times as subsequences in the data sequence, sorted by occurrence count and then
lexicographically. If ``uniquePositions`` is ``true``, then only the first
(shortest) word is retained for words occurring at the same set of positions in
the data sequence. Finally, only the first at most ``maxWords`` resulting words
are returned.  

Parameters
----------
* ``dataSuffixTree`` :  
    a suffix tree representation of the data sequence  
* ``minLength`` :  
    the minimum length for returned words (default 0)  
* ``maxLength`` :  
    the maximum length for returned words, or ``0`` (default) for no limit  
* ``minRelevance`` :  
    the minimum \"relevance\" in the data sequence for returned words (default
    1.0)  
* ``maxWords`` :  
    the maximum number of returned words, or ``0`` (default) for no limit  
* ``prefixUnique`` :  
    if set to ``true``, then only retain the shortest prefix among words
    occurring at the same set of positions in the data sequence (default
    ``false``)  
* ``suffixUnique`` :  
    if set to ``true``, then only retain the shortest suffix among words ending
    at the same set of positions in the data sequence (default ``false``)  
* ``relevance`` :  
    a ``PositionRelevance`` object to compute the *relevance value* for words,
    which defaults to their occurrence count in the data (plus a fraction to
    penalize longer words). This can be customized by passing a ``relevance``
    inherited from ``PositionRelevance`` with an overwritten ``.compute()``
    method.  

Returns
-------
an array of words occurring in the data sequence according to the given
constraints, and sorted by occurrence count (descending) and then
lexicographically.  

Notes
-----  
If two words a and b with |a| < |b| always occur at the same positions in the
data sequence, then a must be a prefix of b, i.e., b = ac for some non-empty
word c. Furthermore, this means that for any word x, the occurrence counts for
the words xa and xb will be the same, and moreover these counts will be based on
exactly the same occurrences in the data sequence. When collecting such
occurrence statistics, it may therefore be useful to use only one (e.g., the
shortest) of the words with identical occurrence positions. This is what
``uniquePositions`` is for.  

Note that in the context of OOM learning, the above means that
``uniquePositions`` makes sense when finding **characteristic** words from the
training data. To achieve the same (avoid exact duplicate data exploitation) for
**indicative** words, one needs to use this function with ``uniquePositions``
for a **reversed** training sequence, and then reverse each word in the
resulting set.  
";

%feature("docstring") tom::transformWeights "
``transformWeights(W, B, covariances=true) -> MatrixXd``  

Return a new weight matrix for ``X``, assuming ``X`` is a solution to the
D(``W``)-weighted WLS problem ``B`` * ``X`` = ``M``.  

Note that the columns of ``X`` can be regarded as coordinate representations for
the columns of ``M`` with respect to a basis given by the columns of ``B``. This
function transforms the given weights for the columns of ``M`` to appropriate
weights for the coordinates in the columns of ``X``. The resulting weight matrix
for ``X`` will therefore be block-diagonal in general, but if ``covariances`` is
set to ``false``, the off-diagonal weights are ignored, resulting in element-
wise weights for ``X``.  

The returned matrix will therefore be  

*   [B^T * D([W]_1) * B, ..., B^T * D([W]_m) * B] of size B.cols() x B.cols *
    M.cols() if ``covariances``  
*   [diag(B^T * D([W]_1) * B), ..., diag(B^T * D([W]_m) * B)] of size B.cols() x
    M.cols() otherwise.  
";

%feature("docstring") tom::getIndicativeSequenceNodes "
``getIndicativeSequenceNodes(reverseDataSuffixTree, minIndCount, maxIndLen) -> std::shared_ptr< std::vector< stree::nidx_t > >``  
";

%feature("docstring") tom::solveGLS "
``solveGLS(A, M, W, transposed=false, method=\"LDLT\") -> MatrixXd``  

Return the D(W1,..., Wm)-weighted least-squares (GLS) solution to the
overdetermined problem ``A`` * ``X`` = ``M`` (or to ``X`` * ``A`` = ``M`` if
``transposed``) using a ``method`` from {\"Cholesky\", \"LDLT\" (default)},
where the block-diagonal symmetric and positive definite weight matrix is given
by ``W`` = [W1,..., Wn], where each ``Wj`` is the full weight matrix for the
column j of ``M``.  

This computes ``X`` that minimizes |``A`` * ``X`` - ``M``|_D(W1,...,Wn) (or
|``X`` * ``A`` - ``M``|_D(W1,...,Wn) if ``transposed``).  

Note that the \"LDLT\" method is essentially the same as \"Cholesky\", but uses
a more robust Cholesky decomposition with pivoting that also avoids taking a
square root. This method is recommended over \"Cholesky\" by Eigen3.  
";

%feature("docstring") tom::solveOLS "
``solveOLS(A, M, transposed=false, method=\"QR\") -> MatrixXd``  

Return the ordinary least-squares (OLS) solution to the problem ``A`` * ``X`` =
``M`` (or if ``transposed`` to ``X`` * ``A`` = ``M``) using a ``method`` from
{\"Cholesky\", \"LDLT\", \"QR\" (default), \"SVD\", \"JacobiSVD\"}.  

The \"Cholesky\" method solves the normal equations using a Cholesky
decomposition. This is the fastest method, but loses most precision and requires
the problem to be overdetermined and ``A`` to have full rank.  

The \"LDLT\" method is essentially the same as \"Cholesky\", but uses a more
robust Cholesky decomposition with pivoting that also avoids taking a square
root. This method is recommended over \"Cholesky\" by Eigen3.  

The \"QR\" method uses a QR decomposition. This is slower than \"Cholesky\", but
gives more precision. The marix ``A`` should have full rank.  

The \"SVD\" uses an SVD decomposition. This is the slowest, but gives best
precision. Also, the matrix ``A`` does not need to have full rank, and in the
case of an underdetermined problem, the least-squares solution with the smallest
norm is returned.  

The \"JacobiSVD\" method is similar to the \"SVD\" method, but uses a different
(slower, but potentially more accurate) svd algorithm.  
";

%feature("docstring") tom::normalizeRows "
``normalizeRows(matrix) -> bool``  

Devide each row of the given ``matrix`` by its sum, i.e., normalize the rows to
have row-sum one.  

Return ``true`` if successful, or ``false`` if a row could not be normalized due
to a zero row-sum.  
";

%feature("docstring") tom::weightedNorm "
``weightedNorm(M, W, squared=false) -> double``  

Return the weighted norm of ``M`` with weights given in ``W``, or the squared
weighted norm if ``squared`` is set to ``true``.  

Depending on the size of ``W``, the given weights are interpreted in different
ways, assuming ``M`` is of size m x n:  

*   if ``W`` is of size zero, then no weights are used and the Frobenius norm
    |M|_F is computed  
*   if ``W`` is of size m+n x 1, then row and column weights [w_r; w_c] = W are
    assumed and |M|_D(w_r w_c^T) is computed  
*   if ``W`` is of size m x n, then element-wise weights are assumed and
    |M|_D(W) is computed  
*   if ``W`` is of size m x mn, then a block-diagonal weight matrix is assumed
    and |M|_D(W1,...,Wn) is computed  
*   if ``W`` is of size mn x mn, then a full weight matrix is assumed and |M|_W
    is computed  
";

%feature("docstring") tom::wordsFromModel "
``wordsFromModel(oom, minLength=0, maxLength=0, minProbability=1e-5, maxWords=0) -> std::shared_ptr< Sequences >``  

Return all words satisfying the given constraints sorted (descending) by their
probability according to the given ``oom``.  

Parameters
----------
* ``oom`` :  
    the ``Oom`` model from which to compute the word probabilities  
* ``minLength`` :  
    the minimum length for returned words (default 0)  
* ``maxLength`` :  
    the maximum length for returned words, or ``0`` (default) for no limit  
* ``minProbability`` :  
    the minimum probability for returned words (default 1e-5)  
* ``maxWords`` :  
    the maximum number of returned words, or ``0`` (default) for no limit  
";

%feature("docstring") tom::rowwiseMean "
``rowwiseMean(matrix, p=1.0) -> VectorXd``  

Return the row-wise generalized mean with exponent ``p`` (default 1) of the
given ``matrix``.  

For ``p`` = 1, 0, -1 this is the arithmetic, geometric and harmonic mean,
respectively.  

Note that for values of ``p`` other than {1, 2k} this requires all matrix
entries to be positive.  
";

// File: _cereal_tom_8h.xml

// File: _core_sequences_8h.xml

// File: _efficiency_sharpening_8h.xml

// File: _estimator_8h.xml

// File: _hmm_8h.xml

// File: _implementations_8h.xml

// File: _linear_algebra_8h.xml

// File: _oom_8h.xml

// File: _policy_8h.xml

// File: _random_8h.xml

// File: _sequence_8h.xml

// File: _stop_condition_8h.xml

// File: _r_b_tree_8h.xml

// File: stree_8h.xml

// File: _s_tree_core_8h.xml

// File: _s_tree_iterators_8h.xml

// File: _s_tree_node_8h.xml

// File: tom_8h.xml

// File: _learner_8h.xml

// File: _oom_tools_8h.xml

// File: _r_e_a_d_m_e_8md.xml

// File: dir_d44c64559bbebec7f509842c48db8b23.xml

// File: dir_12d8770726183bfd67fa25f0fa619395.xml

// File: dir_dcf3794d532a9ac5f87cd327a4d01b07.xml

// File: dir_4a7ba64ee0f8d8f953f2007bf2110ebd.xml

// File: indexpage.xml

