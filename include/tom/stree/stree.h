/**
 * @file   stree.h
 * @brief  Include this file to use the suffix tree library.
 */

#ifndef STREE_H
#define STREE_H

#include <cassert>
#include <cstdint>

#include <vector>
#include <deque>
#include <stack> // for iterative dfs traversal and RBTree
#include <queue> // for iterative bfs traversal

#include <iostream>
#include <sstream>  // for output only
#include <iomanip>  // for output only

/**
 * Collects all the functionality of this suffix tree implementation.
 */
namespace stree {

//Nidx BitFlags:
typedef uint32_t Idx;                   // the index "pointer" type to be used
typedef uint32_t Nidx;

const Nidx VALID  =    2147483648u; // (1 << 31) = 2**31 (highest bit)
const Nidx NODE   =    1073741824u; // (1 << 30) = 2**30
const Nidx COLOR  =     536870912u; // (1 << 29) = 2**29
const Nidx INDEX  =     536870911u; // ~(7 << 29) = 2**29-1 (29 lowest bits)
const Nidx ROOT   =   VALID | NODE; // corresponds to the root node of the suffix tree

#ifdef STREE_STRING_TYPE
typedef STREE_STRING_TYPE String;
#else
typedef std::string String;
#endif // STREE_STRING_TYPE

typedef String::value_type Char;


// forward declarations
class STree;
class STreeNode;
class STreeEdge;
class STreePath;
class STreePos;

class PrefixIterator;
class PostfixIterator;
class DFSIterator;

} // namespace stree

#include "RBTree.h"
#include "STreeCore.h"
#include "STreeNode.h"
#include "STreeIterators.h"

#include "STreeCore.cpp"

#endif // STREE_H
