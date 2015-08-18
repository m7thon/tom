/**
 * @file   EfficiencySharpening.h
 * @author Michael Thon
 *
 * @brief  A reimplementation of the ES algorithm
 */

#ifndef EFFICIENCY_SHARPENING_H
#define EFFICIENCY_SHARPENING_H

#include "tom.h"

namespace tom {

std::shared_ptr<Oom> sharpenEfficiency(const Oom& oom, std::shared_ptr<const stree::STree> rStree, std::shared_ptr<std::vector<stree::nidx_t> > indNodes) throw (std::invalid_argument) {
    if (oom.nU() != 0) throw std::invalid_argument("sharpenEfficiency does not work for IO-OOMs");
  Sequence seq = rStree->sequence().rawSub(0, rStree->sequence().rawSize());
  auto room = oom.reverse();
  MatrixXf* CF_l = room->harvestStates(seq);
  // Note: CF_l is [c_0, ..., c_l] as in the paper, but differently ordered from matlab code

  // (3) Weed ( *** not implemented! (This is a weird thing to do anyway...) ***)

  // (4) Associate reverse states with suffixes (leaves) and sum depth-first through suffix tree
  MatrixXf CF_i = MatrixXf::Zero(oom.dim(), rStree->nInternalNodes());
  for (auto internalSuffixNode = stree::Node(rStree, rStree->deepestInternalSuffixNidx()); internalSuffixNode.isValid(); internalSuffixNode.toSuffix()) {
    CF_i.col(internalSuffixNode.index()) = CF_l->col(seq.length() - internalSuffixNode.depth());
  }
  for (auto it = stree::PostfixIterator(rStree); it.isValid(); it.toNext()) {
    // if (it.isLeaf()) CF_l.col(it.index()) = CF_l.col(it.headIndex());
    if (!(it.isLeaf())) {
      stree::Node child = it.child();
      while (child.isValid()) {
        CF_i.col(it.index()) += child.isLeaf() ? CF_l->col(child.index()) : CF_i.col(child.index());
        child.toSibling();
      }
    }
  }

  // Start initializing the newly learnt oom
  std::shared_ptr<Oom> newOom(new Oom());
	newOom->setSize(oom.dim(), seq.nO(), 0);
  newOom->sig() = RowVectorXd::Ones(oom.dim());

  // (5) Obtain normalized argument-value-pair matrices CF and CFz
  //     and compute tau operators
  MatrixXd CF = MatrixXd::Zero(oom.dim(), indNodes->size());
  for (int i = 0; i < indNodes->size(); ++i) {
    stree::Node indNode = stree::Node(rStree, indNodes->at(i));
    CF.col(i) = indNode.isLeaf() ? CF_l->col(indNode.index()).cast<double>() : CF_i.col(indNode.index()).cast<double>();
  }
  RowVectorXd colSums = CF.colwise().sum().cwiseSqrt();
  for (int i = 0; i < CF.cols(); ++i)
    if (colSums(i) != 0) CF.col(i) /= colSums(i);
  CF = pinv(CF);

  for (Symbol o = 0; o < seq.nO(); ++o) {
    MatrixXd CFz = MatrixXd::Zero(oom.dim(), indNodes->size());
    for (int i = 0; i < indNodes->size(); ++i) {
      stree::Position pos(rStree); // root
      pos.toSymbol(o); pos.toSequence(stree::Node(rStree, indNodes->at(i)).sequence());
      stree::Node node = pos.edge();
      if (node.isValid()) {
        CFz.col(i) = node.isLeaf() ? CF_l->col(node.index()).cast<double>() : CF_i.col(node.index()).cast<double>();
      }
      if (colSums(i) != 0) CFz.col(i) /= colSums(i);
    }
    newOom->tau(o) = (CFz * CF).cast<double>();
  }
  delete CF_l;

  newOom->w0() = newOom->stationaryState();
  newOom->init();

  // (6) Set stabilization parameters of new Oom to those of given one
  newOom->minPrediction_ = oom.minPrediction_;
  newOom->maxPredictionError_ = oom.maxPredictionError_;
  newOom->maxSetback(oom.maxSetback());
	newOom->epsilonZero_ = oom.epsilonZero_;
	newOom->reset();
  return newOom;
}


} // namespace tom

#endif // EFFICIENCY_SHARPENING_H
