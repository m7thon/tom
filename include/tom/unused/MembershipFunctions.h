#ifndef __MEMBERSHIP_FUNCTIONS_H__
#define __MEMBERSHIP_FUNCTIONS_H__

#include "../tom.h"

#include "RandTom.h"
#include <iostream>
#include <sstream>
#include <string>
#include <Eigen/Core>
#include <cmath>
#include <algorithm>

/** 
 * This class implements membership functions that are either piecewise linear of simply a single dirac delta function.
 */
class LinearMembershipFunction {
public:
  /** 
	 * constructs an uninitialized (!) LinearMembershipFunction object
	 */
	LinearMembershipFunction();

	/** 
	 * constructs a dirac delta membership function at \a x
	 */
	LinearMembershipFunction(double x);

	/** 
	 * setup a LinearMembershipFunction from the given \a points.
	 */
	template <typename D>
	void setPoints(const Eigen::MatrixBase<D>& points);

	/** 
	 * initialize the LinearMembershipFunction from given nP and points.
	 */
	void init();

	/** 
	 * evaluates this function at \a x
	 */
	inline double value(double x) const;

	/** 
	 * the mean of this membership function viewed as a distribution
	 */
	inline double mean() const;

	/** 
	 * sample from this membership function viewed as a distribution
	 */
	inline double sample() const;

	/** 
	 * returns the inner product of this membership function with another membership function \a f. I.e., if we denote this membership function by \a g, then \f$\langle f,g\rangle = \int f(x)g(x)\,dx\f$ is calculated.
	 */
	double operator*(const LinearMembershipFunction& f) const;
	
/** @name IO-functions */ //@{
	/** 
	 * read the parameters from the given input stream and initialize. The format must correspond to what the output functions produce.
	 */
	std::istream& operator<<(std::istream& is);

	/**
	 * write the parameters to the given output stream.
	 */
	std::ostream& operator>>(std::ostream& os) const;

	/**
	 * read the parameters from the given string. The format must correspond to what the output functions produce.
	 */
	void from_string(const std::string& str);

	/**
	 * output the parameters as a string.
	 */
	std::string to_string() const;
	
  /** 
   * write the continuous OOM parameters to std::cout
   */
	void show() const;
//@}

private:
	bool dirac; /**< true if this is a dirac delta function */
	double w; /**< the integral of this membership function */
	double m; /**< the mean of this membership function viewed as a distribution */
	int nP; /**< the number of pieces that this piecewise defined function consists of */
	Eigen::ArrayXXd P; /**< the data points that define this function. The data points of the form \f$(x, f(x))\f$ are sorted according to the \a x-coordinate. The entries of \a P have the following meaning:
												- \a P(i,0) is the \a x-coordinate of the i-th data point
												- \a P(i,1) is the value \f$f(x)\f$ of the i-th data point
												- \a P(i,2) is the slope of the i-th piece (i.e., \f$(P(i+1,1) - P(i,1)) / (P(i+1,0) - P(i,0))\f$, or 0, if this expression is undefined)
												- \a P(i,3) is the integral of the i-th piece */
}; // class LinearMembershipFunction

#ifdef SWIG
%template(setPoints) LinearMembershipFunction::setPoints<MatrixMd>;
#endif


#ifndef SWIG
/**
 * write the parameters of the LinearMembershipFunction \a f to the given output stream.
 */
std::ostream& operator<<(std::ostream& os, const LinearMembershipFunction& f);
/** 
 * read the parameters of the LinearMembershipFunction \a f from the given input stream and initialize. The format must correspond to what the output functions produce.
 */
std::istream& operator>>(std::istream& is, LinearMembershipFunction& f);
#endif // SWIG


#ifndef SWIG
// ##################################################################
//                         IMPLEMENTATION
// ##################################################################

LinearMembershipFunction::LinearMembershipFunction() {}

LinearMembershipFunction::LinearMembershipFunction(double x) {
	nP = 0;
	P.resize(1,4);
	P(0,0) = x;
	P(0,1) = 1;
	P(0,2) = 0;
	P(0,3) = 0;
	m = 0;
	w = 1;
}	

template <typename D>
void LinearMembershipFunction::setPoints(const Eigen::MatrixBase<D>& points) {
	nP = points.rows() - 1;
	P = MatrixXd::Zero(nP + 1, 4);
	P.block(0, 0, nP + 1, 2) = points; 
	init();
}

void LinearMembershipFunction::init() {
	// Calculate slopes:
	P(nP,2) = 0;
	for (int n = 0; n < nP; ++n)
		if (P(n+1,0) - P(n,0) == 0) P(n,2) = 0;
		else P(n,2) = (P(n+1,1) - P(n,1)) / (P(n+1,0) - P(n,0));
	// Calculate weights:
	P(nP,3) = 0;
	w = 0;
	for (int n = 0; n < nP; ++n)
		P(n, 3) = P(n,1) * (P(n+1,0) - P(n,0)) + 0.5 * (P(n+1,0) - P(n,0)) * (P(n+1,1) - P(n,1));
	// = (P(n+1,0) - P(n,0)) * (P(n,1) - P(n,2) * P(n,0))
	// + 0.5 * (P(n+1,0) * P(n+1,0) - P(n,0) * P(n,0)) * P(n,2);
	w = P.col(3).sum();
	if (w == 0) { //assume single dirac delta function for now
		dirac = true;
		w = 1;
		m = P(0,0);
	}
	else {
		// Normalize:
		P.block(0, 1, nP+1, 2) /= w;
		w = 1;
		// Calculate m:
		dirac = false;
		m = 0;
		for (int n = 0; n < nP; ++n)
			m += 0.5 * (P(n+1,0) * P(n+1,0) - P(n,0) * P(n,0)) * (P(n,1) - P(n,2) * P(n,0))
				+ 1.0/3 * (P(n+1,0) * P(n+1,0) * P(n+1,0) - P(n,0) * P(n,0) * P(n,0)) * P(n,2);
		m /= w;
	}
}

inline double LinearMembershipFunction::value(double x) const {
	int n = 0;
	if ((x < P(0,0)) or (x > P(nP,0))) return 0;
	if (x == P(nP,0)) return P(nP,1) / w;
	while (x >= P(n+1,0)) n++;
	return (P(n,1) + (x - P(n,0)) * P(n,2)) / w;
}

inline double LinearMembershipFunction::mean() const {
	return m;
}

inline double LinearMembershipFunction::sample() const {
	if (dirac) return P(0,0);
	int n = RandTom::sample(P.col(3));
	if (n == nP) return P(nP,0);
	if (P(n,2) == 0) return P(n,0) + RandTom::r.Fixed() * (P(n+1,0) - P(n,0));
	if (P(n,2) > 0) return P(n,0) - P(n,1) / P(n,2) + std::sqrt(P(n,1) * P(n,1) / (P(n,2) * P(n,2)) + 2 * RandTom::r.Fixed() * P(n,3) / P(n,2));
	if (P(n,2) < 0) return P(n,0) - P(n,1) / P(n,2) - std::sqrt(P(n,1) * P(n,1) / (P(n,2) * P(n,2)) + 2 * RandTom::r.Fixed() * P(n,3) / P(n,2));
	return 0;
}

double LinearMembershipFunction::operator*(const LinearMembershipFunction& f) const {
	if (dirac) return P(0,1)/w * f.value(P(0,0));
	if (f.dirac) return f.P(0,1)/f.w * value(f.P(0,0));		
	int i=0, j=0;
	double x = std::max(P(i,0), f.P(j,0));
	double x_ = x;
	double x_max = std::min(P(nP,0), f.P(f.nP,0));
	double res = 0;
	while (x < x_max) {
		while ((P(i+1, 0) <= x)) i++;
		while ((f.P(j+1, 0) <= x)) j++;
		double x_ = std::min(P(i+1,0), f.P(j+1, 0));
		double y0 = P(i,1) - P(i,2) * P(i,0), m0 = P(i,2);
		double y1 = f.P(j,1) - f.P(j,2) * f.P(j,0), m1 = f.P(j,2);
		res += 1.0/3 * (std::pow(x_,3) - std::pow(x,3)) * m0 * m1
			+ 1.0/2 * (std::pow(x_,2) - std::pow(x,2)) * ( m1 * y0 + m0 * y1 )
			+ (x_ - x) * y0 * y1;
		x = x_;
	}
	return res;
}

std::ostream& LinearMembershipFunction::operator>>(std::ostream &os) const {
	os << "LINEAR_MEMBERSHIP_FUNCTION" << std::endl;
	os << "nP: " << nP << std::endl << P.block(0,0,nP+1,2) << std::endl;
	return os;
}

std::istream& LinearMembershipFunction::operator<<(std::istream& is) {
	std::string dummy;
	is >> dummy; // LINEAR_MEMBERSHIP_FUNCTION
	is >> dummy >> nP;
	P.resize(nP+1,4);
	for (int n = 0; n < nP+1; ++n)
		is >> P(n,0) >> P(n,1);
	init();
	return is;
}

void LinearMembershipFunction::from_string(const std::string& str) {
	std::stringstream iss(str);
	*this << iss;
}

std::string LinearMembershipFunction::to_string() const {
	std::stringstream oss;
	*this >> oss;
	return oss.str();
}

void LinearMembershipFunction::show() const {
	std::cout << *this;
}

std::ostream& operator<<(std::ostream& os, const LinearMembershipFunction& f) {
	return f >> os;
}

std::istream& operator>>(std::istream& is, LinearMembershipFunction& f) {
	return f << is;
}

#endif // SWIG

#endif // __MEMBERSHIP_FUNCTIONS_H__
