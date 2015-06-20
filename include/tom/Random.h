/**
 * @file   Random.h
 * @author Michael Thon
 *
 * @brief  This file provides the \a Random class for random number generation.
 */

#ifndef RANDOM_H
#define RANDOM_H

namespace tom {

/** This class provides basic random number generation functionality. */
class Random {
public:
	/** create a \a Random object initialized with a random seed */
	Random() : engine(std::random_device{}()) {}
	/** create a \a Random object initialized with the given \a seed */
	Random(unsigned int seed) : engine(seed) {}
  /** randomly seed the random number generator and return the seed */
	unsigned int seed() {
		std::random_device rd{};
		unsigned int seedVal = rd();
		engine.seed(seedVal);
		return seedVal;
	}
	/** seed the random number generator with the given \a seedVal */
	void seed(unsigned int seedVal) { engine.seed(seedVal); }
	/** return a random real value in the range [0, 1) */
	double random() { return real_dist(engine); }
	/** return a matrix of size \a m x \a n with uniformly random entries from [0,1) **/
	Eigen::MatrixXd random(int m, int n) {
		Eigen::MatrixXd M(m,n);
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < m; ++i)
				M(i,j) = random();
		return M;
	}
	/** return a random integer in the range [0, \a n) */
	unsigned int integer(unsigned int n) {
		using dist_t = std::uniform_int_distribution<unsigned int>;
		using param_t = dist_t::param_type;
		return int_dist(engine, param_t{0,n-1});
	}
	/** return an integer in the range [0, n), where n is the size of the array \a probs, distributed according to the discrete distribution described by \a probs\. Note that the probabilities in \a probs must sum to 1. */
	template <typename Derived>
	unsigned int sample(const Eigen::MatrixBase<Derived>& prob) {
		double dice = random() * prob.sum();
		double cumul = 0;
		for (unsigned int i = 0; i < prob.size(); i++) {
			cumul += prob(i);
			if (dice < cumul)
				return i;
		}
		return integer(prob.size());
	};
private:
	std::uniform_int_distribution<unsigned int> int_dist = std::uniform_int_distribution<unsigned int>{0,1};
	std::uniform_real_distribution<double> real_dist = std::uniform_real_distribution<double>{0,1};
	std::mt19937 engine; // Mersenne twister MT19937
};

SWIGCODE(%feature("autodoc", "sample(ArrayMd const & prob) -> unsigned int") Random::sample<ArrayMd >;)
SWIGCODE(%feature("docstring") Random::sample<ArrayMd > "This is some additional documentation";)
SWIGCODE(%template(sample) Random::sample<ArrayMd >;)

} // namespace tom

#endif // RANDOM_H
