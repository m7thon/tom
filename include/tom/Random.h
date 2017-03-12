#ifndef RANDOM_H
#define RANDOM_H

#include "tom.h"

namespace tom {

/** This class provides basic random number generation functionality. */
class Random {
public:
	/** Create a \c Random object initialized with a random seed. */
	Random() : engine(std::random_device{}()) {}

	/** Create a \c Random object initialized with the given \c seed. */
	Random(unsigned int seed) : engine(seed) {}

	/** Randomly seed the random number generator and return the seed. */
	unsigned int seed() {
		std::random_device rd{};
		unsigned int seedVal = rd();
		engine.seed(seedVal);
		return seedVal;
	}

	/** Seed the random number generator with the given \c seedValue. */
	void seed(unsigned int seedValue) { engine.seed(seedValue); }

	/** Return a double value sampled uniformly from in the range [0, 1). */
	double random() const { return const_cast<Random&>(*this).real_dist(const_cast<Random&>(*this).engine); }

	/** Return a matrix of size \c m x \c n with uniformly random entries from [0, 1). */
	Eigen::MatrixXd random(int m, int n) const {
		Eigen::MatrixXd M(m,n);
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < m; ++i)
				M(i,j) = random();
		return M;
	}

	/** Return a non-negative integer sampled uniformly from the set {0, ..., `n`-1}. */
	unsigned int integer(unsigned int n) const {
		if (n <= 1) return 0;
		using dist_t = std::uniform_int_distribution<unsigned int>;
		using param_t = dist_t::param_type;
		return const_cast<Random&>(*this).int_dist(const_cast<Random&>(*this).engine, param_t{0,n-1});
	}

	/** Return a non-negative integer from the set {0, ..., n-1}, where n is the size of the given array \c probArray, distributed according to the discrete distribution described by the \c probArray\. Note that the probabilities in \c probArray are assumed to sum to 1. */
	template <typename Derived>
	unsigned int sample(const Eigen::DenseBase<Derived>& probArray) const {
		double dice = random() * probArray.sum();
		double cumul = 0;
		for (unsigned int i = 0; i < probArray.size(); i++) {
			cumul += probArray.coeff(i);
			if (dice < cumul)
				return i;
		}
		return integer((unsigned int) probArray.size());
	};

private:
	std::uniform_int_distribution<unsigned int> int_dist = std::uniform_int_distribution<unsigned int>{0,1};
	std::uniform_real_distribution<double> real_dist = std::uniform_real_distribution<double>{0,1};
	std::mt19937 engine; // Mersenne twister MT19937
};

SWIGCODE(%template(sample) Random::sample<ArrayMd >;)

} // namespace tom

#endif // RANDOM_H
