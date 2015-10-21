#ifndef POMDP_TOOLS_H
#define POMDP_TOOLS_H

#include "tom.h"

namespace tom {

class Policy {
	friend class cereal::access;
public:
	unsigned int nU_;
	double exploration_;
	Policy(unsigned int nU = 0, double exploration = 1) : nU_(nU), exploration_(exploration) {};
	int u(const Eigen::VectorXd& w, Random& r) const;
	Eigen::VectorXd p(const Eigen::VectorXd& w) const;
	void addPlane(int u, const std::vector<int>& indices, const std::vector<double> vals) {
		planes_.push_back(Plane(u, indices, vals));
		if (u >= nU_) nU_ = u + 1;
	}
	INSERT_JSON_IO_FUNCTIONS()
private:
	class Plane {
		friend class cereal::access;
	public:
		Plane() {}
		Plane(int u, const std::vector<int>& indices, const std::vector<double> vals) :
			u_(u), indices_(indices), vals_(vals) {}
		int u() const { return u_; }
		bool applicable(const Eigen::VectorXd& w, const double epsilon = 1e-7) const {
			int currentIndex = 0;
			for (int i = 0; i < w.size(); ++i) {
				if (w(i) > epsilon) {
					// we need to check if the index i is contained in indices_[currentIndex:]
					if (currentIndex == indices_.size()) return false;
					while (i > indices_[currentIndex]) {
						currentIndex++;
						if (currentIndex == indices_.size()) return false;
					}
					if (i != indices_[currentIndex]) return false;
					currentIndex++;
				}
			}
			return true;
		}
		double value(const Eigen::VectorXd& w) const {
			double ret = 0;
			for (int i = 0; i < indices_.size(); ++i) {
				ret += vals_[i] * w(indices_[i]);
			}
			return ret;
		}
	private:
		template<class Archive>
		void serialize(Archive& ar) { CEREALIZE(ar,  u_,  u); CEREALIZE(ar,  indices_,  indices); CEREALIZE(ar,  vals_,  vals); }

		int u_;
		std::vector<int> indices_;
		std::vector<double> vals_;
	};

	template<class Archive>
	void save(Archive& ar) const {
		const std::string type = "POLICY";
		CEREALIZE(ar, type, Type);
		CEREALIZE(ar,  nU_,  nU); CEREALIZE(ar,  exploration_,  exploration);
		CEREALIZE(ar,  planes_,  planes);
	}
	template<class Archive>
	void load(Archive& ar) {
		std::string type = "POLICY";
        CEREALIZE(ar, type, Type);
		CEREALIZE(ar,  nU_,  nU); CEREALIZE(ar,  exploration_,  exploration);
		CEREALIZE(ar,  planes_,  planes);
	}

	std::vector<Plane> planes_;
};

inline int Policy::u(const Eigen::VectorXd& w, Random& r) const {
	if (nU_ == 0) return 0;
	if (r.random() < exploration_) return r.integer(nU_);
	int act = r.integer(nU_);
	double reward = -1e99;
	for (auto p = planes_.begin(); p != planes_.end(); ++p) {
		if (p->applicable(w)) {
			double p_reward = p->value(w);
			if (p_reward > reward) { reward = p_reward; act = p->u(); }
		}
	}
	return act;
}

inline Eigen::VectorXd Policy::p(const Eigen::VectorXd& w) const {
	if (nU_ == 0) return Eigen::VectorXd::Zero(0);
	int act = -1;
	double reward = -1e99;
	for (auto p = planes_.begin(); p != planes_.end(); ++p) {
		if (p->applicable(w)) {
			double p_reward = p->value(w);
			if (p_reward > reward) { reward = p_reward; act = p->u(); }
		}
	}
	if (act == -1) return Eigen::VectorXd::Constant(nU_, double(1) / double(nU_));
	Eigen::VectorXd p = Eigen::VectorXd::Constant(nU_, exploration_ * double(1) / double(nU_));
	p(act) += (1 - exploration_);
	return p;
}


} //namespace tom

#endif // POMDP_TOOLS_H
