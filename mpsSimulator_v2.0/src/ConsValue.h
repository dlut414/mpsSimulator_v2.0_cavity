/*
*/
#pragma once
#include <vector>
#include <Eigen/Dense>
#include "header.h"
#define _USE_MATH_DEFINES
#include <cmath>

namespace SIM {

	template <typename real, enum Dim dim>
	class ConsValue {
		typedef Vec3<real> vec;
		typedef Eigen::Matrix<real, 3, 1> vec3;
		typedef Eigen::Matrix<real, 5, 1> vec5;
		typedef Eigen::Matrix<real, 6, 1> vec6;
		typedef Eigen::Matrix<real, 7, 1> vec7;
		typedef Eigen::Matrix<real, 8, 1> vec8;
		typedef Eigen::Matrix<real, 9, 1> vec9;
	public:
		ConsValue() {}
		~ConsValue() {}

		inline const real ww(const real& r, const real& re) const {
			if (r >= re) {
				return 0.;
			}
			else {
				//return re / r - 1.;
				return pow((1 - r / re), 2);
			}
		}
		//inline const real wPnd(const real& r) const {
		//	if (r >= r0) {
		//		return 0.;
		//	}
		//	else {
		//		return r0 / r - 1.;
		//	}
		//}
		//inline const real wPnd(const real& r, const real& re) const {
		//	if (r >= re) {
		//		return 0.;
		//	}
		//	else {
		//		return re / r - 1.;
		//	}
		//}
		inline const real w1(const real& r) const {
			if (r >= r0) {
				return 0.;
			}
			else {
				return pow((1 - r / r0), 2);
			}
		}
		inline const real w1(const real& r, const real& re) const {
			if (r >= re) {
				return 0.;
			}
			else {
				return pow((1 - r / re), 2);
			}
		}

		inline const real w2(const real& r) const {
			if (r >= r0) {
				return 0.;
			}
			else {
				return pow((1 - r / r0), 2);
			}
		}
		inline const real w3(const real& r) const {
			if (r >= r0) return 0.;
			else {
				return pow((1. - r / r0), 2);
			}
		}
		inline const real w4_bSpline(const real& r) const {
			if (r >= r0) return 0.;
			const real q = r / r0;
			if (q <= 0.5) return (1. - 6 * q*q + 6 * q*q*q);
			else return (2.*pow((1. - q), 3));
		}

	public:
		real dp;
		real k;
		real r0;
		real beta;	//beta: free surface
		real n0;
		real lambda;
		unsigned pn0;
		real eps;
		real eps_mat;
		
	public:
		void init(const real& _k, const real& _beta) {
			k = _k; beta = _beta; r0 = k* dp;
			eps = std::numeric_limits<real>::epsilon();
			eps_mat = std::numeric_limits<real>::epsilon();

			initConst();
		}

	private:
		void initConst() {
			std::vector<vec> v;
			if (dim == TWOD) {
				v.clear();
				for (int i = 0; i < k * 2 + 1; i++) {
					for (int j = 0; j < k * 2 + 1; j++) {
						vec p = vec(real(i), real(j), 0.);
						v.push_back(p);
					}
				}
				//cal n0
				//cal lambda
				n0 = 0.;
				lambda = 0.;
				pn0 = 0;
				for (unsigned i = 0; i < v.size(); i++) {
					real n = 0.;
					real lam = 0.;
					unsigned pn = 0;
					for (unsigned j = 0; j < v.size(); j++) {
						//if (j == i) continue;
						real w = w1((v[i] - v[j]).mag(), k);
						n += w;
						lam += (v[i] - v[j]).mag2() * w;
						if (w > 0.) pn++;
					}
					lam = lam / n;
					n0 = n>n0 ? n : n0;
					lambda = lam > lambda ? lam : lambda;
					pn0 = pn > pn0 ? pn : pn0;
				}
			}
			else if (dim == THREED) {
				v.clear();
				for (int i = 0; i < k * 2 + 1; i++) {
					for (int j = 0; j < k * 2 + 1; j++) {
						for (int k = 0; k < k * 2 + 1; k++) {
							vec p = vec(real(i), real(j), real(k));
							v.push_back(p);
						}
					}
				}
				//cal n0
				//cal lambda
				n0 = 0.;
				lambda = 0.;
				pn0 = 0;
				for (unsigned i = 0; i < v.size(); i++) {
					real n = 0.;
					real lam = 0.;
					unsigned pn = 0;
					for (unsigned j = 0; j < v.size(); j++) {
						//if (j == i) continue;
						real w = w1((v[i] - v[j]).mag(), k);
						n += w;
						lam += (v[i] - v[j]).mag2() * w;
						if (w > 0.) pn++;
					}
					lam = lam / n;
					n0 = n>n0 ? n : n0;
					lambda = lam>lambda ? lam : lambda;
					pn0 = pn > pn0 ? pn : pn0;
				}
			}
			std::cout << " n0: " << n0 << std::endl;
			std::cout << " lambda: " << lambda << std::endl;
			std::cout << " pn0: " << pn0 << std::endl;
			v.clear();
		}
	};

}