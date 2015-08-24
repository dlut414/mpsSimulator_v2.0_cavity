/*
*/
#pragma once
#include <list>
#include "header.h"
#include "Parameter.h"
#include "Particle.h"
//#include "BgPart.h"

namespace SIM {

	template <class real, enum Dim dim>
	class Shifter {
		typedef Vec3<real> vec;
		typedef Mat3<real> mat;
		typedef Eigen::Matrix<real, dim, dim> matEi;
	public:
		Shifter() {}
		~Shifter() {}

		template <class R, enum Dim D, class Der>
		void shiftNearest2d(Particle<R, D, Der>* const part, const Parameter<R>& para) const {
			std::vector<vec> dp(part->pos);
			std::vector<vec> tmp(part->np, vec(0.));
			for (int loop = 0; loop < 1; loop++) {
#if OMP
#pragma omp parallel for
#endif
				for (int p = 0; p < int(part->np); p++) {
					if (part->type[p] != FLUID) continue;
					const auto inf = std::numeric_limits<real>::max();
					std::list<int> knn(8, 0);
					std::list<real> knd(8, inf);
					//auto gc = vec(0., 0., 0.);
					const auto c = part->cell->iCoord(part->pos[p]);
					for (auto k = -1; k <= 1; k++) {
						for (auto j = -1; j <= 1; j++) {
							for (auto i = -1; i <= 1; i++) {
								const auto ne = c + iVec3(i, j, k);
								const auto key = part->cell->hash(ne);
								for (auto m = 0; m < part->cell->linkList[key].size(); m++) {
									const auto q = part->cell->linkList[key][m];
									if (q == p) continue;
									const auto dr = part->pos[q] - part->pos[p];
									const auto dr2 = dr.mag2();
									if (dr2 > part->r0*part->r0) continue;
									std::list<int>::iterator it1 = knn.begin();
									std::list<real>::iterator it2 = knd.begin();
									for (; it1 != knn.end(); ++it1, ++it2) {
										if (dr2 < *it2) {
											knn.insert(it1, q);
											knn.pop_back();
											knd.insert(it2, dr2);
											knd.pop_back();
											break;
										}
									}
								}
							}
						}
					}
					auto dpq = vec(0., 0., 0.);
					for (std::list<int>::iterator it = knn.begin(); it != knn.end(); ++it) {
						const auto dr = part->pos[*it] - part->pos[p];
						const auto dr1 = dr.mag();
						dpq -= w_spring(dr1, part->dp*1.5)* dr.norm();
					}
					//if (part->isFs(p)) {
					//	gc = gc.norm();
					//	dpq = dpq - 1.*(dpq*gc)*gc;
					//}
					part->pos[p] += para.umax*para.dt* dpq;
				}
			}
			/*explicit Euler*/
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				const auto tmp = part->pos[p];
				part->pos[p] = dp[p];
				dp[p] = tmp;
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				if (part->isFs(p)) continue;
				tmp[p] = part->derived().func_mls_a_upwind(part->vel2, p, dp[p]);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] = dp[p];
				part->vel2[p] = tmp[p];
			}
		}

		template <class R, enum Dim D, class Der>
		void shiftPnd(Particle<R, D, Der>* part, Parameter<R>& para) const {
			std::vector<vec> dn(part->np, vec(0.));
			for (unsigned p = 0; p < part->pnd.size(); p++) {
				dn[p] = part->grad_suzuki_pnd(part->pnd, p);
			}
			real nn = 0.;
			for (unsigned p = 0; p < dn.size(); p++) {
				nn += dn[p] * dn[p];
			}
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if (part->type[p] != FLUID) continue;
				vec gc = vec(0., 0., 0.);
				const iVec3 c = part->cell->iCoord(part->pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const unsigned key = part->cell->hash(ne);
							for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
								const unsigned q = part->cell->linkList[key][m];
								if (q == p) continue;
								const vec	dr = part->pos[p] - part->pos[q];
								const real  dr1 = dr.mag();
								if (dr1 > part->r0) continue;
								gc += ( /*part->w1(dr1)*/1. / dr1)* dr;
							}
						}
					}
				}
				const real corr = part->nbd[p] > 0 ? real(part->pn0) / real(part->pn[p]) : 1.;
				vec dp = (part->n0 - part->pnd[p]*corr) / nn * dn[p];
				if ((part->type[p] == FLUID) && part->isFs(p)) {
					gc = gc.norm();
					//dp = dp - (dp*gc)*gc;
				}
				dp = 20.*dp;
				part->vel1[p] += part->grad_suzuki(part->vel1, p) * dp;
				part->pos[p] += dp;
			}
		}

		template <class R, enum Dim D, class Der>
		void shiftPbf(Particle<R, D, Der>* part, Parameter<R>& para) const {
			std::vector<vec> dp(part->pos.size());
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if (part->type[p] != FLUID) continue;
				vec dpq = vec(0., 0., 0.);
				vec gc = vec(0., 0., 0.);
				real lam = 0.;
				const iVec3 c = part->cell->iCoord(part->pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const unsigned key = part->cell->hash(ne);
							for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
								const unsigned q = part->cell->linkList[key][m];
								if (q == p) continue;
								const vec	dr = part->pos[q] - part->pos[p];
								const real	dr2 = dr.mag2();
								const real  dr1 = sqrt(dr2);
								if (dr1 > part->r0) continue;
								const vec cp = (-part->r0 / (dr2 * dr1)) * dr;
								lam += cp*cp;
								dpq -= cp;

								gc -= part->w1(dr1) * dr;
							}
						}
					}
				}
				lam += dpq* dpq;
				lam = (part->n0 - part->pnd[p]) / (lam + 1.e-3);
				if ((part->type[p] == FLUID) && (part->pn[p] < part->pn0 * para.beta)) {
					gc = gc.norm();
					dpq = dpq - (dpq*gc)*gc;
				}
				dp[p] = 0.001*lam * dpq;
			}
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] += part->grad(part->vel1, p) * dp[p];
			}
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] += dp[p];
				part->vel1[p] = part->vel2[p];
			}
		}

		template <class R, enum Dim D, class Der>
		void shiftXu(Particle<R, D, Der>* const part, const Parameter<R>& para) const {
			std::vector<vec> dp(part->pos.size());
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if ((part->type[p] != FLUID)) continue;
				real pa = 0.;
				unsigned count = 0;
				vec dpq = vec(0., 0., 0.);
				vec gc = vec(0., 0., 0.);
				const iVec3 c = part->cell->iCoord(part->pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const unsigned key = part->cell->hash(ne);
							for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
								const unsigned q = part->cell->linkList[key][m];
								if (q == p) continue;
								const vec	dr = part->pos[q] - part->pos[p];
								const real	dr2 = dr.mag2();
								const real  dr1 = sqrt(dr2);
								if (dr1 > part->r0) continue;
								pa += dr1;
								count++;
								dpq -= (1. / dr2) * (dr / dr1);
								gc -= /*w2(dr1, part->r0)**/ (dr / dr1);
							}
						}
					}
				}
				count = (count == 0) ? 1 : count;
				pa = pa / count;

				if (part->isFs(p)) {
					gc = gc.norm();
					dpq = dpq - 1.*(dpq*gc)*gc;
				}
				dp[p] = 0.05* para.umax * para.dt *pa*pa * dpq;
				//dp[p] = 0.05* para.cfl*para.dp *pa*pa * dpq;
			}
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if ((part->type[p] != FLUID)) continue;
				//part->vel2[p] += dp[p] * (part->derived().grad(part->vel1, p));
				//if ((part->vel2[p]).mag() > 5.) part->vel2[p] = 5. * (part->vel2[p]).norm();
			}
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if ((part->type[p] != FLUID)) continue;
				//part->vel1[p] = part->vel2[p];
				part->pos[p] += dp[p];
			}
		}

		template <class R, enum Dim D, class Der>
		void shiftSpring(Particle<R, D, Der>* const part, const Parameter<R>& para) const {
			std::vector<vec> dp(part->np, vec(0.));
			std::vector<vec> tmp1(part->np, vec(0.));
			std::vector<vec> tmp2(part->np, vec(0.));
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				auto dpq = vec(0., 0., 0.);
				auto gc = vec(0., 0., 0.);
				const auto c = part->cell->iCoord(part->pos[p]);
				for (auto k = -1; k <= 1; k++) {
					for (auto j = -1; j <= 1; j++) {
						for (auto i = -1; i <= 1; i++) {
							const auto ne = c + iVec3(i, j, k);
							const auto key = part->cell->hash(ne);
							for (auto m = 0; m < part->cell->linkList[key].size(); m++) {
								const auto q = part->cell->linkList[key][m];
								if (q == p) continue;
								const auto dr = part->pos[q] - part->pos[p];
								const auto dr1 = dr.mag();
								const auto re = part->dp*para.alpha;
								//gc -= w2(dr1, part->r0)* (dr / dr1);
								if (dr1 > re) continue;
								const auto ww = w2(dr1, re);
								dpq -= ww * (dr / dr1);
							}
						}
					}
				}
				//if (part->isFs(p)) {
				//	gc = gc.norm();
				//	dpq = dpq - 1.*(dpq*gc)*gc;
				//}
				dp[p] = para.c* para.umax* para.dt* dpq;
			}
			/*explicit Euler*/
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				if (part->isFs(p)) continue;
				dp[p] = part->pos[p] + dp[p];
				tmp1[p] = part->derived().func_mls_a_upwind(part->vel1, p, dp[p]);
				tmp2[p] = part->derived().func_mls_a_upwind(part->vel2, p, dp[p]);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] = dp[p];
				part->vel1[p] = tmp1[p];
				part->vel2[p] = tmp2[p];
			}
			/*corrected Euler*/
			/*second order Ronge-Kuta*/
			/*fourth order Ronge-Kuta*/
		}

		template <class R, enum Dim D, class Der>
		void shiftLs(Particle<R, D, Der>* const part, const Parameter<R>& para) const {
			std::vector<int> nearb(part->pos.size(), 0);
			std::vector<vec> dp(part->pos.size());
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if ((part->type[p] != FLUID)) continue;
				vec dpq = vec(0., 0., 0.);
				vec gc = vec(0., 0., 0.);
				mat m3 = mat(vec(0., 0., 0.), vec(0., 0., 0.), vec(0., 0., 0.));
				const iVec3 c = part->cell->iCoord(part->pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const unsigned key = part->cell->hash(ne);
							for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
								const unsigned q = part->cell->linkList[key][m];
								if (q == p) continue;
								const vec	dr = part->pos[q] - part->pos[p];
								const real	dr2 = dr.mag2();
								const real  dr1 = sqrt(dr2);
								if (dr1 > part->r0) continue;
								const real  w = part->w1(dr1);
								const vec	npq = dr / dr1;
								m3.x += w * npq.x * npq;
								m3.y += w * npq.y * npq;
								m3.z += w * npq.z * npq;
								gc += (/*part->w1(dr1)*/1. / dr1)* dr;
								if (part->type[q] != FLUID) nearb[p]++;
							}
						}
					}
				}
				m3 = (dim / part->n0) * m3;
				gc = gc.norm();

				matEi mi;
				switch (dim) {
				case 2:
					mi << m3.x.x, m3.x.z, m3.z.x, m3.z.z; break;
				case 3:
					mi << m3.x.x, m3.x.y, m3.x.z, m3.y.x, m3.y.y, m3.y.z, m3.z.x, m3.z.y, m3.z.z; break;
				}
				Eigen::SelfAdjointEigenSolver<matEi> sol(mi);
				Eigen::Matrix<real, dim, 1> eig = sol.eigenvalues();
				Eigen::Matrix<real, dim, dim> eigvs = sol.eigenvectors();
				switch (dim) {
				case 2:
					for (int d = 0; d < dim; d++) {
						dpq.x += eig[d] * eigvs.col(d)[0];
						dpq.z += eig[d] * eigvs.col(d)[1];
					}
					break;
				case 3:
					for (int d = 0; d < dim; d++) {
						dpq.x += eig[d] * eigvs.col(d)[0];
						dpq.y += eig[d] * eigvs.col(d)[1];
						dpq.z += eig[d] * eigvs.col(d)[2];
					}
					break;
				}
				if (dpq * gc > 0) dpq = -1.*dpq;

				if (part->isFs(p)) {
					dpq = dpq - 1.*(dpq*gc)*gc;
				}
				dp[p] = para.dp* /*para.umax * para.dt*/para.dp*para.cfl * dpq;
			}
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if ((part->type[p] != FLUID)) continue;
				//part->vel2[p] += dp[p] * ( part->grad_suzuki(part->vel1, p) );
			}
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if ((part->type[p] != FLUID)) continue;
				//part->vel1[p] = part->vel2[p];
				part->pos[p] += dp[p];
			}
		}

		template <class R, enum Dim D, class Der>
		void shiftPndLs(Particle<R, D, Der>* const part, const Parameter<R>& para) const {
			std::vector<int> nearb(part->pos.size(), 0);
			std::vector<vec> dp(part->pos.size());
			for (unsigned p = 0; p < part->pos.size(); p++) {
				vec dpq = -1.*part->grad_suzuki_dm(part->pnd, p);
				dp[p] = 0.05*para.dp* /*para.umax * para.dt*/para.dp*para.cfl * dpq;
			}
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if ((part->type[p] != FLUID)) continue;
				//part->vel2[p] += dp[p] * ( part->grad_suzuki(part->vel1, p) );
			}
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if ((part->type[p] != FLUID) || part->isFs(p)) continue;
				//part->vel1[p] = part->vel2[p];
				part->pos[p] += dp[p];
			}
		}

		template <class R, enum Dim D, class Der>
		void shiftCo(Particle<R, D, Der>* part, Parameter<R>& para) const {
			std::vector<vec> dp(part->pos.size());
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if (part->type[p] != FLUID) continue;
				vec dpq = vec(0., 0., 0.);
				const iVec3 c = part->cell->iCoord(part->pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) { 
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const unsigned key = part->cell->hash(ne);
							for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
								const unsigned q = part->cell->linkList[key][m];
								if (q == p) continue;
								const vec	dr = part->pos[q] - part->pos[p];
								const vec	dv = part->vel1[q] - part->vel1[p];
								const real	dr2 = dr.mag2();
								const real  dr1 = sqrt(dr2);
								if ((dr1 < 0.5* part->dp) && (dr*dv < 0.)) {
									dpq += (1. / dr2) * (dv* dr) * dr;
								}
							}
						}
					}
				}
				dp[p] = 0.1 * para.dt * dpq;
			}

			for (unsigned p = 0; p < part->pos.size(); p++) {
				if ( (part->type[p] != FLUID) ) continue;
				part->vel1[p] += part->grad(part->vel1, p) * dp[p];
				part->pos[p] += dp[p];
			}
		}

		void shiftRep() const {}

	private:
		inline const real w_spline(const real& r, const real& r0) const {
			/*cubic spline*/
			const real q = r / r0;
			if (q <= 0.5) return (1. - 6 * q*q + 6 * q*q*q);
			else if (q <= 1.) return (2.*pow((1. - q), 3));
			else return 0.;
			/*inverse r2*/
			//return pow(r0 / r, 3);
		}
		inline const real w2(const real& r, const real& r0) const {
			if (r >= r0) {
				return 0.;
			}
			else {
				//return r0 / r - 1.;
				return pow((1 - r / r0), 2);
			}
		}
		inline const real w3(const real& r, const real& r0) const {
			/*parabolic*/
			const real q = r / r0;
			if (q < 1.) return 2.5*(1. - q*q);
			else return 0.;
		}
		inline const real w_spring(const real&r, const real& r0) const {
			if (r >= r0) return 0.;
			return 1. - r / r0;
		}

	};

}