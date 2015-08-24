/*
*/
#pragma once
#include "header.h"
#include "Particle.h"

namespace SIM {

	template <typename real, enum Dim dim>
	class Particle_2d_e : public Particle< real, dim, Particle_2d_e<real, dim> > {
		typedef Vec3<real> vec;
		typedef Mat3<real> mat;
		typedef Eigen::Matrix<real, 1, 3> vec13;
		typedef Eigen::Matrix<real, 5, 1> vec5;
		typedef Eigen::Matrix<real, 5, 5> mat55;
		typedef Eigen::Matrix<real, 5, 3> mat53;
		typedef Eigen::Matrix<real, dim, dim> matEi;
	public:
		Particle_2d_e() : Particle() {}
		~Particle_2d_e() {}

		/*poly2d_e*/
		const real func(const std::vector<real>& phi, const unsigned& p) const {
			vec5  vv = vec5::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec5	npq = poly2d_e(dr);
							vv += w * phi[q] * npq;
						}
					}
				}
			}
			vec5 a = invMat[p] * vv;
			vec5 func = poly2d_e_0;
			return func.dot(a);
		}

		const vec func(const std::vector<vec>& u, const unsigned& p) const {
			mat53 vv = mat53::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec5	npq = poly2d_e(dr);
							vv.block<5, 1>(0, 0) += w * u[q].x * npq;
							vv.block<5, 1>(0, 1) += w * u[q].y * npq;
							vv.block<5, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}

			mat53 a = invMat[p] * vv;
			vec5 func = poly2d_e_0;
			vec13 ret = func.transpose() * a;
			return vec(ret[0], ret[1], ret[2]);
		}

		const vec grad(const std::vector<real>& phi, const unsigned& p) const {
			vec5  vv = vec5::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec5	npq = poly2d_e(dr);
							vv += w * phi[q] * npq;
						}
					}
				}
			}
			vec5  a = invMat[p] * vv;
			vec5 px = poly2d_e_px_0;
			vec5 pz = poly2d_e_pz_0;
			return vec(a.dot(px), 0., a.dot(pz));
		}

		const mat grad(const std::vector<vec>& u, const unsigned& p) const {
			mat53 vv = mat53::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec5	npq = poly2d_e(dr);
							vv.block<5, 1>(0, 0) += w * u[q].x * npq;
							vv.block<5, 1>(0, 1) += w * u[q].y * npq;
							vv.block<5, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			mat53 a = invMat[p] * vv;
			vec5 px = poly2d_e_px_0;
			vec5 pz = poly2d_e_pz_0;
			return mat(vec(px.dot(a.block<5, 1>(0, 0)), 0., px.dot(a.block<5, 1>(0, 2))),
				vec(0.),
				vec(pz.dot(a.block<5, 1>(0, 0)), 0., pz.dot(a.block<5, 1>(0, 2))));
		}

		const real div(const std::vector<vec>& u, const unsigned& p) const {
			mat53 vv = mat53::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec5	npq = poly2d_e(dr);
							vv.block<5, 1>(0, 0) += w * u[q].x * npq;
							vv.block<5, 1>(0, 1) += w * u[q].y * npq;
							vv.block<5, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			mat53 a = invMat[p] * vv;
			vec5 px = poly2d_e_px_0;
			vec5 pz = poly2d_e_pz_0;
			return a.block<5, 1>(0, 0).dot(px) + 0. + a.block<5, 1>(0, 2).dot(pz);
		}

		const real lap(const std::vector<real>& phi, const unsigned& p) const {
			vec5 vv = vec5::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif;
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec5	npq = poly2d_e(dr);
							vv += w * phi[q] * npq;
						}
					}
				}
			}
			vec5 a = invMat[p] * vv;
			vec5 lap = poly2d_e_lap_0;
			return lap.dot(a);
		}

		const vec lap(const std::vector<vec>& u, const unsigned& p) const {
			mat53 vv = mat53::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec5	npq = poly2d_e(dr);
							vv.block<5, 1>(0, 0) += w * u[q].x * npq;
							vv.block<5, 1>(0, 1) += w * u[q].y * npq;
							vv.block<5, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			mat53 a = invMat[p] * vv;
			vec5 lapxz = poly2d_e_lap_0;
			return vec(a.block<5, 1>(0, 0).dot(lapxz), 0., a.block<5, 1>(0, 2).dot(lapxz));
		}

		const real func(const std::vector<real>& phi, const vec& p) const {
			mat55 mm = mat55::Zero();
			vec5  vv = vec5::Zero();
			const iVec3 c = cell->iCoord(p);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(q)) continue;
#endif
							const vec	dr = pos[q] - p;
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec5	npq = poly2d_e(dr);
							vv += w * phi[q] * npq;
							mm.block<1, 5>(0, 0) += w * npq[0] * npq;
							mm.block<1, 5>(1, 0) += w * npq[1] * npq;
							mm.block<1, 5>(2, 0) += w * npq[2] * npq;
							mm.block<1, 5>(3, 0) += w * npq[3] * npq;
							mm.block<1, 5>(4, 0) += w * npq[4] * npq;
						}
					}
				}
			}
			mat55 inv = mat55::Zero();
			if (abs(mm.determinant()) < eps_mat) {
				real mm_ = mm(0, 0);
				if (abs(mm_) < eps_mat) {
					inv = mat55::Zero();
				}
				else inv(0, 0) = 1. / mm_;
			}
			else {
				inv = mm.inverse();
			}
			vec5 a = inv* vv;
			vec5 func = poly2d_e_0;
			return func.dot(a);
		}

		const vec func(const std::vector<vec>& u, const vec& p) const {
			mat55 mm = mat55::Zero();
			mat53 vv = mat53::Zero();
			const iVec3 c = cell->iCoord(p);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(q)) continue;
#endif
							const vec	dr = pos[q] - p;
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec5	npq = poly2d_e(dr);
							vv.block<5, 1>(0, 0) += w * u[q].x * npq;
							vv.block<5, 1>(0, 1) += w * u[q].y * npq;
							vv.block<5, 1>(0, 2) += w * u[q].z * npq;
							mm.block<1, 5>(0, 0) += w * npq[0] * npq;
							mm.block<1, 5>(1, 0) += w * npq[1] * npq;
							mm.block<1, 5>(2, 0) += w * npq[2] * npq;
							mm.block<1, 5>(3, 0) += w * npq[3] * npq;
							mm.block<1, 5>(4, 0) += w * npq[4] * npq;
						}
					}
				}
			}
			mat55 inv = mat55::Zero();
			if (abs(mm.determinant()) < eps_mat) {
				real mm_ = mm(0, 0);
				if (abs(mm_) < eps_mat) {
					inv = mat55::Zero();
				}
				else inv(0, 0) = 1. / mm_;
			}
			else {
				inv = mm.inverse();
			}
			mat53 a = inv* vv;
			vec5 func = poly2d_e_0;
			vec13 ret = func.transpose() * a;
			return vec(ret[0], ret[1], ret[2]);
		}

		void initInvMat() {
			invMat.clear();
			for (int p = 0; p < int(np); p++) {
				invMat.push_back(mat55());
			}
		}
		void updateInvMat() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(np); p++) {
				mat55 mm = mat55::Zero();
				const iVec3 c = cell->iCoord(pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const unsigned key = cell->hash(ne);
							for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
								const unsigned q = cell->linkList[key][m];
#if BD_OPT
								if (bdOpt(p, q)) continue;
#endif
								const vec	dr = pos[q] - pos[p];
								const real	dr1 = dr.mag();
								if (dr1 > r0) continue;
								const real  w = w3(dr1);
								const vec5	npq = poly2d_e(dr);
								mm.block<1, 5>(0, 0) += w * npq[0] * npq;
								mm.block<1, 5>(1, 0) += w * npq[1] * npq;
								mm.block<1, 5>(2, 0) += w * npq[2] * npq;
								mm.block<1, 5>(3, 0) += w * npq[3] * npq;
								mm.block<1, 5>(4, 0) += w * npq[4] * npq;
							}
						}
					}
				}
				if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
					std::cout << mm.determinant() << std::endl;
#endif
					real mm_ = mm(0, 0);
					if (abs(mm_) < eps_mat) {
						invMat[p] = mat55::Zero();
						continue;
					}
					invMat[p] = mat55::Zero();
					invMat[p](0, 0) = 1. / mm_;
					continue;
				}
				invMat[p] = mm.inverse();
			}
		}

		void updateDiver() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(np); p++) {
				diver[p] = div_mls_poly2d_e(part->vel1, p);
			}
		}

	public:
		std::vector<mat55> invMat;

	};

}