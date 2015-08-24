/*
*/
#pragma once
#include "header.h"
#include "Particle.h"

namespace SIM {

	template <typename real, enum Dim dim>
	class Particle_2d_c : public Particle< real, dim, Particle_2d_c<real, dim> > {
		typedef Vec3<real> vec;
		typedef Mat3<real> mat;
		typedef Eigen::Matrix<real, 1, 3> vec13;
		typedef Eigen::Matrix<real, 5, 1> vec5;
		typedef Eigen::Matrix<real, 6, 1> vec6;
		typedef Eigen::Matrix<real, 7, 1> vec7;
		typedef Eigen::Matrix<real, 8, 1> vec8;
		typedef Eigen::Matrix<real, 9, 1> vec9;
		typedef Eigen::Matrix<real, 5, 5> mat55;
		typedef Eigen::Matrix<real, 6, 6> mat66;
		typedef Eigen::Matrix<real, 7, 7> mat77;
		typedef Eigen::Matrix<real, 8, 8> mat88;
		typedef Eigen::Matrix<real, 9, 9> mat99;
		typedef Eigen::Matrix<real, 5, 3> mat53;
		typedef Eigen::Matrix<real, 6, 3> mat63;
		typedef Eigen::Matrix<real, 7, 3> mat73;
		typedef Eigen::Matrix<real, 8, 3> mat83;
		typedef Eigen::Matrix<real, 9, 3> mat93;
		typedef Eigen::Matrix<real, dim, dim> matEi;
	public:
		Particle_2d_c() : Particle() {}
		~Particle_2d_c() {}

		/*poly2d_d*/
		const real func(const std::vector<real>& phi, const unsigned& p) const {
			vec7  vv = vec7::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							//if (type[q] == BD2) continue;
							//if (team[p] != team[q]) continue;
							//if (type[q] == BD1 && isFs(q)) continue;
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec7	npq = poly2d_c(dr);
							vv += w * phi[q] * npq;
						}
					}
				}
			}
			vec7 a = invMat[p] * vv;
			vec7 func = poly2d_c_0;
			return func.dot(a);
		}

		const vec func(const std::vector<vec>& u, const unsigned& p) const {
			mat73 vv = mat73::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							//if (type[q] == BD2) continue;
							//if (team[p] != team[q]) continue;
							//if (type[q] == BD1 && isFs(q)) continue;
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec7	npq = poly2d_c(dr);
							vv.block<7, 1>(0, 0) += w * u[q].x * npq;
							vv.block<7, 1>(0, 1) += w * u[q].y * npq;
							vv.block<7, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}

			mat73 a = invMat[p] * vv;
			vec7 func = poly2d_c_0;
			vec13 ret = func.transpose() * a;
			return vec(ret[0], ret[1], ret[2]);
		}

		const vec grad(const std::vector<real>& phi, const unsigned& p) const {
			vec7  vv = vec7::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							//if (type[q] == BD2) continue;
							//if (team[p] != team[q]) continue;
							//if (type[q] == BD1 && isFs(q)) continue;
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec7	npq = poly2d_c(dr);
							vv += w * phi[q] * npq;
						}
					}
				}
			}
			vec7  a = invMat[p] * vv;
			vec7 px = poly2d_c_px_0;
			vec7 pz = poly2d_c_pz_0;
			return vec(a.dot(px), 0., a.dot(pz));
		}

		const mat grad(const std::vector<vec>& u, const unsigned& p) const {
			mat73 vv = mat73::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							//if (type[q] == BD2) continue;
							//if (team[p] != team[q]) continue;
							//if (type[q] == BD1 && isFs(q)) continue;
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec7	npq = poly2d_c(dr);
							vv.block<7, 1>(0, 0) += w * u[q].x * npq;
							vv.block<7, 1>(0, 1) += w * u[q].y * npq;
							vv.block<7, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			mat73 a = invMat[p] * vv;
			vec7 px = poly2d_c_px_0;
			vec7 pz = poly2d_c_pz_0;
			return mat(vec(px.dot(a.block<7, 1>(0, 0)), 0., px.dot(a.block<7, 1>(0, 2))),
				vec(0.),
				vec(pz.dot(a.block<7, 1>(0, 0)), 0., pz.dot(a.block<7, 1>(0, 2))));
		}

		const real div(const std::vector<vec>& u, const unsigned& p) const {
			mat73 vv = mat73::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							//if (type[q] == BD2) continue;
							//if (team[p] != team[q]) continue;
							//if (type[q] == BD1 && isFs(q)) continue;
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec7	npq = poly2d_c(dr);
							vv.block<7, 1>(0, 0) += w * u[q].x * npq;
							vv.block<7, 1>(0, 1) += w * u[q].y * npq;
							vv.block<7, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			mat73 a = invMat[p] * vv;
			vec7 px = poly2d_c_px_0;
			vec7 pz = poly2d_c_pz_0;
			return a.block<7, 1>(0, 0).dot(px) + 0. + a.block<7, 1>(0, 2).dot(pz);
		}

		const vec lap(const std::vector<vec>& u, const unsigned& p) const {
			mat73 vv = mat73::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							//if (type[q] == BD2) continue;
							//if (team[p] != team[q]) continue;
							//if (type[q] == BD1 && isFs(q)) continue;
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec7	npq = poly2d_c(dr);
							vv.block<7, 1>(0, 0) += w * u[q].x * npq;
							vv.block<7, 1>(0, 1) += w * u[q].y * npq;
							vv.block<7, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			mat73 a = invMat[p] * vv;
			vec7 lapxz = poly2d_c_lap_0;
			return vec(a.block<7, 1>(0, 0).dot(lapxz), 0., a.block<7, 1>(0, 2).dot(lapxz));
		}

		const real func(const std::vector<real>& phi, const vec& p) const {
			mat77 mm = mat77::Zero();
			vec9  vv = vec9::Zero();
			const iVec3 c = cell->iCoord(p);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							//if (type[q] == BD2) continue;
							//if (team[p] != team[q]) continue;
							//if (type[q] == BD1 && isFs(q)) continue;
							const vec	dr = pos[q] - p;
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec7	npq = poly2d_c(dr);
							vv += w * phi[q] * npq;
							mm.block<1, 7>(0, 0) += w * npq[0] * npq;
							mm.block<1, 7>(1, 0) += w * npq[1] * npq;
							mm.block<1, 7>(2, 0) += w * npq[2] * npq;
							mm.block<1, 7>(3, 0) += w * npq[3] * npq;
							mm.block<1, 7>(4, 0) += w * npq[4] * npq;
							mm.block<1, 7>(5, 0) += w * npq[5] * npq;
							mm.block<1, 7>(6, 0) += w * npq[6] * npq;
						}
					}
				}
			}
			mat77 inv = mat77::Zero();
			if (abs(mm.determinant()) < part->eps_mat) {
				mat55 mm_ = mm.block<5, 5>(0, 0);
				if (abs(mm_.determinant()) < part->eps_mat) {
					real mm__ = mm_(0, 0);
					if (abs(mm__) < part->eps_mat) {
						invMat[p] = mat77::Zero();
						continue;
					}
					else inv(0, 0) = 1. / mm__;
				}
				else {
					inv.block<5, 5>(0, 0) = mm_.inverse();
				}
			}
			else {
				inv = mm.inverse();
			}
			vec7 a = inv* vv;
			vec7 func = poly2d_c_0;
			return func.dot(a);
		}

		const vec func(const std::vector<vec>& u, const vec& p) const {
			mat77 mm = mat77::Zero();
			mat73 vv = mat73::Zero();
			const iVec3 c = cell->iCoord(p);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							//if (type[q] == BD2) continue;
							//if (team[p] != team[q]) continue;
							//if (type[q] == BD1 && isFs(q)) continue;
							const vec	dr = pos[q] - p;
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vec7	npq = poly2d_c(dr);
							vv.block<7, 1>(0, 0) += w * u[q].x * npq;
							vv.block<7, 1>(0, 1) += w * u[q].y * npq;
							vv.block<7, 1>(0, 2) += w * u[q].z * npq;
							mm.block<1, 7>(0, 0) += w * npq[0] * npq;
							mm.block<1, 7>(1, 0) += w * npq[1] * npq;
							mm.block<1, 7>(2, 0) += w * npq[2] * npq;
							mm.block<1, 7>(3, 0) += w * npq[3] * npq;
							mm.block<1, 7>(4, 0) += w * npq[4] * npq;
							mm.block<1, 7>(5, 0) += w * npq[5] * npq;
							mm.block<1, 7>(6, 0) += w * npq[6] * npq;
						}
					}
				}
			}
			mat77 inv = mat77::Zero();
			if (abs(mm.determinant()) < eps_mat) {
				mat55 mm_ = mm.block<5, 5>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					real mm__ = mm_(0, 0);
					if (abs(mm__) < eps_mat) {
						inv = mat77::Zero();
					}
					else inv(0, 0) = 1. / mm__;
				}
				else {
					inv.block<5, 5>(0, 0) = mm_.inverse();
				}
			}
			else {
				inv = mm.inverse();
			}
			mat73 a = inv* vv;
			vec7 func = poly2d_c_0;
			vec13 ret = func.transpose() * a;
			return vec(ret[0], ret[1], ret[2]);
		}

		void initInvMat() {
			invMat.clear();
			for (int p = 0; p < int(np); p++) {
				invMat.push_back(mat77());
			}
		}
		void updateInvMat() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(np); p++) {
				mat77 mm = mat77::Zero();
				const iVec3 c = cell->iCoord(pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const unsigned key = cell->hash(ne);
							for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
								const unsigned q = cell->linkList[key][m];
								//if (type[q] == BD2) continue;
								//if (team[p] != team[q]) continue;
								//if (type[q] == BD1 && isFs(q)) continue;
								const vec	dr = pos[q] - pos[p];
								const real	dr1 = dr.mag();
								if (dr1 > r0) continue;
								const real  w = w3(dr1);
								const vec7	npq = poly2d_c(dr);
								mm.block<1, 7>(0, 0) += w * npq[0] * npq;
								mm.block<1, 7>(1, 0) += w * npq[1] * npq;
								mm.block<1, 7>(2, 0) += w * npq[2] * npq;
								mm.block<1, 7>(3, 0) += w * npq[3] * npq;
								mm.block<1, 7>(4, 0) += w * npq[4] * npq;
								mm.block<1, 7>(5, 0) += w * npq[5] * npq;
								mm.block<1, 7>(6, 0) += w * npq[6] * npq;
							}
						}
					}
				}
				if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
					std::cout << mm.determinant() << std::endl;
#endif
					mat55 mm_ = mm.block<5, 5>(0, 0);
					if (abs(mm_.determinant()) < eps_mat) {
						real mm__ = mm_(0, 0);
						if (abs(mm__) < eps_mat) {
							invMat[p] = mat77::Zero();
							continue;
						}
						invMat[p] = mat77::Zero();
						invMat[p](0, 0) = 1. / mm__;
						continue;
					}
					invMat[p] = mat77::Zero();
					invMat[p].block<5, 5>(0, 0) = mm_.inverse();
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
				diver[p] = div_mls_poly2d_c(part->vel1, p);
			}
		}

	public:
		std::vector<mat77> invMat;

	};

}