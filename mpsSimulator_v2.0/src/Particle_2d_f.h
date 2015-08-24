/*
*/
#pragma once
#include "header.h"
#include "Particle.h"

namespace SIM {

	template <typename real, enum Dim dim>
	class Particle_2d_f : public Particle< real, dim, Particle_2d_f<real, dim> > {
		typedef Vec3<real> vec;
		typedef Mat3<real> mat;
		typedef Eigen::Matrix<real, 1, 3> vec13;
		typedef Eigen::Matrix<real, 8, 1> vec8;
		typedef Eigen::Matrix<real, 4, 4> mat44;
		typedef Eigen::Matrix<real, 8, 8> mat88;
		typedef Eigen::Matrix<real, 8, 3> mat83;
		typedef Eigen::Matrix<real, dim, dim> matEi;
	public:
		Particle_2d_f() {}
		~Particle_2d_f() {}

		/*poly2d_d*/
		const real func(const std::vector<real>& phi, const unsigned& p) const {
			vec8  vv = vec8::Zero();
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
							const vec8	npq = poly2d_f(dr);
							vv += w * phi[q] * npq;
						}
					}
				}
			}
			vec8 a = invMat[p] * vv;
			vec8 func = poly2d_f_0;
			return func.dot(a);
		}

		const vec func(const std::vector<vec>& u, const unsigned& p) const {
			mat83 vv = mat83::Zero();
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
							const vec8	npq = poly2d_f(dr);
							vv.block<8, 1>(0, 0) += w * u[q].x * npq;
							vv.block<8, 1>(0, 1) += w * u[q].y * npq;
							vv.block<8, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}

			mat83 a = invMat[p] * vv;
			vec8 func = poly2d_f_0;
			vec13 ret = func.transpose() * a;
			return vec(ret[0], ret[1], ret[2]);
		}

		const vec grad(const std::vector<real>& phi, const unsigned& p) const {
			vec8  vv = vec8::Zero();
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
							const vec8	npq = poly2d_f(dr);
							vv += w * phi[q] * npq;
						}
					}
				}
			}
			vec8  a = invMat[p] * vv;
			vec8 px = poly2d_f_px_0;
			vec8 pz = poly2d_f_pz_0;
			return vec(a.dot(px), 0., a.dot(pz));
		}

		const mat grad(const std::vector<vec>& u, const unsigned& p) const {
			mat83 vv = mat83::Zero();
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
							const vec8	npq = poly2d_f(dr);
							vv.block<8, 1>(0, 0) += w * u[q].x * npq;
							vv.block<8, 1>(0, 1) += w * u[q].y * npq;
							vv.block<8, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			mat83 a = invMat[p] * vv;
			vec8 px = poly2d_f_px_0;
			vec8 pz = poly2d_f_pz_0;
			return mat(vec(px.dot(a.block<8, 1>(0, 0)), 0., px.dot(a.block<8, 1>(0, 2))),
						vec(0.),
						vec(pz.dot(a.block<8, 1>(0, 0)), 0., pz.dot(a.block<8, 1>(0, 2))));
		}

		const real div(const std::vector<vec>& u, const unsigned& p) const {
			mat83 vv = mat83::Zero();
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
							const vec8	npq = poly2d_f(dr);
							vv.block<8, 1>(0, 0) += w * u[q].x * npq;
							vv.block<8, 1>(0, 1) += w * u[q].y * npq;
							vv.block<8, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			mat83 a = invMat[p] * vv;
			vec8 px = poly2d_f_px_0;
			vec8 pz = poly2d_f_pz_0;
			return a.block<8, 1>(0, 0).dot(px) + 0. + a.block<8, 1>(0, 2).dot(pz);
		}

		const vec lap(const std::vector<vec>& u, const unsigned& p) const {
			mat83 vv = mat83::Zero();
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
							const vec8	npq = poly2d_f(dr);
							vv.block<8, 1>(0, 0) += w * u[q].x * npq;
							vv.block<8, 1>(0, 1) += w * u[q].y * npq;
							vv.block<8, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			mat83 a = invMat[p] * vv;
			vec8 lapxz = poly2d_f_lap_0;
			return vec(a.block<8, 1>(0, 0).dot(lapxz), 0., a.block<8, 1>(0, 2).dot(lapxz));
		}

		const real func(const std::vector<real>& phi, const vec& p) const {
			mat88 mm = mat88::Zero();
			vec8  vv = vec8::Zero();
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
							const vec8	npq = poly2d_f(dr);
							vv += w * phi[q] * npq;
							mm.block<1, 8>(0, 0) += w * npq[0] * npq;
							mm.block<1, 8>(1, 0) += w * npq[1] * npq;
							mm.block<1, 8>(2, 0) += w * npq[2] * npq;
							mm.block<1, 8>(3, 0) += w * npq[3] * npq;
							mm.block<1, 8>(4, 0) += w * npq[4] * npq;
							mm.block<1, 8>(5, 0) += w * npq[5] * npq;
							mm.block<1, 8>(6, 0) += w * npq[6] * npq;
							mm.block<1, 8>(7, 0) += w * npq[7] * npq;
						}
					}
				}
			}
			mat88 inv = mat88::Zero();
			if (abs(mm.determinant()) < eps_mat) {
				mat44 mm_ = mm.block<4, 4>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					invMat[p] = mat88::Zero();
					continue;
				}
				else {
					inv.block<4, 4>(0, 0) = mm_.inverse();
				}
			}
			else {
				inv = mm.inverse();
			}
			vec8 a = inv* vv;
			vec8 func = poly2d_f_0;
			return func.dot(a);
		}

		const vec func(const std::vector<vec>& u, const vec& p) const {
			mat88 mm = mat88::Zero();
			mat83 vv = mat83::Zero();
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
							const vec8	npq = poly2d_f(dr);
							vv.block<8, 1>(0, 0) += w * u[q].x * npq;
							vv.block<8, 1>(0, 1) += w * u[q].y * npq;
							vv.block<8, 1>(0, 2) += w * u[q].z * npq;
							mm.block<1, 8>(0, 0) += w * npq[0] * npq;
							mm.block<1, 8>(1, 0) += w * npq[1] * npq;
							mm.block<1, 8>(2, 0) += w * npq[2] * npq;
							mm.block<1, 8>(3, 0) += w * npq[3] * npq;
							mm.block<1, 8>(4, 0) += w * npq[4] * npq;
							mm.block<1, 8>(5, 0) += w * npq[5] * npq;
							mm.block<1, 8>(6, 0) += w * npq[6] * npq;
							mm.block<1, 8>(7, 0) += w * npq[7] * npq;
						}
					}
				}
			}
			mat88 inv = mat88::Zero();
			if (abs(mm.determinant()) < eps_mat) {
				mat44 mm_ = mm.block<4, 4>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					inv = mat88::Zero();
				}
				else {
					inv.block<4, 4>(0, 0) = mm_.inverse().block<4, 4>(0, 0);
				}
			}
			else {
				inv = mm.inverse();
			}
			mat83 a = inv* vv;
			vec8 func = poly2d_f_0;
			vec13 ret = func.transpose() * a;
			return vec(ret[0], ret[1], ret[2]);
		}

		void initInvMat() {
			invMat.clear();
			for (int p = 0; p < int(np); p++) {
				invMat.push_back(mat88());
			}
		}
		void updateInvMat() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(np); p++) {
				mat88 mm = mat88::Zero();
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
								const vec8	npq = poly2d_f(dr);
								mm.block<1, 8>(0, 0) += w * npq[0] * npq;
								mm.block<1, 8>(1, 0) += w * npq[1] * npq;
								mm.block<1, 8>(2, 0) += w * npq[2] * npq;
								mm.block<1, 8>(3, 0) += w * npq[3] * npq;
								mm.block<1, 8>(4, 0) += w * npq[4] * npq;
								mm.block<1, 8>(5, 0) += w * npq[5] * npq;
								mm.block<1, 8>(6, 0) += w * npq[6] * npq;
								mm.block<1, 8>(7, 0) += w * npq[7] * npq;
							}
						}
					}
				}
				if (abs(mm.determinant()) < eps_mat) {
					mat44 mm_ = mm.block<4, 4>(0, 0);
					if (abs(mm_.determinant()) < eps_mat) {
						invMat[p] = mat88::Zero();
						continue;
					}
					invMat[p] = mat88::Zero();
					invMat[p].block<4, 4>(0, 0) = mm_.inverse().block<4, 4>(0, 0);
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
				diver[p] = div_mls_poly2d_f(vel1, p);
			}
		}

	public:
		std::vector<mat88> invMat;

	};

}