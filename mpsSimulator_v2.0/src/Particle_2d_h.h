/*
*/
#pragma once
#include "header.h"
#include "Particle.h"

namespace SIM {

	template <typename real, enum Dim dim>
	class Particle_2d_h : public Particle< real, dim, Particle_2d_h<real, dim> > {
		typedef Vec3<real> vec;
		typedef Mat3<real> mat;
		typedef Eigen::Matrix<real, 1, 3> vec13;
		typedef Eigen::Matrix<real, 6, 1> vecp;
		typedef Eigen::Matrix<real, 6, 3> matp3;
		typedef Eigen::Matrix<real, 3, 3> mat33;
		typedef Eigen::Matrix<real, 6, 6> matpp;
		typedef Eigen::Matrix<real, dim, dim> matEi;
	public:
		Particle_2d_h() : Particle() {}
		~Particle_2d_h() {}

		inline const vecp poly(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret << 1.,
				s.x,
				s.z,
				s.x*s.z,
				s.x*s.x,
				s.z*s.z;
			return ret;
		}
		inline const vecp poly_px(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret << 0.,
				1. / varrho[0],
				0.,
				s.z / varrho[0],
				2.*s.x / varrho[0],
				0.;
			return ret;
		}
		inline const vecp poly_pz(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret << 0.,
				0.,
				1. / varrho[0],
				s.x / varrho[0],
				0.,
				2.*s.z / varrho[0];
			return ret;
		}
		inline const vecp poly_lap(const vec& v) const {
			vecp ret;
			vec s = omega[0] * (v);
			ret << 0.,
				0.,
				0.,
				0.,
				2. / (varrho[0] * varrho[0]),
				2. / (varrho[0] * varrho[0]);
			return ret;
		}

		const real func(const std::vector<real>& phi, const unsigned& p) const {
			vecp  vv = vecp::Zero();
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
							const vecp	npq = poly(dr);
							vv += w * phi[q] * npq;
						}
					}
				}
			}
			vecp a = invMat[p] * vv;
			vecp func = poly_0;
			return func.dot(a);
		}

		const vec func(const std::vector<vec>& u, const unsigned& p) const {
			matp3 vv = matp3::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							vv.block<6, 1>(0, 0) += w * u[q].x * npq;
							vv.block<6, 1>(0, 1) += w * u[q].y * npq;
							vv.block<6, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}

			auto a = invMat[p] * vv;
			auto func = poly_0;
			auto ret = func.transpose() * a;
			return vec(ret[0], ret[1], ret[2]);
		}

		const vec grad(const std::vector<real>& phi, const unsigned& p) const {
			vecp  vv = vecp::Zero();
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
							const vecp	npq = poly(dr);
							vv += w * phi[q] * npq;
						}
					}
				}
			}
			vecp  a = invMat[p] * vv;
			vecp px = poly_px_0;
			vecp pz = poly_pz_0;
			return vec(px.dot(a), 0., pz.dot(a));
		}

		const mat grad(const std::vector<vec>& u, const unsigned& p) const {
			matp3 vv = matp3::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							vv.block<6, 1>(0, 0) += w * u[q].x * npq;
							vv.block<6, 1>(0, 1) += w * u[q].y * npq;
							vv.block<6, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			auto a = invMat[p] * vv;
			auto px = poly_px_0;
			auto pz = poly_pz_0;
			return mat( vec(px.dot(a.block<6, 1>(0, 0)), 0., px.dot(a.block<6, 1>(0, 2))),
						vec(0.),
						vec(pz.dot(a.block<6, 1>(0, 0)), 0., pz.dot(a.block<6, 1>(0, 2))) );
		}

		const real div(const std::vector<vec>& u, const unsigned& p) const {
			matp3 vv = matp3::Zero();
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
							const vecp	npq = poly(dr);
							vv.block<6, 1>(0, 0) += w * u[q].x * npq;
							vv.block<6, 1>(0, 1) += w * u[q].y * npq;
							vv.block<6, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			matp3 a = invMat[p] * vv;
			vecp px = poly_px_0;
			vecp pz = poly_pz_0;
			return a.block<6, 1>(0, 0).dot(px) + 0. + a.block<6, 1>(0, 2).dot(pz);
		}

		const real lap(const std::vector<real>& phi, const unsigned& p) const {
			vecp vv = vecp::Zero();
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
							const vecp	npq = poly(dr);
							vv += w * phi[q] * npq;
						}
					}
				}
			}
			vecp a = invMat[p] * vv;
			vecp lap = poly_lap_0;
			return lap.dot(a);
		}

		const vec lap(const std::vector<vec>& u, const unsigned& p) const {
			matp3 vv = matp3::Zero();
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
							const vecp	npq = poly(dr);
							vv.block<6, 1>(0, 0) += w * u[q].x * npq;
							vv.block<6, 1>(0, 1) += w * u[q].y * npq;
							vv.block<6, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			matp3 a = invMat[p] * vv;
			vecp lap = poly_lap_0;
			return vec(a.block<6, 1>(0, 0).dot(lap), 0., a.block<6, 1>(0, 2).dot(lap));
		}

		const vec rot(const std::vector<vec>& u, const unsigned& p) const {
			matp3 vv = matp3::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							vv.block<6, 1>(0, 0) += w * u[q].x * npq;
							vv.block<6, 1>(0, 1) += w * u[q].y * npq;
							vv.block<6, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			auto a = invMat[p] * vv;
			auto px = poly_px_0;
			auto pz = poly_pz_0;
			return vec(0., px.dot(a.block<6, 1>(0, 2)) - pz.dot(a.block<6, 1>(0, 0)), 0.);
		}

		const real func(const std::vector<real>& phi, const vec& p) const {
			matpp mm = matpp::Zero();
			vecp  vv = vecp::Zero();
			const iVec3 c = cell->iCoord(p);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							//if (type[q] == BD2) continue;
							//if (team[p] != team[q]) continue;
							if (type[q] == BD1 && isFs(q)) continue;
#endif
							const vec	dr = pos[q] - p;
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vecp	npq = poly(dr);
							vv += w * phi[q] * npq;
							mm.block<1, 6>(0, 0) += w * npq[0] * npq;
							mm.block<1, 6>(1, 0) += w * npq[1] * npq;
							mm.block<1, 6>(2, 0) += w * npq[2] * npq;
							mm.block<1, 6>(3, 0) += w * npq[3] * npq;
							mm.block<1, 6>(4, 0) += w * npq[4] * npq;
							mm.block<1, 6>(5, 0) += w * npq[5] * npq;
						}
					}
				}
			}
			matpp inv = matpp::Zero();
			if (abs(mm.determinant()) < eps_mat) {
				mat33 mm_ = mm.block<3, 3>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					real mm__ = mm(0, 0);
					if (abs(mm__) < eps_mat) {
						inv = matpp::Zero();
					}
					else inv(0, 0) = 1. / mm__;
				}
				else {
					inv.block<3, 3>(0, 0) = mm_.inverse();
				}
			}
			else {
				inv = mm.inverse();
			}
			vecp a = inv* vv;
			vecp func = poly_0;
			return func.dot(a);
		}

		const vec func(const std::vector<vec>& u, const vec& p) const {
			matpp mm = matpp::Zero();
			matp3 vv = matp3::Zero();
			const auto c = cell->iCoord(p);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							//if (type[q] == BD2) continue;
							//if (team[p] != team[q]) continue;
							if (type[q] == BD1 && isFs(q)) continue;
#endif
							const auto dr = pos[q] - p;
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							vv.block<6, 1>(0, 0) += w * u[q].x * npq;
							vv.block<6, 1>(0, 1) += w * u[q].y * npq;
							vv.block<6, 1>(0, 2) += w * u[q].z * npq;
							mm.block<1, 6>(0, 0) += w * npq[0] * npq;
							mm.block<1, 6>(1, 0) += w * npq[1] * npq;
							mm.block<1, 6>(2, 0) += w * npq[2] * npq;
							mm.block<1, 6>(3, 0) += w * npq[3] * npq;
							mm.block<1, 6>(4, 0) += w * npq[4] * npq;
							mm.block<1, 6>(5, 0) += w * npq[5] * npq;
						}
					}
				}
			}
			matpp inv = matpp::Zero();
			if (abs(mm.determinant()) < eps_mat) {
				mat33 mm_ = mm.block<3, 3>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					real mm__ = mm(0, 0);
					if (abs(mm__) < eps_mat) {
						inv = matpp::Zero();
					}
					else inv(0, 0) = 1. / mm__;
				}
				else {
					inv.block<3, 3>(0, 0) = mm_.inverse();
				}
			}
			else {
				inv = mm.inverse();
			}
			auto a = inv* vv;
			auto func = poly_0;
			auto ret = func.transpose() * a;
			return vec(ret[0], ret[1], ret[2]);
		}

		const vec func_nobd2(const std::vector<vec>& u, const vec& p) const {
			matpp mm = matpp::Zero();
			matp3 vv = matp3::Zero();
			const iVec3 c = cell->iCoord(p);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (type[q] == BD2) continue;
							//if (team[p] != team[q]) continue;
							if (type[q] == BD1 && isFs(q)) continue;
#endif
							const vec	dr = pos[q] - p;
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vecp	npq = poly2d_h(dr);
							vv.block<6, 1>(0, 0) += w * u[q].x * npq;
							vv.block<6, 1>(0, 1) += w * u[q].y * npq;
							vv.block<6, 1>(0, 2) += w * u[q].z * npq;
							mm.block<1, 6>(0, 0) += w * npq[0] * npq;
							mm.block<1, 6>(1, 0) += w * npq[1] * npq;
							mm.block<1, 6>(2, 0) += w * npq[2] * npq;
							mm.block<1, 6>(3, 0) += w * npq[3] * npq;
							mm.block<1, 6>(4, 0) += w * npq[4] * npq;
							mm.block<1, 6>(5, 0) += w * npq[5] * npq;
						}
					}
				}
			}
			matpp inv = matpp::Zero();
			if (abs(mm.determinant()) < eps_mat) {
				mat33 mm_ = mm.block<3, 3>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					real mm__ = mm(0, 0);
					if (abs(mm__) < eps_mat) {
						inv = matpp::Zero();
					}
					else inv(0, 0) = 1. / mm__;
				}
				else {
					inv.block<3, 3>(0, 0) = mm_.inverse();
				}
			}
			else {
				inv = mm.inverse();
			}
			matp3 a = inv* vv;
			vecp func = poly_0;
			vec13 ret = func.transpose() * a;
			return vec(ret[0], ret[1], ret[2]);
		}

		const mat grad_upwind(const std::vector<vec>& u, const unsigned& p) const {
			matpp mm = matpp::Zero();
			matp3 vv = matp3::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto up = u[p];
							if (up* dr > 0.) continue;
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							vv.block<6, 1>(0, 0) += w * u[q].x * npq;
							vv.block<6, 1>(0, 1) += w * u[q].y * npq;
							vv.block<6, 1>(0, 2) += w * u[q].z * npq;
							mm.block<1, 6>(0, 0) += w * npq[0] * npq;
							mm.block<1, 6>(1, 0) += w * npq[1] * npq;
							mm.block<1, 6>(2, 0) += w * npq[2] * npq;
							mm.block<1, 6>(3, 0) += w * npq[3] * npq;
							mm.block<1, 6>(4, 0) += w * npq[4] * npq;
							mm.block<1, 6>(5, 0) += w * npq[5] * npq;
						}
					}
				}
			}
			matpp inv = matpp::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << mm.determinant() << std::endl;
#endif
				mat33 mm_ = mm.block<3, 3>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					real mm__ = mm(0, 0);
					if (abs(mm__) < eps_mat) {
						inv = matpp::Zero();
					}
					else inv(0, 0) = 1. / mm__;
				}
				else {
					inv.block<3, 3>(0, 0) = mm_.inverse();
				}
			}
			else {
				inv = mm.inverse();
			}
			auto a = inv * vv;
			auto px = poly_px_0;
			auto pz = poly_pz_0;
			return mat(	vec(px.dot(a.block<6, 1>(0, 0)), 0., px.dot(a.block<6, 1>(0, 2))),
						vec(0.),
						vec(pz.dot(a.block<6, 1>(0, 0)), 0., pz.dot(a.block<6, 1>(0, 2))) );
		}

		void updateInvMat() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(np); p++) {
				if (type[p] == BD2) continue;
				matpp mm = matpp::Zero();
				const auto c = cell->iCoord(pos[p]);
				for (auto k = -1; k <= 1; k++) {
					for (auto j = -1; j <= 1; j++) {
						for (auto i = -1; i <= 1; i++) {
							const auto ne = c + iVec3(i, j, k);
							const auto key = cell->hash(ne);
							for (auto m = 0; m < cell->linkList[key].size(); m++) {
								const auto q = cell->linkList[key][m];
#if BD_OPT
								if (bdOpt(p, q)) continue;
#endif
								const auto dr = pos[q] - pos[p];
								const auto dr1 = dr.mag();
								if (dr1 > r0) continue;
								const auto w = w3(dr1);
								const auto npq = poly(dr);
								mm.block<1, 6>(0, 0) += w * npq[0] * npq;
								mm.block<1, 6>(1, 0) += w * npq[1] * npq;
								mm.block<1, 6>(2, 0) += w * npq[2] * npq;
								mm.block<1, 6>(3, 0) += w * npq[3] * npq;
								mm.block<1, 6>(4, 0) += w * npq[4] * npq;
								mm.block<1, 6>(5, 0) += w * npq[5] * npq;
							}
						}
					}
				}
				invMat[p] = matpp::Zero();
				if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
					std::cout << mm.determinant() << std::endl;
#endif
					auto mm_ = mm.block<3, 3>(0, 0);
					if (abs(mm_.determinant()) < eps_mat) {
						auto mm__ = mm(0, 0);
						if (abs(mm__) < eps_mat) {
							invMat[p] = matpp::Zero();
							continue;
						}
						invMat[p](0, 0) = 1. / mm__;
						continue;
					}
					invMat[p].block<3, 3>(0, 0) = mm_.inverse();
					continue;
				}
				invMat[p] = mm.inverse();
			}
		}

		void init2d_x() {
			invMat.clear();
			for (int p = 0; p < int(np); p++) {
				invMat.push_back(matpp());
			}

			varrho.clear();
			varrho.push_back(1.*dp);
			poly_0 = poly(vec(0.));
			poly_px_0 = poly_px(vec(0.));
			poly_pz_0 = poly_pz(vec(0.));
			poly_lap_0 = poly_lap(vec(0.));
		}
	public:
		std::vector<matpp> invMat;

		std::vector<real> varrho;
		vecp poly_0;
		vecp poly_px_0;
		vecp poly_pz_0;
		vecp poly_lap_0;
	};

}