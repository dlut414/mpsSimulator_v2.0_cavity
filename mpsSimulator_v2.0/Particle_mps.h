/*
*/
#pragma once
#include "header.h"
#include "Particle.h"

namespace SIM {

	template <typename real, enum Dim dim>
	class Particle_mps : public Particle< real, dim, Particle_mps<real, dim> > {
		typedef Vec3<real> vec;
		typedef Mat3<real> mat;
		typedef Eigen::Matrix<real, 1, 3> vec13;
		typedef Eigen::Matrix<real, 5, 1> vecp;
		typedef Eigen::Matrix<real, 5, 3> matp3;
		typedef Eigen::Matrix<real, 3, 3> mat33;
		typedef Eigen::Matrix<real, 5, 5> matpp;
	public:
		Particle_mps() : Particle() {}
		~Particle_mps() {}

		inline const vecp poly(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
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
			ret <<
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
			ret <<
				0.,
				1. / varrho[0],
				s.x / varrho[0],
				0.,
				2.*s.z / varrho[0];
			return ret;
		}
		inline const vecp poly_pxx(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
				0.,
				0.,
				0.,
				2. / (varrho[0]*varrho[0]),
				0.;
			return ret;
		}
		inline const vecp poly_pzz(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
				0.,
				0.,
				0.,
				0.,
				2. / (varrho[0]*varrho[0]);
			return ret;
		}
		inline const vecp poly_pxz(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
				0.,
				0.,
				1. / (varrho[0]*varrho[0]),
				0.,
				0.;
			return ret;
		}

		const real func(const std::vector<real>& phi, const unsigned& p) const {
			return phi[p];
		}

		const vec func(const std::vector<vec>& u, const unsigned& p) const {
			return u[p];
		}

		const vec grad_hat(const std::vector<real>& phi, const unsigned& p) const {
			auto ret = vec(0.);
			auto hat = phi[p];
			const auto c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							if (type[q] == BD2) continue;
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr1 = (pos[q] - pos[p]).mag();
							if (dr1 > r0) continue;
							if (phi[q] < hat) hat = phi[q];
						}
					}
				}
			}
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
							if (type[q] == BD2 || p == q) continue;
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr2 = dr.mag2();
							const auto dr1 = sqrt(dr2);
							if (dr1 > r0) continue;
							ret += (phi[q] - hat) * w1(dr1) / dr2 * dr;
						}
					}
				}
			}
			ret = dim / n0 * ret;
			return ret;
		}

		const vec grad(const std::vector<real>& phi, const unsigned& p) const {
			auto ret = vec(0.);
			const auto c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
							if (type[q] == BD2 || p == q) continue;
#if BD_OPT
		if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr2 = dr.mag2();
							const auto dr1 = sqrt(dr2);
							if (dr1 > r0) continue;
							ret += (phi[q] - phi[p]) * w1(dr1) / dr2 * dr;
						}
					}
				}
			}
			ret = dim / n0 * ret;
			return ret;
		}

		const mat grad(const std::vector<vec>& u, const unsigned& p) const {
			auto ret = mat(vec(0., 0., 0.), vec(0., 0., 0.), vec(0., 0., 0.));
			const auto c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
							if (type[q] == BD2 || p == q) continue;
#if BD_OPT
		if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr2 = dr.mag2();
							const auto dr1 = sqrt(dr2);
							if (dr1 > r0) continue;
							const auto tmp = w1(dr1) / dr2 * dr;
							ret.x += (u[q].x - u[p].x) * tmp;
							ret.y += (u[q].y - u[p].y) * tmp;
							ret.z += (u[q].z - u[p].z) * tmp;
						}
					}
				}
			}
			ret = dim / n0 * ret;
			return ret.trans();
		}

		const real div(const std::vector<vec>& u, const unsigned& p) const {
			auto ret = 0.;
			const auto c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
							if (type[q] == BD2 || p == q) continue;
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr2 = dr.mag2();
							const auto dr1 = sqrt(dr2);
							if (dr1 > r0) continue;
							const auto dv = u[q] - u[p];
							ret += w1(dr1) / dr2 * (dv * dr);
						}
					}
				}
			}
			ret = (2. * dim / n0) * ret;
			return ret;
		}

		const real lap(const std::vector<real>& phi, const unsigned& p) const {
			auto ret = 0.;
			const auto c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
							if (p == q) continue;
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							ret += (phi[q] - phi[p]) * w1(dr1);
						}
					}
				}
			}
			ret = 2. * dim / (n0 * lambda) * ret;
			return ret;
		}

		const vec lap(const std::vector<vec>& u, const unsigned& p) const {
			auto ret = vec(0.);
			const auto c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
							if (p == q) continue;
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							ret += (u[q] - u[p]) * w1(dr1);
						}
					}
				}
			}
			ret = 2. * dim / (n0 * lambda) * ret;
			return ret;
		}

		const vec rot(const std::vector<vec>& u, const unsigned& p) const {
			auto mm = grad(u, p);
			return vec(0., mm.x.z - mm.z.x, 0.);
		}

		const real func(const std::vector<real>& phi, const vec& p) const {
			auto ret = real(0.);
			auto ww = real(0.);
			const auto c = cell->iCoord(p);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							if (type[q] == BD2) continue;
#if BD_OPT
							if (bdOpt(q)) continue;
#endif
							const auto dr1 = (pos[q] - p).mag();
							if (dr1 > r0) continue;
							const auto w = w1(dr1);
							ww += w;
							ret += w * phi[q];
						}
					}
				}
			}
			if (abs(ww) < eps) ww = 1.;
			return ret / ww;
		}

		const vec func(const std::vector<vec>& u, const vec& p) const {
			auto ret = vec(0.);
			auto ww = real(0.);
			const auto c = cell->iCoord(p);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							if (type[q] == BD2) continue;
#if BD_OPT
							if (bdOpt(q)) continue;
#endif
							const auto dr1 = (pos[q] - p).mag();
							if (dr1 > r0) continue;
							const auto w = w1(dr1);
							ww += w;
							ret += w * u[q];
						}
					}
				}
			}
			if (abs(ww) < eps) ww = 1.;
			return ret / ww;
		}

		const vec func_mafl(const std::vector<vec>& u, const unsigned& p, const vec& p_new) const {
			const auto re = 1.5* dp;
			const auto dx = 1.5* dp;
			const auto p_i = pos[p];
			const auto dmove = p_new - p_i;
			const auto up = dmove.norm();
			vec pLocal[6];

			for (int i = -3; i <= 2; i++) {
				pLocal[i + 3] = p_i - i*dx*up;
			}

			vec ret[6];
			ret[3] = u[p];
			for (int fp = 1; fp <= 4; fp++) {
				if (fp == 3) continue;
				ret[fp] = vec(0.);
				auto ww = real(0.);
				auto c = cell->iCoord(pLocal[fp]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const auto ne = c + iVec3(i, j, k);
							const auto key = cell->hash(ne);
							for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
								const unsigned q = cell->linkList[key][m];
								if (type[q] == BD2) continue;
#if BD_OPT
								if (bdOpt(q)) continue;
#endif
								const auto dr1 = (pos[q] - pLocal[fp]).mag();
								const auto dr1_m1 = (pos[q] - pLocal[fp-1]).mag();
								const auto dr1_p1 = (pos[q] - pLocal[fp+1]).mag();
								if (dr1 > re) continue;
								if (dr1 > dr1_m1 || dr1 > dr1_p1) continue;
								const auto w = w1(dr1);
								ww += w;
								ret[fp] += w * u[q];
							}
						}
					}
				}
				if (abs(ww) < eps) ww = 1.;
				ret[fp] = ret[fp] / ww;
			}
			return ret[3] - (dmove.mag()/dx)* (0.125* ret[1] - 0.875* ret[2] + 0.375* ret[3] + 0.375* ret[4]);
		}

		const vec func_mafl_mmt(const std::vector<vec>& u, const unsigned& p, const vec& p_new) const {
			const auto re = 1.5* dp;
			const auto dx = 1.5* dp;
			const auto p_i = pos[p];
			const auto dmove = p_new - p_i;
			const auto up = dmove.norm();
			vec pLocal[6];

			for (int i = -3; i <= 2; i++) {
				pLocal[i + 3] = p_i - i*dx*up;
			}

			vec ret[6];
			ret[3] = u[p];
			for (int fp = 1; fp <= 4; fp++) {
				if (fp == 3) continue;
				ret[fp] = vec(0.);
				auto ww = real(0.);
				auto c = cell->iCoord(pLocal[fp]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const auto ne = c + iVec3(i, j, k);
							const auto key = cell->hash(ne);
							for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
								const unsigned q = cell->linkList[key][m];
								if (type[q] == BD2) continue;
#if BD_OPT
								if (bdOpt(q)) continue;
#endif
								const auto dr1 = (pos[q] - pLocal[fp]).mag();
								const auto dr1_m1 = (pos[q] - pLocal[fp - 1]).mag();
								const auto dr1_p1 = (pos[q] - pLocal[fp + 1]).mag();
								if (dr1 > re) continue;
								if (dr1 > dr1_m1 || dr1 > dr1_p1) continue;
								const auto w = w1(dr1);
								ww += w;
								ret[fp] += w * u[q];
							}
						}
					}
				}
				if (abs(ww) < eps) continue;
				ret[fp] = ret[fp] / ww;
			}
			auto ret_mmt = ret[3] - (dmove.mag() / dx)* (0.125* ret[1] - 0.875* ret[2] + 0.375* ret[3] + 0.375* ret[4]);
			auto ret_min = ret[1];
			auto ret_max = ret[1];
			for (int fp = 1; fp <= 4; fp++) {
				if (ret[fp].mag2() < ret_min.mag2()) ret_min = ret[fp];
				if (ret[fp].mag2() > ret_max.mag2()) ret_max = ret[fp];
			}
			if (ret_mmt.mag2() < ret_min.mag2()) ret_mmt = ret_min* ret_mmt.norm();
			if (ret_mmt.mag2() > ret_max.mag2()) ret_mmt = ret_max* ret_mmt.norm();
			return ret_mmt;
		}

		const vec func_mls_a(const std::vector<vec>& u, const unsigned& p, const vec& p_new) const {
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
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							mm.block<1, 5>(0, 0) += w * npq[0] * npq;
							mm.block<1, 5>(1, 0) += w * npq[1] * npq;
							mm.block<1, 5>(2, 0) += w * npq[2] * npq;
							mm.block<1, 5>(3, 0) += w * npq[3] * npq;
							mm.block<1, 5>(4, 0) += w * npq[4] * npq;
							vv.block<5, 1>(0, 0) += w * (u[q].x - u[p].x) * npq;
							vv.block<5, 1>(0, 1) += w * (u[q].y - u[p].y) * npq;
							vv.block<5, 1>(0, 2) += w * (u[q].z - u[p].z) * npq;
						}
					}
				}
			}
			matpp inv = matpp::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					inv = matpp::Zero();
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();
			
			const auto a = inv * vv;
			const auto ax = a.block<5, 1>(0, 0);
			const auto az = a.block<5, 1>(0, 2);
			const auto px = poly_px_0;
			const auto pz = poly_pz_0;
			const auto pxx = poly_pxx_0;
			const auto pzz = poly_pzz_0;
			const auto pxz = poly_pxz_0;
			const auto dp = p_new - pos[p];
			auto ret = u[p];
			ret.x += dp.x* (px.dot(ax)) + dp.z* (pz.dot(ax))
				+ 0.5* dp.x* dp.x* (pxx.dot(ax)) + dp.x* dp.z* (pxz.dot(ax)) + 0.5* dp.z* dp.z* (pzz.dot(ax));
			ret.z += dp.x* (px.dot(az)) + dp.z* (pz.dot(az))
				+ 0.5* dp.x* dp.x* (pxx.dot(az)) + dp.x* dp.z* (pxz.dot(az)) + 0.5* dp.z* dp.z* (pzz.dot(az));
			return ret;
		}

		const vec func_mls_a_upwind(const std::vector<vec>& u, const unsigned& p, const vec& p_new) const {
			const auto dp = p_new - pos[p];
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
							if (dr*dp < 0) continue;
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							mm.block<1, 5>(0, 0) += w * npq[0] * npq;
							mm.block<1, 5>(1, 0) += w * npq[1] * npq;
							mm.block<1, 5>(2, 0) += w * npq[2] * npq;
							mm.block<1, 5>(3, 0) += w * npq[3] * npq;
							mm.block<1, 5>(4, 0) += w * npq[4] * npq;
							vv.block<5, 1>(0, 0) += w * (u[q].x - u[p].x) * npq;
							vv.block<5, 1>(0, 1) += w * (u[q].y - u[p].y) * npq;
							vv.block<5, 1>(0, 2) += w * (u[q].z - u[p].z) * npq;
						}
					}
				}
			}
			matpp inv = matpp::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					inv = matpp::Zero();
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto ax = a.block<5, 1>(0, 0);
			const auto az = a.block<5, 1>(0, 2);
			const auto px = poly_px_0;
			const auto pz = poly_pz_0;
			const auto pxx = poly_pxx_0;
			const auto pzz = poly_pzz_0;
			const auto pxz = poly_pxz_0;
			auto ret = u[p];
			ret.x += dp.x* (px.dot(ax)) + dp.z* (pz.dot(ax))
				+ 0.5* dp.x* dp.x* (pxx.dot(ax)) + dp.x* dp.z* (pxz.dot(ax)) + 0.5* dp.z* dp.z* (pzz.dot(ax));
			ret.z += dp.x* (px.dot(az)) + dp.z* (pz.dot(az))
				+ 0.5* dp.x* dp.x* (pxx.dot(az)) + dp.x* dp.z* (pxz.dot(az)) + 0.5* dp.z* dp.z* (pzz.dot(az));
			return ret;
		}

		const vec func_mls_b(const std::vector<vec>& u, const unsigned& p, const vec& p_new) const {
			matpp mm = matpp::Zero();
			matp3 vv = matp3::Zero();
			const auto c = cell->iCoord(p_new);
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
							const auto dr = pos[q] - p_new;
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							mm.block<1, 5>(0, 0) += w * npq[0] * npq;
							mm.block<1, 5>(1, 0) += w * npq[1] * npq;
							mm.block<1, 5>(2, 0) += w * npq[2] * npq;
							mm.block<1, 5>(3, 0) += w * npq[3] * npq;
							mm.block<1, 5>(4, 0) += w * npq[4] * npq;
							vv.block<5, 1>(0, 0) += w * u[q].x * npq;
							vv.block<5, 1>(0, 1) += w * u[q].y * npq;
							vv.block<5, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			matpp inv = matpp::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					inv = matpp::Zero();
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto ret = poly_0.transpose()* a;
			return vec(ret[0], ret[1], ret[2]);
		}

		const vec func_mls_b_upwind(const std::vector<vec>& u, const unsigned& p, const vec& p_new) const {
			const auto dp = p_new - pos[p];
			matpp mm = matpp::Zero();
			matp3 vv = matp3::Zero();
			const auto c = cell->iCoord(p_new);
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
							const auto dr = pos[q] - p_new;
							if (dr*dp < 0) continue;
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							mm.block<1, 5>(0, 0) += w * npq[0] * npq;
							mm.block<1, 5>(1, 0) += w * npq[1] * npq;
							mm.block<1, 5>(2, 0) += w * npq[2] * npq;
							mm.block<1, 5>(3, 0) += w * npq[3] * npq;
							mm.block<1, 5>(4, 0) += w * npq[4] * npq;
							vv.block<5, 1>(0, 0) += w * u[q].x * npq;
							vv.block<5, 1>(0, 1) += w * u[q].y * npq;
							vv.block<5, 1>(0, 2) += w * u[q].z * npq;
						}
					}
				}
			}
			matpp inv = matpp::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					inv = matpp::Zero();
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto ret = poly_0.transpose()* a;
			return vec(ret[0], ret[1], ret[2]);
		}

		void updateDiver() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(np); p++) {
				diver[p] = div(part->vel1, p);
			}
		}

		void init2d_x() {
			varrho.clear();
			varrho.push_back(1.*dp);
			poly_0 = poly(vec(0.));
			poly_px_0 = poly_px(vec(0.));
			poly_pz_0 = poly_pz(vec(0.));
			poly_pxx_0 = poly_pxx(vec(0.));
			poly_pzz_0 = poly_pzz(vec(0.));
			poly_pxz_0 = poly_pxz(vec(0.));
		}

		//const vec grad_suzuki(const std::vector<real>& phi, const unsigned& p) const {
		//	mat	m3 = mat(vec(0.), vec(0.), vec(0.));
		//	vec v3 = vec(0.);
		//	const iVec3 c = cell->iCoord(pos[p]);
		//	for (int k = -1; k <= 1; k++) {
		//		for (int j = -1; j <= 1; j++) {
		//			for (int i = -1; i <= 1; i++) {
		//				const iVec3 ne = c + iVec3(i, j, k);
		//				const unsigned key = cell->hash(ne);
		//				for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
		//					const unsigned q = cell->linkList[key][m];
		//					if (q == p || type[q] == BD2) continue;
		//					//if (team[p] != team[q]) continue;
		//					//if (type[q]==BD1 && isFs(q)) continue;
		//					//if (type[p]==FLUID && isFs(q) == 2) continue;
		//					//if (type[p]==FLUID && type[q]==FLUID && pn_nb[q] < 0.1*pn0) continue;
		//					const vec	dr = pos[q] - pos[p];
		//					const real	dr1 = dr.mag();
		//					if (dr1 > r0) continue;
		//					const real  df = phi[q] - phi[p];
		//					const real  w = w1(dr1);
		//					const vec	npq = dr/dr1;
		//					m3.x += w * npq.x * npq;
		//					m3.y += w * npq.y * npq;
		//					m3.z += w * npq.z * npq;
		//					v3 += w * (df / dr1) * npq;
		//				}
		//			}
		//		}
		//	}
		//	matEi pinv;
		//	switch (dim) {
		//	case 2:
		//		pinv << m3.x.x, m3.x.z, m3.z.x, m3.z.z; break;
		//	case 3:
		//		pinv << m3.x.x, m3.x.y, m3.x.z, m3.y.x, m3.y.y, m3.y.z, m3.z.x, m3.z.y, m3.z.z; break;
		//	}
		//	if (abs(pinv.determinant()) < eps_mat) return grad_hat(phi, p); //vec(0.);
		//	matEi inv = pinv.inverse();
		//	switch (dim) {
		//	case 2:
		//		m3.x.x = inv(0, 0);	m3.x.y = 0.; m3.x.z = inv(0, 1);
		//		m3.y.x = 0.;		m3.y.y = 0.; m3.y.z = 0.;
		//		m3.z.x = inv(1, 0); m3.z.y = 0.; m3.z.z = inv(1, 1);
		//		break;
		//	case 3:
		//		m3.x.x = inv(0, 0); m3.x.y = inv(0, 1); m3.x.z = inv(0, 2);
		//		m3.y.x = inv(1, 0); m3.y.y = inv(1, 1); m3.y.z = inv(1, 2);
		//		m3.z.x = inv(2, 0); m3.z.y = inv(2, 1); m3.z.z = inv(2, 2);
		//		break;
		//	}
		//	return m3 * v3;
		//}

		//const mat grad_suzuki(const std::vector<vec>& u, const unsigned& p) const {
		//	mat	m3 = mat(vec(0.), vec(0.), vec(0.));
		//	mat v3 = mat(vec(0.), vec(0.), vec(0.));
		//	const iVec3 c = cell->iCoord(pos[p]);
		//	for (int k = -1; k <= 1; k++) {
		//		for (int j = -1; j <= 1; j++) {
		//			for (int i = -1; i <= 1; i++) {
		//				const iVec3 ne = c + iVec3(i, j, k);
		//				const unsigned key = cell->hash(ne);
		//				for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
		//					const unsigned q = cell->linkList[key][m];
		//					if (q == p || type[q] == BD2) continue;
		//					if (team[p] != team[q]) continue;
		//					if (type[q] == BD1 && isFs(q)) continue;
		//					//if (type[p]==FLUID && isFs(q)==2) continue;
		//					//if (type[p]==FLUID && type[q]==FLUID && pn_nb[q] < 0.1*pn0) continue;
		//					const vec	dr = pos[q] - pos[p];
		//					const vec  df = u[q] - u[p];
		//					const real	dr1 = dr.mag();
		//					if (dr1 > r0) continue;
		//					const real  w = w1(dr1);
		//					const vec	npq = (1. / dr1) * dr;
		//					m3.x += w * npq.x * npq;
		//					m3.y += w * npq.y * npq;
		//					m3.z += w * npq.z * npq;
		//					v3.x += w * (df.x / (dr1)) * npq;
		//					v3.y += w * (df.y / (dr1)) * npq;
		//					v3.z += w * (df.z / (dr1)) * npq;
		//				}
		//			}
		//		}
		//	}
		//	matEi pinv;
		//	switch (dim) {
		//	case 2:
		//		pinv << m3.x.x, m3.x.z, m3.z.x, m3.z.z; break;
		//	case 3:
		//		pinv << m3.x.x, m3.x.y, m3.x.z, m3.y.x, m3.y.y, m3.y.z, m3.z.x, m3.z.y, m3.z.z; break;
		//	}
		//	if (abs(pinv.determinant()) < eps_mat) return grad(u, p); //mat(0.);
		//	matEi inv = pinv.inverse();
		//	switch (dim) {
		//	case 2:
		//		m3.x.x = inv(0, 0);	m3.x.y = 0.; m3.x.z = inv(0, 1);
		//		m3.y.x = 0.;		m3.y.y = 0.; m3.y.z = 0.;
		//		m3.z.x = inv(1, 0); m3.z.y = 0.; m3.z.z = inv(1, 1);
		//		break;
		//	case 3:
		//		m3.x.x = inv(0, 0); m3.x.y = inv(0, 1); m3.x.z = inv(0, 2);
		//		m3.y.x = inv(1, 0); m3.y.y = inv(1, 1); m3.y.z = inv(1, 2);
		//		m3.z.x = inv(2, 0); m3.z.y = inv(2, 1); m3.z.z = inv(2, 2);
		//		break;
		//	}
		//	return m3 * (v3.trans());
		//}

		//const mat grad_suzuki_nb(const std::vector<vec>& u, const unsigned& p) const {
		//	mat	m3 = mat(vec(0.), vec(0.), vec(0.));
		//	mat v3 = mat(vec(0.), vec(0.), vec(0.));
		//	const iVec3 c = cell->iCoord(pos[p]);
		//	for (int k = -1; k <= 1; k++) {
		//		for (int j = -1; j <= 1; j++) {
		//			for (int i = -1; i <= 1; i++) {
		//				const iVec3 ne = c + iVec3(i, j, k);
		//				const unsigned key = cell->hash(ne);
		//				for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
		//					const unsigned q = cell->linkList[key][m];
		//					if (q == p || type[q] != FLUID) continue;
		//					const vec	dr = pos[q] - pos[p];
		//					const vec  df = u[q] - u[p];
		//					const real	dr1 = dr.mag();
		//					const real  w = w1(dr1);
		//					const vec	npq = (1. / dr1) * dr;
		//					m3.x += w * npq.x * npq;
		//					m3.y += w * npq.y * npq;
		//					m3.z += w * npq.z * npq;
		//					v3.x += w * npq * (df.x / (dr1));
		//					v3.y += w * npq * (df.y / (dr1));
		//					v3.z += w * npq * (df.z / (dr1));
		//				}
		//			}
		//		}
		//	}
		//	matEi pinv;
		//	switch (dim) {
		//	case 2:
		//		pinv << m3.x.x, m3.x.z, m3.z.x, m3.z.z; break;
		//	case 3:
		//		pinv << m3.x.x, m3.x.y, m3.x.z, m3.y.x, m3.y.y, m3.y.z, m3.z.x, m3.z.y, m3.z.z; break;
		//	}
		//	if (abs(pinv.determinant()) < eps_mat) return grad_nb(u, p);
		//	matEi inv = pinv.inverse();
		//	switch (dim) {
		//	case 2:
		//		m3.x.x = inv(0, 0);	m3.x.y = 0.; m3.x.z = inv(0, 1);
		//		m3.y.x = 0.;		m3.y.y = 0.; m3.y.z = 0.;
		//		m3.z.x = inv(1, 0); m3.z.y = 0.; m3.z.z = inv(1, 1);
		//		break;
		//	case 3:
		//		m3.x.x = inv(0, 0); m3.x.y = inv(0, 1); m3.x.z = inv(0, 2);
		//		m3.y.x = inv(1, 0); m3.y.y = inv(1, 1); m3.y.z = inv(1, 2);
		//		m3.z.x = inv(2, 0); m3.z.y = inv(2, 1); m3.z.z = inv(2, 2);
		//		break;
		//	}
		//	return m3 * (v3.trans());
		//}

		//const vec grad_suzuki_pnd(const std::vector<real>& phi, const unsigned& p) const {
		//	mat	m3 = mat(vec(0.), vec(0.), vec(0.));
		//	vec v3 = vec(0.);
		//	const iVec3 c = cell->iCoord(pos[p]);
		//	for (int k = -1; k <= 1; k++) {
		//		for (int j = -1; j <= 1; j++) {
		//			for (int i = -1; i <= 1; i++) {
		//				const iVec3 ne = c + iVec3(i, j, k);
		//				const unsigned key = cell->hash(ne);
		//				for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
		//					const unsigned q = cell->linkList[key][m];
		//					if (q == p) continue;
		//					const vec	dr = pos[q] - pos[p];
		//					const real	dr1 = dr.mag();
		//					if (dr1 > r0) continue;
		//					const real corp = nbd[p] > 0 ? real(pn0) / real(pn[p]) : 1.;
		//					const real corq = nbd[q] > 0 ? real(pn0) / real(pn[q]) : 1.;
		//					const real  df = phi[q]*corq - phi[p]*corp;
		//					const real  w = w1(dr1);
		//					const vec	npq = dr / dr1;
		//					m3.x += w * npq.x * npq;
		//					m3.y += w * npq.y * npq;
		//					m3.z += w * npq.z * npq;
		//					v3 += w * (df / dr1) * npq;
		//				}
		//			}
		//		}
		//	}
		//	matEi pinv;
		//	switch (dim) {
		//	case 2:
		//		pinv << m3.x.x, m3.x.z, m3.z.x, m3.z.z; break;
		//	case 3:
		//		pinv << m3.x.x, m3.x.y, m3.x.z, m3.y.x, m3.y.y, m3.y.z, m3.z.x, m3.z.y, m3.z.z; break;
		//	}
		//	if (abs(pinv.determinant()) < eps_mat) return grad_hat(phi, p); //vec(0.);
		//	matEi inv = pinv.inverse();
		//	switch (dim) {
		//	case 2:
		//		m3.x.x = inv(0, 0);	m3.x.y = 0.; m3.x.z = inv(0, 1);
		//		m3.y.x = 0.;		m3.y.y = 0.; m3.y.z = 0.;
		//		m3.z.x = inv(1, 0); m3.z.y = 0.; m3.z.z = inv(1, 1);
		//		break;
		//	case 3:
		//		m3.x.x = inv(0, 0); m3.x.y = inv(0, 1); m3.x.z = inv(0, 2);
		//		m3.y.x = inv(1, 0); m3.y.y = inv(1, 1); m3.y.z = inv(1, 2);
		//		m3.z.x = inv(2, 0); m3.z.y = inv(2, 1); m3.z.z = inv(2, 2);
		//		break;
		//	}
		//	return m3 * v3;
		//}

		//const vec grad_nb(const std::vector<real>& phi, const unsigned& p) const {
		//	vec ret = vec(0.);
		//	const iVec3 c = cell->iCoord(pos[p]);
		//	for (int k = -1; k <= 1; k++) {
		//		for (int j = -1; j <= 1; j++) {
		//			for (int i = -1; i <= 1; i++) {
		//				const iVec3 ne = c + iVec3(i, j, k);
		//				const unsigned key = cell->hash(ne);
		//				for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
		//					const unsigned q = cell->linkList[key][m];
		//					if (q == p || type[q] != FLUID) continue;
		//					const vec	dr = pos[q] - pos[p];
		//					const real	dr2 = dr.mag2();
		//					ret += (phi[q] - phi[p]) * w1(sqrt(dr2)) / dr2 * dr;
		//				}
		//			}
		//		}
		//	}
		//	ret = dim / n0 * ret;
		//	return ret;
		//}

		//const mat grad_nb(const std::vector<vec>& u, const unsigned& p) const {
		//	mat ret = mat(vec(0., 0., 0.), vec(0., 0., 0.), vec(0., 0., 0.));
		//	const iVec3 c = cell->iCoord(pos[p]);
		//	for (int k = -1; k <= 1; k++) {
		//		for (int j = -1; j <= 1; j++) {
		//			for (int i = -1; i <= 1; i++) {
		//				const iVec3 ne = c + iVec3(i, j, k);
		//				const unsigned key = cell->hash(ne);
		//				for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
		//					const unsigned q = cell->linkList[key][m];
		//					if (q == p || type[q] != FLUID) continue;
		//					const vec	dr = pos[q] - pos[p];
		//					const real	dr2 = dr.mag2();
		//					const vec	tmp = w1(sqrt(dr2)) / dr2 * dr;
		//					ret.x += (u[q].x - u[p].x) * tmp;
		//					ret.y += (u[q].y - u[p].y) * tmp;
		//					ret.z += (u[q].z - u[p].z) * tmp;
		//				}
		//			}
		//		}
		//	}
		//	ret = dim / n0 * ret;
		//	return ret.trans();
		//}

		//const real div_suzuki(const std::vector<vec>& u, const unsigned& p) const {
		//	mat m3 = grad_suzuki(u, p);
		//	return m3.x.x + m3.y.y + m3.z.z;
		//}

		//const real lap_nb(const std::vector<real>& phi, const unsigned& p) const {
		//	real ret = 0.;
		//	const iVec3 c = cell->iCoord(pos[p]);
		//	for (int k = -1; k <= 1; k++) {
		//		for (int j = -1; j <= 1; j++) {
		//			for (int i = -1; i <= 1; i++) {
		//				const iVec3 ne = c + iVec3(i, j, k);
		//				const unsigned key = cell->hash(ne);
		//				for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
		//					const unsigned q = cell->linkList[key][m];
		//					if (q == p || type[q] != FLUID) continue;
		//					const vec	dr = pos[q] - pos[p];
		//					ret += (phi[q] - phi[p]) * w1(dr.mag());
		//				}
		//			}
		//		}
		//	}
		//	ret = 2. * dim / (n0 * lambda) * ret;
		//	return ret;
		//}

		//const vec lap_nb(const std::vector<vec>& u, const unsigned& p) const {
		//	vec ret = vec(0.);
		//	const iVec3 c = cell->iCoord(pos[p]);
		//	for (int k = -1; k <= 1; k++) {
		//		for (int j = -1; j <= 1; j++) {
		//			for (int i = -1; i <= 1; i++) {
		//				const iVec3 ne = c + iVec3(i, j, k);
		//				const unsigned key = cell->hash(ne);
		//				for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
		//					const unsigned q = cell->linkList[key][m];
		//					if (q == p || type[q] != FLUID) continue;
		//					const vec	dr = pos[q] - pos[p];
		//					ret += (u[q] - u[p]) * w1(dr.mag());
		//				}
		//			}
		//		}
		//	}
		//	ret = 2. * dim / (n0 * lambda) * ret;
		//	return ret;
		//}

	public:
		std::vector<real> varrho;
		vecp poly_0;
		vecp poly_px_0;
		vecp poly_pz_0;
		vecp poly_pxx_0;
		vecp poly_pzz_0;
		vecp poly_pxz_0;
	};

}