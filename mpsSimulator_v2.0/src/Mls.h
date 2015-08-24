/*
*/
#pragma once
#include "Simulator.h"

namespace SIM {

	template <typename real, enum Dim dim>
	class Mls : public Simulator< real, dim, Mls<real, dim> > {
		typedef Vec3<real> vec;
		typedef Eigen::Matrix<real, 5, 1> vec5;
		typedef Eigen::Matrix<real, 5, 5> mat55;
		typedef Eigen::Triplet<real> tpl;
	public:
		Mls() {}
		~Mls() {}

		void step() {
			//convect1();
			part->updateTeam();
			calPnd();
			makeFs();
			makeMat();
			makeSource();
			solvMat();
			//convect2();
			makeFs();
			//collision();
			//shift();
			check();
		}

	private:
		void convect1() {
			/*standard*/
			/*
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] = part->vel1[p] + para.dt * (para.g + para.niu * part->lap(part->vel1, p));
			}
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] += para.dt * part->vel2[p];
			}
			*/
			//-------------------------------------------------------------------------//
			/*Euler*/
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] = part->vel1[p] + para.dt * (para.g + para.niu * part->lap(part->vel1, p));
			}
		}

		void convect2() {
			/*standard*/
			/*
			std::vector<vec> dash(part->np);
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				dash[p] = -para.dt / para.rho * part->grad_mls_poly2d_a(part->pres, p);
			}
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				part->vel1[p] = part->vel2[p] = part->vel2[p] + dash[p];
				part->pos[p] += para.dt * dash[p];
			}
			*/
			//------------------------------------------------------------------------//
			/*Euler*/
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] += -(para.dt / para.rho) * part->grad_suzuki(part->pres, p);
			}
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] += 0.5* para.dt * (part->vel1[p] + part->vel2[p]);
				part->vel1[p] = part->vel2[p];
			}
			//------------------------------------------------------------------------//
			/*fourth-order Ronge Kuta*/
			/*
			std::vector<vec> x0(part->np);
			std::vector<vec> k2(part->np);
			std::vector<vec> k3(part->np);
			std::vector<vec> k4(part->np);
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			x0[p] = part->pos[p];
			part->pos[p] += 0.5* para.dt* part->vel1[p];
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			k2[p] = part->vel2[p] + 0.5* para.dt* (-1./para.rho)* part->grad_suzuki(part->pres, p);
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			part->pos[p] = x0[p] + 0.5* para.dt* k2[p];
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			k3[p] = part->vel2[p] + 0.5* para.dt* (-1. / para.rho)* part->grad_suzuki(part->pres, p);
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			part->pos[p] = x0[p] + para.dt* k3[p];
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			k4[p] = part->vel2[p] + para.dt* (-1. / para.rho)* part->grad_suzuki(part->pres, p);
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			part->pos[p] = x0[p] + (para.dt / 6.)* (part->vel1[p] + 2.*k2[p] + 2.*k3[p] + k4[p]);
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			part->vel2[p] += -(para.dt / para.rho) * part->grad_suzuki(part->pres, p);
			part->vel1[p] = part->vel2[p];
			}
			*/
			//------------------------------------------------------------------------//
		}

		void makeMat() {
			std::vector<tpl> coef;
			coef.clear();
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					coef.push_back(tpl(p, p, 1.));
					continue;
				}
				if (part->isFs(p)) {
					coef.push_back(tpl(p, p, 1.));
					continue;
				}
				mat55 m5 = mat55::Zero();
				//std::vector<unsigned> used;
				const iVec3 c = part->cell->iCoord(part->pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const unsigned key = part->cell->hash(ne);
							//for (unsigned us = 0; us < used.size(); us++) { if (key == used[us]) std::cout << "used!!!!!!!!!!" << std::endl; }
							//used.push_back(key);
							for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
								const unsigned q = part->cell->linkList[key][m];
								if (part->type[q] == BD2) continue;
								if (part->team[p] != part->team[q]) continue;
								if (part->type[q] == BD1 && part->isFs(q)) continue;
								const vec	dr = part->pos[q] - part->pos[p];
								const real	dr1 = dr.mag();
								if (dr1 > part->r0) continue;
								const real  w = part->w3(dr1);
								const vec5	npq = part->poly2d_a(dr);
								m5.block<1, 5>(0, 0) += w * npq[0] * npq;
								m5.block<1, 5>(1, 0) += w * npq[1] * npq;
								m5.block<1, 5>(2, 0) += w * npq[2] * npq;
								m5.block<1, 5>(3, 0) += w * npq[3] * npq;
								m5.block<1, 5>(4, 0) += w * npq[4] * npq;
							}
						}
					}
				}
				if (abs(m5.determinant()) < part->eps_mat) {
					//std::cout << m5.determinant() << " " << p << std::endl;
					coef.push_back(tpl(p, p, 1.));
					continue;
				}
				const mat55 inv = m5.inverse();

				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const unsigned key = part->cell->hash(ne);
							for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
								const unsigned q = part->cell->linkList[key][m];
								if (part->type[q] == BD2) continue;
								if (part->team[p] != part->team[q]) continue;
								if (part->type[q] == BD1 && part->isFs(q)) continue;
								const vec	dr = part->pos[q] - part->pos[p];
								const real	dr1 = dr.mag();
								if (dr1 > part->r0) continue;
								const real  w = part->w3(dr1);
								const vec5	npq = part->poly2d_a(dr);
								const vec5  a = inv* (w* npq);
								real pq = a.dot( part->poly2d_a_lap(vec(0.)) );
								coef.push_back(tpl(p, q, pq));
							}
						}
					}
				}
			}
			mSol->a.setFromTriplets(coef.begin(), coef.end());
		}

		void makeSource() {
			for (unsigned p = 0; p < part->np; p++) {
				if ((part->type[p] == BD2) || part->isFs(p)) {
					//mSol->b[p] = 0.;
					mSol->b[p] = pow(part->pos[p].x, 2) + pow(part->pos[p].z, 2) + 5000.;
					continue;
				}
				mSol->b[p] = 4.;
				/*
				const real coef = para.rho / para.dt;
				const real corr = part->nbd[p] > 0 ? real(part->pn0) / real(part->pn[p]) : 1.;
				//mSol->b[p] = coef * part->div_suzuki(part->vel2, p);
				//mSol->b[p] = coef * (part->n0-part->pnd[p]) / (part->n0* para.dt);
				//mSol->b[p] = coef * ( 1.* part->div_suzuki(part->vel2, p) + 0.03* (part->n0 - part->pnd[p]*corr) / (part->n0* para.dt) );
				const real divTerm = part->div_mls_poly2d_a(part->vel1, p);
				const real pndTerm = (part->n0 - part->pnd[p] * corr) / (part->n0);
				mSol->b[p] = coef* (part->div_mls_poly2d_a(part->vel2, p) + abs(pndTerm)* divTerm + abs(divTerm)* pndTerm);
				//if (part->type[p]==FLUID&&part->pn_nb[p] < 0.1*part->pn0) mSol->b[p] = coef * (part->n0 - part->pnd[p]) / (part->n0* para.dt);
				*/
			}
		}
	};

}