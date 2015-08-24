/*
*/
#pragma once
#include "Simulator.h"
#include "Particle_mps.h"

namespace SIM {

	template <typename real, enum Dim dim>
	class Mps : public Simulator < real, dim, Mps<real, dim> > {
		typedef Vec3<real> vec;
		typedef Eigen::Triplet<real> tpl;
	public:
		Mps() {}
		~Mps() {}

		void step() {
			//shift();
			calCell(); //update cells
			//insertRand();
			convect1();
			calPnd();
			//makeFs();
			makeDirichlet();
			makeMat();
			makeSource();
			//bvpSource();
			solvMat();
			convect2();
			//makeFs();
			collision();
			calForVis();
			check();
			//profileOut();
			//sensorOut();
			////bvpAvgError();
			//bvpMaxError();
			//lapMaxError();
			//gradMaxError();
			//surfCol();
			//pthOrderVelSpatialFilter();
		}

		void convect1() {
			/*standard*/
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] = part->vel1[p] + para.dt * (para.g + para.niu * part->lap(part->vel1, p));
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] += para.dt * part->vel2[p];
			}
			//-------------------------------------------------------------------------//
			/*Euler*/
//#if OMP
//#pragma omp parallel for
//#endif
//			for (int p = 0; p < int(part->np); p++) {
//				if (part->type[p] != FLUID) continue;
//				part->vel2[p] = part->vel1[p] + para.dt * (para.g + para.niu * part->lap(part->vel1, p));
//			}
		}

		void convect2() {
			/*standard*/
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->dash[p] = -para.dt / para.rho * part->grad_hat(part->pres, p);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel1[p] = part->vel2[p] = part->vel2[p] + part->dash[p];
				part->pos[p] += para.dt * part->dash[p];
			}
			//------------------------------------------------------------------------//
			/*Euler*/
//#if OMP
//#pragma omp parallel for
//#endif
//			for (int p = 0; p < int(part->np); p++) {
//				if (part->type[p] != FLUID) continue;
//				part->vel2[p] += -(para.dt / para.rho) * part->grad(part->pres, p);
//			}
//#if OMP
//#pragma omp parallel for
//#endif
//			for (int p = 0; p < int(part->np); p++) {
//				if (part->type[p] != FLUID) continue;
//				part->pos[p] += 0.5* para.dt * (part->vel1[p] + part->vel2[p]);
//				part->vel1[p] = part->vel2[p];
//			}
		}

		void makeMat() {
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
				real pp = 0.;
				//std::vector<unsigned> used;
				const auto c = part->cell->iCoord(part->pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const auto ne = c + iVec3(i, j, k);
							const auto key = part->cell->hash(ne);
							//for (unsigned us = 0; us < used.size(); us++) { if (key == used[us]) std::cout << "used!!!!!!!!!!" << std::endl; }
							//used.push_back(key);
							for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
								const auto q = part->cell->linkList[key][m];
								if (part->type[q] == BD2) continue;
#if BD_OPT
								if (part->bdOpt(p, q)) continue;
#endif
								const auto dr = part->pos[q] - part->pos[p];
								const auto dr1 = dr.mag();
								if (dr1 > part->r0) continue;
								if (q == p) continue;
								const auto w = part->w1(dr1);
								pp -= w;
								if (part->isFs(q)) continue;
								coef.push_back(tpl(p, q, w));
							}
						}
					}
				}
				pp -= para.rho * part->lambda * (1./(para.rho*40.*40.)) * part->n0 / (2 * dim * para.dt * para.dt);
				if (pp == 0.) pp = 1.;
				coef.push_back(tpl(p, p, pp));
			}
			mSol->a.setFromTriplets(coef.begin(), coef.end());
		}

		void makeSource() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = 0.;
					continue;
				}
				if (part->isFs(p)) {
					mSol->b[p] = 0.;
					continue;
				}
				const auto coef = (part->lambda* part->n0* para.rho) / (2. * dim* para.dt);
				mSol->b[p] = coef * (1.* part->div(part->vel2, p) + 0.02* (part->n0 - part->pnd[p]) / (part->n0* para.dt));
				//mSol->b[p] = coef * part->div(part->vel2, p);
			}
		}

		void init_() {
			part = new Particle_mps<real, dim>();
			part->clean();
			*part << "Geo.in";
			part->init(para.k, para.beta);
			part->buildCell();
			part->b2b();
			part->b2norm();
			//part->updateTeam();
			part->init2d_x();
		}

	public:
		Particle_mps<real, dim>*  part;

	private:
		std::vector<tpl> coef;
	};

}