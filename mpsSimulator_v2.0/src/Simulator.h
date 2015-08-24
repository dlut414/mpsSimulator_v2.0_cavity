/*
*/
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include "header.h"
#include "Parameter.h"
#include "Particle.h"
#include "MatSolver.h"
#include "Shifter.h"
#include "Bvp.h"
#include "Sensor.h"

namespace SIM {

	template <typename real, enum Dim dim, typename Derived>
	class Simulator {
		typedef Vec3<real> vec;
		typedef Mat3<real> mat;
		typedef Eigen::Matrix<real, dim, dim> matEi;
		typedef Eigen::Triplet<real> tpl;
	public:
		Simulator() { timeStep = 0; }
		~Simulator() {}

		Derived& derived() { return *static_cast<Derived*>(this); }
		const Derived& derived() const { return *static_cast<const Derived*>(this); }

		void operator >> (const std::string& str) const {
			saveData(str);
		}
		void operator << (const std::string& str) {
			std::ifstream file(str);
			if (!file.is_open()) std::cout << " no file Para. found ! " << std::endl;
			file >> para.k >> para.rho >> para.niu >> para.dtMax >> para.cfl >> para.tt >> para.eps >> para.beta >> para.alpha >> para.c;
			std::cout << para.k << " " << para.rho << " " << para.niu << std::endl;
			std::cout << para.dtMax << " " << para.cfl << " " << para.tt << " " << para.eps << " " << para.beta << std::endl;
			std::cout << para.alpha << " " << para.c << std::endl;
			std::cout << " reading Para. done " << std::endl;
			file.close();
		}

		void init() {
			derived().init_();
			mSol = new MatSolver<real, dim>(unsigned(derived().part->np), para.eps);
			bvp = new Bvp<real>(sinu_b);
			std::cout << " part num " << derived().part->np << std::endl;
			sen << "Sensor.in";
			real tmp = cfl();
			para.dt = tmp < para.dtMax ? tmp : para.dtMax;
			timeStep = int(derived().part->ct / para.dt);
		}

		void mainLoop() {
			auto* const part = derived().part;
			while (part->ct <= para.tt) {
				std::cout << " step ----------------------------------> " << timeStep << std::endl;
				real tmp = cfl();
				para.dt = tmp < para.dtMax ? tmp : para.dtMax;
				part->updateCell();
				derived().step();
				part->ct += para.dt;	timeStep++;
				std::cout << " time --------> " << part->ct << std::endl;
				std::cout << " dt ----------> " << para.dt << std::endl;
			}
		}

		real stepGL() {
			auto* const part = derived().part;
			if (part->ct > para.tt) {
				saveData();
				std::exit(0);
			}
			std::cout << " step ----------------------------------> " << timeStep << std::endl;
			real tmp = cfl();
			para.dt = tmp < para.dtMax ? tmp : para.dtMax;
			part->updateCell();
			derived().step();
			part->ct += para.dt;	timeStep++;
			std::cout << " time --------> " << part->ct << std::endl;
			std::cout << " dt ----------> " << para.dt << std::endl;
			return part->ct;
		}

		void sensorOut() {
			static int i = 0;
			std::ostringstream convert;
			convert << i++;
			sen.writeVect(derived().part);
			sen >> convert.str();
		}

		void profileOut() {
			auto* const part = derived().part;
			static std::string pf = "profile";
			sen.writeScal(derived().part);
			sen.profile(rt, pf);
		}

		void profileOut_avgVel2() {
			auto* const part = derived().part;
			static std::string pf = "profile";
			auto sum = 0.;
			auto count = 0;
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] == FLUID || part->type[p] == BD1) {
					sum += part->vel2[p].mag2();
					count++;
				}
			}
			sum = sum / count;
			std::ofstream file("./out/" + pf + ".out", std::ofstream::app);
			file << std::setprecision(6) << std::scientific << timeStep*para.cfl << " "
				<< std::setprecision(6) << std::scientific << sum
				<< std::endl;
			file.close();
			std::cout << " writing profile. done " << std::endl;
		}

		void saveData() const {
			static int i = 0;
			std::ostringstream convert;
			convert << i++;
			*(derived().part) >> ("./out/" + convert.str() + ".out");
		}
		void saveData(const std::string& str) const {
			*(derived().part) >> ("./out/" + str + ".out");
		}

		void fina() {}

	public:
		Parameter<real> para;
		MatSolver<real, dim>* mSol;
		Shifter<real, dim> shi;
		Sensor<real> sen;

	protected:
		void step() {}
		void convect() {}
		void visTerm_e() {}
		void visTerm_i() {}
		void presTerm_e() {}
		void presTerm_i() {}

		void makeDirichlet_p_op() {
			auto* const part = derived().part;
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] != BD1) continue;
				part->fs[p] = 1;
				break;
			}
		}

		void makeDirichlet_p_avg() {
			auto* const part = derived().part;
			Eigen::SparseMatrix<real> d(part->np, part->np);
			std::vector<tpl> coef;
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] != BD1) continue;
				for (int q = 0; q<int(part->np); q++) {
					if (part->type[q] == BD2) continue;
					coef.push_back(tpl(p, q, 1.));
				}
				break;
			}
			d.setFromTriplets(coef.begin(), coef.end());
			mSol->a = mSol->a + d;
		}

		void makeDirchlet_v() {}

		void makeNeumann_p() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] != BD1) continue;
				//part->neumann[p] = para.rho* (para.niu* part->lap(part->vel2, p) + para.g)* part->bdnorm.at(p);
				part->neumann[p] = 0.;
			}
		}

		void solvMat_p() {
			mSol->biCg();
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->pres[p] = mSol->x[p];
				if (part->pres[p] < -1.e5) part->pres[p] = -1.e5;
				if (part->pres[p] > 1.e5) part->pres[p] = 1.e5;
			}
		}

		void solvMat_phi() {
			auto* const part = derived().part;
			mSol->ccBiCg_augment(part->type);
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->phi[p] = mSol->x[p];
				if (part->phi[p] < -1.e5) part->phi[p] = -1.e5;
				if (part->phi[p] > 1.e5) part->phi[p] = 1.e5;
			}
		}

		void solvMat_v() {
			mSol->biCg_v();
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				switch (dim) {
				case TWOD:
					part->vel2[p].x = mSol->u[dim*p];
					part->vel2[p].z = mSol->u[dim*p + 1];
					break;
				case THREED:
					part->vel2[p].x = mSol->u[dim*p];
					part->vel2[p].y = mSol->u[dim*p + 1];
					part->vel2[p].z = mSol->u[dim*p + 2];
					break;
				default:
					break;
				}
			}
		}

		real cfl() {
			real umax = 0.;
			const auto* const part = derived().part;
			for (unsigned p = 0; p < part->np; p++) {
				real tmp = part->vel1[p].mag();
				if (tmp > umax) umax = tmp;
			}
			para.umax = umax;
			return para.cfl * part->dp / umax;
		}

		void calBdNoSlip() {
			derived().part->bdNoSlip();
		}
		void bdSetZero() {
			derived().part->bdSetZero();
		}

		void calCell() {
			derived().part->updateCell();
		}

		void calPnd() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->pnd[p] = part->cPnd(p);
				//part->pn[p] = part->cPn(p);
				//part->nbd[p] = part->cNbd(p);
			}
		}

		void calInvMat() {
			derived().part->updateInvMat();
		}

		void shift() {
			//shi.shiftPnd(part, para);
			//shi.shiftXu(derived().part, para);
			shi.shiftSpring(derived().part, para);
			//shi.shiftNearest2d(derived().part, para);
			//shi.shiftPndLs(part, para);
			//shi.shiftLs(part, para);
			//shi.shiftCo(part, para);
			//shi.shiftPbf(part, para);
		}

		void makeFs() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) part->fs[p] = part->_isFs(p);
			//derived().part->updateTeam();
		}

		void pthOrderPresSpatialFilter() {
			static auto counter = 0;
			if (counter++ % 1 != 0) return;
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] == BD2 || part->isFs(p)) continue;
				part->pres[p] = part->func(part->pres, p);
			}
		}
		void pthOrderVelSpatialFilter() {
			static auto counter = 0;
			if (counter++ % 30 != 0) return;
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] != FLUID || part->isFs(p)) continue;
				part->vel2[p] = part->func(part->vel1, p);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] != FLUID || part->isFs(p)) continue;
				part->vel1[p] = part->vel2[p];
			}
		}
		void surfCol() {
			auto* const part = derived().part;
			std::vector<vec> cor(part->np, vec(0.));
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID || !part->isFs(p)) continue;
				const auto c = part->cell->iCoord(part->pos[p]);
				for (auto k = -1; k <= 1; k++) {
					for (auto j = -1; j <= 1; j++) {
						for (auto i = -1; i <= 1; i++) {
							const auto ne = c + iVec3(i, j, k);
							const auto key = part->cell->hash(ne);
							for (auto m = 0; m < part->cell->linkList[key].size(); m++) {
								const auto q = part->cell->linkList[key][m];
								if (!part->isFs(q)) continue;
								if (q == p) continue;
								const auto dr = part->pos[q] - part->pos[p];
								const auto dv = part->vel1[q] - part->vel1[p];
								const auto dr1 = dr.mag();
								if (dr1 > part->r0) continue;
								const auto pro = dv*dr / dr1;
								if (abs(pro) > 0.8*para.umax) {
									cor[p] += 0.5* abs(pro)* dv.norm();
								}
							}
						}
					}
				}
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				part->vel1[p] += cor[p];
				part->vel2[p] = part->vel1[p];
			}
		}
		void collision() {
			auto* const part = derived().part;
			std::vector<vec> cor(part->np, vec(0.));
#if OMP
#pragma omp parallel for
#endif
			for(int p = 0; p<int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				const auto c = part->cell->iCoord(part->pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const auto ne = c + iVec3(i, j, k);
							const auto key = part->cell->hash(ne);
							for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
								const auto q = part->cell->linkList[key][m];
								if (p == q) continue;
								const auto dr = part->pos[q] - part->pos[p];
								const auto dr1 = dr.mag();
								if (dr1 < 0.5*part->dp) {
									const auto dv = part->vel1[q] - part->vel1[p];
									const auto tmp = dr*dv;
									if (tmp < 0.) {
										const auto coef = (0.5 / (dr1*dr1)*tmp*(1. + 0.2));
										cor[p] += coef* dr;
									}
								}
							}
						}
					}
				}
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel1[p] += cor[p];
			}
		}

		void damping() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID || !part->isFs(p)) continue;
				const auto n = part->vel1[p].norm();
				const auto gd = part->grad(part->vel1, p)* n;
				part->vel2[p] -= 0.01*part->dp* gd.mag()* n;
				part->vel1[p] = part->vel2[p];
			}
		}

		void calForVis() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				part->vort[p] = part->rot(part->vel1, p);
			}
		}

		void check() const {
			const auto* const part = derived().part;
			real dis = std::numeric_limits<real>::max();
			real vel = std::numeric_limits<real>::min();
			real phi = std::numeric_limits<real>::min();
			//real premi = std::numeric_limits<real>::max();
			//real prema = std::numeric_limits<real>::min();
			unsigned iv = 0, id = 0;
			for (unsigned p = 0; p < part->np; p++) {
				//const iVec3 c = part->cell->iCoord(part->pos[p]);
				//for (int k = -1; k <= 1; k++) {
				//	for (int j = -1; j <= 1; j++) {
				//		for (int i = -1; i <= 1; i++) {
				//			const iVec3 ne = c + iVec3(i, j, k);
				//			const unsigned key = part->cell->hash(ne);
				//			for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
				//				const unsigned q = part->cell->linkList[key][m];
				//				if (q == p) continue;
				//				const real dr1 = (part->pos[q] - part->pos[p]).mag();
				//				dis = dr1 < dis ? dr1 : dis;
				//			}
				//		}
				//	}
				//}
				const real v = part->vel1[p].mag();
				if (v > vel) {
					vel = v;
					iv = p;
				}
				real d = part->phi[p];
				if (abs(d) > phi) {
					phi = abs(d);
					id = p;
				}
				//if (part->pres[p] < premi) premi = part->pres[p];
				//if (part->pres[p] > prema) prema = part->pres[p];
				//if (part->isFs(p)) part->pres[p] = 0.;
			}
			//std::cout << " min dis: " << dis / part->dp << std::endl;
			std::cout << " max vel: " << vel << " --- id: " << iv << std::endl;
			std::cout << " max phi: " << phi << " --- id: " << id << std::endl;
			//std::cout << " max pres: " << prema << " --- min pres: " << premi << std::endl;
		}

		void bvpSource() {
			const auto* const part = derived().part;
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = 0.;
					continue;
				}
				if (part->isFs(p)) {
					mSol->b[p] = bvp->func(part->pos[p]);
					continue;
				}
				mSol->b[p] = bvp->lap(part->pos[p]);
			}
		}

		void bvpAvgError() {
			const auto* const part = derived().part;
			int n = 0;
			real err = 0.;
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID || part->isFs(p)) continue;
				real ext = bvp->func(part->pos[p]);
				err += abs(part->pres[p] - ext);
				n++;
			}
			std::cout << " bvp --- avg Error: " << err / n << std::endl;
			std::ofstream file("./out/out.txt", std::ofstream::app);
			file << part->dp << " " << err / n << std::endl;
			file.close();
		}

		void bvpMaxError() {
			const auto* const part = derived().part;
			real err = 0.;
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID || part->isFs(p)) continue;
				real ext = bvp->func(part->pos[p]);
				real tmp = abs(part->pres[p] - ext);
				if (tmp > err) err = tmp;
			}
			std::cout << " bvp --- max Error: " << err << std::endl;
			std::ofstream file("./out/out.txt", std::ofstream::app);
			file << part->dp << " " << err << std::endl;
			file.close();
		}

		void gradMaxError() {
			auto* const part = derived().part;
			for (unsigned p = 0; p < part->np; p++) {
				part->pres[p] = bvp->func(part->pos[p]);
			}
			real err = 0.;
			for (unsigned p = 0; p < part->np; p++) {
				vec ext = bvp->grad(part->pos[p]);
				real tmp = (part->grad(part->pres, p) - ext).mag();
				if (tmp > err) err = tmp;
			}
			std::cout << " |grad| --- max Error: " << err << std::endl;
			std::ofstream file("./out/out.txt", std::ofstream::app);
			file << part->dp << " " << err << std::endl;
			file.close();
		}

		void lapMaxError() {
			auto* const part = derived().part;
			for (unsigned p = 0; p < part->np; p++) {
				part->pres[p] = bvp->func(part->pos[p]);
			}
			real err = 0.;
			for (unsigned p = 0; p < part->np; p++) {
				real ext = bvp->lap(part->pos[p]);
				real tmp = abs(part->lap(part->pres, p) - ext);
				if (tmp > err) err = tmp;
			}
			std::cout << " lap --- max Error: " << err << std::endl;
			std::ofstream file("./out/out.txt", std::ofstream::app);
			file << part->dp << " " << err << std::endl;
			file.close();
		}

		void insertRand() {
			auto* const part = derived().part;
			real coef = 0.1;
			std::default_random_engine gen;
			std::normal_distribution<real> dis(0., 0.5);
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				vec dr = coef* part->dp* vec(dis(gen), 0., dis(gen));
				part->pos[p] += dr;
			}
		}

	protected:
		int timeStep;
		Timer tim;
		Bvp<real>* bvp;

	};

}