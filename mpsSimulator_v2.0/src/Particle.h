/*
*/
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "header.h"
#include "LinkCell.h"
#include "ConsValue.h"

namespace SIM {
	
	template <class real, enum Dim dim, class Derived>
	class Particle : public ConsValue<real, dim> {
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
		Particle() {}
		~Particle() {}

		Derived& derived() { return *static_cast<Derived*>(this); }
		const Derived& derived() const { return *static_cast<const Derived*>(this); }

		void clean() {
			type.clear(); pos.clear(); vel1.clear(); vel2.clear(); vel_m1.clear();
			pnd.clear(); pres.clear(); pn.clear(); nbd.clear(); fs.clear();
			team.clear(); phi.clear(); vort.clear(); dash.clear(); norm.clear();
		}
		void operator >> (std::string str) const {
			std::ofstream file(str, std::ofstream::out);
			file << ct << std::endl;
			file << dp << std::endl;
			file << np << " " << bd1 << " " << bd2 << std::endl;
			for (unsigned p = 0; p < np; p++) {
				file << type[p] << " " << pos[p].x << " " << pos[p].y << " " << pos[p].z << " "
					<< vel1[p].x << " " << vel1[p].y << " " << vel1[p].z << std::endl;
			}
			std::cout << " writing Geo. done " << std::endl;
			file.close();
		}
		void operator << (std::string str) {
			int n;	 int t;		vec p;		vec	v;
			std::ifstream file(str);
			if (!file.is_open()) std::cout << " no file Geo. found ! " << std::endl;
			file >> ct >> dp >> np >> bd1 >> bd2;
			n = np;
			while (n-- > 0) {
				file >> t >> p.x >> p.y >> p.z >> v.x >> v.y >> v.z;
				addPart(pType(t), p, v); 
			}
			file.close();
			std::cout << " reading Geo. done " << std::endl;
		}

		void addPart (const pType& t, const vec& p, const vec& v) {
			type.push_back(t);	pos.push_back(p);
			vel1.push_back(v);	vel2.push_back(v); vel_m1.push_back(v);
			pnd.push_back(0.);	pres.push_back(0.);
			pn.push_back(0);	nbd.push_back(0);
			fs.push_back(0);	team.push_back(0);	
			phi.push_back(0.); vort.push_back(vec(0.));
			dash.push_back(vec(0.));
			norm.push_back(vec(0.));
		}

		const real cPnd(const unsigned& p) const {
			real ret = 0.;
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							//if (q == p) continue;
							//if (type[q]==FLUID && (team[p]!=team[q])) continue;
							const real dr1 = (pos[q] - pos[p]).mag();
							ret += w1(dr1);
						}
					}
				}
			}
			return ret;
		}

		const real cPnd(const vec& p) const {
			real ret = 0.;
			const iVec3 c = cell->iCoord(p);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							const real dr1 = (pos[q] - p).mag();
							ret += w1(dr1);
						}
					}
				}
			}
			return ret;
		}

		const real cPn(const unsigned& p) const {
			real ret = 0;
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							if (q == p) continue;
							const real	dr1 = (pos[q] - pos[p]).mag();
							//if (dr1 < r0) ret += (isFs(q) ? 0.7 : 1.);
							if (dr1 < r0) ret += 1.;
						}
					}
				}
			}
			return ret;
		}

		const unsigned cNbd(const unsigned& p) const {
			unsigned ret = 0;
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							const real	dr1 = (pos[q] - pos[p]).mag();
							if (q == p || dr1 > r0) continue;
							if (isFs(q)) ret++;
						}
					}
				}
			}
			return ret;
		}

		void buildCell() {
			BBox<real> b = BBox<real>();
			for (unsigned i = 0; i < pos.size(); i++) {
				b += pos[i];
			}
			b.Expand(0.1);
			cell = new LinkCell<real>(b, r0);
			updateCell();
		}
		inline void updateCell() {
			cell->update(pos);
		}

		void updateTeam() {
			for (unsigned p = 0; p < np; p++) team[p] = 0;
			auto t = 1;
			for (unsigned p = 0; p < np; p++) {
				if (type[p] != FLUID) continue;
				dfs(p, t++);
			}
		}

		void dfs(const unsigned& p, const int& t) {
			if (type[p] == BD2) return;
			if (team[p] != 0) return;
			if (type[p] == BD1) { team[p] = t; return; }
			team[p] = t;
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							const real	dr1 = (pos[q] - pos[p]).mag();
							if (q == p || dr1 > 2.*dp) continue;
							dfs(q, t);
						}
					}
				}
			}
			return;
		}

		void b2b() {
			bbMap.clear();
			for (unsigned p = 0; p < pos.size(); p++) {
				if (type[p] != BD2) continue;
				real tmpdr = std::numeric_limits<real>::infinity();
				unsigned tmpbb = 0;
				const iVec3 c = cell->iCoord(pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const unsigned key = cell->hash(ne);
							for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
								const unsigned q = cell->linkList[key][m];
								if (q == p || type[q] != BD1) continue;
								const real dr1 = (pos[q] - pos[p]).mag();
								if (dr1 < tmpdr) {
									tmpdr = dr1;
									tmpbb = q;
								}
							}
						}
					}
				}
				bbMap[p] = tmpbb;
			}
		}

		void b2norm() {
			bdnorm.clear();
			for (const auto& p : bbMap) {
				const auto n = pos[p.second] - pos[p.first];
				bdnorm[p.second] = n;
			}
			for (const auto& p : bbMap) {
				const auto n = pos[p.second] - pos[p.first];
				const auto tmp = bdnorm.at(p.second);
				if(tmp.mag() < n.mag())bdnorm[p.second] = n;
			}
			for (const auto& p : bbMap) {
				bdnorm[p.second] = bdnorm.at(p.second).norm();
			}
		}

		void b2neumann() {
			neumann.clear();
			for (unsigned p = 0; p < np; p++) {
				if (type[p] != BD1) continue;
				neumann[p] = 0.;
			}
		}

		inline const int isFs(const unsigned& p) const {
			return fs[p];
			//return _isFs(p);
		}

		inline const int isFs(const vec& p) const {
			return (cPnd(p) < beta * n0);
		}

		const int _isFs(const unsigned& p) {
			//return (pnd[p] < beta * n0);

			/*tamai surface detection*/
			if (type[p] == BD2) return 0;
			if (pnd[p] > (beta * n0)) return 0;
			mat m3 = mat(vec(0., 0., 0.), vec(0., 0., 0.), vec(0., 0., 0.));
			//vec gc = vec(0., 0., 0.);
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							if (q == p) continue;
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w2(dr1);
							const vec	npq = dr / dr1;
							m3.x += w * npq.x * npq;
							m3.y += w * npq.y * npq;
							m3.z += w * npq.z * npq;
							//gc += w* npq;
						}
					}
				}
			}
			m3 = (dim / n0) * m3;
			//gc = gc.norm();
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
			real eigmin = eig[0];
			Eigen::Matrix<real, dim, 1> eigv = eigvs.col(0);

			if (eigmin <= 0.2) return 2;
			if (eigmin > 0.8) return 0;
			vec neigv1, neigv2;
			switch (dim) {
			case 2:
				neigv1.x = eigv[0]; neigv1.z = eigv[1]; break;
			case 3:
				neigv1.x = eigv[0]; neigv1.y = eigv[1]; neigv1.z = eigv[2]; break;
			}
			neigv2 = -neigv1;
			//if (neigv * gc > 0.) neigv = -1.*neigv;
			//gc = -1. * gc;
			const real root2 = 1.415;
			int flag1 = 1, flag2 = 1;
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							if (q == p) continue;
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 < root2 * dp) {
								if ((dr / dr1) * neigv1 >(root2 / 2.)) flag1 = 0;
								if ((dr / dr1) * neigv2 >(root2 / 2.)) flag2 = 0;
							}
							else {
								if ((pos[p] + dp * neigv1 - pos[q]).mag() < dp) flag1 = 0;
								if ((pos[p] + dp * neigv2 - pos[q]).mag() < dp) flag2 = 0;
							}
						}
					}
				}
			}
			if (flag1) {
				norm[p] = neigv1;
				return 1;
			}
			if (flag2) {
				norm[p] = neigv2;
				return 1;
			}
			return 0;
		}

		inline const int bdOpt(const unsigned& p, const unsigned& q) const {
			//if (type[q] == BD2) return 1;
			//if (team[p] != team[q]) return 1;
			if (type[q] == BD1 && isFs(q)) return 1;
			//if (type[q] == BD2 && isFs(bbMap.at(q))) return 1;
			//if (type[p] == FLUID && type[q] == BD2 && isFs(bbMap.at(q))) return 1;
			return 0;
		}
		inline const int bdOpt(const unsigned& q) const {
			//if (type[q] == BD2) return 1;
			//if (team[p] != team[q]) return 1;
			if (type[q] == BD1 && isFs(q)) return 1;
			//if (type[q] == BD2 && isFs(bbMap.at(q))) return 1;
			//if (type[p] == FLUID && type[q] == BD2 && isFs(bbMap.at(q))) return 1;
			return 0;
		}
		inline const int bdSlip(const unsigned& q) const {
			///slip
			if (type[q] != FLUID) return 1;
			return 0;
		}

		void bdNoSlip() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(np); p++) {
				if (type[p] != BD2) continue;
				const auto mid = bbMap[p];
				const auto mirror = 2.* pos[mid] - pos[p];
				const auto vmir = derived().func_nobd2(vel1, mirror);
				///axis symmetric
				//const auto norm = (pos[mid] - pos[p]).norm();
				//const auto vnor = (vmir* norm)* norm;
				//vel1[p] = vel2[p] = 2.* vnor - vmir;
				///point symmetric
				vel1[p] = vel2[p] = -vmir;
			}
		}
		void bdSetZero() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(np); p++) {
				if (type[p] != BD2) continue;
				vel1[p] = vel2[p] = 0.;
			}
		}

	public:
		real ct;
		unsigned np, bd1, bd2;
		std::vector<vec> pos;
		std::vector<vec> vel1;
		std::vector<vec> vel2;
		std::vector<vec> vel_m1;
		std::vector<vec> dash;
		std::vector<vec> norm;

		std::vector<real>	pnd;
		std::vector<real>	pres;
		std::vector<pType>	type;
		std::vector<real>	pn;
		std::vector<unsigned> nbd;
		std::vector<int> fs;
		std::vector<int> team;
		std::vector<real> phi;
		std::vector<vec> vort;
		std::map<unsigned, vec> bdnorm;
		std::map<unsigned, real> neumann;
		std::map<unsigned, unsigned> bbMap;

		LinkCell<real>*		cell;

	private:
	};

}