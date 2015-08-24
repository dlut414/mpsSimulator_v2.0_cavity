/*
*/
#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "Particle.h"

namespace SIM {

	template <class R, enum Dim D>
	class Mls_2d_d;

	template <typename real>
	class Sensor {
		typedef Vec3<real> vec;
	public:
		Sensor() {}
		~Sensor() {}

		void operator >> (const std::string& str) const {
			std::ofstream file("./out/s" + str + ".out", std::ofstream::out);
			for (unsigned s = 0; s < pos.size(); s++) {
				file << std::setprecision(6) << std::scientific << pos[s].x << " "
					<< std::setprecision(6) << std::scientific << pos[s].y << " "
					<< std::setprecision(6) << std::scientific << pos[s].z << " "
					<< std::setprecision(6) << std::scientific << vect[s].x << " "
					<< std::setprecision(6) << std::scientific << vect[s].y << " "
					<< vect[s].z
					<< std::endl;
			}
			file.close();
			std::cout << " writing Sensor. done " << std::endl;
		}
		void operator << (const std::string& str) {
			std::ifstream file(str);
			if (!file.is_open()) std::cout << " no file Sensor. found ! " << std::endl;
			pos.clear();	vect.clear();	scal.clear();
			while (file.good()) {
				vec p;
				file >> p.x >> p.y >> p.z;
				pos.push_back(p);
				vect.push_back(vec(0.));
				scal.push_back(0.);
			}
			file.close();
			std::cout << " reading Sensor. done " << std::endl;
		}
		void profile(const real& time, const std::string& str) const {
			std::ofstream file("./out/s" + str + ".out", std::ofstream::app);
			for (unsigned s = 0; s < pos.size(); s++) {
				file << std::setprecision(6) << std::scientific << time << " "
					<< std::setprecision(6) << std::scientific << scal[s]
					<< std::endl;
			}
			file.close();
			std::cout << " writing profile. done " << std::endl;
		}

		template <class R, enum Dim D, class Der>
		void writeVect(const Particle<R, D, Der>* const part) {
			const auto& pt = part->derived();
			for (unsigned s = 0; s < pos.size(); s++) {
				vect[s] = pt.func(pt.vel1, pos[s]);
			}
		}
		template <class R, enum Dim D, class Der>
		void writeScal(const Particle<R, D, Der>* const part) {
			const auto& pt = part->derived();
			for (unsigned s = 0; s < pos.size(); s++) {
				scal[s] = pt.func(pt.pres, pos[s]);
			}
		}

	public:
		std::vector<vec> pos;
		std::vector<vec> vect;
		std::vector<real> scal;
	};

}