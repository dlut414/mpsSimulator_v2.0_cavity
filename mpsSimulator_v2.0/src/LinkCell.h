/*
*/
#pragma once
#include <iostream>
#include <vector>
#include <cassert>
#include "header.h"

namespace SIM {

	template <typename real>
	class LinkCell {
		typedef Vec3<real> vec;
	public:
		LinkCell(const BBox<real>& b, const real c) : box(b), cSize(c) { init(); }
		~LinkCell() { fina(); }

		inline const unsigned hash(const iVec3& c) const {
			//return unsigned((541 * c.x + 79 * c.y + 31 * c.z + 100000) % cNum);
			return unsigned((c.x + dv.x* c.y + sheet* c.z + 100000) % cNum);
		}
		inline const unsigned hash(const vec& p) const {
			return hash( iCoord(p) );
		}
		inline const iVec3 iCoord(const vec& p) const {
			return iVec3( int(p.x / cSize), int(p.y / cSize), int(p.z / cSize) );
		}

		void update(const std::vector<vec>& pos) {
			for (unsigned i = 0; i < linkList.size(); i++) {
				linkList[i].clear();
			}
			for (unsigned i = 0; i < pos.size(); i++) {
				linkList[hash(pos[i])].push_back(i);
			}
		}

	public:
		BBox<real> box;
		std::vector< std::vector<unsigned> > linkList;

	private:
		void init() {
			vec v = vec(box.pMax - box.pMin);
			dv.x = unsigned(v.x / cSize) + 1; dv.x = dv.x > 3 ? dv.x : 3;
			dv.y = unsigned(v.y / cSize) + 1; dv.y = dv.y > 3 ? dv.y : 3;
			dv.z = unsigned(v.z / cSize) + 1; dv.z = dv.z > 3 ? dv.z : 3;
			cNum = dv.x * dv.y * dv.z;
			sheet = dv.x * dv.y;

			if (cNum >= linkList.max_size()) std::cout << " linkList overflow ! " << std::endl;
			linkList.clear();
			linkList = std::vector< std::vector<unsigned> >(cNum);
			std::cout << " cell num: " << cNum << std::endl;
		}
		void fina() {}

	private:
		iVec3 dv;
		unsigned sheet;
		unsigned cNum;
		real cSize;
	};

}