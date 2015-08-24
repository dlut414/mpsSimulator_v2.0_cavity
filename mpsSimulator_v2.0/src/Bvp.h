/*
*/
#pragma once
#include "header.h"

enum bvpType {
	poly_a,
	poly_b,
	sinu_a,
	sinu_b,
	expo_a,
	expo_b,
};

template <class real>
class Bvp {
public:
	typedef Vec3<real> vec;
	Bvp(const bvpType& _type) : type(_type) {
	}
	~Bvp() {}

	const real func(const vec& v) const {
		switch (type) {
		case poly_a:
			return bvpFunc_poly_a(v);
		case poly_b:
			return bvpFunc_poly_b(v);
		case sinu_a:
			return bvpFunc_sinu_a(v);
		case sinu_b:
			return bvpFunc_sinu_b(v);
		case expo_a:
			return bvpFunc_expo_a(v);
		case expo_b:
			return bvpFunc_expo_b(v);
		default:
			std::cout << " !!! no type found !!! " << std::endl;
			return 0.;
		}
	}
	const real lap(const vec& v) const {
		switch (type) {
		case poly_a:
			return bvpLap_poly_a(v);
		case poly_b:
			return bvpLap_poly_b(v);
		case sinu_a:
			return bvpLap_sinu_a(v);
		case sinu_b:
			return bvpLap_sinu_b(v);
		case expo_a:
			return bvpLap_expo_a(v);
		case expo_b:
			return bvpLap_expo_b(v);
		default:
			std::cout << " !!! no type found !!! " << std::endl;
			return 0.;
		}
	}
	const vec grad(const vec& v) const {
		switch (type) {
		case poly_a:
			return bvpGrad_poly_a(v);
		case poly_b:
			return bvpGrad_poly_b(v);
		case sinu_a:
			return bvpGrad_sinu_a(v);
		case sinu_b:
			return bvpGrad_sinu_b(v);
		case expo_a:
			return bvpGrad_expo_a(v);
		case expo_b:
			return bvpGrad_expo_b(v);
		default:
			std::cout << " !!! no type found !!! " << std::endl;
			return 0.;
		}
	}

public:
	bvpType type;

private:
	/*poly_a*/
	const real bvpFunc_poly_a(const vec& v) const {
		return pow(v.x, 2) + pow(v.z, 2);
	}
	const real bvpLap_poly_a(const vec& v) const {
		return 4.;
	}
	const vec bvpGrad_poly_a(const vec& v) const {
		return vec(2.*v.x, 0., 2.*v.z);
	}

	/*poly_b*/
	const real bvpFunc_poly_b(const vec& v) const {
		return pow(v.x, 2) + pow(v.z, 2) + v.x*v.z;
	}
	const real bvpLap_poly_b(const vec& v) const {
		return 4.;
	}
	const vec bvpGrad_poly_b(const vec& v) const {
		return vec(2.*v.x + v.z, 0., 2.*v.z + v.x);
	}

	/*sinu_a*/
	const real bvpFunc_sinu_a(const vec& v) const {
		return (cos(M_PI*v.x) + cos(M_PI*v.z));
	}
	const real bvpLap_sinu_a(const vec& v) const {
		return (-M_PI*M_PI *(cos(M_PI*v.x) + cos(M_PI*v.z)));
	}
	const vec bvpGrad_sinu_a(const vec& v) const {
		return vec(-M_PI*sin(M_PI* v.x), 0., -M_PI*sin(M_PI* v.z));
	}

	/*sinu_b*/
	const real bvpFunc_sinu_b(const vec& v) const {
		return (cos(M_PI*v.x) + cos(M_PI*v.z) + cos(M_PI*v.x*v.z));
	}
	const real bvpLap_sinu_b(const vec& v) const {
		return (-M_PI*M_PI *(cos(M_PI*v.x) + cos(M_PI*v.z) + cos(M_PI*v.x*v.z)*(v.x*v.x + v.z*v.z)));
	}
	const vec bvpGrad_sinu_b(const vec& v) const {
		return vec(-M_PI*sin(M_PI* v.x) - M_PI*v.z*sin(M_PI*v.x*v.z), 0., -M_PI*sin(M_PI* v.z) - M_PI*v.x*sin(M_PI*v.x*v.z));
	}

	/*expo_a*/
	const real bvpFunc_expo_a(const vec& v) const {
		return (exp(M_PI*v.x) + exp(M_PI*v.z));
	}
	const real bvpLap_expo_a(const vec& v) const {
		return (M_PI*M_PI *(exp(M_PI*v.x) + exp(M_PI*v.z)));
	}
	const vec bvpGrad_expo_a(const vec& v) const {
		return vec(M_PI*exp(M_PI* v.x), 0., M_PI*exp(M_PI* v.z));
	}

	/*expo_b*/
	const real bvpFunc_expo_b(const vec& v) const {
		return (exp(M_PI*v.x) + exp(M_PI*v.z) + exp(M_PI*v.x*v.z));
	}
	const real bvpLap_expo_b(const vec& v) const {
		return (M_PI*M_PI *(exp(M_PI*v.x) + exp(M_PI*v.z) + (v.z*v.z + v.x*v.x)*exp(M_PI*v.x*v.z)));
	}
	const vec bvpGrad_expo_b(const vec& v) const {
		return vec(M_PI*(exp(M_PI* v.x) + v.z*exp(M_PI*v.x*v.z)), 0., M_PI*(exp(M_PI* v.z) + v.x*exp(M_PI*v.x*v.z)));
	}

};

