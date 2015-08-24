/*
*/
#pragma once
#include "header.h"

namespace SIM {

	template <typename real>
	class Parameter {
		typedef Vec3<real> vec;
	public:
		Parameter() { init(); }
		~Parameter() {}

	public:
		real k;		//radius k*dp
		real rho;	//density
		real niu;	//kinetic viscosity
		real dtMax;	//max time step
		real cfl;	//cfl num
		real tt;	//total time
		real eps;	//eps
		real beta;	//beta: free surface
		real alpha;	//particle shifting
		real c;		//particle shifting
		
		vec g;		//gravity
		real dt;	//current time step
		real umax;	//current umax
	private:
		void init() {
			g = vec(0., 0., 0.); dt = 1.;
		}
	};

}