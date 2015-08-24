/*
*/
#pragma once
#include "iVec3.h"
#include "Vec3.h"
#include "Mat3.h"
#include "BBox.h"
#include "Timer.h"

#define DEBUG 1
#define NTHREAD 4
#define OMP	1
#define BD_OPT 0

enum pType {
	FLUID = 0,
	BD1 = 1,
	BD2 = 2,
};
enum Dim {
	TWOD = 2,
	THREED = 3,
};
enum dType {
	SCALAR,
	VECTOR,
};
enum baseType {
	A,
	B,
	C,
	D,
	E,
};