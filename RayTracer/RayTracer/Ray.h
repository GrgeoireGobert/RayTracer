#pragma once
#include "Vector.h"

class Ray
{
public:

	Vector origine;
	Vector direction;

	Ray(const Vector& Origin, const Vector& Dir);
};

