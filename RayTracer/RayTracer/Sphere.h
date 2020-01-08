#pragma once
#include "Vector.h"
#include "Ray.h"

class Sphere
{
public:

	Vector centre;
	double rayon;

	Sphere(const Vector& Centre, const double radius);
	bool intersect(Ray& rayon);
};

