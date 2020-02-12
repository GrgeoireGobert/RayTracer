#pragma once
#include "Object.h"
#include "Vector.h"
#include "Ray.h"

class Triangle : public Object
{
public:

	Vector A;
	Vector B;
	Vector C;
	Vector n;

	Triangle(const Vector& ptA,const Vector& ptB,const Vector& ptC,const Vector& color, const double albed);
	virtual intersection_details intersect(Ray& rayon);
};

