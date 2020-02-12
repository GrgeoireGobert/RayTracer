#pragma once
#include "Vector.h"
#include "Ray.h"
#include "Object.h"



//
// Classe qui décrit une sphere
//

class Sphere : public Object
{
public:

	Vector centre;
	double rayon;

	Sphere(const Vector& Centre, const double radius, const Vector& Color, const double albed,bool miror, bool transp,double ind_sphere, double intensity);
	virtual intersection_details intersect(Ray& rayon);
};
