#pragma once
#include "Vector.h"

//
// Classe qui décrit un rayon lumineux
//

class Ray
{
public:

	Vector origine;
	Vector direction;

	Ray(const Vector& Origin, const Vector& Dir);
};

