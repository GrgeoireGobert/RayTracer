#include "Ray.h"

//
// Les m�thodes
//

Ray::Ray(const Vector& Origin, const Vector& Dir)
{
	origine = Origin;
	direction = Dir;
	// On normalise le vecteur directeur
	direction.normalize();
}