#include "Sphere.h"
#include <math.h>

//
// Les m�thodes
//

// Constructeur
Sphere::Sphere(const Vector& Centre, const double radius)
{
	centre = Centre;
	rayon = radius;
}
// Booleen qui indique si il y a une intersection
bool Sphere::intersect(Ray& ray)
{
	double a = ray.direction.norme();
	double b = 2 * ray.direction.dot(ray.origine - centre);
	double c = (ray.origine - this->centre).norme2() - rayon * rayon;
	// Discriminant
	double Delta = b * b - 4 * a * c;
	// Si delta n�gatif, pas d'intersection
	if (Delta < 0.0)
	{
		return false;
	}
	// Sinon, on regarde si la plus grand des solutions est n�gative
	else if((-b+sqrt(Delta))/(2*a)<0.0)
	{
		return false;
	}
	else
	{
		return true;
	}

}