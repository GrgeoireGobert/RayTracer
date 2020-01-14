#include "Sphere.h"
#include <math.h>

//
// Les méthodes
//

// Constructeur
Sphere::Sphere(const Vector& Centre, const double radius, const Vector& Color, const double albed, bool miror)
{
	centre = Centre;
	rayon = radius;
	couleur = Color;
	albedo = albed;
	miroir = miror;
}

/////////////////////////////
// Structures contenant les infos sur une
// possible intersection
/////////////////////////////
intersection_details Sphere::intersect(Ray& ray)
{
	double a = ray.direction.norme();
	double b = 2 * ray.direction.dot(ray.origine - centre);
	double c = (ray.origine - this->centre).norme2() - rayon * rayon;
	// Discriminant
	double Delta = b * b - 4 * a * c;
	intersection_details infos;
	// Si delta négatif, pas d'intersection
	if (Delta < 0.0)
	{
		return infos;
	}
	// Sinon, on regarde si la plus grand des solutions est négative
	else if((-b+sqrt(Delta))/(2*a)<0.0)
	{
		return infos;
	}
	// Sinon, alors on a une intersection
	else
	{
		// Intersection la plus proche
		double t_proche = (-b + sqrt(Delta)) / (2 * a);
		if ((-b - sqrt(Delta)) / (2 * a) > 0.0)
		{
			t_proche = (-b - sqrt(Delta)) / (2 * a);
		}

		// Details sur l'intersection
		infos.intersection = true;
		infos.t = t_proche;
		infos.inter_Pos = ray.origine + t_proche * ray.direction; // Position
		infos.inter_Color = couleur; // Couleur
		infos.inter_Norm = infos.inter_Pos-centre; // Normale
		infos.inter_Norm.normalize();
		infos.albedo = albedo; // Albedo
		infos.miroir = miroir;

		return infos;
	}

}