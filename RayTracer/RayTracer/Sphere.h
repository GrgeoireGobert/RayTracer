#pragma once
#include "Vector.h"
#include "Ray.h"

// Structure renvoyee lors d'une intersection
struct  intersection_details
{
	bool intersection = false; // Il y a une intersection ?
	double t = 0.0; // pour savoir le plus proche
	Vector inter_Pos; // Position monde de l'intersection
	Vector inter_Norm; // Normale au point d'intersection
	Vector inter_Color; // Couleur de l'objet
	double albedo = 1.0; // Albedo
	bool miroir = false; // Miroir
	bool transparent = false; // Transparent
	double indice_sphere = 1.0; // SI transparent, indice du milieu
	double intensite = 0.0; // Si source
	double rayon = 0.0; // SI transparent, indice du milieu
};

//
// Classe qui d�crit une sphere
//

class Sphere
{
public:

	Vector centre;
	double rayon;
	Vector couleur;
	double albedo;
	bool miroir;
	bool transparent;
	double indice_sphere;
	double intensite;

	Sphere(const Vector& Centre, const double radius, const Vector& Color, const double albed,bool miror, bool transp,double ind_sphere, double intensity);
	intersection_details intersect(Ray& rayon);
};
