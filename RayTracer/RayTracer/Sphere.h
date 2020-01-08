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
};

//
// Classe qui décrit une sphere
//

class Sphere
{
public:

	Vector centre;
	double rayon;
	Vector couleur;
	double albedo;

	Sphere(const Vector& Centre, const double radius, const Vector& Color, const double albed);
	intersection_details intersect(Ray& rayon);
};
