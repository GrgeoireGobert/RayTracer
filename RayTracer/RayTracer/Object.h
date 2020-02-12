#pragma once
#include "Vector.h"
#include "Ray.h"

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

class Object
{
public:
	Vector couleur;
	double albedo;
	bool miroir;
	bool transparent;
	double indice_sphere;
	double intensite=0.0;


	virtual intersection_details intersect(Ray& rayon)
	{
		intersection_details renvoi;
		return renvoi;
	}
};

