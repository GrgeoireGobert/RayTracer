#pragma once
#include "Sphere.h"
#include "Vector.h"
#include "Ray.h"
#include <vector>

//
// Structure pour une lumiere
//
struct light_source
{
	Vector light_Pos; // Position
	double intensite=0.0; // Intensite
};

class Scene
{
public:

	std::vector<Sphere*> liste_objets;
	std::vector<light_source*> liste_sources;

	Scene();
	void add_sphere(Sphere* sphere);
	void add_light(light_source* light);
	intersection_details intersection(Ray& ray);
	Vector getColorMiroir(Ray& ray, intersection_details infos_intersection, int n);
	Vector getColorTransparent(Ray& ray, intersection_details infos_intersection, double indice_milieu_incident, int n);
	Vector getColorLambert(intersection_details infos_intersection);
	Vector getDirect(Ray& ray, intersection_details infos_intersection, double indice_milieu_incident, int n_miroir, int n_transp);
	Vector getIndirect(intersection_details infos_intersection, int n_rebonds);
};

