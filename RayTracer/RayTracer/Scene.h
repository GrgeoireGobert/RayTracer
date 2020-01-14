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
	double intensite; // Intensite
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
	Vector getColorLambert(intersection_details infos_intersection);
};

