#include "Scene.h"


// Constructeur par défaut
Scene::Scene()
{

}
// Ajout sphere
void Scene::add_sphere(Sphere* sphere)
{
	liste_objets.push_back(sphere);
}
// Ajout lumiere
void Scene::add_light(light_source* light)
{
	liste_sources.push_back(light);
}
// Renvoie les infos sur l'intersection la plus proche
intersection_details Scene::intersection(Ray& ray)
{
	intersection_details infos;

	for (int i=0; i<liste_objets.size();i++)
	{
		Sphere* sphere = liste_objets.at(i);
		intersection_details infos_inter_sphere = sphere->intersect(ray);
		if (infos_inter_sphere.intersection==true)
		{
			if (infos.intersection == false)
			{
				infos = infos_inter_sphere;
			}
			else if (infos_inter_sphere.t < infos.t)
			{
				infos = infos_inter_sphere;
			}
		}
	}

	return infos;
}
