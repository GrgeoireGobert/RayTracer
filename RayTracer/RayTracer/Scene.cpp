#include "Scene.h"
#include <iostream>

// Constructeur par défaut
Scene::Scene()
{

}

/////////////////////////////
// Ajout sphere
/////////////////////////////
void Scene::add_sphere(Sphere* sphere)
{
	liste_objets.push_back(sphere);
}

/////////////////////////////
// Ajout lumiere
/////////////////////////////
void Scene::add_light(light_source* light)
{
	liste_sources.push_back(light);
}

/////////////////////////////
// Renvoie les infos sur l'intersection la plus proche
/////////////////////////////
intersection_details Scene::intersection(Ray& ray)
{
	// A priori, pas d'intersection
	intersection_details infos;
	// Parcours des objets de la scene
	for (int i=0; i<liste_objets.size();i++)
	{
		Sphere* sphere = liste_objets.at(i);
		intersection_details infos_inter_sphere = sphere->intersect(ray);
		// Si il y a intersection avec la sphere
		if (infos_inter_sphere.intersection==true)
		{
			// Si il n'y avait pas encore eu d'intersection trouvee
			if (infos.intersection == false)
			{
				infos = infos_inter_sphere;
			}
			// Sinon, si l'intersection est plus proche
			else if (infos_inter_sphere.t < infos.t)
			{
				infos = infos_inter_sphere;
			}
		}
	}

	return infos;
}

/////////////////////////////
// Determine la couleur dans le cas d'un materiau diffus
// infos_intersection est la structure décrivant l'intersection
/////////////////////////////
Vector Scene::getColorLambert(intersection_details infos_intersection)
{
	Vector light_Pos = this->liste_sources.at(0)->light_Pos;
	double light_power = this->liste_sources.at(0)->intensite;
	// On verifie si le point est dans l'ombre ou pas
	Vector to_Light = light_Pos - infos_intersection.inter_Pos;
	to_Light.normalize();
	double epsilon = 0.01;
	//// Rayon pour la visibilite
	Ray ray_ps(infos_intersection.inter_Pos + epsilon * infos_intersection.inter_Norm, to_Light);
	intersection_details infos_eclairage = this->intersection(ray_ps);
	
	// On verifie, si on a une intersection, si celle-ci se situe avant ou après la lumière !
	double dist_to_inter = (infos_eclairage.inter_Pos - infos_intersection.inter_Pos).norme2();
	double dist_to_light = to_Light.norme2();
	// Si le point est éclairé
	if (infos_eclairage.intersection == false || (infos_eclairage.intersection == true && dist_to_light<dist_to_inter))
	{
		//Calcul de la couleur pour le modèle Lambertien (diffus)
		double cos_theta = infos_intersection.inter_Norm.dot(to_Light);
		if (cos_theta < 0.0) { cos_theta = 0.0; }
		double facteur = light_power * infos_intersection.albedo * cos_theta / 3.1415 / ((light_Pos - infos_intersection.inter_Pos).norme2());
		Vector RGB = infos_intersection.inter_Color;
		RGB = facteur * RGB;

		return RGB;
	}
	// Sinon, le point est dans l'ombre
	else
	{
		return Vector(0.0, 0.0, 0.0);
	}

}

/////////////////////////////
// Determine la couleur dans le cas d'un materiau reflechissant
// ray est le rayon incident
// infos_intersection est la structure décrivant l'intersection
// n est la variable de fin de récursion
/////////////////////////////
Vector Scene::getColorMiroir(Ray& ray, intersection_details infos_intersection, int n)
{
	// Fin de la recursion
	if (n <= 0)
	{
		return Vector(0.0, 0.0, 0.0);
	}
	else
	{
		// Calcul du rayon réfléchi
		double epsilon = 0.001;
		Vector dir_reflexion = ray.direction - 2 * ray.direction.dot(infos_intersection.inter_Norm) * infos_intersection.inter_Norm;
		Ray new_ray(infos_intersection.inter_Pos+epsilon* infos_intersection.inter_Norm, dir_reflexion);
		// Lancer du rayon reflechi
		intersection_details infos = this->intersection(new_ray);
		// Renvoi de la couleur
		return this->getColor(new_ray, infos, 1.0, n - 1, 5);
	}
}

/////////////////////////////
// Determine la couleur dans le cas d'un materiau transparent
// ray est le rayon incident
// infos_intersection est la structure décrivant l'intersection
// indice_milieu_incident donne l'indice du milieu incident
// n est la variable de fin de récursion
/////////////////////////////
Vector Scene::getColorTransparent(Ray& ray, intersection_details infos_intersection, double indice_milieu_incident, int n)
{
	// Fin de la recursion
	if (n <= 0)
	{
		return Vector(0.0, 0.0, 0.0);
	}
	else
	{
		double indice_milieu_objet = infos_intersection.indice_sphere;
		double disc = 1 - pow(indice_milieu_incident / indice_milieu_objet, 2) * (1 - pow(ray.direction.dot(infos_intersection.inter_Norm), 2));
		// Cas de la reflexion totale
		if (disc < 0)
		{
			return this->getColorMiroir(ray, infos_intersection, n);
		}
		// Cas de la refraction
		else
		{
			//Si on rentre dans l'objet (air->objet)
			if (ray.direction.dot(infos_intersection.inter_Norm) < 0)
			{
				double epsilon = 0.001;
				Vector refr_dir = (indice_milieu_incident / indice_milieu_objet) * ray.direction - ((indice_milieu_incident / indice_milieu_objet) * ray.direction.dot(infos_intersection.inter_Norm) + sqrt(disc)) * infos_intersection.inter_Norm;
				Ray refr_ray(infos_intersection.inter_Pos - epsilon * infos_intersection.inter_Norm, refr_dir);
				intersection_details refr_ray_infos = this->intersection(refr_ray);

				return this->getColorTransparent(refr_ray, refr_ray_infos, infos_intersection.indice_sphere, n - 1);
			}
			//Si on quitte l'objet (objet->air)
			{
				//Inverser la normale pour les calculs
				infos_intersection.inter_Norm = (-1.0)*infos_intersection.inter_Norm;
				// Calcul
				double epsilon = 0.001;
				Vector refr_dir = (indice_milieu_objet/indice_milieu_incident) * ray.direction - ((indice_milieu_objet / indice_milieu_incident) * ray.direction.dot(infos_intersection.inter_Norm) + sqrt(disc)) * infos_intersection.inter_Norm;
				Ray refr_ray(infos_intersection.inter_Pos - epsilon * infos_intersection.inter_Norm, refr_dir);
				intersection_details refr_ray_infos = this->intersection(refr_ray);
				// Renvoi de la couleur
				return this->getColor(refr_ray, refr_ray_infos, 1.0, 5, 5);
			}
			
		}
	}
}

/////////////////////////////
// Fait la distinction de cas entre les effets physiques
// ray est le rayon incident
// infos_intersection est la structure décrivant l'intersection
// indice_milieu_incident donne l'indice du milieu incident
// n_miroir est la variable de fin de récursion pour l'effet miroir
// n_miroir est la variable de fin de récursion pour l'effet transparent
/////////////////////////////
Vector Scene::getColor(Ray& ray, intersection_details infos_intersection, double indice_milieu_incident, int n_miroir, int n_transp)
{
	if (infos_intersection.intersection == false)
	{
		return Vector(0.0, 0.0, 0.0);
	}
	//  Intersection
	else
	{
		if (infos_intersection.miroir == true)
		{
			return this->getColorMiroir(ray, infos_intersection, n_miroir);
		}
		else if (infos_intersection.transparent == true)
		{
			return this->getColorTransparent(ray, infos_intersection, indice_milieu_incident, n_transp);
		}
		else
		{
			return this->getColorLambert(infos_intersection);
		}
	}
}