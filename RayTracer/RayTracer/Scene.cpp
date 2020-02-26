#include "Scene.h"
#include <iostream>
#include <random>
#include <cmath>


std::default_random_engine engine;
std::uniform_real_distribution<double> distrib(0, 1);
double pi = 3.14159265359;

// Constructeur par défaut
Scene::Scene()
{

}

/////////////////////////////
// Ajout sphere
/////////////////////////////
void Scene::add_objet(Object* object)
{
	liste_objets.push_back(object);
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
		Object* sphere = liste_objets.at(i);
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
	//On parcourt les spheres
	for (int i = 0; i < liste_objets.size(); i++)
	{
		//On verifie si c'est une source de lumière
		if (liste_objets.at(i)->intensite > 0.0)
		{
			Sphere* source = dynamic_cast<Sphere*>(liste_objets.at(i));
			//On prend un point au hasard sur cette source
			//Construction repère local
			Vector e_z =  infos_intersection.inter_Pos + (-1.0) * source->centre;
			e_z.normalize();
			Vector e_x(0.0, 0.0, 0.0);
			if (std::abs(e_z.x) <= std::abs(e_z.y) && std::abs(e_z.x) < std::abs(e_z.z))
			{
				e_x.x = 0.0;
				e_x.y = (-1.0) * e_z.z;
				e_x.z = (1.0) * e_z.y;
			}
			else if (std::abs(e_z.y) <= std::abs(e_z.x) && std::abs(e_z.y) < std::abs(e_z.z))
			{
				e_x.x = (-1.0) * e_z.z;
				e_x.y = 0.0;
				e_x.z = (1.0) * e_z.x;
			}
			else
			{
				e_x.x = (-1.0) * e_z.y;
				e_x.y = (1.0) * e_z.x;
				e_x.z = 0.0;
			}
			e_x.normalize();
			Vector e_y = e_z.prod_vect(e_x);

			//On génère une direction aléatoire
			double r1 = distrib(engine);
			double r2 = distrib(engine);

			double x = cos(2 * pi * r1) * sqrt(1 - r2);
			double y = sin(2 * pi * r1) * sqrt(1 - r2);
			double z = sqrt(r2);

			Vector dir = x * e_x + y * e_y + z * e_z;
			dir.normalize();

			//On a un point sur la surface de la source:
			Vector pt_on_source = source->centre + source->rayon * dir;
			
			//On calcule sa contribution
			double light_power = source->intensite;
			Vector to_Light_unnormalized = pt_on_source - infos_intersection.inter_Pos;
			Vector to_Light = pt_on_source - infos_intersection.inter_Pos;
			to_Light.normalize();
			double epsilon = 0.01;
			//// Rayon pour la visibilite
			Ray ray_ps(infos_intersection.inter_Pos + epsilon * infos_intersection.inter_Norm, to_Light);
			intersection_details infos_eclairage = this->intersection(ray_ps);
			// On verifie, si on a une intersection, si celle-ci se situe avant ou après la lumière !
			double dist_to_inter = (infos_eclairage.inter_Pos - infos_intersection.inter_Pos).norme2();
			double dist_to_light = to_Light_unnormalized.norme2();
			// Si le point est éclairé
			if (infos_eclairage.intersection == false || (infos_eclairage.intersection == true && infos_eclairage.intensite>0.0))
			{
				//Calcul de la couleur pour le modèle Lambertien (diffus)
				double cos_theta = infos_intersection.inter_Norm.dot(to_Light);
				if (cos_theta < 0.0) { cos_theta = 0.0; }

				Vector n_source = pt_on_source - source->centre;
				n_source.normalize();
				Vector to_inter = infos_intersection.inter_Pos - pt_on_source;
				to_inter.normalize();
				double cos_theta_prime = n_source.dot(to_inter);
				if (cos_theta_prime < 0.0) { cos_theta_prime = 0.0; }

				Vector OP = (infos_intersection.inter_Pos  + (-1.0) * source->centre);
				OP.normalize();
				double cos_alpha = n_source.dot(OP);

				double facteur = light_power * infos_intersection.albedo /pi * cos_theta * cos_theta_prime / (dist_to_light) / (cos_alpha/pi);
				// Le 4*\pi*R*R s'annule entre l'intensité et la densité de proba
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

	}
}

/////////////////////////////
// Determine la couleur dans le cas d'un materiau partiellement diffus/speculaire (Phong)
// infos_intersection est la structure décrivant l'intersection
/////////////////////////////
Vector Scene::getColorPhong(intersection_details infos_intersection, double rho_diff, int n)
{
	double rho_spec = 1.0 - rho_diff;

	//On parcourt les spheres
	for (int i = 0; i < liste_objets.size(); i++)
	{
		//On verifie si c'est une source de lumière
		if (liste_objets.at(i)->intensite > 0.0)
		{
			Sphere* source = dynamic_cast<Sphere*>(liste_objets.at(i));
			//On prend un point au hasard sur cette source
			//Construction repère local
			Vector e_z = infos_intersection.inter_Pos + (-1.0) * source->centre;
			e_z.normalize();
			Vector e_x(0.0, 0.0, 0.0);
			if (std::abs(e_z.x) <= std::abs(e_z.y) && std::abs(e_z.x) < std::abs(e_z.z))
			{
				e_x.x = 0.0;
				e_x.y = (-1.0) * e_z.z;
				e_x.z = (1.0) * e_z.y;
			}
			else if (std::abs(e_z.y) <= std::abs(e_z.x) && std::abs(e_z.y) < std::abs(e_z.z))
			{
				e_x.x = (-1.0) * e_z.z;
				e_x.y = 0.0;
				e_x.z = (1.0) * e_z.x;
			}
			else
			{
				e_x.x = (-1.0) * e_z.y;
				e_x.y = (1.0) * e_z.x;
				e_x.z = 0.0;
			}
			e_x.normalize();
			Vector e_y = e_z.prod_vect(e_x);

			//On génère une direction aléatoire
			double r1 = distrib(engine);
			double r2 = distrib(engine);

			double x = cos(2 * pi * r1) * sqrt(1 - r2);
			double y = sin(2 * pi * r1) * sqrt(1 - r2);
			double z = sqrt(r2);

			Vector dir = x * e_x + y * e_y + z * e_z;
			dir.normalize();

			//On a un point sur la surface de la source:
			Vector pt_on_source = source->centre + source->rayon * dir;

			//On calcule sa contribution
			double light_power = source->intensite;
			Vector to_Light_unnormalized = pt_on_source - infos_intersection.inter_Pos;
			Vector to_Light = pt_on_source - infos_intersection.inter_Pos;
			to_Light.normalize();
			double epsilon = 0.01;
			//// Rayon pour la visibilite
			Ray ray_ps(infos_intersection.inter_Pos + epsilon * infos_intersection.inter_Norm, to_Light);
			intersection_details infos_eclairage = this->intersection(ray_ps);
			// On verifie, si on a une intersection, si celle-ci se situe avant ou après la lumière !
			double dist_to_inter = (infos_eclairage.inter_Pos - infos_intersection.inter_Pos).norme2();
			double dist_to_light = to_Light_unnormalized.norme2();
			// Si le point est éclairé
			if (infos_eclairage.intersection == false || (infos_eclairage.intersection == true && infos_eclairage.intensite > 0.0))
			{
				//Calcul de la couleur pour le modèle Lambertien (diffus)
				double cos_theta = infos_intersection.inter_Norm.dot(to_Light);
				if (cos_theta < 0.0) { cos_theta = 0.0; }

				Vector n_source = pt_on_source - source->centre;
				n_source.normalize();
				Vector to_inter = infos_intersection.inter_Pos - pt_on_source;
				to_inter.normalize();
				double cos_theta_prime = n_source.dot(to_inter);
				if (cos_theta_prime < 0.0) { cos_theta_prime = 0.0; }

				Vector OP = (infos_intersection.inter_Pos + (-1.0) * source->centre);
				OP.normalize();
				double cos_alpha = n_source.dot(OP);

				double facteur=0.0;
				//On genere un nombre aleatoire :
				if (distrib(engine) <= rho_diff)
				{
					//Partie diffuse
					facteur = light_power * rho_diff / pi * cos_theta * cos_theta_prime / (dist_to_light) / (cos_alpha / pi);
					// Le 4*\pi*R*R s'annule entre l'intensité et la densité de proba

					//Normalisation
					facteur = facteur / rho_diff;
				}
				else
				{
					//Partie speculaire
					Vector perfect_dir = to_inter - 2 * to_inter.dot(infos_intersection.inter_Norm) * infos_intersection.inter_Norm;
					perfect_dir.normalize();
					double cos_theta_sec = (perfect_dir).dot((-1.0) * infos_intersection.ray_dir);
					if (cos_theta_sec < 0.0) { cos_theta_sec = 0.0; }


					facteur = light_power * rho_spec*(n+2) / (2*pi) * pow(cos_theta_sec, n) * cos_theta_prime / (dist_to_light) / (cos_alpha / pi);
					//Normalisation
					facteur = facteur / rho_spec;
				}


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
		return this->getDirect(new_ray, infos, 1.0, n - 1, 5);
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
				return this->getDirect(refr_ray, refr_ray_infos, 1.0, 5, 5);
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
Vector Scene::getDirect(Ray& ray, intersection_details infos_intersection, double indice_milieu_incident, int n_miroir, int n_transp)
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
		else if (infos_intersection.n_phong > 0)
		{
			return this->getColorPhong(infos_intersection, infos_intersection.rho_diff_phong, infos_intersection.n_phong);
		}
		else
		{
			return this->getColorLambert(infos_intersection);
		}
	}
}



/////////////////////////////
// Calcule l'eclairage indirect
// infos_intersection est la structure décrivant l'intersection
/////////////////////////////
Vector Scene::getIndirect(intersection_details infos_intersection,int n_rebonds)
{
	//On regarde l'energie incidente indirecte
	if (n_rebonds > 0)
	{
		//Construction repère local
		Vector e_z = infos_intersection.inter_Norm;
		Vector e_x(0.0,0.0, 0.0);
		if (std::abs(e_z.x) <= std::abs(e_z.y) && std::abs(e_z.x) < std::abs(e_z.z))
		{
			e_x.x = 0.0;
			e_x.y = (-1.0) * e_z.z;
			e_x.z = (1.0) * e_z.y;
		}
		else if (std::abs(e_z.y) <= std::abs(e_z.x) && std::abs(e_z.y) < std::abs(e_z.z))
		{
			e_x.x = (-1.0) * e_z.z;
			e_x.y = 0.0;
			e_x.z = (1.0) * e_z.x;
		}
		else
		{
			e_x.x = (-1.0) * e_z.y;
			e_x.y = (1.0) * e_z.x;
			e_x.z = 0.0;
		}
		e_x.normalize();
		Vector e_y= e_z.prod_vect(e_x);
		//Stockage de la couleur renvoyée
		Vector accumulated_color(0, 0, 0);

		//On génère une direction aléatoire
		double r1 = distrib(engine);
		double r2 = distrib(engine);
		
		double x = cos(2 * pi * r1) * sqrt(1 - r2);
		double y= sin(2 * pi * r1) * sqrt(1 - r2);
		double z = sqrt(r2);

		Vector dir = x * e_x + y * e_y + z * e_z;

		//On regarde l'eclairage indirect dans cette direction
		double epsilon = 0.01;
		Ray ray_dir(infos_intersection.inter_Pos+epsilon* infos_intersection.inter_Norm, dir);
		intersection_details inter_dir = this->intersection(ray_dir);
		// Si on intersecte autre chose qu'une source
		if (inter_dir.intersection == true && inter_dir.intensite==0.0)
		{
			accumulated_color = accumulated_color + (1 / pi) * infos_intersection.albedo * getIndirect(inter_dir, n_rebonds - 1);
			// Le cos de la PDF s'annule avec celui de l'angle solide projeté
		}
		//Renvoi couleur
		return accumulated_color;
	}
	//Sinon, on renvoie l'éclairage direct
	else
	{
		return this->getColorLambert(infos_intersection);
	}
	
}


