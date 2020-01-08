#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>
#include "Vector.h"
#include "Sphere.h"
#include "Ray.h"
#include "Scene.h"
#include <math.h>


int main()
{
	//Dimensions de l'image 
	int W = 512;
	int H = 512;

	//Decription de la scene
	Scene scene;
	Vector camera_centre(0.0, 0.0, 55.0); // Centre de la caméra
	//// Creation des spheres
	Sphere sphere1(Vector(0.0, 0.0, 0.0), 10.0, Vector(255.0, 100.0, 50.0), 1.0); // Objet sphere (centre,rayon,couleur,albedo
	Sphere sphere2(Vector(20.0, 0.0, 0.0), 3.0, Vector(0.0, 200.0, 50.0), 1.0); // Objet sphere (centre,rayon,couleur,albedo
	Sphere sphere3(Vector(-20.0, 0.0, 0.0), 5.0, Vector(0.0, 0.0, 150.0), 1.0); // Objet sphere (centre,rayon,couleur,albedo
	//// Ajout des spheres a la scene
	scene.add_sphere(&sphere1); 
	scene.add_sphere(&sphere2);
	scene.add_sphere(&sphere3);
	//// Caracs camera
	double fov = 60; // Angle de vue en degres
	double d = W / (2 * tan(fov * 3.1415 / 360.0)); // distance centre_camera -> plan_image
	//// Caracs source
	Vector light_Pos(-10.0, 20.0, 40.0); // Position source lumineuse
	double light_power = 1000000.0; // Intensité de la source lumineuse


	//Création de l'image sous forme d'un vecteur
	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			//Calcul du rayon
			Vector X_ij(j - W / 2.0, -i + H / 2.0, -d);
			Vector U = X_ij - camera_centre;
			U.normalize();
			Ray ray_ij(camera_centre, U);

			// On determine si il y a intersection
			intersection_details infos_intersection = scene.intersection(ray_ij);
			// Si il y a une intersection
			if (infos_intersection.intersection)
			{
				//Calcul de la couleur pour le modèle Lambertien (diffus)
				Vector to_Light = light_Pos - infos_intersection.inter_Pos;
				to_Light.normalize();
				double cos_theta = infos_intersection.inter_Norm.dot(to_Light);
				if (cos_theta < 0.0) { cos_theta = 0.0; }
				double facteur = light_power * infos_intersection.albedo *cos_theta / 3.1415 / ((light_Pos - infos_intersection.inter_Pos).norme2());
				Vector RGB = infos_intersection.inter_Color;
				RGB = facteur * RGB;

				// Correction Gamma
				RGB.x = pow(RGB.x, 0.45);
				RGB.y = pow(RGB.y, 0.45);
				RGB.z = pow(RGB.z, 0.45);

				// Clamping
				if (RGB.x > 255) { RGB.x = 255; }
				if (RGB.y > 255) { RGB.y = 255; }
				if (RGB.z > 255) { RGB.z = 255; }

				// Ecriture de la couleur
				image[(i * W + j) * 3 + 0] = (int)RGB.x;
				image[(i * W + j) * 3 + 1] = (int)RGB.y;
				image[(i * W + j) * 3 + 2] = (int)RGB.z;

			}
			// Sinon
			else
			{
				image[(i * W + j) * 3 + 0] = 0;
				image[(i * W + j) * 3 + 1] = 0;
				image[(i * W + j) * 3 + 2] = 0;
			}

		}
	}
	//Enregistrement de l'image sur le disque dur
	stbi_write_png("image.png", W, H, 3, &image[0], 0);
	std::cout << "Image calculee" << std::endl;

	return 0;
}


