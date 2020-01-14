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
	Sphere sphere1(Vector(0.0, 0.0, 0.0), 15.0, Vector(105.0, 105.0, 0.0), 1.0,true,false,1.33); // Objet sphere (centre,rayon,couleur,albedo,miroir,transparent,indice)
	Sphere sphere2(Vector(0.0, 0.0, -1000.0), 940.0, Vector(0.0, 255.0, 0.0), 1.0, false,false,1.0); // Objet sphere (centre,rayon,couleur,albedo,miroir,transparent,indice)
	Sphere sphere3(Vector(0.0, 1000.0, 0.0), 940.0, Vector(255.0, 0.0, 0.0), 1.0, false, false, 1.0); // Objet sphere (centre,rayon,couleur,albedo,miroir,transparent,indice)
	////
	Sphere sphere4(Vector(10.0, 15.0, 10.0), 4.0, Vector(0.0, 255.0, 255.0), 1.0, false, true, 2.33); // Objet sphere (centre,rayon,couleur,albedo,miroir,transparent,indice)
	////
	Sphere sphere5(Vector(0.0, -1000.0, 0.0), 940.0, Vector(0.0, 0, 255.0), 1.0, false, false, 1.0); // Objet sphere (centre,rayon,couleur,albedo,miroir,transparent,indice)
	Sphere sphere6(Vector(0.0, 0.0, 1000.0), 940.0, Vector(255.0, 255.0, 0.0), 1.0, false, false, 1.0); // Objet sphere (centre,rayon,couleur,albedo,miroir,transparent,indice)
	Sphere sphere7(Vector(1000.0, 0.0, 0.0), 940.0, Vector(255.0, 0.0, 255.0), 1.0, false, false, 1.0); // Objet sphere (centre,rayon,couleur,albedo,miroir,transparent,indice)
	Sphere sphere8(Vector(-1000.0, 0.0, 0.0), 940.0, Vector(255.0, 255.0, 255.0), 1.0, false, false, 1.0); // Objet sphere (centre,rayon,couleur,albedo,miroir,transparent,indice)
	//// Ajout des spheres a la scene
	scene.add_sphere(&sphere1); 
	scene.add_sphere(&sphere2);
	scene.add_sphere(&sphere3);
	scene.add_sphere(&sphere4);
	scene.add_sphere(&sphere5);
	scene.add_sphere(&sphere6);
	scene.add_sphere(&sphere7);
	scene.add_sphere(&sphere8);
	//// Caracs camera
	double fov = 60; // Angle de vue en degres
	double d = H / (2 * tan(fov * 3.1415 / 360.0)); // distance centre_camera -> plan_image
	//// Caracs source
	Vector light_Position(-10.0,0.0, 40.0); // Position source lumineuse
	double light_Power = 2000000.0; // Intensité de la source lumineuse
	light_source main_light;
	main_light.light_Pos = light_Position;
	main_light.intensite = light_Power;
	scene.add_light(&main_light);


	//Création de l'image sous forme d'un vecteur
	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			//Calcul du rayon
			Vector X_ij(-i + W / 2.0 , -j + H / 2.0 ,-d);
			//Vector U = X_ij - camera_centre;
			Vector U = X_ij;
			U.normalize();
			Ray ray_ij(camera_centre, U);

			// On determine si il y a intersection
			intersection_details infos_intersection = scene.intersection(ray_ij);
			// Si il y a une intersection
			if (infos_intersection.intersection)
			{
				// Variable pour stocker la couleur du pixel
				Vector RGB;
				RGB = scene.getColor(ray_ij, infos_intersection, 1.0, 5, 5);
				
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
			// Sinon, pas d'intersection
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


