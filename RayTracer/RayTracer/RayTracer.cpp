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
#include <math.h>

int main()
{
	//Dimensions de l'image 
	int W = 512;
	int H = 512;

	//Decription de la scene
	Vector camera_centre(0.0, 0.0, 55.0); // Centre de la caméra
	Sphere sphere1(Vector(0.0, 0.0, 0.0), 10.0); // Objet scene
	double fov = 60; // Angle de vue en degres
	double d = W / (2 * tan(fov * 3.1415 / 360.0)); // distance centre_camera -> plan_image


	//Création de l'image sous forme d'un vecteur
	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			//Calcul du rayon
			Vector X_ij(i - W / 2.0, j - H / 2.0, -d);
			Vector U = X_ij - camera_centre;
			U.normalize();
			Ray ray_ij(camera_centre, U);

			// On determine si il y a intersection
			if (sphere1.intersect(ray_ij))
			{
				image[(i * W + j) * 3 + 0] = 255;
				image[(i * W + j) * 3 + 1] = 255;
				image[(i * W + j) * 3 + 2] = 255;
			}
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
	std::cout << "Image calculée" << std::endl;

	return 0;
}


