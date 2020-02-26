#include "Triangle.h"
#include "iostream"


//Constructeur de classe
Triangle::Triangle(const Vector& ptA, const Vector& ptB, const Vector& ptC, const Vector& color, const double albed, bool miror, bool transp, double ind_sphere, double intensity, int N_phong, double Rho_diff_phong)
{
	A = ptA;
	B = ptB;
	C = ptC;

	Vector AB = B + (-1.0) * A;
	Vector AC = C + (-1.0) * A;
	n = AB.prod_vect(AC);

	albedo = albed; miroir = miror;
	transparent = transp;
	indice_sphere = ind_sphere;
	intensite = intensity;
	n_phong = N_phong;
	rho_diff_phong = Rho_diff_phong;
}


//Intersection rayon triangle
intersection_details Triangle::intersect(Ray& rayon)
{
	intersection_details infos;
	//Intersection rayon-plan

	//Cas ou le rayon et le plan sont parallèles
	if(rayon.direction.dot(n)==0.0)
	{
		//Pt dans triangle
		if ((A + (-1.0) * rayon.origine).dot(n) == 0.0)
		{
			// Details sur l'intersection
			infos.intersection = true;
			infos.t = 0.0;
			infos.inter_Pos = rayon.origine; // Position
			infos.inter_Color = couleur; // Couleur
			infos.inter_Norm = n; // Normale
			infos.albedo = albedo; // Albedo
			infos.miroir = miroir;
			infos.transparent = transparent;
			infos.indice_sphere = indice_sphere;
			infos.intensite = intensite;
			infos.ray_dir = rayon.direction;
			infos.n_phong = n_phong;
			infos.rho_diff_phong = rho_diff_phong;

		}
		//Pt hors triangle
		else
		{
			return infos;
		}
	}
	//Sinon
	else
	{
		double num = (A + (-1.0) * rayon.origine).dot(n);
		double den = rayon.direction.dot(n);
		double t = num / den;

		//Il y a instersection devant
		if (t > 0.0)
		{
			//On verifie si l'intersection se situe dans le triangle
			double a = (B + (-1.0) * A).norme2();
			double b = (C + (-1.0) * A).dot(B + (-1.0) * A);
			double c = (B + (-1.0) * A).dot(C + (-1.0) * A);
			double d = (C + (-1.0) * A).norme2();

			double det =( a * d) - (b * c);

			double y1 = (rayon.origine + t * rayon.direction + (-1.0) * A).dot(B + (-1.0) * A);
			double y2 = (rayon.origine + t * rayon.direction + (-1.0) * A).dot(C + (-1.0) * A);

			double beta = ((d * y1) - (b * y2)) / det;
			double gamma = ((a * y2) - (c * y1)) / det;
			double alpha = 1.0 - beta - gamma;


			//Si coords barycentriques OK : dans le triangle
			if (alpha >= 0.0 && beta >= 0.0 && gamma >= 0.0 && alpha <= 1.0 && beta <= 1.0 && gamma <= 1.0)
			{
				// Details sur l'intersection
				infos.intersection = true;
				infos.t = t;
				infos.inter_Pos = alpha*A+beta*B+gamma*C; // Position
				infos.inter_Color = couleur; // Couleur
				infos.inter_Norm = n; // Normale
				infos.albedo = albedo; // Albedo
				infos.miroir = miroir;
				infos.transparent = transparent;
				infos.indice_sphere = indice_sphere;
				infos.intensite = intensite;
				infos.ray_dir = rayon.direction;
				infos.n_phong = n_phong;
				infos.rho_diff_phong = rho_diff_phong;


			}
			// Sinon : hors du triangle 
			else
			{
				return infos;
			}

		}
		//Il y a intersection derriere
		else
		{
			return infos;
		}
	}
}