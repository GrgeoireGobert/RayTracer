#pragma once
#include <vector>
#include <string>
#include <map>
#include "Vector.h"
#include "Object.h"

class BBox
{
public:
	double xmin, xmax, ymin, ymax, zmin, zmax;
	Vector centre;

	BBox()
	{
		BBox(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	}
	BBox(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)
	{
		xmin = x_min;
		xmax = x_max;
		ymin = y_min;
		ymax = y_max;
		zmin = z_min;
		zmax = z_max;
		centre.x = 0.5 * (xmin + xmax);
		centre.y = 0.5 * (ymin + ymax);
		centre.z = 0.5 * (zmin + zmax);
	}

	bool is_inside(Vector& point)
	{
		if (point.x <= xmax && point.x >= xmin)
		{
			if (point.y <= ymax && point.y >=ymin)
			{
				if (point.z <= zmax && point.z >= zmin)
				{
					return true;
				}
			}
		}
		return false;
	}

	bool intersection(Ray& rayon)
	{
		Vector normal(0.0, 0.0, 0.0);
		Vector A(0.0, 0.0, 0.0);
		double tx_min, tx_max, ty_min, ty_max, tz_min, tz_max;
		bool parallele_x = false;
		bool parallele_y = false;
		bool parallele_z = false;

		double denom;

		//Plans (x,y)
		normal = Vector(0.0, 0.0, 1.0);
		//min
		A = Vector(centre.x, centre.y, zmin);
		denom = rayon.direction.dot(normal);
		if (denom == 0) { parallele_z = true; }
		else { tz_min = ((rayon.origine + (-1.0) * A).dot(normal)) / denom; }
		//max
		A = Vector(centre.x, centre.y, zmax);
		denom = rayon.direction.dot(normal);
		if (denom == 0) { parallele_z = true; }
		else { tz_max = ((rayon.origine + (-1.0) * A).dot(normal)) / denom; }


		//Plans (x,z)
		normal = Vector(0.0, 1.0, 0.0);
		//min
		A = Vector(centre.x, ymin, centre.z);
		denom = rayon.direction.dot(normal);
		if (denom == 0) { parallele_y = true; }
		else { ty_min = ((rayon.origine + (-1.0) * A).dot(normal)) / denom; }
		//max
		A = Vector(centre.x, ymax, centre.z);
		denom = rayon.direction.dot(normal);
		if (denom == 0) { parallele_y = true; }
		else { ty_max = ((rayon.origine + (-1.0) * A).dot(normal)) / denom; }


		//Plans (y,z)
		normal = Vector(1.0, 0.0, 0.0);
		//min
		A = Vector(xmin, centre.y, centre.z);
		denom = rayon.direction.dot(normal);
		if (denom == 0) { parallele_x = true; }
		else { tx_min = ((rayon.origine + (-1.0) * A).dot(normal)) / denom; }
		//max
		A = Vector(xmax, centre.y, centre.z);
		denom = rayon.direction.dot(normal);
		if (denom == 0) { parallele_x = true; }
		else { tx_max = ((rayon.origine + (-1.0) * A).dot(normal)) / denom; }


		//Cas parallele
		if (parallele_x == true)
		{
			if (rayon.origine.y >= ymin && rayon.origine.y <= ymax && rayon.origine.z >= zmin && rayon.origine.z <= ymax)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		if (parallele_y == true)
		{
			if (rayon.origine.x >= xmin && rayon.origine.x <= xmax && rayon.origine.z >= zmin && rayon.origine.z <= ymax)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		if (parallele_z == true)
		{
			if (rayon.origine.y >= ymin && rayon.origine.y <= ymax && rayon.origine.x >= xmin && rayon.origine.x <= xmax)
			{
				return true;
			}
			else
			{
				return false;
			}
		}

		//Cas general
		//On reordonne
		if (tx_min > tx_max) { double save = tx_max; tx_max = tx_min; tx_min = save; }
		if (ty_min > ty_max) { double save = ty_max; ty_max = ty_min; ty_min = save; }
		if (tz_min > tz_max) { double save = tz_max; tz_max = tz_min; tz_min = save; }

		double max_des_min = tx_min;
		if (ty_min > max_des_min) { max_des_min = ty_min; };
		if (tz_min > max_des_min) { max_des_min = tz_min; };

		double min_des_max = tx_max;
		if (ty_max < min_des_max) { min_des_max = ty_max; }
		if (tz_max < min_des_max) { min_des_max = tz_max; }

		return (max_des_min < min_des_max);
	}
};


//Structure node pour l'arbre
struct node
{
	BBox bbox;
	int begin=0;
	int end=0;
	node* left_leaf=nullptr;
	node* right_leaf=nullptr;
};


class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk) {
	};
	int vtxi, vtxj, vtxk;
	int uvi, uvj, uvk;
	int ni, nj, nk;
	int faceGroup;
};

class Geometry : public Object {
public:
	Geometry() {};
	Geometry(const char* obj, double scaling, const Vector& offset) {
		readOBJ(obj);
		for (int i = 0; i < vertices.size(); i++) {
			vertices[i] = vertices[i] * scaling + offset;
		}

		/// AFFICHAGE BOUNDING BOX
		double xmin, xmax, ymin, ymax, zmin, zmax;
		xmin = vertices.at(0).x;
		xmax = vertices.at(0).x;
		ymin = vertices.at(0).y;
		ymax = vertices.at(0).y;
		zmin = vertices.at(0).z;
		zmax = vertices.at(0).z;
		for (int i = 0; i < vertices.size(); i++)
		{
			if (vertices.at(i).x < xmin) { xmin = vertices.at(i).x; }
			if (vertices.at(i).x > xmax) { xmax = vertices.at(i).x; }
			if (vertices.at(i).y < ymin) { ymin = vertices.at(i).y; }
			if (vertices.at(i).y > ymax) { ymax = vertices.at(i).y; }
			if (vertices.at(i).z < zmin) { zmin = vertices.at(i).z; }
			if (vertices.at(i).z > zmax) { zmax = vertices.at(i).z; }
		}
		std::cout << "------------------" << std::endl;
		std::cout << vertices.size() << " sommets" << std::endl;
		std::cout << indices.size() << " triangles" << std::endl;
		std::cout << " " << std::endl;
		std::cout << "x : " << xmin << " --> " << xmax << std::endl;
		std::cout << "y : " << ymin << " --> " << ymax << std::endl;
		std::cout << "z : " << zmin << " --> " << zmax << std::endl;
		std::cout << "------------------" << std::endl;

		//Création de la BBox
		bounding_box = BBox(xmin, xmax, ymin, ymax, zmin, zmax);
		//Création de la racine
		root.bbox = bounding_box;
		root.begin = 0;
		root.end = indices.size();
		//Construction de la hierachie
		build_hierarchy(&root);
		std::cout << "BVH calcul termine" << std::endl;


	}

	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");

		std::map<std::string, int> groupNames;
		int curGroup = -1;
		while (!feof(f)) {

			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				if (groupNames.find(std::string(grp)) != groupNames.end()) {
					curGroup = groupNames[std::string(grp)];
				}
				else {
					curGroup = groupNames.size();
					groupNames[std::string(grp)] = curGroup;
				}
			}
			if (line[0] == 'm' && line[1] == 't' && line[2] == 'l') {
				sscanf(line, "mtllib %[^\n]\n", matfile);
			}
			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;
				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec.x, &vec.z, &vec.y, &col.x, &col.y, &col.z) == 6) {
					vertices.push_back(vec);
					vertexcolors.push_back(col);
				}
				else {
					sscanf(line, "v %lf %lf %lf\n", &vec.x, &vec.z, &vec.y);  // helmet
																				 //vec.z = -vec.z; //car2
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf_s(line, "vn %lf %lf %lf\n", &vec.x, &vec.z, &vec.y); //girl
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec.x, &vec.y);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;

				char* consumedline = line + 1;
				int offset;
				t.faceGroup = curGroup;
				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;

					indices.push_back(t);
				}
				else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					}
					else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						}
						else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}


				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.faceGroup = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					}
					else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						}
						else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							}
							else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								}
								else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}


		}
		fclose(f);

		
	}

	void add_texture(const char* filename) {

		int width, height, channels;
		stbi_set_flip_vertically_on_load(true);
		unsigned char* image = stbi_load(filename, &width, &height, &channels, STBI_rgb);
		textures.push_back(image);
		w.push_back(width);
		h.push_back(height);
		std::cout << "Texture OK" << std::endl;
	}

	std::vector<double> compute_bounding_box(std::vector<int> tri_inds)
	{
		if (tri_inds.size() == 0)
		{
			return { 0,0,0,0,0,0 };
		}
		double xmin, ymin, zmin, xmax, ymax, zmax;

		xmin = vertices.at(indices.at(tri_inds.at(0)).vtxi).x;
		xmax = vertices.at(indices.at(tri_inds.at(0)).vtxi).x;
		ymin = vertices.at(indices.at(tri_inds.at(0)).vtxi).y;
		ymax = vertices.at(indices.at(tri_inds.at(0)).vtxi).y;
		zmin = vertices.at(indices.at(tri_inds.at(0)).vtxi).z;
		zmax = vertices.at(indices.at(tri_inds.at(0)).vtxi).z;

		//Parcourt des traingles
		for (int i = 0; i < tri_inds.size(); i++)
		{
			//Sommet i
			double x = vertices.at(indices.at(tri_inds.at(i)).vtxi).x;
			double y = vertices.at(indices.at(tri_inds.at(i)).vtxi).y;
			double z = vertices.at(indices.at(tri_inds.at(i)).vtxi).z;

			if (x < xmin) { xmin = x; }
			if (x > xmax) { xmax = x; }
			if (y < ymin) { ymin = y; }
			if (y > ymax) { ymax = y; }
			if (z < zmin) { zmin = z; }
			if (z > zmax) { zmax = z; }

			//Sommet j
			x = vertices.at(indices.at(tri_inds.at(i)).vtxj).x;
			y = vertices.at(indices.at(tri_inds.at(i)).vtxj).y;
			z = vertices.at(indices.at(tri_inds.at(i)).vtxj).z;

			if (x < xmin) { xmin = x; }
			if (x > xmax) { xmax = x; }
			if (y < ymin) { ymin = y; }
			if (y > ymax) { ymax = y; }
			if (z < zmin) { zmin = z; }
			if (z > zmax) { zmax = z; }

			//Sommet k
			x = vertices.at(indices.at(tri_inds.at(i)).vtxk).x;
			y = vertices.at(indices.at(tri_inds.at(i)).vtxk).y;
			z = vertices.at(indices.at(tri_inds.at(i)).vtxk).z;

			if (x < xmin) { xmin = x; }
			if (x > xmax) { xmax = x; }
			if (y < ymin) { ymin = y; }
			if (y > ymax) { ymax = y; }
			if (z < zmin) { zmin = z; }
			if (z > zmax) { zmax = z; }
		}

		return { xmin,xmax,ymin,ymax,zmin,zmax };
	}

	//Construction Bounding Volumes Hierarchy de manière récursive
	void build_hierarchy(node* tree_node)
	{
		//condition d'arret
		if ((tree_node->end - tree_node->begin) > 10) //Le nombre de triangles max dans chaque feuille 
		{
			//On cherche la plus grande direction
			double x_length = tree_node->bbox.xmax - tree_node->bbox.xmin;
			double y_length = tree_node->bbox.ymax - tree_node->bbox.ymin;
			double z_length = tree_node->bbox.zmax - tree_node->bbox.zmin;

			std::vector<int> left_node_indices;
			std::vector<int> right_node_indices;

			std::vector<double> left_node_bounds;
			std::vector<double> right_node_bounds;

			//Si plus grande direction : x
			if (x_length > y_length&& x_length > z_length)
			{
				double xmoy = 0.5 * (tree_node->bbox.xmax + tree_node->bbox.xmin);

				//On parcourt les triangles du pere
				for (int i = tree_node->begin; i < tree_node->end; i++)
				{
					//Calcul du centre du ieme triangle
					TriangleIndices& triangle = indices.at(i);
					Vector centre = vertices.at(triangle.vtxi);
					centre = centre + vertices.at(triangle.vtxj);
					centre = centre + vertices.at(triangle.vtxk);
					centre = (1.0 / 3.0) * centre;

					//On check si le centre du triangle est dans le fils de gauche
					if (centre.x < xmoy)
					{
						left_node_indices.push_back(i);
					}
					//Sinon, il est dans le fils de droite
					else
					{
						right_node_indices.push_back(i);
					}
				}
			}
			//Si plus grande direction : y
			else if (y_length > x_length&& y_length > z_length)
			{
				double ymoy = 0.5 * (tree_node->bbox.ymax + tree_node->bbox.ymin);
				//On parcourt les triangles du pere
				for (int i = tree_node->begin; i < tree_node->end; i++)
				{
					//Calcul du centre du ieme triangle
					TriangleIndices& triangle = indices.at(i);
					Vector centre = vertices.at(triangle.vtxi);
					centre = centre + vertices.at(triangle.vtxj);
					centre = centre + vertices.at(triangle.vtxk);
					centre = (1.0 / 3.0) * centre;

					//On check si le centre du triangle est dans le fils de gauche
					if (centre.y < ymoy)
					{
						left_node_indices.push_back(i);
					}
					//Sinon, il est dans le fils de droite
					else
					{
						right_node_indices.push_back(i);
					}
				}
			}
			//Si plus grande direction : z
			else
			{
				double zmoy = 0.5 * (tree_node->bbox.zmax + tree_node->bbox.zmin);
				//On parcourt les triangles du pere
				for (int i = tree_node->begin; i < tree_node->end; i++)
				{
					//Calcul du centre du ieme triangle
					TriangleIndices& triangle = indices.at(i);
					Vector centre = vertices.at(triangle.vtxi);
					centre = centre + vertices.at(triangle.vtxj);
					centre = centre + vertices.at(triangle.vtxk);
					centre = (1.0 / 3.0) * centre;

					//On check si le centre du triangle est dans le fils de gauche
					if (centre.z < zmoy)
					{
						left_node_indices.push_back(i);
					}
					//Sinon, il est dans le fils de droite
					else
					{
						right_node_indices.push_back(i);
					}
				}
			}
			//On va calculer la bounding box de gauche et celle de droite
			left_node_bounds = compute_bounding_box(left_node_indices);
			right_node_bounds = compute_bounding_box(right_node_indices);
			//On réordonne la liste
			std::vector<TriangleIndices> indices_copy(indices);
			//std::vector<TriangleIndices> indices_copy;
			//for (int i = 0; i < indices.size(); i++) { indices_copy.push_back(indices.at(i)); }
			//D'abord ceux de gauche
			for (int k = 0; k < left_node_indices.size(); k++)
			{
				indices.at(tree_node->begin + k) = indices_copy.at(left_node_indices.at(k));
			}
			//Puis ceux de droite
			for (int k = 0; k < right_node_indices.size(); k++)
			{
				indices.at(tree_node->begin + left_node_indices.size() + k) = indices_copy.at(right_node_indices.at(k));
			}


			tree_node->left_leaf = new node;
			tree_node->left_leaf->bbox = BBox(left_node_bounds.at(0), left_node_bounds.at(1), left_node_bounds.at(2), left_node_bounds.at(3), left_node_bounds.at(4), left_node_bounds.at(5));
			tree_node->left_leaf->begin = tree_node->begin;
			tree_node->left_leaf->end = tree_node->begin + left_node_indices.size();

			tree_node->right_leaf = new node;
			tree_node->right_leaf->bbox = BBox(right_node_bounds.at(0), right_node_bounds.at(1), right_node_bounds.at(2), right_node_bounds.at(3), right_node_bounds.at(4), right_node_bounds.at(5));
			tree_node->right_leaf->begin = tree_node->begin + left_node_indices.size();
			tree_node->right_leaf->end = tree_node->begin + left_node_indices.size()+ right_node_indices.size();


			if (left_node_indices.size() && right_node_indices.size() > 0) { build_hierarchy(tree_node->left_leaf); }
			if (left_node_indices.size() && right_node_indices.size() > 0) { build_hierarchy(tree_node->right_leaf); }
		}
	}


	std::vector<int> check_intersection(Ray& rayon, node& tree_node)
	{
		//Si on a une intersection
		if (tree_node.bbox.intersection(rayon) == true)
		{
			if (tree_node.left_leaf != nullptr && tree_node.left_leaf->bbox.intersection(rayon) == true && tree_node.right_leaf != nullptr && tree_node.right_leaf->bbox.intersection(rayon) == true)
			{
				std::vector<int> res_left = check_intersection(rayon, *(tree_node.left_leaf));
				std::vector<int> res_right = check_intersection(rayon, *(tree_node.right_leaf));
				return { res_left.at(0),res_right.at(1) };
				//return { tree_node.left_leaf->begin,tree_node.right_leaf->end };
			}
			else if (tree_node.left_leaf!= nullptr && tree_node.left_leaf->bbox.intersection(rayon) == true)
			{
				return check_intersection(rayon, *(tree_node.left_leaf));
				//return { tree_node.left_leaf->begin,tree_node.left_leaf->end };
			}
			else if (tree_node.right_leaf != nullptr &&  tree_node.right_leaf->bbox.intersection(rayon) == true)
			{
				return check_intersection(rayon, *(tree_node.right_leaf));
				//return { tree_node.right_leaf->begin,tree_node.right_leaf->end };
			}
			else
			{
				return { tree_node.begin,tree_node.end };
			}
		}
		else
		{
			return { 0,0 };
		}
	}


	// Intersection rayon-mesh
	virtual intersection_details intersect(Ray& rayon)
	{

		//On teste l'intersection les BBoxs
		
		
		std::vector<int> boundaries = check_intersection(rayon, root);

		int left_bound = boundaries.at(0);
		int right_bound = boundaries.at(1);

		

		/*
		//Si pas d'intersection
		if (bounding_box.intersection(rayon)==false)
		{
			intersection_details renvoi;
			return renvoi;
		}
		int left_bound = 0;
		int right_bound = indices.size();
		*/

		//Sinon, on continue
		std::vector<intersection_details> all_intersections;

		std::vector<int> Tri_inds;
		//Parcourt les triangles
		for (int t = left_bound; t < right_bound; t++)
		{
			//creation du triangle

			//A MODIFIER POUR CHOISIR LA BRDF VOULUE
			Triangle tri(vertices.at(indices.at(t).vtxi), vertices.at(indices.at(t).vtxj), vertices.at(indices.at(t).vtxk), Vector(200.0, 200.0, 200.0), 1.0,false, false, 1.0, 0.0, 0, 0.0);
			
			// Inversion de la normale si besoin !!!
			tri.n = (-1.0) * tri.n;

			//On checke l'intersection
			intersection_details infos = tri.intersect(rayon);
			if (infos.intersection == true)
			{ 
				all_intersections.push_back(infos); 
				Tri_inds.push_back(t);
			}
		}

		//Si pas d'interscetion
		if (all_intersections.size() == 0)
		{
			intersection_details renvoi;
			return renvoi;
		}
		//Sinon
		else
		{
			intersection_details renvoi;
			renvoi = all_intersections.at(0);
			int final_tri_ind=Tri_inds.at(0);
			for (int i = 1; i < all_intersections.size(); i++)
			{
				if (all_intersections.at(i).t < renvoi.t) 
				{ 
					renvoi = all_intersections.at(i); 
					final_tri_ind = Tri_inds.at(i);
				}
			}

			/////////////
			// TEXTURE //
			/////////////
			
			//On ne gere que le cas avec 1 texture
			if (textures.size() == 1)
			{
				int tri_ind = final_tri_ind;
				// Les 3 sommets:
				Vector A = vertices.at(indices.at(tri_ind).vtxi);
				Vector B = vertices.at(indices.at(tri_ind).vtxj);
				Vector C = vertices.at(indices.at(tri_ind).vtxk);
				Triangle tri(A, B, C, Vector(0, 0, 0), 0.5, false, false, 1.0, 0.0, 0, 0.0);
				// Les uv des 3 sommets:
				Vector UV_A = uvs.at(indices.at(tri_ind).uvi);
				Vector UV_B = uvs.at(indices.at(tri_ind).uvj);
				Vector UV_C = uvs.at(indices.at(tri_ind).uvk);
				// Les coords barycentriques du pt d'intersection:
				Vector P = renvoi.inter_Pos;
				// Calcul des aires
				double ABC = 0.5 * ((B - A).prod_vect(C - A)).dot(tri.n);
				double PBC = 0.5 * ((B - P).prod_vect(C - P)).dot(tri.n);
				double APC = 0.5 * ((P - A).prod_vect(C - A)).dot(tri.n);
				double ABP = 0.5 * ((B - A).prod_vect(P - A)).dot(tri.n);
				// Coords barycentriques
				double a = PBC / ABC;
				double b = APC / ABC;
				double c = ABP / ABC;
				// Coords UV du pt d'intersection
				Vector UV_P = a * UV_A + b * UV_B + c * UV_C;
				Vector local_Norm = a * normals.at(indices.at(tri_ind).ni) + b * normals.at(indices.at(tri_ind).nj) + c * normals.at(indices.at(tri_ind).nk);

				//Cas où la texture se répète
				UV_P.x = UV_P.x - (int)(UV_P.x);
				UV_P.y = UV_P.y - (int)(UV_P.y);
				UV_P.z = UV_P.z - (int)(UV_P.z);

				//Image texture
				unsigned char* img_tex = textures.at(0);

				// Valeurs pixels
				auto pixel_R = (int)img_tex[(int)(UV_P.y * h.at(0)) * w.at(0) * 3 + (int)(UV_P.x * w.at(0)) * 3];
				auto pixel_V = (int)img_tex[(int)(UV_P.y * h.at(0)) * w.at(0) * 3 + (int)(UV_P.x * w.at(0)) * 3 + 1];
				auto pixel_B = (int)img_tex[(int)(UV_P.y * h.at(0)) * w.at(0) * 3 + (int)(UV_P.x * w.at(0)) * 3 + 2];

				//Modification couleur et normale renvoyées
				renvoi.inter_Color = Vector((double)pixel_R, (double)pixel_V, (double)pixel_B);
				renvoi.inter_Norm = local_Norm;
			}

			return renvoi;
		}


	}


	std::vector<TriangleIndices> indices;
	node root;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs; // Vector en 3D mais on n'utilise que 2 composantes
	std::vector<Vector> vertexcolors;
	BBox bounding_box;

	std::vector<unsigned char*> textures;
	std::vector<int> w, h;
};
