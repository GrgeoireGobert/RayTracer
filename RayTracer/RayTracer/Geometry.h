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

	bool intersection(Ray& rayon)
	{
		//Plans (x,y):
		Vector n(0.0, 0.0, 1.0);
		Vector pt_plan_min(centre.x, centre.y, zmin);
		double tz_min;
		double den = rayon.direction.dot(n);
		if (den == 0) { tz_min = 10000000000000.0; }
		else 
		{
			double num = (pt_plan_min + (-1.0) * rayon.origine).dot(n);
			tz_min = num/den;
		}
		Vector pt_plan_max(centre.x, centre.y, zmax);
		double tz_max;
		den = rayon.direction.dot(n);
		if (den == 0) { tz_max = 10000000000000.0; }
		else
		{
			double num = (pt_plan_max + (-1.0) * rayon.origine).dot(n);
			tz_max = num / den;
		}
		//Plans (x,z):
		n= Vector(0.0, 1.0, 0.0);
		pt_plan_min=Vector(centre.x, ymin, centre.z);
		double ty_min;
		den = rayon.direction.dot(n);
		if (den == 0) { ty_min = 10000000000000.0; }
		else
		{
			double num = (pt_plan_min + (-1.0) * rayon.origine).dot(n);
			ty_min = num / den;
		}
		pt_plan_max= Vector(centre.x, ymax, centre.z);
		double ty_max;
		den = rayon.direction.dot(n);
		if (den == 0) { ty_max = 10000000000000.0; }
		else
		{
			double num = (pt_plan_max + (-1.0) * rayon.origine).dot(n);
			ty_max = num / den;
		}
		//Plans (y,z):
		n= Vector(1.0, 0.0, 0.0);
		pt_plan_min= Vector(xmin, centre.y, centre.z);
		double tx_min;
		den = rayon.direction.dot(n);
		if (den == 0) { tx_min = 10000000000000.0; }
		else
		{
			double num = (pt_plan_min + (-1.0) * rayon.origine).dot(n);
			tx_min = num / den;
		}
		pt_plan_max= Vector(xmax, centre.y, centre.z);
		double tx_max;
		den = rayon.direction.dot(n);
		if (den == 0) { tx_max = 10000000000000.0; }
		else
		{
			double num = (pt_plan_max +(-1.0)* rayon.origine).dot(n);
			tx_max = num / den;
		}

		//
		double max_des_min = tx_min;
		if (ty_min > max_des_min) { max_des_min = ty_min; }
		if (tz_min > max_des_min) { max_des_min = tz_min; }
		double min_des_max = tx_max;
		if (ty_max < min_des_max) { min_des_max = ty_max; }
		if (tz_max < min_des_max) { min_des_max = tz_max; }

		return (max_des_min > min_des_max);
	}
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

		bounding_box = BBox(xmin, xmax, ymin, ymax, zmin, zmax);
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

		textures.resize(textures.size() + 1);
		w.resize(w.size() + 1);
		h.resize(h.size() + 1);

		FILE* f;
		f = fopen(filename, "rb");
		unsigned char info[54];
		fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

		w[w.size() - 1] = *(int*)&info[18]; // extract image height and width from header
		h[h.size() - 1] = *(int*)&info[22];

		int size = 3 * w[w.size() - 1] * h[h.size() - 1];
		textures[textures.size() - 1].resize(size); // allocate 3 bytes per pixel
		fread(&textures[textures.size() - 1][0], sizeof(unsigned char), size, f); // read the rest of the data at once
		fclose(f);

		for (int i = 0; i < size; i += 3) {
			std::swap(textures[textures.size() - 1][i], textures[textures.size() - 1][i + 2]);
		}
	}


	virtual intersection_details intersect(Ray& rayon)
	{

		//On teste d'abord l'intersection avec la BBox
		if (bounding_box.intersection(rayon)==false)
		{
			intersection_details renvoi;
			return renvoi;
		}


		std::vector<intersection_details> all_intersections;

		//Parcourt les triangles
		for (int t = 0; t < indices.size(); t++)
		{
			//creation du triangle
			Triangle tri(vertices.at(indices.at(t).vtxi), vertices.at(indices.at(t).vtxj), vertices.at(indices.at(t).vtxk), Vector(255.0, 255.0, 255.0), 1.0);
			
			// Inversion de la normale si besoin !!!
			tri.n = (-1.0) * tri.n;

			//On checke l'intersection
			intersection_details infos = tri.intersect(rayon);
			if (infos.intersection == true) { all_intersections.push_back(infos); }
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
			for (int i = 1; i < all_intersections.size(); i++)
			{
				if (all_intersections.at(i).t < renvoi.t) { renvoi = all_intersections.at(i); }
			}
			return renvoi;
		}


	}


	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs; // Vector en 3D mais on n'utilise que 2 composantes
	std::vector<Vector> vertexcolors;
	BBox bounding_box;

	std::vector<std::vector<unsigned char> > textures;
	std::vector<int> w, h;
};
