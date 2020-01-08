#pragma once

//
// Classe qui d�finit un vecteur � 3 composantes
// Avec des op�rations (+,*,/,normalisation,norme2,...)
// 

class Vector
{
public :

	double x;
	double y;
	double z;

	Vector();
	Vector(double X, double Y, double Z);
	double norme();
	double norme2();
	void normalize();
	double dot(const Vector& B);
	Vector prod_vect(const Vector& B);
};


///Les op�rateurs externes
Vector operator+(const Vector& A, const Vector& B);
Vector operator-(const Vector& A, const Vector& B);
Vector operator*(const Vector& A, const double k);
Vector operator*(const double k, const Vector& A);
Vector operator/(const Vector& A, const double k);

