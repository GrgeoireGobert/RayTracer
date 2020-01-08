#include "Vector.h"
#include <assert.h>
#include <math.h>

///
/// Les méthodes
///

// Constructeur par défaut
Vector::Vector()
{
	x = 0.0;
	y = 0.0;
	z = 0.0;
}
// Constructeur
Vector::Vector(double X, double Y, double Z)
{
	x = X;
	y = Y;
	z = Z;
}
// Norme au carré
double Vector::norme2()
{
	double norme2 = 0.0;
	norme2 += x * x;
	norme2 += y * y;
	norme2 += z * z;

	return norme2;
}
// Norme
double Vector::norme()
{
	double norme = 0.0;
	norme += x * x;
	norme += y * y;
	norme += z * z;

	return sqrt(norme);
}
// Normalisation
void Vector::normalize()
{
	double norme = this->norme();
	assert(norme>0.0);
	x = x / norme;
	y = y / norme;
	z = z / norme;
}
// Produit scalaire
double Vector::dot(const Vector& B)
{
	double result = 0.0;
	result += this->x * B.x;
	result += this->y * B.y;
	result += this->z * B.z;

	return result;
}
// Produit vectoriel
Vector Vector::prod_vect(const Vector& B)
{
	Vector Result(0.0, 0.0, 0.0);
	Result.x = this->y * B.z - this->z * B.y;
	Result.y = this->z * B.x - this->x * B.z;
	Result.z = this->x * B.y - this->y * B.x;

	return Result;
}


///
/// Les opérateurs externes
///

// Addition
Vector operator+(const Vector& A, const Vector& B)
{
	Vector Result(A.x+B.x, A.y+B.y, A.z+B.z);
	return Result;
}
// Soustraction
Vector operator-(const Vector& A, const Vector& B)
{
	Vector Result(A.x - B.x, A.y - B.y, A.z - B.z);
	return Result;
}
// Multiplication
Vector operator*(const Vector& A, const double k)
{
	Vector Result(k * A.x, k * A.y, k * A.z);
	return Result;
}
// Multiplication
Vector operator*(const double k,const Vector& A)
{
	Vector Result(k * A.x, k * A.y, k * A.z);
	return Result;
}
// Division
Vector operator/(const Vector& A, const double k)
{
	assert(k != 0.0);
	Vector Result(A.x/k, A.y/k, A.z/k);
	return Result;
}