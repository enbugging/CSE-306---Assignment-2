#ifndef VECTOR_H
#define VECTOR_H

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#include <cmath>

class Vector {
public:
	explicit Vector(double x = 0, double y = 0) {
		data[0] = x;
		data[1] = y;
	}
	double norm2() const {
		return data[0] * data[0] + data[1] * data[1];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
	}

	Vector operator+=(const Vector& b) {
		data[0] += b[0];
		data[1] += b[1];
		return *this;
	}

	Vector operator-=(const Vector& b) {
		data[0] -= b[0];
		data[1] -= b[1];
		return *this;
	}

	double operator[](int i) const { return data[i]; };
	double& operator[](int i) { return data[i]; };
	double data[2];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1]);
}
Vector operator-(const Vector& a) {
	return Vector(-a[0], -a[1]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1];
}

#define M_PI 3.14159265358979323846

double sqr(double x) { return x * x; }

#endif