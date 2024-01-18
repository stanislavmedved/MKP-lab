#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#define PI 3.14159265358

const double e = 0.6566970091;
const double eps = 0.00001;
const double a = 7690000;
const double mu = 2.2297143 * pow(10, 13);
const double M = 0.33 * pow(10, 24);

double timeFind() {
	return (2 * PI * sqrt(a * a * a / mu));
}
double nFind() {
	return 2 * PI / timeFind();
}
double radFind(double theta, double p) {
	return p / (1 + e * cos(theta));
}
double speedRad(double p, double theta) {
	return sqrt(mu / p) * e * sin(theta);
}
double speedN(double p, double theta) {
	return sqrt(mu / p) * (1 + e * cos(theta));
}
double fullSpeed(double p, double theta) {
	double rad = speedRad(p, theta);
	double n = speedN(p, theta);
	return sqrt(rad * rad + n * n);
}
double excentricToTrue(double E) {
	if (atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)) * 2 > 0) {
		return atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)) * 2;
	}
	else {
		return atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)) * 2 + 2 * PI;
	}
}
double iterationMethod(double ENext, double ENow, double M) {
	if (fabs(ENow - ENext) < eps) {
		return ENext;
	}
	else {
		return iterationMethod(e * sin(ENext) + M, ENext, M);
	}
}
double halfDivisionMethod(double a, double b, double M) {
	if (fabs(b - a) < 2 * eps || fabs(((a + b) / 2) - e * sin((a + b / 2) - M) < eps) < eps) {
		return b + a / 2;
	}
	else if ((a - e * sin(a) - M) * (((a + b) / 2) - e * sin((a + b) / 2) - M) < 0) {
		return halfDivisionMethod(a, (a + b) / 2, M);
	}
	else {
		return halfDivisionMethod((a + b) / 2, b, M);
	}
}
double goldenRatioMethod(double a, double b, double X, double M) {
	if (fabs(b - a) < 2 * eps || fabs(((a + b) / 2) - e * sin((a + b) / 2) - M) < eps) {
		return a + (b - a) / X;
	}
	else if ((a - e * sin(a) - M) * (((a + b) / 2) - e * sin((a + b) / 2) - M) < 0) {
		return goldenRatioMethod(a, a + (b - a) / X, X, M);
	}
	else {
		return goldenRatioMethod(a + (b - a) / X, b, X, M);
	}
}
double newtonMethod(double dif, double E, double M) {
	if (dif < eps) {
		return E - dif;
	}
	else {
		return newtonMethod((E - e * sin(E) - M) / ((eps - e * sin(E + eps) - e * sin(E)) / eps),
			E - (E - e * sin(E) - M) / ((eps - e * sin(E + eps) - e * sin(E)) / eps), M);
	}
}

int main() {
	double p = a * (1 - e * e);
	double n = nFind();
	int time = timeFind();
	for (int i = 0; i < time; i += 60) {
		double teta = excentricToTrue(
			iterationMethod(e * sin(i * 2 * PI / time) + (i * 2 * PI / time), i * 2 * PI / time, i * 2 * PI / time)
		);
		std::cout << std::fixed << radFind(teta, p) << '\t' << speedRad(p, teta) << '\t' << speedN(p, teta) << '\t' << fullSpeed(p, teta) << std::endl;
	}

	std::cout << std::endl << timeFind() / 60 << std::endl;
	return 0;
}