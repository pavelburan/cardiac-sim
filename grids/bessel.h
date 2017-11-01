#ifndef BESSEL_H
#define BESSEL_H

#define BIG  1e30

double bessel_J0(double x);
double bessel_J1(double x);
double bessel_J(int n, double x);

double bessel_Y0(double x);
double bessel_Y1(double x);
double bessel_Y(int n, double x);

double bessel_K0(double x);
double bessel_K1(double x);
double bessel_K(int n, double x);

double bessel_I0(double x);
double bessel_I1(double x);
double bessel_I(int n, double x);

#endif //BESSEL_H
