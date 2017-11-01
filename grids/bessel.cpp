#include "bessel.h"
#include <math.h>

double bessel_J0(double x) {
	const double
	p1=1.0, p2=-0.1098628627E-2, p3=0.2734510407E-4,
	p4=-0.2073370639E-5, p5= 0.2093887211E-6,
	q1=-0.1562499995E-1, q2= 0.1430488765E-3, q3=-0.6911147651E-5,
	q4= 0.7621095161E-6, q5=-0.9349451520E-7,
	r1= 57568490574.0, r2=-13362590354.0, r3=651619640.7,
	r4=-11214424.18, r5= 77392.33017, r6=-184.9052456,
	s1= 57568490411.0, s2=1029532985.0, s3=9494680.718,
	s4= 59272.64853, s5=267.8532712, s6=1.0;
	
	double ax, fr, fs, z, fp, fq, xx, y, tmp;
	
	if (x==0.0) return 1.0;
	ax = fabs(x);
	if (ax < 8.0) {
		y = x*x;
		fr = r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))));
		fs = s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6))));
		tmp = fr/fs;
	}
	else {
		z = 8./ax;
		y = z*z;
		xx = ax-0.785398164;
		fp = p1+y*(p2+y*(p3+y*(p4+y*p5)));
		fq = q1+y*(q2+y*(q3+y*(q4+y*q5)));
		tmp = sqrt(0.636619772/ax)*(fp*cos(xx)-z*fq*sin(xx));
	}
	return tmp;
}

double bessel_J1(double x) {
	const double  
	p1=1.0, p2=0.183105E-2, p3=-0.3516396496E-4, p4=0.2457520174E-5,
	p5=-0.240337019E-6,  p6=0.636619772,
	q1= 0.04687499995, q2=-0.2002690873E-3, q3=0.8449199096E-5,
	q4=-0.88228987E-6, q5= 0.105787412E-6,
	r1= 72362614232.0, r2=-7895059235.0, r3=242396853.1,
	r4=-2972611.439,   r5=15704.48260,  r6=-30.16036606,
	s1=144725228442.0, s2=2300535178.0, s3=18583304.74,
	s4=99447.43394,    s5=376.9991397,  s6=1.0;
	
	double ax, fr, fs, y, z, fp, fq, xx, tmp;
	
	ax = fabs(x);
	if (ax < 8.0) {
		y = x*x;
		fr = r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))));
		fs = s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6))));
		tmp = x*(fr/fs);
	}
	else {
		z = 8.0/ax;
		y = z*z;
		xx = ax-2.35619491;
		fp = p1+y*(p2+y*(p3+y*(p4+y*p5)));
		fq = q1+y*(q2+y*(q3+y*(q4+y*q5)));
		tmp = sqrt(p6/ax)*(cos(xx)*fp-z*sin(xx)*fq);
		if (x < 0.0) tmp = -tmp;

	}
	return tmp;
}

double bessel_J(int n, double x) {
	const int IACC = 40; 
	const double BIGNO = 1e10,  BIGNI = 1e-10;
	
	double tox, bjm, bj, bjP, sum, tmp;
	int j, jsum, m;
	
	if (n == 0) return bessel_J0(x);
	if (n == 1) return bessel_J1(x);
	if (x == 0.0) return 0.0;
	
	tox = 2.0/x;
	if (x > 1.0*n) {
		bjm = bessel_J0(x);
		bj  = bessel_J1(x);
		for (j=1; j<n; j++) {
			bjP = j*tox*bj-bjm;
			bjm = bj;
			bj  = bjP;
		}
		return bj;
	}
	else {
		m = (int) (2*((n+floor(sqrt(1.0*(IACC*n))))/2));
		tmp = 0.0;
		jsum = 0;
		sum = 0.0;
		bjP = 0.0;
		bj  = 1.0;
		for (j=m; j>0; j--) {
			bjm = j*tox*bj-bjP;
			bjP = bj;
			bj  = bjm;
			if (fabs(bj) > BIGNO) {
            bj  = bj*BIGNI;
            bjP = bjP*BIGNI;
            tmp = tmp*BIGNI;
            sum = sum*BIGNI;
			}
			if (jsum != 0)  sum += bj;
			jsum = 1-jsum;
			if (j == n)  tmp = bjP;
		}
		sum = 2.0*sum-bj;
		return (tmp/sum);
	}
}

double bessel_Y0(double x) {
	const double
	p1= 1.0, p2=-0.1098628627E-2, p3=0.2734510407E-4,
	p4=-0.2073370639E-5, p5= 0.2093887211E-6,
	q1=-0.1562499995E-1, q2= 0.1430488765E-3, q3=-0.6911147651E-5,
	q4= 0.7621095161E-6, q5=-0.9349451520E-7,
	r1=-2957821389.0, r2=7062834065.0, r3=-512359803.6,
	r4= 10879881.29,  r5=-86327.92757, r6=228.4622733,
	s1= 40076544269.0, s2=745249964.8, s3=7189466.438,
	s4= 47447.26470,   s5=226.1030244, s6=1.0;
	
	double fs, fr, z, fp, fq, xx, y;
	
	if (x == 0.0) return -1e30;
	if (x < 8.0) {
		y = x*x;
		fr = r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))));
		fs = s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6))));
		return (fr/fs+0.636619772*bessel_J0(x)*log(x));
	}
	else {
		z = 8.0/x;
		y = z*z;
		xx = x-0.785398164;
		fp = p1+y*(p2+y*(p3+y*(p4+y*p5)));
		fq = q1+y*(q2+y*(q3+y*(q4+y*q5)));
		return (sqrt(0.636619772/x)*(fp*sin(xx)+z*fq*cos(xx)));
	}
}

double bessel_Y1(double x) {
	const double
	p1= 1.0, p2=0.183105E-2, p3=-0.3516396496E-4,
	p4= 0.2457520174E-5, p5=-0.240337019E-6,
	q1= 0.04687499995, q2=-0.2002690873E-3, q3=0.8449199096E-5,
	q4=-0.88228987E-6, q5= 0.105787412E-6,
	r1=-0.4900604943E13, r2= 0.1275274390E13, r3=-0.5153438139E11,
	r4= 0.7349264551E9,  r5=-0.4237922726E7,  r6= 0.8511937935E4,
	s1= 0.2499580570E14, s2= 0.4244419664E12, s3= 0.3733650367E10,
	s4= 0.2245904002E8,  s5= 0.1020426050E6,  s6= 0.3549632885E3, s7=1.0;
	
	double fr, fs, z, fp, fq, xx, y;
	
	if (x == 0.0) return -1e30;
	if (x < 8.0) {
		y = x*x;
		fr = r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))));
		fs = s1+y*(s2+y*(s3+y*(s4+y*(s5+y*(s6+y*s7)))));
		return (x*(fr/fs)+0.636619772*(bessel_J1(x)*log(x)-1.0/x));
	}
	else {
		z = 8./x;
		y = z*z;
		xx = x-2.356194491;
		fp = p1+y*(p2+y*(p3+y*(p4+y*p5)));
		fq = q1+y*(q2+y*(q3+y*(q4+y*q5)));
		return (sqrt(0.636619772/x)*(sin(xx)*fp+z*cos(xx)*fq));
	}
}

double bessel_Y(int n, double x) {
	double tox, by, bym, byp;
	int j;
	
	if (n == 0) return bessel_Y0(x);
	if (n == 1) return bessel_Y1(x);
	if (x == 0.0) return -1e30;
	tox = 2.0/x;
	by  = bessel_Y1(x);
	bym = bessel_Y0(x);
	for (j=1; j<n; j++) {
		byp = j*tox*by-bym;
		bym = by;
		by  = byp;
	};
	return by;
}

double bessel_I0(double x)  {
	const double
	p1=1.0, p2=3.5156229, p3=3.0899424, p4=1.2067429,
	p5=0.2659732, p6=0.360768e-1, p7=0.45813e-2,
	q1=0.39894228, q2=0.1328592e-1, q3=0.225319e-2, q4=-0.157565e-2,
	q5=0.916281e-2, q6=-0.2057706e-1, q7=0.2635537e-1,
	q8=-0.1647633e-1, q9=0.392377e-2;
	
	double ax, bx, y;
	
	if (fabs(x) < 3.75) {
		y=pow((x/3.75),2);
		return (p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
	}
	else {
		ax=fabs(x);
		y=3.75/ax;
		bx=exp(ax)/sqrt(ax);
		ax=q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9)))))));
		return ax*bx;
	}
}

double bessel_I1(double x)  {
	const double
	p1=0.5, p2=0.87890594, p3=0.51498869, p4=0.15084934,
	p5=0.2658733e-1, p6=0.301532e-2, p7=0.32411e-3,
	q1=0.39894228, q2=-0.3988024e-1, q3=-0.362018e-2,
	q4=0.163801e-2,q5=-0.1031555e-1, q6= 0.2282967e-1,
	q7=-0.2895312e-1, q8=0.1787654e-1, q9=-0.420059e-2;
	
	double ax, bx, y;
	
	if (fabs(x) < 3.75) {
		y=pow((x/3.75),2);
		return (x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))));
	}
	else {
		ax=fabs(x);
		y=3.75/ax;
		bx=exp(ax)/sqrt(ax);
		ax=q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9)))))));
		return ax*bx;
	}
}

double bessel_I(int n, double x) {
	const int IACC = 40; 
	const double BIGNO = 1e10, BIGNI = 1e-10;
	double tox, bim, bi, bip, bsi;
	int j, m;
	
	if (n==0)  return (bessel_I0(x));
	if (n==1)  return (bessel_I1(x));
	if (x==0.0) return 0.0;
	
	tox = 2.0/x;
	bip = 0.0;
	bi  = 1.0;
	bsi = 0.0;
	m = (int) (2*((n+floor(sqrt(IACC*n)))));
	for (j = m; j>0; j--) {
		bim = bip+j*tox*bi;
		bip = bi;
		bi  = bim;
		if (fabs(bi) > BIGNO) {
			bi  = bi*BIGNI;
			bip = bip*BIGNI;
			bsi = bsi*BIGNI;
		}
		if (j==n)  bsi = bip;
	}
	return (bsi*bessel_I0(x)/bi);
}

double bessel_K0(double x)  {
	const double
	p1=-0.57721566, p2= 0.42278420, p3=0.23069756, p4=0.3488590e-1,
	p5= 0.262698e-2, p6=0.10750e-3, p7=0.74e-5,
	q1= 1.25331414, q2=-0.7832358e-1, q3=0.2189568e-1, q4=-0.1062446e-1,
	q5= 0.587872e-2, q6=-0.251540e-2, q7=0.53208e-3;
	
	double tmp, ax, y;
	
	if (x == 0.0) return BIG;  //arbitrary big value
	if (x <= 2.0) {
		y=x*x/4.0;
		ax=-log(x/2.0)*bessel_I0(x);
		tmp = ax+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
		return tmp;
	}
	else {
		y=2.0/x;
		ax=exp(-x)/sqrt(x);
		tmp = ax*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
		return tmp;
	}
}

double bessel_K1(double x)  {
	const double
	p1=1.0, p2=0.15443144, p3=-0.67278579, p4=-0.18156897,
	p5=-0.1919402e-1, p6=-0.110404e-2, p7=-0.4686e-4,
	q1=1.25331414, q2=0.23498619, q3=-0.3655620e-1, q4=0.1504268e-1,
	q5=-0.780353e-2, q6=0.325614e-2, q7=-0.68245e-3;
	
	double ax, y;
	
	if (x == 0.0) return BIG;
	if (x <= 2.0) {
		y=x*x/4.0;
		ax=log(x/2.0)*bessel_I1(x);
		return (ax+(1.0/x)*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))));
	}
	else {
		y=2.0/x;
		ax=exp(-x)/sqrt(x);
		return (ax*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7)))))));
	}
}

double bessel_K(int n, double x) {
	double tox, bk, bkm, bkp;
	int j;
	if (n == 0) return bessel_K0(x);
	if (n == 1) return bessel_K1(x);
	if (x == 0.0) return BIG;  //arbitrary big value
	tox = 2.0/x;
	bk  = bessel_K1(x);
	bkm = bessel_K0(x);
	for (j=1; j<n; j++) {
		bkp = bkm + j*tox*bk;
		bkm = bk;
		bk  = bkp;
	}
	return bk;
}
