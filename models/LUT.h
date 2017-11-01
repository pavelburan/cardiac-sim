#ifndef LUT_H
#define LUT_H
#include <math.h>
#include <stdlib.h>

//#define LUTERROR
//#define LUTSAVE
#define LUTORIGINAL

class LUT_function
{
public:
	//Konstruktor
	LUT_function(double (*f)(double), double xmin, double xmax, int N):xmin(xmin),xmax(xmax),dxinv(0.0),lut_array(new double[N]),lut(NULL)
	{
		double dx = (xmax-xmin)/(N-1);
		dxinv = 1/dx;
		lut = lut_array - int(xmin*dxinv+0.5);
		for (int i=0; i<N; i++)
		{ 
			lut_array[i] = f(xmin+i*dx);
		}
		forigin = f;
	}
	//Destruktor
	~LUT_function(){
		delete lut_array;
	}
		
	double operator()(double x) const
	{ 
		#ifdef LUTERROR
		if(x < xmin || x > xmax){
			std::cerr<<x<<" < "<<xmin<<" oder "<<x<<" > "<<xmax<<std::endl;
			exit(1);
		}
		#endif
		
		#ifdef LUTORIGINAL
		return forigin(x);
		#elif LUTSAVE
		return (x>xmax && x<xmin) ? forigin(x) : lut[int(x*dxinv+0.5)];	
		#else
		return lut[int(x*dxinv+0.5)];
		#endif
	}
	
private:
	double xmin;
	double xmax;
	double (*forigin)(double);
	double dxinv;
	double *lut_array;
	double *lut; 
};

#define LUT_EXP_N 20000000
#define LUT_EXP_XMIN -700.0
#define LUT_EXP_XMAX 700.0
#define LUT_EXP_DX ((LUT_EXP_XMAX-LUT_EXP_XMIN)/LUT_EXP_N)
#define LUT_EXP_DXINV (1/LUT_EXP_DX)
class LUT_exp
{
public:
	//Konstruktor
	LUT_exp():lut_array(new double[LUT_EXP_N]),lut(NULL)
	{
		lut = lut_array - int(LUT_EXP_XMIN*LUT_EXP_DXINV+0.5);
		for (int i=0; i<LUT_EXP_N; i++)
		{ 
			lut_array[i] = exp(LUT_EXP_XMIN+i*LUT_EXP_DX);
		}
	}
	//Destruktor
	~LUT_exp(){
		delete lut_array;
	}
		
	double operator()(double x) const
	{ 
		#ifdef LUTERROR
		if(x < LUT_EXP_XMIN || x > LUT_EXP_XMAX){
			std::cerr<<x<<" < "<<LUT_EXP_XMIN<<" oder "<<x<<" > "<<LUT_EXP_XMAX<<std::endl;
			exit(1);
		}
		#endif
		
		#ifdef LUTORIGINAL
		return exp(x);
		#elif LUTSAVE
		if(x < LUT_EXP_XMIN)
			return 0.0;
		else if( x < LUT_EXP_XMAX)
			return lut[int(x*LUT_EXP_DXINV+0.5)];
		else
			return std::numeric_limits<double>::infinity();
		#else
		return lut[int(x*LUT_EXP_DXINV+0.5)];
		#endif
	} 
private:
	double *lut_array;
	double *lut; 
};

#define LUT_TANH_N 10000000
#define LUT_TANH_XMIN -15.0
#define LUT_TANH_XMAX 15.0
#define LUT_TANH_DX ((LUT_TANH_XMAX-LUT_TANH_XMIN)/LUT_TANH_N)
#define LUT_TANH_DXINV (1/LUT_TANH_DX)
class LUT_tanh
{
public:
	//Konstruktor
	LUT_tanh():lut_array(new double[LUT_TANH_N]),lut(NULL)
	{
		lut = lut_array - int(LUT_TANH_XMIN*LUT_TANH_DXINV+0.5);
		for (int i=0; i<LUT_TANH_N; i++)
		{ 
			lut_array[i] = tanh(LUT_TANH_XMIN+i*LUT_TANH_DX);
		}
	}
	//Destruktor
	~LUT_tanh(){
		delete lut_array;
	}
		
	double operator()(double x) const
	{ 
		#ifdef LUTERROR
		if(x < LUT_TANH_XMIN || x > LUT_TANH_XMAX){
			std::cerr<<x<<" < "<<LUT_TANH_XMIN<<" oder "<<x<<" > "<<LUT_TANH_XMAX<<std::endl;
			exit(1);
		}
		#endif
		
		#ifdef LUTORIGINAL
		return tanh(x);
		#elif LUTSAVE
		if(x < LUT_TANH_XMIN)
			return -1.0;
		else if( x < LUT_TANH_XMAX)
			return lut[int(x*LUT_TANH_DXINV+0.5)];
		else
			return 1.0;
		#else
		return lut[int(x*LUT_TANH_DXINV+0.5)];
		#endif
	} 
private:
	double *lut_array;
	double *lut; 
};

#endif //LUT_H

