#include "model_luorudy.h"
#include "../configuration.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
//omp//#include <omp.h>

Model_LuoRudy::Model_LuoRudy(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix, double VmThresh, double dVmdtThresh, const double LuoRudy_Params[], double VmMin, double VmMax)
											:Model(grid,configFileName,keyPrefix,8,6,1,VmThresh,dVmdtThresh,VmMin,VmMax),E_Na(p[0]),K_0(p[1]),E_K(p[2]),G_K(p[3]),E_K1(p[4]),G_K1(p[5]),E_Kp(p[6]),G_Kp(p[7]),G_Na(p[8]),G_Si(p[9]),c_h(p[10]),c_j(p[11]),c_m(p[12]),c_d(p[13]),c_f(p[14]),c_X(p[15])
{
	E_Na = 54.4;
	K_0 = 5.4;
	E_K = -77.0;
	G_K = 0.282;
	E_K1 = -87.25;
	G_K1 = 0.6047;
	E_Kp = -87.25;
	G_Kp = 0.0183;
	G_Na = 23.0;
	G_Si = 0.09;
	c_h = 1.0;
	c_j = 1.0;
	c_m = 1.0;
	c_d = 1.0;
	c_f = 1.0;
	c_X = 1.0;
	
	if(LuoRudy_Params != NULL)
		for(int i=0;i<16;i++)
			p[i] = (LuoRudy_Params[i] >= 0.0) ? LuoRudy_Params[i] : p[i];
	
	restingState[0] = -84.3801107371;
	restingState[1] = 0.982660523699656;
	restingState[2] = 0.989108212766685;
	restingState[3] = 0.00171338077730188;
	restingState[4] = 0.00302126301779861;
	restingState[5] = 0.999967936476325;
	restingState[6] = 0.0417603108167287;
	restingState[7] = 0.00017948816388306;
	
	refractoryState[0] = -27.316947957774961;
	refractoryState[1] = 0.0;
	refractoryState[2] = 0.0;
	refractoryState[3] = .88539327531610046;
	refractoryState[4] = 0.64108519691731369;
	refractoryState[5] = 0.48045551606739073;
	refractoryState[6] = 0.46867351095665494;
	refractoryState[7] = 0.004488599864854453;
	
	excitedState[0] = 0.0;
	excitedState[1] = 0.013466785134096645*0.5 + restingState[1]*0.5;
	excitedState[2] = 0.0080369236341538974*0.5 + restingState[2]*0.5;
	excitedState[3] = 0.16770627638045082*0.5 + restingState[3]*0.5;
	excitedState[4] = 0.22938134304378971*0.5 + restingState[4]*0.5;
	excitedState[5] = 0.63725773300847988*0.5 + restingState[5]*0.5;
	excitedState[6] = 0.39515751930580828*0.5 + restingState[6]*0.5;
	excitedState[7] = 0.003234453122370621*0.5 + restingState[7]*0.5;

}

void Model_LuoRudy::readModelParams(){
	Model::readModelParams();
	std::stringstream paramName;
	for(int i=0;i<16;i++){
		paramName.str("");
		paramName << "Param_" << i;
		cfg.readInto(p[i], paramName.str());
	}
}

__attribute__((always_inline)) void Model_LuoRudy::getf_Vm(const double &Vm, double &dVmdt, const double &h, const double &j, const double &m, const double &d, const double &f, const double &X, const double &Ca) const{
	double alpha, beta;
	double I_Na, I_si, I_K, I_K1, I_Kp, I_b;
	double E_si, X_i, K_p;
		
	I_Na = G_Na*m*m*m*h*j*(Vm-E_Na);
		
	E_si = 7.7 - 13.0287*log(Ca);
	I_si = G_Si*d*f*(Vm - E_si);
		
	X_i = (Vm >= -100.0) ? 2.837*(exp(0.04*(Vm+77.0)) - 1.0) / ((Vm+77.0)*exp(0.04*(Vm+35.0))) : 1.0;
	I_K = G_K*X*X_i*(Vm-E_K);
		
	alpha = 1.02/(1.0 + exp(0.2385*(Vm-E_K1-59.215)));
	beta = (0.49124*exp(0.08032*(Vm-E_K1+5.476)) + exp(0.06175*(Vm-E_K1-594.31))) / (1.0 + exp(-0.5143*(Vm-E_K1+4.753)));
	I_K1 = G_K1*alpha/(alpha+beta)*(Vm-E_K1);
		
	K_p = 1.0/(1.0 + exp((7.488-Vm)/5.98));
	I_Kp = 0.0183*K_p*(Vm-E_Kp);
		
	I_b = 0.03921*(Vm+59.87);
		 
	dVmdt = -(I_Na + I_si + I_K + I_K1 + I_Kp + I_b);
}

__attribute__((always_inline)) void Model_LuoRudy::getab_HH(const double &Vm, double &a_h, double &b_h, double &a_j, double &b_j, double &a_m, double &b_m, double &a_d, double &b_d, double &a_f, double &b_f, double &a_X, double &b_X) const{
	double alpha, beta;
		
	alpha = (Vm >= -40.0) ? 0.0 : 0.135*exp(-(Vm+80.0)/6.8); 
	beta = (Vm >= -40.0) ? 1.0/(0.13*(1.0+exp(-(Vm+10.66)/11.1))) : 3.56*exp(0.079*Vm) + 3.1e5*exp(0.35*Vm);
	a_h = -(alpha+beta)*c_h;
	b_h = alpha*c_h;	
		
	alpha = (Vm >= -40.00) ? 0.0 : (Vm+37.78) * (-1.2714e5*exp(0.2444*Vm) - 3.474e-5*exp(-0.04391*Vm)) / (1.0 + exp(0.311*(Vm+79.23)));
	beta = (Vm >= -39.826) ? 0.3*exp(-2.535e-7*Vm) / (1.0 + exp(-0.1*(Vm+32.0))) : 0.1212*exp(-0.01052*Vm) / (1.0 + exp(-0.1378*(Vm+40.14)));
	a_j = -(alpha+beta)*c_j;
	b_j = alpha*c_j;
		
	alpha = 0.32*(Vm + 47.13) / (1.0 - exp(-0.1*(Vm+47.13)));
	beta = 0.08*exp(-Vm/11.0);
	a_m = -(alpha+beta)*c_m;
	b_m = alpha*c_m;
		
	alpha = 0.095*exp(-0.01*(Vm-5.0)) / (1.0 + exp(-0.072*(Vm-5.0)));
	beta = 0.07*exp(-0.017*(Vm+44.0)) / (1.0 + exp(0.05*(Vm+44.0)));
	a_d = -(alpha+beta)*c_d;
	b_d = alpha*c_d;	
	
	alpha = 0.012*exp(-0.008*(Vm+28.0)) / (1.0 + exp(0.15*(Vm+28.0)));
	beta = 0.0065*exp(-0.02*(Vm+30.0)) / (1.0 + exp(-0.2*(Vm+30.0)));
	a_f = -(alpha+beta)*c_f;
	b_f = alpha*c_f;	
		
	alpha = 0.0005*exp(0.083*(Vm+50.0)) / (1.0 + exp(0.057*(Vm+50.0)));
	beta = 0.0013*exp(-0.06*(Vm+20.0)) / (1.0 + exp(-0.04*(Vm+20.0)));
	a_X = -(alpha+beta)*c_X;
	b_X = alpha*c_X;
}

__attribute__((always_inline)) void Model_LuoRudy::getf_NL(const double &Vm, const double &d, const double &f, const double &Ca, double  &dCadt) const{
	double I_si, E_si;
		
	E_si = 7.7 - 13.0287*log(Ca);
	I_si = G_Si*d*f*(Vm - E_si);
	
	dCadt = -1e-4*I_si + 0.07*(1e-4 - Ca);
}

void Model_LuoRudy::f(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double *restrict f_vars, double t, int n) const{
	double a_h, a_j, a_m, a_d, a_f, a_X;
	double b_h, b_j, b_m, b_d, b_f, b_X;
	getf_Vm(Vm[0], f_Vm[0], vars[0], vars[n], vars[2*n], vars[3*n], vars[4*n], vars[5*n], vars[6*n]);
	getab_HH(Vm[0], a_h, b_h, a_j, b_j, a_m, b_m, a_d, b_d, a_f, b_f, a_X, b_X);
	f_vars[0] = a_h*vars[0]+b_h;
	f_vars[n] = a_j*vars[n]+b_j;
	f_vars[2*n] = a_m*vars[2*n]+b_m;
	f_vars[3*n] = a_d*vars[3*n]+b_d;
	f_vars[4*n] = a_f*vars[4*n]+b_f;
	f_vars[5*n] = a_X*vars[5*n]+b_X;
	getf_NL(Vm[0], vars[3*n], vars[4*n], vars[6*n], f_vars[6*n]);
}

void Model_LuoRudy::ab(const double *restrict Vm, const double *restrict vars, double *restrict a_Vm, double *restrict b_Vm, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	getf_Vm(Vm[0], b_Vm[0], vars[0], vars[n], vars[2*n], vars[3*n], vars[4*n], vars[5*n], vars[6*n]);
	a_Vm[0] = 0;
	getab_HH(Vm[0], a_vars[0], b_vars[0], a_vars[n], b_vars[n], a_vars[2*n], b_vars[2*n], a_vars[3*n], b_vars[3*n], a_vars[4*n], b_vars[4*n], a_vars[5*n], b_vars[5*n]);
	getf_NL(Vm[0], vars[3*n], vars[4*n], vars[6*n], b_vars[6*n]);
	a_vars[6*n] = 0;
}

void Model_LuoRudy::f_Vm(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double t, int n) const{
	getf_Vm(Vm[0], f_Vm[0], vars[0], vars[n], vars[2*n], vars[3*n], vars[4*n], vars[5*n], vars[6*n]);
}

void Model_LuoRudy::ab_HH(const double *restrict Vm, const double *restrict vars, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	getab_HH(Vm[0], a_vars[0], b_vars[0], a_vars[n], b_vars[n], a_vars[2*n], b_vars[2*n], a_vars[3*n], b_vars[3*n], a_vars[4*n], b_vars[4*n], a_vars[5*n], b_vars[5*n]);
}

void Model_LuoRudy::f_NL(const double *restrict Vm, const double *restrict vars, double *f_vars, double t, int n) const{
	getf_NL(Vm[0], vars[3*n], vars[4*n], vars[6*n], f_vars[6*n]);
}

void Model_LuoRudy::f_vec(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double *restrict f_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		double a_h, a_j, a_m, a_d, a_f, a_X;
		double b_h, b_j, b_m, b_d, b_f, b_X;
		getf_Vm(Vm[i], f_Vm[i], vars[i], vars[n+i], vars[2*n+i], vars[3*n+i], vars[4*n+i], vars[5*n+i], vars[6*n+i]);
		getab_HH(Vm[i], a_h, b_h, a_j, b_j, a_m, b_m, a_d, b_d, a_f, b_f, a_X, b_X);
		f_vars[i] = a_h*vars[i]+b_h;
		f_vars[n+i] = a_j*vars[n+i]+b_j;
		f_vars[2*n+i] = a_m*vars[2*n+i]+b_m;
		f_vars[3*n+i] = a_d*vars[3*n+i]+b_d;
		f_vars[4*n+i] = a_f*vars[4*n+i]+b_f;
		f_vars[5*n+i] = a_X*vars[5*n+i]+b_X;
		getf_NL(Vm[i], vars[3*n+i], vars[4*n+i], vars[6*n+i], f_vars[6*n+i]);
	}
}

void Model_LuoRudy::ab_vec(const double *restrict Vm, const double *restrict vars, double *restrict a_Vm, double *restrict b_Vm, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		getf_Vm(Vm[i], b_Vm[i], vars[i], vars[n+i], vars[2*n+i], vars[3*n+i], vars[4*n+i], vars[5*n+i], vars[6*n+i]);
		a_Vm[i] = 0;
		getab_HH(Vm[i], a_vars[i], b_vars[i], a_vars[n+i], b_vars[n+i], a_vars[2*n+i], b_vars[2*n+i], a_vars[3*n+i], b_vars[3*n+i], a_vars[4*n+i], b_vars[4*n+i], a_vars[5*n+i], b_vars[5*n+i]);
		getf_NL(Vm[i], vars[3*n+i], vars[4*n+i], vars[6*n+i], b_vars[6*n+i]);
		a_vars[6*n+i] = 0;
	}
}

void Model_LuoRudy::f_Vm_vec(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		getf_Vm(Vm[i], f_Vm[i], vars[i], vars[n+i], vars[2*n+i], vars[3*n+i], vars[4*n+i], vars[5*n+i], vars[6*n+i]);
	}
}

void Model_LuoRudy::ab_HH_vec(const double *restrict Vm, const double *restrict vars, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		getab_HH(Vm[i], a_vars[i], b_vars[i], a_vars[n+i], b_vars[n+i], a_vars[2*n+i], b_vars[2*n+i], a_vars[3*n+i], b_vars[3*n+i], a_vars[4*n+i], b_vars[4*n+i], a_vars[5*n+i], b_vars[5*n+i]);
	}
}

void Model_LuoRudy::f_NL_vec(const double *restrict Vm, const double *restrict vars, double *f_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		getf_NL(Vm[i], vars[3*n+i], vars[4*n+i], vars[6*n+i], f_vars[6*n+i]);
	}
}
