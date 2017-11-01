#include "model_mitchschaeff.h"
#include "../configuration.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
//omp//#include <omp.h>

Model_MitchSchaeff::Model_MitchSchaeff(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix, double VmThresh, double dVmdtThresh, const double MitchSchaeff_Params[], double VmMin, double VmMax)
						:Model(grid,configFileName,keyPrefix,2,1,0,VmThresh,dVmdtThresh,VmMin,VmMax),u_gate(p[0]),
						 u_min(p[1]),u_max(p[2]),tau_in(p[3]),tau_out(p[4]),tau_open(p[5]),tau_close(p[6]),u_range(0.0),i_u_range(0.0),i_tau_in(0.0),i_tau_out(0.0),i_tau_open(0.0),i_tau_close(0.0)
{
	u_gate = 0.13;
	u_min = -70.0;
	u_max = 30.0;
	tau_in = 0.3;
	tau_out = 6.0;
	tau_open = 150.0;
	tau_close = 120.0;
	
	u_range = u_max - u_min;
	i_u_range = 1.0 / u_range;
	i_tau_in = 1.0 / tau_in;
	i_tau_out = 1.0 / tau_out;
	i_tau_open = 1.0 / tau_open;
	i_tau_close = 1.0 / tau_close;
	
	if(MitchSchaeff_Params != NULL)
		for(int i=0;i<7;i++)
			p[i] = MitchSchaeff_Params[i];
	
	restingState[0] = 0.0;
	restingState[0] = u_min + restingState[0]*u_range;
	restingState[1] = 1.0;

	excitedState[0] = 0.9;
	excitedState[0] = u_min + excitedState[0]*u_range;
	excitedState[1] = 0.9;
	
}


void Model_MitchSchaeff::readModelParams(){
	std::stringstream paramName;
	for(int i=0;i<7;i++){
		paramName.str("");
		paramName << "Param_" << i;
		cfg.readInto(p[i], paramName.str());
	}
	u_range = u_max - u_min;
	i_u_range = 1.0 / u_range;
	i_tau_in = 1.0 / tau_in;
	i_tau_out = 1.0 / tau_out;
	i_tau_open = 1.0 / tau_open;
	i_tau_close = 1.0 / tau_close;
}

inline void Model_MitchSchaeff::getf_u(const double &u, double &dVmdt, const double &h) const{
	double J_in, J_out;
		
	J_in = -h*u*u*(1.0-u)*i_tau_in; 
	J_out = u*i_tau_out;
		
	dVmdt = -(J_in+J_out);
}

inline void Model_MitchSchaeff::getab_HH(const double &u, double &a_h, double &b_h) const{	
	a_h = (u < u_gate) ? -i_tau_open : -i_tau_close;
	b_h = (u < u_gate) ? i_tau_open : 0.0;
}

void Model_MitchSchaeff::f(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double *restrict f_vars, double t, int n) const{
	double a_h, b_h;
	
	double u = (Vm[0]-u_min) * i_u_range;
	getf_u(u, f_Vm[0], vars[0]);
	f_Vm[0] *= u_range;
	
	getab_HH(Vm[0], a_h, b_h);
	f_vars[0] = a_h*vars[0]+b_h;
}

void Model_MitchSchaeff::ab(const double *restrict Vm, const double *restrict vars, double *restrict a_Vm, double *restrict b_Vm, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	double u = (Vm[0]-u_min) * i_u_range;
	getf_u(u, b_Vm[0], vars[0]);
	a_Vm[0] = 0;
	b_Vm[0] *= u_range;
	
	getab_HH(u, a_vars[0], b_vars[0]);
}

void Model_MitchSchaeff::f_Vm(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double t, int n) const{
	double u = (Vm[0]-u_min) * i_u_range;
	getf_u(u, f_Vm[0], vars[0]);
	f_Vm[0] *= u_range;
}

void Model_MitchSchaeff::ab_HH(const double *restrict Vm, const double *restrict vars, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	double u = (Vm[0]-u_min) * i_u_range;
	getab_HH(u, a_vars[0], b_vars[0]);
}

void Model_MitchSchaeff::f_NL(const double *restrict Vm, const double *restrict vars, double *restrict f_vars, double t, int n) const{
}

void Model_MitchSchaeff::f_vec(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double *restrict f_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		double a_h, b_h;
	
		double u = (Vm[i]-u_min) * i_u_range;
			
		getf_u(u, f_Vm[i], vars[i]);
		f_Vm[i] *= u_range;
	
		getab_HH(Vm[i], a_h, b_h);
		f_vars[i] = a_h*vars[i]+b_h;
	}
}

void Model_MitchSchaeff::ab_vec(const double *restrict Vm, const double *restrict vars, double *restrict a_Vm, double *restrict b_Vm, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		double u = (Vm[i]-u_min) * i_u_range;
		getf_u(u, b_Vm[i], vars[i]);
		a_Vm[i] = 0;
		b_Vm[i] *= u_range;
	
		getab_HH(u, a_vars[i], b_vars[i]);
	}
}

void Model_MitchSchaeff::f_Vm_vec(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		double u = (Vm[i]-u_min) * i_u_range;
		getf_u(u, f_Vm[i], vars[i]);
		f_Vm[i] *= u_range;
	}
}

void Model_MitchSchaeff::ab_HH_vec(const double *restrict Vm, const double *restrict vars, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		double u = (Vm[i]-u_min) * i_u_range;
		getab_HH(u, a_vars[i], b_vars[i]);
	}
}

void Model_MitchSchaeff::f_NL_vec(const double *restrict Vm, const double *restrict vars, double *restrict f_vars, double t, int n) const{
}
