#include "model_fentonkarma.h"
#include "../configuration.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
//omp//#include <omp.h>

Model_FentonKarma::Model_FentonKarma(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix, int modNum, double VmThresh, double dVmdtThresh, const double FentonKarma_Params[], double VmMin, double VmMax)
						:Model(grid,configFileName,keyPrefix,3,2,0,VmThresh,dVmdtThresh,VmMin,VmMax),modNum(0),
						u_c(p[0]),u_v(p[1]),u_csi(p[2]),u_min(p[3]),u_max(p[4]),tau_d(p[5]),tau_o(p[6]),tau_r(p[7]),tau_si(p[8]),tau_v_p(p[9]),tau_v1_m(p[10]),tau_v2_m(p[11]),
						tau_w_p(p[12]),tau_w_m(p[13]),k(p[14]),u_range(0.0),i_u_range(0.0),i_tau_v_p(0.0),i_tau_v1_m(0.0),i_tau_v2_m(0.0),i_tau_w_p(0.0),i_tau_w_m(0.0),i_tau_d(0.0),i_tau_o(0.0),i_tau_r(0.0),i_tau_si(0.0)
{
	setModel(modNum);
	
	if(FentonKarma_Params != NULL)
		for(int i=0;i<15;i++)
			p[i] = FentonKarma_Params[i];
	
	restingState[0] = 0.0;
	restingState[0] = restingState[0]*u_range + u_min;
	restingState[1] = 1.0;
	restingState[2] = 1.0;

	excitedState[0] = 0.9;
	excitedState[0] = excitedState[0]*u_range + u_min;
	excitedState[1] = 0.9;
	excitedState[2] = 0.9;
	
}

void Model_FentonKarma::setModel(int newModNum){
#ifdef _ERROR_
	if(newModNum < 0 || newModNum >4){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Modellnumber must be between 0 and 0!"<<std::endl;
		exit(1);
	}
#endif
	modNum = newModNum;
	if(modNum == 0){
		//EPI
		u_c = 0.13;
		u_v = 0.04;
		u_csi = 0.861;
		tau_d = 1.0 / 1.75;
		tau_o = 12.5;
		tau_r = 33.83;
		tau_si = 29.0;
		tau_v_p = 7.99;
		tau_v1_m = 9.8;
		tau_v2_m = 312.5;
		tau_w_p = 870.0;
		tau_w_m = 41.0;
		k = 10.0;
	}
	u_min = -84.0;
	u_max = 1.7;
	u_range = u_max - u_min;
	i_u_range = 1.0 / u_range;
	i_tau_v1_m = 1.0 / tau_v1_m;
	i_tau_v2_m = 1.0 / tau_v2_m ;
	i_tau_v_p = 1.0 / tau_v_p;
	i_tau_w_p = 1.0 / tau_w_p;
	i_tau_w_m = 1.0 / tau_w_m ;
	i_tau_d = 1.0 / tau_d;
	i_tau_o = 1.0 / tau_o;
	i_tau_r = 1.0 / tau_r;
	i_tau_si = 1.0 / tau_si;
}

void Model_FentonKarma::readModelParams(){
	Model::readModelParams();
	cfg.readInto(modNum, "modNum");
	setModel(modNum);
	std::stringstream paramName;
	for(int i=0;i<15;i++){
		paramName.str("");
		paramName << "Param_" << i;
		cfg.readInto(p[i], paramName.str());
	}
	u_range = u_max - u_min;
	i_u_range = 1.0 / u_range;
	i_tau_v1_m = 1.0 / tau_v1_m;
	i_tau_v2_m = 1.0 / tau_v2_m ;
	i_tau_v_p = 1.0 / tau_v_p;
	i_tau_w_p = 1.0 / tau_w_p;
	i_tau_w_m = 1.0 / tau_w_m ;
	i_tau_d = 1.0 / tau_d;
	i_tau_o = 1.0 / tau_o;
	i_tau_r = 1.0 / tau_r;
	i_tau_si = 1.0 / tau_si;
	
	restingState[0] = 0.0;
	restingState[0] = restingState[0]*u_range + u_min;
	restingState[1] = 1.0;
	restingState[2] = 1.0;

	excitedState[0] = 0.9;
	excitedState[0] = excitedState[0]*u_range + u_min;
	excitedState[1] = 0.9;
	excitedState[2] = 0.9;
}

inline void Model_FentonKarma::getf_u(const double &u, double &dVmdt, const double &v, const double &w) const{
	double J_fi, J_so, J_si;
	
	J_fi = (u >= u_c) ? -v*i_tau_d*(1.0-u)*(u-u_c) : 0.0;
	J_so = (u_c >= u) ? u*i_tau_o : i_tau_r;
	J_si = -0.5*i_tau_si*w*(1.0 + tanh(k*(u-u_csi)));
			
	dVmdt = -(J_fi + J_so + J_si);
}

inline void Model_FentonKarma::getab_HH(const double &u, double &a_v, double &b_v, double &a_w, double &b_w) const{	
	double i_tau_v_m = (u >= u_v) ? i_tau_v1_m : i_tau_v2_m; 
	a_v = (u_c >= u) ? -i_tau_v_m : -i_tau_v_p;
	b_v = (u_c >= u) ? i_tau_v_m : 0.0;
	a_w = (u_c >= u) ? -i_tau_w_m : -i_tau_w_p;
	b_w = (u_c >= u) ? i_tau_w_m : 0.0;
}

void Model_FentonKarma::f(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double *restrict f_vars, double t, int n) const{
	double a_v, b_v, a_w, b_w;
	
	double u = (Vm[0]-u_min) * i_u_range;
	getf_u(u, f_Vm[0], vars[0], vars[n]);
	f_Vm[0] *= u_range;
	
	getab_HH(Vm[0], a_v, b_v, a_w, b_w);
	f_vars[0] = a_v*vars[0]+b_v;
	f_vars[n] = a_w*vars[n]+b_w;
}

void Model_FentonKarma::ab(const double *restrict Vm, const double *restrict vars, double *restrict a_Vm, double *restrict b_Vm, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	double u = (Vm[0]-u_min) * i_u_range;
	getf_u(u, b_Vm[0], vars[0], vars[n]);
	a_Vm[0] = 0;
	b_Vm[0] *= u_range;
	
	getab_HH(u, a_vars[0], b_vars[0], a_vars[n], b_vars[n]);
}

void Model_FentonKarma::f_Vm(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double t, int n) const{
	double u = (Vm[0]-u_min) * i_u_range;
	getf_u(u, f_Vm[0], vars[0], vars[n]);
	f_Vm[0] *= u_range;
}

void Model_FentonKarma::ab_HH(const double *restrict Vm, const double *restrict vars, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	double u = (Vm[0]-u_min) * i_u_range;
	getab_HH(u, a_vars[0], b_vars[0], a_vars[n], b_vars[n]);
}

void Model_FentonKarma::f_NL(const double *restrict Vm, const double *restrict vars, double *restrict f_vars, double t, int n) const{
}

void Model_FentonKarma::f_vec(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double *restrict f_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		double a_v, b_v, a_w, b_w;
	
		double u = (Vm[i]-u_min) * i_u_range;		
		getf_u(u, f_Vm[i], vars[i], vars[n+i]);
		f_Vm[i] *= u_range;
	
		getab_HH(Vm[i], a_v, b_v, a_w, b_w);
		f_vars[i] = a_v*vars[i]+b_v;
		f_vars[n+i] = a_w*vars[n+i]+b_w;
	}
}

void Model_FentonKarma::ab_vec(const double *restrict Vm, const double *restrict vars, double *restrict a_Vm, double *restrict b_Vm, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		double u = (Vm[i]-u_min) * i_u_range;
		getf_u(u, b_Vm[i], vars[i], vars[n+i]);
		a_Vm[i] = 0;
		b_Vm[i] *= u_range;
	
		getab_HH(u, a_vars[i], b_vars[i], a_vars[n+i], b_vars[n+i]);
	}
}

void Model_FentonKarma::f_Vm_vec(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		double u = (Vm[i]-u_min) * i_u_range;
		getf_u(u, f_Vm[i], vars[i], vars[n+i]);
		f_Vm[i] *= u_range;
	}
}

void Model_FentonKarma::ab_HH_vec(const double *restrict Vm, const double *restrict vars, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		double u = (Vm[i]-u_min) * i_u_range;
		getab_HH(u, a_vars[i], b_vars[i], a_vars[n+i], b_vars[n+i]);
	}
}

void Model_FentonKarma::f_NL_vec(const double *restrict Vm, const double *restrict vars, double *restrict f_vars, double t, int n) const{
}
