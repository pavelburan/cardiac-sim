#include "model_bocf.h"
#include "../configuration.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
//omp//#include <omp.h>

Model_BOCF::Model_BOCF(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix, int modNum, double VmThresh, double dVmdtThresh, const double BOCF_Params[], double VmMin, double VmMax)
						:Model(grid,configFileName,keyPrefix,4,3,0,VmThresh,dVmdtThresh,VmMin,VmMax),modNum(0),u_o(p[0]),
						 u_u(p[1]),u_min(p[2]),u_max(p[3]),theta_v(p[4]),theta_w(p[5]),theta_m_v(p[6]),theta_o(p[7]),tau_m_v1(p[8]),tau_m_v2(p[9]),tau_p_v(p[10]),tau_m_w1(p[11]),tau_m_w2(p[12]),
						 k_m_w(p[13]),u_m_w(p[14]),tau_p_w(p[15]),tau_fi(p[16]),tau_o1(p[17]),tau_o2(p[18]),tau_so1(p[19]),tau_so2(p[20]),k_so(p[21]),u_so(p[22]),tau_s1(p[23]),
						 tau_s2(p[24]),k_s(p[25]),u_s(p[26]),tau_si(p[27]),tau_winf(p[28]),w_infstar(p[29]),u_range(0.0),i_u_range(0.0),i_tau_p_v(0.0),i_tau_p_w(0.0),i_tau_fi(0.0),i_tau_si(0.0),i_tau_winf(0.0)
{
	setModel(modNum);
	
	if(BOCF_Params != NULL)
		for(int i=0;i<30;i++)
			p[i] = BOCF_Params[i];
	
	restingState[0] = 0.0;
	restingState[0] = restingState[0]*u_range + u_min;
	restingState[1] = 1.0;
	restingState[2] = 1.0;
	restingState[3] = 0.0;

	excitedState[0] = 1.387;
	excitedState[0] = excitedState[0]*u_range + u_min;
	excitedState[1] = 0.215;
	excitedState[2] = 0.747;
	excitedState[3] = 0.312;
}

void Model_BOCF::setModel(int newModNum){
#ifdef _ERROR_
	if(newModNum < 0 || newModNum >4){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Modellnumber must be between 0 and 4!"<<std::endl;
		exit(1);
	}
#endif
	modNum = newModNum;
	if(modNum == 0){
		//EPI
		u_o = 0.0;
		u_u = 1.55;
		theta_v = 0.3;
		theta_w = 0.13;
		theta_m_v = 0.006;
		theta_o = 0.006;
		tau_m_v1 = 60.0;
		tau_m_v2 = 1150;
		tau_p_v = 1.4506;
		tau_m_w1 = 60.0;
		tau_m_w2 = 15.0;
		k_m_w = 65.0;
		u_m_w = 0.03;
		tau_p_w = 200.0;
		tau_fi = 0.11;
		tau_o1 = 400.0;
		tau_o2 = 6.0;
		tau_so1 = 30.0181;
		tau_so2 = 0.9957;
		k_so = 2.0458;
		u_so = 0.65;
		tau_s1 = 2.7342;
		tau_s2 = 16.0;
		k_s = 2.0994;
		u_s = 0.9087;
		tau_si = 1.8875;
		tau_winf = 0.07;
		w_infstar = 0.94;
	}
	else if(modNum == 1){
		//ENDO
		u_o = 0.0;
		u_u = 1.56;
		theta_v = 0.3;
		theta_w = 0.13;
		theta_m_v = 0.2;
		theta_o = 0.006;
		tau_m_v1 = 75.0;
		tau_m_v2 = 10.0;
		tau_p_v = 1.4506;
		tau_m_w1 = 6.0;
		tau_m_w2 = 140.0;
		k_m_w = 200.0;
		u_m_w = 0.016;
		tau_p_w = 280.0;
		tau_fi = 0.1;
		tau_o1 = 470.0;
		tau_o2 = 6.0;
		tau_so1 = 40.0;
		tau_so2 = 1.2;
		k_so = 2.0;
		u_so = 0.65;
		tau_s1 = 2.7342;
		tau_s2 = 2.0;
		k_s = 2.0994;
		u_s = 0.9087;
		tau_si = 2.9013;
		tau_winf = 0.0273;
		w_infstar = 0.78;
		
	}
	else if(modNum == 2){
		//M
		u_o = 0.0;
		u_u = 1.61;
		theta_v = 0.3;
		theta_w = 0.13;
		theta_m_v = 0.1;
		theta_o = 0.005;
		tau_m_v1 = 80.0;
		tau_m_v2 = 1.4506;
		tau_p_v = 1.4506;
		tau_m_w1 = 70.0;
		tau_m_w2 = 8.0;
		k_m_w = 200.0;
		u_m_w = 0.016;
		tau_p_w = 280.0;
		tau_fi = 0.078;
		tau_o1 = 410.0;
		tau_o2 = 7.0;
		tau_so1 = 91.0;
		tau_so2 = 0.8;
		k_so = 2.1;
		u_so = 0.6;
		tau_s1 = 2.7342;
		tau_s2 = 4.0;
		k_s = 2.0994;
		u_s = 0.9087;
		tau_si = 3.3849;
		tau_winf = 0.01;
		w_infstar = 0.5;
		
	}
	else if(modNum == 3){
		//PB
		u_o = 0.0;
		u_u = 1.45;
		theta_v = 0.35;
		theta_w = 0.13;
		theta_m_v = 0.175;
		theta_o = 0.006;
		tau_m_v1 = 10.0;
		tau_m_v2 = 1150.0;
		tau_p_v = 1.4506;
		tau_m_w1 = 140.0;
		tau_m_w2 = 6.25;
		k_m_w = 65.0;
		u_m_w = 0.015;
		tau_p_w = 326.0;
		tau_fi = 0.105;
		tau_o1 = 400.0;
		tau_o2 = 6.0;
		tau_so1 = 30.0181;
		tau_so2 = 0.9957;
		k_so = 2.0458;
		u_so = 0.65;
		tau_s1 = 2.7342;
		tau_s2 = 16.0;
		k_s = 2.0994;
		u_s = 0.9087;
		tau_si = 1.8875;
		tau_winf = 0.175;
		w_infstar = 0.9;
	}
	else if(modNum == 4){
		//TNNP
		u_o = 0.0;
		u_u = 1.58;
		theta_v = 0.3;
		theta_w = 0.015;
		theta_m_v = 0.015;
		theta_o = 0.006;
		tau_m_v1 = 60.0;
		tau_m_v2 = 1150.0;
		tau_p_v = 1.4506;
		tau_m_w1 = 70.0;
		tau_m_w2 = 20.0;
		k_m_w = 65.0;
		u_m_w = 0.03;
		tau_p_w = 280.0;
		tau_fi = 0.11;
		tau_o1 = 6.0;
		tau_o2 = 6.0;
		tau_so1 = 43;
		tau_so2 = 0.2;
		k_so = 2.0;
		u_so = 0.65;
		tau_s1 = 2.7342;
		tau_s2 = 3.0;
		k_s = 2.0994;
		u_s = 0.9087;
		tau_si = 2.8723;
		tau_winf = 0.07;
		w_infstar = 0.94;
	}
	u_min = -84.0;
	u_max = 1.7;
	u_range = u_max - u_min;
	i_u_range = 1.0 / u_range;
	i_tau_p_v = 1.0/tau_p_v;
	i_tau_p_w = 1.0/tau_p_w;
	i_tau_fi = 1.0/tau_fi;
	i_tau_si = 1.0/tau_si;
	i_tau_winf = 1.0/tau_winf;
}

void Model_BOCF::readModelParams(){
	Model::readModelParams();
	cfg.readInto(modNum, "modNum");
	setModel(modNum);
	std::stringstream paramName;
	for(int i=0;i<30;i++){
		paramName.str("");
		paramName << "Param_" << i;
		cfg.readInto(p[i], paramName.str());
	}
}

inline void Model_BOCF::getf_u(const double &u, double &dVmdt, const double &v, const double &w, const double &s) const{
	double h_v, h_w, h_o;
	double i_tau_so, i_tau_o;
	double J_fi, J_so, J_si;

	h_v = (u >= theta_v)? 1. : 0.;
	h_w = (u >= theta_w)? 1. : 0.;
	h_o = (u >= theta_o)? 1. : 0.;
		
	i_tau_so = 1. / (tau_so1 + (tau_so2-tau_so1)*(1.+tanh(k_so*(u-u_so)))*.5);
	i_tau_o = 1. / (h_o*tau_o2 + (1.-h_o)*tau_o1);
		
	J_fi = -h_v*v*(u-theta_v)*(u_u-u)*i_tau_fi; 
	J_so = h_w*i_tau_so + (1.-h_w)*(u-u_o)*i_tau_o;
	J_si = -h_w*w*s*i_tau_si;
		
	dVmdt = -(J_fi+J_so+J_si);
}

inline void Model_BOCF::getab_HH(const double &u, double &a_v, double &b_v, double &a_w, double &b_w, double &a_s, double &b_s) const{
	double h_m_v, h_v, h_w;
	double i_tau_m_v, i_tau_m_w, i_tau_s;
	double v_inf, w_inf;
		
	h_m_v = (u >= theta_m_v)? 1. : 0.;
	h_v = (u >= theta_v)? 1. : 0.;
	h_w = (u >= theta_w)? 1. : 0.;
		
	i_tau_m_v = 1. / (h_m_v*tau_m_v2 + (1.-h_m_v)*tau_m_v1);
	i_tau_m_w = 1. / (tau_m_w1 + (tau_m_w2-tau_m_w1)*(1.+tanh(k_m_w*(u-u_m_w)))*.5);
	i_tau_s = 1. / (h_w*tau_s2 + (1.-h_w)*tau_s1);
		
	v_inf = (1.-h_m_v);
	w_inf = h_w*w_infstar + (1.-h_w)*(1.-u*i_tau_winf); //todo w=o
	a_v = -h_v*i_tau_p_v - (1.-h_v)*i_tau_m_v;
	b_v = (1.-h_v)*v_inf*i_tau_m_v;
	a_w = -h_w*i_tau_p_w - (1.-h_w)*i_tau_m_w;
	b_w = (1.-h_w)*w_inf*i_tau_m_w;
	a_s = -1*i_tau_s;
	b_s = (1.+tanh(k_s*(u-u_s)))*.5*i_tau_s;
}

void Model_BOCF::f(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double *restrict f_vars, double t, int n) const{
	double a_v, a_w, a_s;
	double b_v, b_w, b_s;
	
	double u = (Vm[0]-u_min) * i_u_range;
	getf_u(u, f_Vm[0], vars[0], vars[n], vars[2*n]);
	f_Vm[0] *= u_range;
	
	getab_HH(Vm[0], a_v, b_v, a_w, b_w, a_s, b_s);
	f_vars[0] = a_v*vars[0]+b_v;
	f_vars[n] = a_w*vars[n]+b_w;
	f_vars[2*n] = a_s*vars[2*n]+b_s;
}

void Model_BOCF::ab(const double *restrict Vm, const double *restrict vars, double *restrict a_Vm, double *restrict b_Vm, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	double u = (Vm[0]-u_min) * i_u_range;
	getf_u(u, b_Vm[0], vars[0], vars[n], vars[2*n]);
	a_Vm[0] = 0;
	b_Vm[0] *= u_range;
	
	getab_HH(u, a_vars[0], b_vars[0], a_vars[n], b_vars[n], a_vars[2*n], b_vars[2*n]);
}

void Model_BOCF::f_Vm(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double t, int n) const{
	double u = (Vm[0]-u_min) * i_u_range;
	getf_u(u, f_Vm[0], vars[0], vars[n], vars[2*n]);
	f_Vm[0] *= u_range;
}

void Model_BOCF::ab_HH(const double *restrict Vm, const double *restrict vars, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	double u = (Vm[0]-u_min) * i_u_range;
	getab_HH(u, a_vars[0], b_vars[0], a_vars[n], b_vars[n], a_vars[2*n], b_vars[2*n]);
}

void Model_BOCF::f_NL(const double *restrict Vm, const double *restrict vars, double *restrict f_vars, double t, int n) const{
}

void Model_BOCF::f_vec(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double *restrict f_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		double a_v, a_w, a_s;
		double b_v, b_w, b_s;
	
		double u = (Vm[i]-u_min) * i_u_range;
		getf_u(u, f_Vm[i], vars[i], vars[n+i], vars[2*n+i]);
		f_Vm[i] *= u_range;
	
		getab_HH(Vm[i], a_v, b_v, a_w, b_w, a_s, b_s);
		f_vars[i] = a_v*vars[i]+b_v;
		f_vars[n+i] = a_w*vars[n+i]+b_w;
		f_vars[2*n+i] = a_s*vars[2*n+i]+b_s;
	}
}

void Model_BOCF::ab_vec(const double *restrict Vm, const double *restrict vars, double *restrict a_Vm, double *restrict b_Vm, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		double u = (Vm[i]-u_min) * i_u_range;
		getf_u(u, b_Vm[i], vars[i], vars[n+i], vars[2*n+i]);
		a_Vm[i] = 0;
		b_Vm[i] *= u_range;
	
		getab_HH(u, a_vars[i], b_vars[i], a_vars[n+i], b_vars[n+i], a_vars[2*n+i], b_vars[2*n+i]);
	}
}

void Model_BOCF::f_Vm_vec(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		double u = (Vm[i]-u_min) * i_u_range;
		getf_u(u, f_Vm[i], vars[i], vars[n+i], vars[2*n+i]);
		f_Vm[i] *= u_range;
	}
}

void Model_BOCF::ab_HH_vec(const double *restrict Vm, const double *restrict vars, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		double u = (Vm[0]-u_min) * i_u_range;
		getab_HH(u, a_vars[i], b_vars[i], a_vars[n+i], b_vars[n+i], a_vars[2*n+i], b_vars[2*n+i]);
	}
}

void Model_BOCF::f_NL_vec(const double *restrict Vm, const double *restrict vars, double *restrict f_vars, double t, int n) const{
}
