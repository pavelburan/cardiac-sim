#include "model_tnnp2006b.h"
#include "../configuration.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
//omp//#include <omp.h>

Model_TNNP2006b::Model_TNNP2006b(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix, int modNum, double VmThresh, double dVmdtThresh, const double TNNP2006b_Params[], double VmMin, double VmMax)
											:Model(grid,configFileName,keyPrefix,9,8,0,VmThresh,dVmdtThresh,VmMin,VmMax),modNum(modNum),
											R(p[0]),T(p[1]),F(p[2]),p_KNa(p[3]),K_pCa(p[4]),k_NaCa(p[5]),K_mNai(p[6]),K_mCa(p[7]),I_NaCa_alpha(p[8]),I_NaCa_gamma(p[9]),
											P_NaK(p[10]),K_mK(p[11]),K_mNa(p[12]),k_sat(p[13]),Na_o(p[14]),Na_i(p[15]),Ca_o(p[16]),Ca_i(p[17]),K_o(p[18]),K_i(p[19]),
											G_Na(p[20]),G_to(p[21]),G_Kr(p[22]),G_Cal(p[23]),G_Ks(p[24]),G_K1(p[25]),G_pCa(p[26]),G_pK(p[27]),G_bNa(p[28]),G_bCa(p[29]),c_j(p[30]),c_f(p[31]),
											E_Na(0.0),E_K(0.0),E_Ks(0.0),E_Ca(0.0),FdivRT(0.0),I_NaCa_Koeff1(0.0),I_NaCa_Koeff2(0.0),I_NaCa_Koeff3(0.0),I_NaK_Koeff1(0.0),
											s_infVm(0.0),s_koeff1(0.0),s_koeff2(0.0),s_koeff3(0.0),s_koeff4(0.0),s_koeff5(0.0)
{		
	setModel(modNum);
	
	if(TNNP2006b_Params != NULL)
		for(int i=0;i<32;i++)
			p[i] = TNNP2006b_Params[i];
	
	restingState[0] = -85.25677562;
	restingState[1] = 0.00170569;
	restingState[2] = 0.74565644;
	restingState[3] = 0.74565644;
	restingState[4] = 0.99990964;
	restingState[5] = 0.99948984;
	restingState[6] = 0.99998937;
	restingState[7] = 0.0032281;
	restingState[8] = 0.0001225;
	
	excitedState[0] = 0.0;
	excitedState[1] = 0.99978524*0.7 + 0.3*0.00170569;
	excitedState[2] = 4.31031549e-12*0.7 + 0.3*0.74565644;
	excitedState[3] = 4.29888885e-12*0.7 + 0.3*0.74565644;
	excitedState[4] = 0.28796308*0.7 + 0.3*0.99990964;
	excitedState[5] = 0.33074071*0.7 + 0.3*0.99948984;
	excitedState[6] = 0.0003797*0.7 + 0.3*0.99998937;
	excitedState[7] = 0.02526923*0.7 + 0.3*0.0032281;
	excitedState[8] = 0.01556574*0.7 + 0.3*0.0001225;
	
	/*excitedState[0] = 0.0;
	excitedState[1] = 0.00170569;
	excitedState[2] = 0.74565644;
	excitedState[3] = 0.74565644;
	excitedState[4] = 0.99990964;
	excitedState[5] = 0.99948984;
	excitedState[6] = 0.99998937;
	excitedState[7] = 0.0032281;
	excitedState[8] = 0.0001225;*/
	
	/*
	excitedState[0] = 0.0;
	excitedState[1] = 0.013466785134096645*0.5 + restingState[1]*0.5;
	excitedState[2] = 0.0080369236341538974*0.5 + restingState[2]*0.5;
	excitedState[3] = 0.16770627638045082*0.5 + restingState[3]*0.5;
	excitedState[4] = 0.22938134304378971*0.5 + restingState[4]*0.5;
	excitedState[5] = 0.63725773300847988*0.5 + restingState[5]*0.5;
	excitedState[6] = 0.39515751930580828*0.5 + restingState[6]*0.5;
	excitedState[7] = 0.003234453122370621*0.5 + restingState[7]*0.5;
	*/
}

void Model_TNNP2006b::setModel(int newModNum){
#ifdef _ERROR_
	if(newModNum < 0 || newModNum >2){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Modellnumber must be between 0 and 2!"<<std::endl;
		exit(1);
	}
#endif 
	//std::cerr<<"setModel"<<std::endl;
	modNum = newModNum;
	R = 8.3143;
	T = 310.0;
	F = 96.4867;
	p_KNa = 0.03;
	K_pCa = 0.0005;
	k_NaCa = 1000;
	K_mNai = 87.5;
	K_mCa = 1.38;
	I_NaCa_alpha = 2.5;
	I_NaCa_gamma = 0.35;
	P_NaK = 2.724;
	K_mK = 1.0;
	K_mNa = 40.0;
	k_sat = 0.1;
	Na_o = 140.0;
	Na_i = 7.67;
	Ca_o = 2.0;
	Ca_i = 0.00007;
	K_o = 5.4;
	K_i = 138.3;
	G_Na = 14.838;
	G_Kr = 0.101;//0.153;
	G_Cal = 0.2786;//3.980e-5;
	G_K1 = 5.405;
	G_pCa = 0.1238;
	G_pK = 0.0293;//0.0146;
	G_bNa = 0.00029;
	G_bCa = 0.000592;
	c_j = 1.0;
	c_f = 1.0;
	if(modNum == 0){
		//EPI
		G_to = 0.294;
		G_Ks = 0.257;//0.392;
		s_infVm = 20.0;
		s_koeff1 = 85.0;
		s_koeff2 = 45.0;
		s_koeff3 = 320.0;
		s_koeff4 = 5.0;
		s_koeff5 = 3.0;
	}
	else if(modNum == 1){
		//M
		G_to = 0.294;
		G_Ks = 0.098;
		s_infVm = 20.0;
		s_koeff1 = 85.0;
		s_koeff2 = 45.0;
		s_koeff3 = 320.0;
		s_koeff4 = 5.0;
		s_koeff5 = 3.0;
	}
	else if(modNum == 2){
		//ENDO
		G_to = 0.073;
		G_Ks = 0.257;//0.392;
		s_infVm = 28.0;
		s_koeff1 = 1000.0;
		s_koeff2 = 67.0;
		s_koeff3 = 1000.0;
		s_koeff4 = 0.0;
		s_koeff5 = 8.0;
	}
	
	E_Na = R*T/F*log(Na_o/Na_i);
	E_Ca = R*T/(F*2.0)*log(Ca_o/Ca_i);
	E_K = R*T/F*log(K_o/K_i);
	E_Ks = R*T/F*log((K_o+p_KNa*Na_o)/(K_i+p_KNa*Na_i));
	
	G_K1 *= sqrt(K_o/5.4);
	
	FdivRT = F/(R*T);
	
	I_NaCa_Koeff1 = k_NaCa / ((K_mNai*K_mNai*K_mNai + Na_o*Na_o*Na_o)*(K_mCa+Ca_o));
	I_NaCa_Koeff2 = Na_i*Na_i*Na_i*Ca_o;
	I_NaCa_Koeff3 = Na_o*Na_o*Na_o*Ca_i*I_NaCa_alpha;
	I_NaK_Koeff1 = P_NaK*K_o*Na_i / ((K_o+K_mK)*(Na_i+K_mNa));
	
}

void Model_TNNP2006b::readModelParams(){
	Model::readModelParams();
	cfg.readInto(modNum, "modNum");
	setModel(modNum);
	std::stringstream paramName;
	for(int i=0;i<32;i++){
		paramName.str("");
		paramName << "Param_" << i;
		cfg.readInto(p[i], paramName.str());
	}
}

__attribute__((always_inline)) void Model_TNNP2006b::getf_Vm(const double &Vm, double &dVmdt, const double &m, const double &h, const double &j, const double &f, const double &f2, const double &s, const double &xs, const double &xr1) const{
	double I_Na,I_to,I_Kr,I_Cal,I_Ks,I_K1,I_NaCa,I_NaK,I_pCa,I_pK,I_bNa,I_bCa;
	double temp;
	double alpha,beta;
		
	I_Na = G_Na*m*m*m*h*j*(Vm-E_Na);
	
	temp = 1.0 / (1.0+exp((20.0-Vm)/6.0)); //rInf
	I_to = G_to*temp*s*(Vm-E_K);
	
	temp = 1.0 / (1.0 + exp((Vm+88.0)/24.0)); //xr2Inf
	I_Kr = G_Kr*xr1*temp*(Vm-E_K);
	
	temp = 1.0 / (1.0+exp((-8.0-Vm)/7.5)); //dInf
	I_Cal = G_Cal*temp*f*f2*(Vm-60.0);
	
	I_Ks = G_Ks*xs*xs*(Vm-E_Ks);
	
	alpha = 0.1 / (1.0 + exp(0.06*(Vm-E_K-200.0)));
	beta = (3.0*exp(0.0002*(Vm-E_K+100.0)) + exp(0.1*(Vm-E_K-10.0))) / (1.0 + exp(-0.5*(Vm-E_K)));
	I_K1 = G_K1*alpha/(alpha+beta)*(Vm-E_K);
	
	temp = exp((I_NaCa_gamma-1.0)*Vm*FdivRT);
	I_NaCa = I_NaCa_Koeff1*( exp(I_NaCa_gamma*Vm*FdivRT)*I_NaCa_Koeff2 - temp*I_NaCa_Koeff3 ) / (1.0 + k_sat*temp); 
	
	I_NaK = I_NaK_Koeff1 / (1.0 + 0.1245*exp(-0.1*Vm*FdivRT) + 0.0353*exp(-1.0*Vm*FdivRT)); 
	
	I_pCa = G_pCa*Ca_i/(K_pCa+Ca_i); //opt
	
	I_pK = G_pK*(Vm-E_K)/(1.0+exp((25.0-Vm)/5.98));
	
	I_bNa = G_bNa*(Vm-E_Na);
	
	I_bCa = G_bCa*(Vm-E_Ca);
		 
	dVmdt = -(I_Na + I_to + I_Kr + I_Cal + I_Ks + I_K1 + I_NaCa + I_NaK + I_pCa + I_pK + I_bNa + I_bCa);
}

__attribute__((always_inline)) void Model_TNNP2006b::getab_HH(const double &Vm, double &a_m, double &b_m, double &a_h, double &b_h, double &a_j, double &b_j, double &a_f, double &b_f, double &a_f2, double &b_f2, double &a_s, double &b_s, double &a_xs, double &b_xs, double &a_xr1, double &b_xr1) const{
	double varInf;
	double alpha, beta, gamma;
	
		/*
		//LuoRudy
		alpha = 0.32*(Vm + 47.13) / (1.0 - exp(-0.1*(Vm+47.13)));
		beta = 0.08*exp(-Vm/11.0);
		a_m = -1.0*(alpha+beta);
		b_m = alpha;
		
		alpha = (Vm >= -40.0) ? 0.0 : 0.135*exp(-(Vm+80.0)/6.8); 
		beta = (Vm >= -40.0) ? 1.0/(0.13*(1.0+exp(-(Vm+10.66)/11.1))) : 3.56*exp(0.079*Vm) + 3.1e5*exp(0.35*Vm);
		a_h = -1.0*(alpha+beta);
		b_h = alpha;	
		
		alpha = (Vm >= -40.00) ? 0.0 : (Vm+37.78) * (-1.2714e5*exp(0.2444*Vm) - 3.474e-5*exp(-0.04391*Vm)) / (1.0 + exp(0.311*(Vm+79.23)));
		beta = (Vm >= -39.826) ? 0.3*exp(-2.535e-7*Vm) / (1.0 + exp(-0.1*(Vm+32.0))) : 0.1212*exp(-0.01052*Vm) / (1.0 + exp(-0.1378*(Vm+40.14)));
		a_j = -1.0*(alpha+beta)*c_j;
		b_j = alpha*c_j;
		*/
		
		alpha = 1.0 + exp((-60.0-Vm)/5.0);
		beta = 1.0 / ( 0.1 / (1.0+exp((Vm+35.0)/5.0)) + 0.1 / (1.0+exp((Vm-50.0)/200.0)) );
		varInf = 1.0 / pow(1.0+exp((-56.86-Vm)/9.03),2);
		a_m = -1.0 * alpha * beta;
		b_m = a_m * (-1.0) * varInf; 

		alpha = (Vm >= -40.0) ? 0.0 : 0.057*exp(-(Vm+80.0)/6.8); 
		beta = (Vm >= -40.0) ? 0.77/(0.13*(1.0+exp(-(Vm+10.66)/11.1))) : 2.7*exp(0.079*Vm) + 3.1e5*exp(0.3485*Vm);
		varInf = 1.0 / pow(1.0+exp((Vm+71.55)/7.43),2);
		a_h = -1.0 * (alpha + beta);
		b_h = a_h * (-1.0) * varInf; 
			
		alpha = (Vm >= -40.0) ? 0.0 : (Vm+37.78) * (-2.5428e4*exp(0.2444*Vm) - 6.948e-6*exp(-0.04391*Vm)) / (1.0 + exp(0.311*(Vm+79.23)));
		beta = (Vm >= -40.0) ? 0.6*exp(0.057*Vm) / (1.0 + exp(-0.1*(Vm+32.0))) : 0.02424*exp(-0.01052*Vm) / (1.0 + exp(-0.1378*(Vm+40.14)));
		//varInf = 1.0 / pow(1.0+exp((Vm+71.55)/7.43),2);
		a_j = -1.0 * (alpha + beta) * c_j;
		b_j = a_j * (-1.0) * varInf;

	alpha = 1102.5*exp(-1.0*pow((Vm+27.0)/15.0,2));
	beta = 200.0 / (1.0 + exp((13.0-Vm)/10.0));
	gamma = 180.0 / (1.0 + exp((Vm+30.0)/10.0)) + 20.0;
	varInf = 1.0 / (1.0+exp((Vm+20.0)/7.0));
	a_f = -1.0 / (alpha + beta + gamma) * ( (Vm > 0.0) ? c_f : 1.0 );
	b_f = a_f * (-1.0) * varInf;
	
	alpha = 600.0*exp(-1.0*pow((Vm+27.0),2)/170.0);
	beta = 7.75 / (1.0 + exp((25.0-Vm)/10.0));
	gamma = 16.0 / (1.0 + exp((Vm+30.0)/10.0));
	varInf = 0.67 / (1.0+exp((Vm+35.0)/7.0)) + 0.33;
	a_f2 = -1.0 / (alpha + beta + gamma);
	b_f2 = a_f2 * (-1.0) * varInf;
	
	varInf = 1.0 / (1.0+exp((Vm+s_infVm)/5.0));
	a_s = -1.0 / ( s_koeff1*exp(-pow(Vm+s_koeff2,2)/s_koeff3) + s_koeff4/(1.0+exp((Vm-20.0)/5.0)) + s_koeff5 );
	b_s = a_s * (-1.0) * varInf;
	
	alpha = 1400.0 / sqrt(1.0+exp((5.0-Vm)/6.0));
	beta = 1.0 / (1.0+exp((Vm-35.0)/15.0));
	varInf = 1.0 / (1.0+exp((-5.0-Vm)/14.0));
	a_xs = -1.0 / (alpha*beta + 80.0);
	b_xs = a_xs * (-1.0) * varInf;
	
	alpha = 1.0/450.0 *(1.0+exp((-45.0-Vm)/10.0));
	beta = 1.0/6.0 * (1.0+exp((Vm+30.0)/11.5));
	varInf = 1.0 / (1.0+exp((-26.0-Vm)/7.0));
	a_xr1 = -1.0 * alpha * beta;
	b_xr1 = a_xr1 * (-1.0) * varInf;
}

void Model_TNNP2006b::f(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double *restrict f_vars, double t, int n) const{
	double a_m, a_h, a_j, a_f, a_f2, a_s, a_xs, a_xr1;
	double b_m, b_h, b_j, b_f, b_f2, b_s, b_xs, b_xr1;
	getf_Vm(Vm[0], f_Vm[0], vars[0], vars[n], vars[2*n], vars[3*n], vars[4*n], vars[5*n], vars[6*n], vars[7*n]);
	getab_HH(Vm[0], a_m, b_m, a_h, b_h, a_j, b_j, a_f, b_f, a_f2, b_f2, a_s, b_s, a_xs, b_xs, a_xr1, b_xr1);
	f_vars[0] = a_m*vars[0]+b_m;
	f_vars[n] = a_h*vars[n]+b_h;
	f_vars[2*n] = a_j*vars[2*n]+b_j;
	f_vars[3*n] = a_f*vars[3*n]+b_f;
	f_vars[4*n] = a_f2*vars[4*n]+b_f2;
	f_vars[5*n] = a_s*vars[5*n]+b_s;
	f_vars[6*n] = a_xs*vars[6*n]+b_xs;
	f_vars[7*n] = a_xr1*vars[7*n]+b_xr1;
}

void Model_TNNP2006b::ab(const double *restrict Vm, const double *restrict vars, double *restrict a_Vm, double *restrict b_Vm, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	getf_Vm(Vm[0], b_Vm[0], vars[0], vars[n], vars[2*n], vars[3*n], vars[4*n], vars[5*n], vars[6*n], vars[7*n]);
	a_Vm[0] = 0;
	getab_HH(Vm[0], a_vars[0], b_vars[0], a_vars[n], b_vars[n], a_vars[2*n], b_vars[2*n], a_vars[3*n], b_vars[3*n], a_vars[4*n], b_vars[4*n], a_vars[5*n], b_vars[5*n], a_vars[6*n], b_vars[6*n], a_vars[7*n], b_vars[7*n]);
}

void Model_TNNP2006b::f_Vm(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double t, int n) const{
	getf_Vm(Vm[0], f_Vm[0], vars[0], vars[n], vars[2*n], vars[3*n], vars[4*n], vars[5*n], vars[6*n], vars[7*n]);
}

void Model_TNNP2006b::ab_HH(const double *restrict Vm, const double *restrict vars, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	getab_HH(Vm[0], a_vars[0], b_vars[0], a_vars[n], b_vars[n], a_vars[2*n], b_vars[2*n], a_vars[3*n], b_vars[3*n], a_vars[4*n], b_vars[4*n], a_vars[5*n], b_vars[5*n], a_vars[6*n], b_vars[6*n], a_vars[7*n], b_vars[7*n]);
}

void Model_TNNP2006b::f_NL(const double *restrict Vm, const double *restrict vars, double *f_vars, double t, int n) const{
}

void Model_TNNP2006b::f_vec(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double *restrict f_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		double a_m, a_h, a_j, a_f, a_f2, a_s, a_xs, a_xr1;
		double b_m, b_h, b_j, b_f, b_f2, b_s, b_xs, b_xr1;
		getf_Vm(Vm[i], f_Vm[i], vars[i], vars[n+i], vars[2*n+i], vars[3*n+i], vars[4*n+i], vars[5*n+i], vars[6*n+i], vars[7*n+i]);
		getab_HH(Vm[i], a_m, b_m, a_h, b_h, a_j, b_j, a_f, b_f, a_f2, b_f2, a_s, b_s, a_xs, b_xs, a_xr1, b_xr1);
		f_vars[i] = a_m*vars[i]+b_m;
		f_vars[n+i] = a_h*vars[n+i]+b_h;
		f_vars[2*n+i] = a_j*vars[2*n+i]+b_j;
		f_vars[3*n+i] = a_f*vars[3*n+i]+b_f;
		f_vars[4*n+i] = a_f2*vars[4*n+i]+b_f2;
		f_vars[5*n+i] = a_s*vars[5*n+i]+b_s;
		f_vars[6*n+i] = a_xs*vars[6*n+i]+b_xs;
		f_vars[7*n+i] = a_xr1*vars[7*n+i]+b_xr1;
	}
}

void Model_TNNP2006b::ab_vec(const double *restrict Vm, const double *restrict vars, double *restrict a_Vm, double *restrict b_Vm, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		getf_Vm(Vm[i], b_Vm[i], vars[i], vars[n+i], vars[2*n+i], vars[3*n+i], vars[4*n+i], vars[5*n+i], vars[6*n+i], vars[7*n+i]);
		a_Vm[i] = 0;
		getab_HH(Vm[i], a_vars[i], b_vars[i], a_vars[n+i], b_vars[n+i], a_vars[2*n+i], b_vars[2*n+i], a_vars[3*n+i], b_vars[3*n+i], a_vars[4*n+i], b_vars[4*n+i], a_vars[5*n+i], b_vars[5*n+i], a_vars[6*n+i], b_vars[6*n+i], a_vars[7*n+i], b_vars[7*n+i]);
	}
}

void Model_TNNP2006b::f_Vm_vec(const double *restrict Vm, const double *restrict vars, double *restrict f_Vm, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		getf_Vm(Vm[i], f_Vm[i], vars[i], vars[n+i], vars[2*n+i], vars[3*n+i], vars[4*n+i], vars[5*n+i], vars[6*n+i], vars[7*n+i]);
	}
}

void Model_TNNP2006b::ab_HH_vec(const double *restrict Vm, const double *restrict vars, double *restrict a_vars, double *restrict b_vars, double t, int n) const{
	#pragma omp simd
	#pragma ivdep
	for(int i=0;i<VECLEN;i++){
		getab_HH(Vm[i], a_vars[i], b_vars[i], a_vars[n+i], b_vars[n+i], a_vars[2*n+i], b_vars[2*n+i], a_vars[3*n+i], b_vars[3*n+i], a_vars[4*n+i], b_vars[4*n+i], a_vars[5*n+i], b_vars[5*n+i], a_vars[6*n+i], b_vars[6*n+i], a_vars[7*n+i], b_vars[7*n+i]);
	}
}

void Model_TNNP2006b::f_NL_vec(const double *restrict Vm, const double *restrict vars, double *f_vars, double t, int n) const{
}
