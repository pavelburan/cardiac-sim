#ifndef MODEL_TNNP2006B_H
#define MODEL_TNNP2006B_H
#include "LUT.h"
#include "model.h"


class Grid;
class Model_TNNP2006b : public Model
{
public:
	//Konstruktoren
	Model_TNNP2006b(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix="", int modNum=0, double VmThresh=-40.0, double dVmdtThresh=10.0, const double TNNP2006b_Params[] = NULL, double VmMin=-150.0, double VmMax=150.0);
	//Destruktor
	~Model_TNNP2006b(){}
private:
	Model_TNNP2006b();
	Model_TNNP2006b(const Model_TNNP2006b &model);
	Model_TNNP2006b &operator=(const Model_TNNP2006b &model);
	
public:		
	//Initialisieren
	void setModel(int newModNum);
	void readModelParams();
	
	//Rechte Seite
private:
	void getf_Vm(const double &Vm, double &dVmdt, const double &m, const double &h, const double &j, const double &f, const double &f2, const double &s, const double &xs, const double &xr1) const;
	void getab_HH(const double &Vm, double &a_m, double &b_m, double &a_h, double &b_h, double &a_j, double &b_j, double &a_f, double &b_f, double &a_f2, double &b_f2, double &a_s, double &b_s, double &a_xs, double &b_xs, double &a_xr1, double &b_xr1) const;
public:
	void f(const double *Vm, const double *vars, double *f_Vm, double *f_vars, double t, int n) const;
	void ab(const double *Vm, const double *vars, double *a_Vm, double *b_Vm, double *a_vars, double *b_vars, double t, int n) const;
	void f_Vm(const double *Vm, const double *vars, double *f_Vm, double t, int n) const;
	void ab_HH(const double *Vm, const double *vars, double *a_vars, double *b_vars, double t, int n) const;
	void f_NL(const double *Vm, const double *vars, double *f_vars, double t, int n) const;
	void f_vec(const double *Vm, const double *vars, double *f_Vm, double *f_vars, double t, int n) const;
	void ab_vec(const double *Vm, const double *vars, double *a_Vm, double *b_Vm, double *a_vars, double *b_vars, double t, int n) const;
	void f_Vm_vec(const double *Vm, const double *vars, double *f_Vm, double t, int n) const;
	void ab_HH_vec(const double *Vm, const double *vars, double *a_vars, double *b_vars, double t, int n) const;
	void f_NL_vec(const double *Vm, const double *vars, double *f_vars, double t, int n) const;
	
	//Zugriffsfunktionen
	double getRestingState(int num)const{return (num >=0 && num <= 8) ? restingState[num] : 0.0;}
	double getRefractoryState(int num)const{return (num >=0 && num <= 8) ? refractoryState[num] : 0.0;} //todo
	double getExcitedState(int num)const{return (num >=0 && num <= 8) ? excitedState[num] : 0.0;}
	
private:
	int modNum; // 0=EPI, 1=M, 2=ENDO
	double restingState[9];
	double refractoryState[9];
	double excitedState[9];
	//Modell Parameter
	double p[32];
	double &R;
	double &T;
	double &F;
	double &p_KNa;
	double &K_pCa;
	double &k_NaCa;
	double &K_mNai;
	double &K_mCa;
	double &I_NaCa_alpha;
	double &I_NaCa_gamma;
	double &P_NaK;
	double &K_mK;
	double &K_mNa;
	double &k_sat;
	double &Na_o;
	double &Na_i;
	double &Ca_o;
	double &Ca_i;
	double &K_o;
	double &K_i;
	double &G_Na;
	double &G_to;
	double &G_Kr;
	double &G_Cal;
	double &G_Ks;
	double &G_K1;
	double &G_pCa;
	double &G_pK;
	double &G_bNa;
	double &G_bCa;
	double &c_j;
	double &c_f;
	
	double E_Na;
	double E_K;
	double E_Ks;
	double E_Ca;
	double FdivRT;
	double I_NaCa_Koeff1;
	double I_NaCa_Koeff2;
	double I_NaCa_Koeff3;
	double I_NaK_Koeff1;
	double s_infVm;
	double s_koeff1;
	double s_koeff2;
	double s_koeff3;
	double s_koeff4;
	double s_koeff5;
};
#endif //MODEL_TNNP2006B_H
