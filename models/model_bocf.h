#ifndef MODEL_BOCF_H
#define MODEL_BOCF_H
#include "LUT.h"
#include "model.h"
class Grid;
class Model_BOCF : public Model
{
public:
	//Konstruktoren
	Model_BOCF(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix="", int modNum=0, double VmThresh=-60.0, double dVmdtThresh=10.0, const double BOCF_Params[] = NULL, double VmMin=-1000.0, double VmMax=1000.0);
	//Destruktor
	~Model_BOCF(){}
private:
	Model_BOCF();
	Model_BOCF(const Model_BOCF &model);
	Model_BOCF &operator=(const Model_BOCF &model);
	
public:	
	//Initialisieren
	void setModel(int newModNum);
	void readModelParams();
private:
	void getf_u(const double &Vm, double &dVmdt, const double &v, const double &w, const double &s) const;
	void getab_HH(const double &Vm, double &a_v, double &b_v, double &a_w, double &b_w, double &a_s, double &b_s) const;
public:
	//Rechte Seite
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
	int getModelNumber() const {return modNum;}
	double getRestingState(int num)const{return (num >=0 && num <= 4) ? restingState[num] : 0.0;}
	double getRefractoryState(int num)const{return (num >=0 && num <= 8) ? refractoryState[num] : 0.0;} //todo
	double getExcitedState(int num)const{return (num >=0 && num <= 4) ? excitedState[num] : 0.0;}

private:
	int modNum; // 0=EPI, 1=ENDO, 2=M, 3=PB, 4=TNNP
	double restingState[4];
	double refractoryState[4];
	double excitedState[4];
	//Modell Parameter
	double p[30];
	double &u_o;
	double &u_u;
	double &u_min;
	double &u_max;
	double &theta_v;
	double &theta_w;
	double &theta_m_v;
	double &theta_o;
	double &tau_m_v1;
	double &tau_m_v2;
	double &tau_p_v;
	double &tau_m_w1;
	double &tau_m_w2;
	double &k_m_w;
	double &u_m_w;
	double &tau_p_w;
	double &tau_fi;
	double &tau_o1;
	double &tau_o2;
	double &tau_so1;
	double &tau_so2;
	double &k_so;
	double &u_so;
	double &tau_s1;
	double &tau_s2;
	double &k_s;
	double &u_s;
	double &tau_si;
	double &tau_winf;	
	double &w_infstar;
	double u_range;
	double i_u_range;
	double i_tau_p_v;
	double i_tau_p_w;
	double i_tau_fi;
	double i_tau_si;
	double i_tau_winf;
};
#endif //MODEL_BOCF_H

