#ifndef MODEL_FENTONKARNA_H
#define MODEL_FENTONKARNA_H
#include "model.h"

class Grid;
class Model_FentonKarma : public Model
{
public:
	//Konstruktoren
	Model_FentonKarma(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix="", int modNum=0, double VmThresh=-60.0, double dVmdtThresh=10.0, const double MitchSchaeff_Params[] = NULL, double VmMin=-150.0, double VmMax=150.0);
	//Destruktor
	~Model_FentonKarma(){}
private:
	Model_FentonKarma();
	Model_FentonKarma(const Model_FentonKarma &model);
	Model_FentonKarma &operator=(const Model_FentonKarma &model);
	
public:	
	//Initialisieren
	void setModel(int newModNum);
	void readModelParams();
private:
	void getf_u(const double &Vm, double &dVmdt, const double &v, const double &w) const;
	void getab_HH(const double &Vm, double &a_v, double &b_v, double &a_w, double &b_w) const;
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
	double getRestingState(int num)const{return (num >=0 && num <= 3) ? restingState[num] : 0.0;}
	double getRefractoryState(int num)const{return (num >=0 && num <= 8) ? refractoryState[num] : 0.0;} //todo
	double getExcitedState(int num)const{return (num >=0 && num <= 3) ? excitedState[num] : 0.0;}

private:
	int modNum;
	double restingState[3];
	double refractoryState[3];
	double excitedState[3];
	//Modell Parameter
	double p[15];
	double &u_c;
	double &u_v;
	double &u_csi;
	double &u_min;
	double &u_max;
	double &tau_d;
	double &tau_o;
	double &tau_r;
	double &tau_si;
	double &tau_v_p;
	double &tau_v1_m;
	double &tau_v2_m;
	double &tau_w_p;
	double &tau_w_m;
	double &k;
	double u_range;
	double i_u_range;
	double i_tau_v_p;
	double i_tau_v1_m;
	double i_tau_v2_m; 
	double i_tau_w_p;
	double i_tau_w_m;
	double i_tau_d;
	double i_tau_o;
	double i_tau_r;
	double i_tau_si;
};
#endif //MODEL_FENTONKARNA_H

