#ifndef MODEL_LUORUDY_H
#define MODEL_LUORUDY_H
#include "LUT.h"
#include "model.h"


class Grid;
class Model_LuoRudy : public Model
{
public:
	//Konstruktoren
	Model_LuoRudy(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix="", double VmThresh=-40.0, double dVmdtThresh=10.0, const double LuoRudy_Params[] = NULL, double VmMin=-150.0, double VmMax=150.0);
	//Destruktor
	~Model_LuoRudy(){}
private:
	Model_LuoRudy();
	Model_LuoRudy(const Model_LuoRudy &model);
	Model_LuoRudy &operator=(const Model_LuoRudy &model);
	
public:		
	//Initialisieren
	void readModelParams();
	
	//Rechte Seite
private:
	void getf_Vm(const double &Vm, double &dVmdt, const double &h, const double &j, const double &m, const double &d, const double &f, const double &X, const double &Ca) const;
	void getab_HH(const double &Vm, double &a_h, double &b_h, double &a_j, double &b_j, double &a_m, double &b_m, double &a_d, double &b_d, double &a_f, double &b_f, double &a_X, double &b_X) const;
	void getf_NL(const double &Vm, const double &d, const double &f, const double &Ca, double  &dCadt) const;
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
	double getRefractoryState(int num)const{return (num >=0 && num <= 8) ? refractoryState[num] : 0.0;}
	double getExcitedState(int num)const{return (num >=0 && num <= 8) ? excitedState[num] : 0.0;}
	
private:
	double restingState[8];
	double refractoryState[8];
	double excitedState[8];
	//Modell Parameter
	double p[16];
	double &E_Na;
	double &K_0;
	double &E_K;
	double &G_K;
	double &E_K1;
	double &G_K1;
	double &E_Kp;
	double &G_Kp;
	double &G_Na;
	double &G_Si;
	double &c_h;
	double &c_j;
	double &c_m;
	double &c_d;
	double &c_f;
	double &c_X;
};
#endif //MODEL_LUORUDY_H
