#ifndef MODEL_MITCHSCHAEFF_H
#define MODEL_MITCHSCHAEFF_H
#include "model.h"

class Grid;
class Model_MitchSchaeff : public Model
{
public:
	//Konstruktoren
	Model_MitchSchaeff(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix="", double VmThresh=-60.0, double dVmdtThresh=10.0, const double MitchSchaeff_Params[] = NULL, double VmMin=-1000.0, double VmMax=1000.0);
	//Destruktor
	~Model_MitchSchaeff(){}
private:
	Model_MitchSchaeff();
	Model_MitchSchaeff(const Model_MitchSchaeff &model);
	Model_MitchSchaeff &operator=(const Model_MitchSchaeff &model);
	
public:	
	//Initialisieren
	void readModelParams();
private:
	void getf_u(const double &Vm, double &dVmdt, const double &h) const;
	void getab_HH(const double &Vm, double &a_h, double &b_h) const;
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
	double getRestingState(int num)const{return (num >=0 && num <= 2) ? restingState[num] : 0.0;}
	double getRefractoryState(int num)const{return (num >=0 && num <= 8) ? refractoryState[num] : 0.0;} //todo
	double getExcitedState(int num)const{return (num >=0 && num <= 2) ? excitedState[num] : 0.0;}

private:
	double restingState[2];
	double refractoryState[2];
	double excitedState[2];
	//Modell Parameter
	double p[7];
	double &u_gate;
	double &u_min;
	double &u_max;
	double &tau_in;
	double &tau_out;
	double &tau_open;
	double &tau_close;
	double u_range;
	double i_u_range;
	double i_tau_in;
	double i_tau_out;
	double i_tau_open;
	double i_tau_close;
	
};
#endif //MODEL_MITCHSCHAEFF_H

