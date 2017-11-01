#ifndef MODEL_H
#define MODEL_H
#include "LUT.h"
#include "../configuration.h"

class Grid;
class Model
{
public:
	//Konstruktoren
	Model(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix="", int nVars=0, int nHH=0, int nNL=0, double VmThresh=0.0, double dVmdtThresh=0.0, double VmMin=0.0, double VmMax=0.0);
	//Destruktor
	virtual ~Model(){};
private:
	Model();
	Model(const Model &model);
	Model &operator=(const Model &model);

public:	
	//Initialisieren
	virtual void readModelParams() = 0;
	
	//Rechte Seite
	virtual void f(const double *Vm, const double *vars, double *f_Vm, double *f_vars, double t, int n) const=0;
	virtual void ab(const double *Vm, const double *vars, double *a_Vm, double *b_Vm, double *a_vars, double *b_vars, double t, int n) const=0;
	virtual void f_Vm(const double *Vm, const double *vars, double *f_Vm, double t, int n) const=0;
	virtual void ab_HH(const double *Vm, const double *vars, double *a_vars, double *b_vars, double t, int n) const=0;
	virtual void f_NL(const double *Vm, const double *vars, double *f_vars, double t, int n) const=0;
	virtual void f_vec(const double *Vm, const double *vars, double *f_Vm, double *f_vars, double t, int n) const=0;
	virtual void ab_vec(const double *Vm, const double *vars, double *a_Vm, double *b_Vm, double *a_vars, double *b_vars, double t, int n) const=0;
	virtual void f_Vm_vec(const double *Vm, const double *vars, double *f_Vm, double t, int n) const=0;
	virtual void ab_HH_vec(const double *Vm, const double *vars, double *a_vars, double *b_vars, double t, int n) const=0;
	virtual void f_NL_vec(const double *Vm, const double *vars, double *f_vars, double t, int n) const=0;
	
	//Zugriffsfunktionen
	int getnVars()const{return nVars;}
	int getnHH()const{return nHH;}
	int getnNL()const{return nNL;}
	double getVmThresh()const{return VmThresh;}
	double getdVmdtThresh()const{return dVmdtThresh;}
	double getVmMin()const{return VmMin;}
	double getVmMax()const{return VmMax;}
	virtual double getRestingState(int num)const = 0;
	virtual double getRefractoryState(int num)const = 0;
	virtual double getExcitedState(int num)const = 0;
	
	//Modelcreator
	static Model* newModel(const std::string& modelType, const Grid& grid, const std::string& configFileName, const std::string& keyPrefix);
	
protected:
//Configuration
	Configuration cfg;
//System
	const Grid& grid;
//Modelparameter	
	const int nVars;
	const int nHH;
	const int nNL;
	double VmThresh;
	double dVmdtThresh;
	double VmMin;
	double VmMax;	
};
#endif //MODEL_H
