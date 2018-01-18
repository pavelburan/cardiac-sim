#ifndef OBS_STIMULUS_H
#define OBS_STIMULUS_H
#include "observer.h"
#include "../configuration.h"

class System;
class Obs_stimulus : public Observer  
{
public:
	//Konstruktor
	Obs_stimulus(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Obs_stimulus(){}
private:
	Obs_stimulus();
	Obs_stimulus(const Obs_stimulus &observer);
	Obs_stimulus &operator=(const Obs_stimulus &observer);
public:
	//Zugriffsmethoden
	bool resumeAble()const{ return true;}
	bool subTimeAble()const{ return true;}

	
	void readParams();
	void init();
	bool observe(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t, double timeStep, bool isResumeAbleStep);
	void finalize(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t){}
	
protected:
	double t0;
	double dt;
	double te;
	double xi0;
	double dxi;
	double eta0;
	double deta;
	double zeta0;
	double dzeta;
	double tLast;
	std::vector<int> posIndices;

};


#endif //OBS_STIMULUS_H
