#ifndef OBS_MINPACEPERIOD_H
#define OBS_MINPACEPERIOD_H
#include "observer.h"
#include "../configuration.h"

class System;
class Obs_minPacePeriod : public Observer  
{
public:
	//Konstruktor
	Obs_minPacePeriod(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Obs_minPacePeriod(){}
private:
	Obs_minPacePeriod();
	Obs_minPacePeriod(const Obs_minPacePeriod &observer);
	Obs_minPacePeriod &operator=(const Obs_minPacePeriod &observer);
public:
	//Zugriffsmethoden
	bool resumeAble()const{ return false;}
	bool subTimeAble()const{ return false;}

	
	void readParams();
	void init();
	bool observe(double *y_prev, double *y, double *dVmdt, double *temp, double t, double timeStep, bool isResumeAbleStep);
	void finalize(double *y_prev, double *y, double *dVmdt, double *temp, double t);
	
protected:
	double startT1;
	double dT1;
	double startT2;
	double dT2;
	double dtCheck;
	double tPuls;
	double tCheck;
	std::vector<double> Ts;
	double T;
};


#endif //OBS_MINPACEPERIOD_H
