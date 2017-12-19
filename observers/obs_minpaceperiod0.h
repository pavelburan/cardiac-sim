#ifndef OBS_MINPACEPERIOD0_H
#define OBS_MINPACEPERIOD0_H
#include "observer.h"
#include "../configuration.h"

class System;
class Obs_minPacePeriod0 : public Observer  
{
public:
	//Konstruktor
	Obs_minPacePeriod0(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Obs_minPacePeriod0(){}
private:
	Obs_minPacePeriod0();
	Obs_minPacePeriod0(const Obs_minPacePeriod0 &observer);
	Obs_minPacePeriod0 &operator=(const Obs_minPacePeriod0 &observer);
public:
	//Zugriffsmethoden
	bool resumeAble()const{ return false;}
	bool subTimeAble()const{ return false;}

	
	void readParams();
	void init();
	bool observe(double *y_prev, double *y, double *dVmdt, double *temp, double t, double timeStep, bool isResumeAbleStep);
	void finalize(double *y_prev, double *y, double *dVmdt, double *temp, double t);
	
protected:
	double startT;
	double dT;
	double dtCheck;
	double tCheck;
	std::vector<double> Ts;
	double T;
};


#endif //OBS_MINPACEPERIOD0_H
