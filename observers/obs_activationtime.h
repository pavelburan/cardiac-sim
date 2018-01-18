#ifndef OBS_ACTIVATIONTIME_H
#define OBS_ACTIVATIONTIME_H
#include "observer.h"
#include "../configuration.h"

class System;
class Obs_activationTime : public Observer  
{
public:
	//Konstruktor
	Obs_activationTime(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Obs_activationTime(){}
private:
	Obs_activationTime();
	Obs_activationTime(const Obs_activationTime &observer);
	Obs_activationTime &operator=(const Obs_activationTime &observer);
public:
	//Zugriffsmethoden
	bool resumeAble()const{ return true;}
	bool subTimeAble()const{ return false;}

	
	void readParams();
	void init();
	bool observe(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t, double timeStep, bool isResumeAbleStep);
	void finalize(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t);
	
protected:
	double ratio;
	double tau;
};


#endif //OBS_ACTIVATIONTIME_H
