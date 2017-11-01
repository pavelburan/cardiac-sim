#ifndef OBS_MINPERIOD_H
#define OBS_MINPERIOD_H
#include "observer.h"
#include "../configuration.h"

class System;
class Obs_minPeriod : public Observer  
{
public:
	//Konstruktor
	Obs_minPeriod(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Obs_minPeriod(){}
private:
	Obs_minPeriod();
	Obs_minPeriod(const Obs_minPeriod &observer);
	Obs_minPeriod &operator=(const Obs_minPeriod &observer);
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
	double pos[3];
	int posIndex;
	std::vector<double> Ts;
	std::vector<double> tSend;
	std::vector<double> tReceive;
	double T;
};


#endif //OBS_MINPERIOD_H
