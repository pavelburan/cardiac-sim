#ifndef OBS_VELOCITY_H
#define OBS_VELOCITY_H
#include "observer.h"
#include "../configuration.h"

class System;
class Obs_velocity : public Observer  
{
public:
	//Konstruktor
	Obs_velocity(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Obs_velocity(){}
private:
	Obs_velocity();
	Obs_velocity(const Obs_velocity &observer);
	Obs_velocity &operator=(const Obs_velocity &observer);
public:
	//Zugriffsmethoden
	bool resumeAble()const{ return false;}
	bool subTimeAble()const{ return false;}

	
	void readParams();
	void init();
	bool observe(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t, double timeStep, bool isResumeAbleStep);
	void finalize(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t);
	
protected:
	double regionBeg[6];
	double regionEnd[6];
	std::vector<int> posIndicesBeg;
	std::vector<int> posIndicesEnd;
	std::vector<double> tBeg;
	std::vector<double> tEnd;
	int nTotal;
	int n;
};


#endif //OBS_VELOCITY_H
