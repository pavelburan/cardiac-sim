#ifndef OBS_TERMINATIONTIME_H
#define OBS_TERMINATIONTIME_H
#include "observer.h"
#include "../configuration.h"

class System;
class Obs_terminationTime : public Observer  
{
public:
	//Konstruktor
	Obs_terminationTime(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Obs_terminationTime(){}
private:
	Obs_terminationTime();
	Obs_terminationTime(const Obs_terminationTime &observer);
	Obs_terminationTime &operator=(const Obs_terminationTime &observer);
public:
	//Zugriffsmethoden
	bool resumeAble()const{ return true;}
	bool subTimeAble()const{ return true;}

	void readParams();
	void init(){};
	bool observe(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t, double timeStep, bool isResumeAbleStep);
	void finalize(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t);
	
private:
	void loadObsPoints(const std::string& fileName);
	void saveObsPoints(const std::string& fileName);
	
protected:
	bool finishAfterTermination;
	double terminationTime;
	double timeAfterTermination;
};


#endif //OBS_TERMINATIONTIME_H
