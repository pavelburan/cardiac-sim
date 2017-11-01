#ifndef OBS_STIMULUSRANDOM_H
#define OBS_STIMULUSRANDOM_H
#include "observer.h"
#include "../configuration.h"

class System;
class Obs_stimulusRandom : public Observer  
{
public:
	//Konstruktor
	Obs_stimulusRandom(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Obs_stimulusRandom(){}
private:
	Obs_stimulusRandom();
	Obs_stimulusRandom(const Obs_stimulusRandom &observer);
	Obs_stimulusRandom &operator=(const Obs_stimulusRandom &observer);
public:
	//Zugriffsmethoden
	bool resumeAble()const{ return false;}
	bool subTimeAble()const{ return false;}

	
	void readParams();
	void init();
	bool observe(double *y_prev, double *y, double *dVmdt, double *temp, double t, double timeStep, bool isResumeAbleStep);
	void finalize(double *y_prev, double *y, double *dVmdt, double *temp, double t){}
	
protected:
	double t0;
	double dxi;
	double deta;
	double dzeta;
	int n;
	std::vector< std::vector<int> > posIndices;

};


#endif //OBS_STIMULUSRANDOM_H
