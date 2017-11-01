#ifndef OBS_OBSERVEPOINTS_H
#define OBS_OBSERVEPOINTS_H
#include "observer.h"
#include "../configuration.h"

class System;
class Obs_observePoints : public Observer  
{
public:
	//Konstruktor
	Obs_observePoints(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Obs_observePoints(){}
private:
	Obs_observePoints();
	Obs_observePoints(const Obs_observePoints &observer);
	Obs_observePoints &operator=(const Obs_observePoints &observer);
public:
	//Zugriffsmethoden
	bool resumeAble()const{ return true;}
	bool subTimeAble()const{ return false;}

	void readParams();
	void init();
	bool observe(double *y_prev, double *y, double *dVmdt, double *temp, double t, double timeStep, bool isResumeAbleStep);
	void finalize(double *y_prev, double *y, double *dVmdt, double *temp, double t);
	
private:
	void loadObsPoints(const std::string& fileName);
	void saveObsPoints(const std::string& fileName);
	
protected:
	int nRandObsPoints;
	std::string fixedObsPoints;
	std::vector<int> posIndices;
	std::vector< std::vector<double> > obsPointsData;
};


#endif //OBS_OBSERVEPOINTS_H
