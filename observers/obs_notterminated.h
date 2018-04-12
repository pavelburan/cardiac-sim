#ifndef OBS_NOTTERMINATED_H
#define OBS_NOTTERMINATED_H
#include "observer.h"
#include "../configuration.h"

class System;
class Obs_notTerminated : public Observer  
{
public:
	//Konstruktor
	Obs_notTerminated(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Obs_notTerminated(){}
private:
	Obs_notTerminated();
	Obs_notTerminated(const Obs_notTerminated &observer);
	Obs_notTerminated &operator=(const Obs_notTerminated &observer);
public:
	//Zugriffsmethoden
	bool resumeAble()const{ return true;}
	bool subTimeAble()const{ return true;}

	void readParams();
	void init();
	bool observe(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t, double timeStep, bool isResumeAbleStep);
	void finalize(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t);
	
private:
	void loadObsPoints(const std::string& fileName);
	void saveObsPoints(const std::string& fileName);
	
protected:
	double tSave;
	std::string tSaveVecRaw;
	std::vector<double> tSaveVec;
	double fracNotTerminated;
	std::vector<double> fracNotTerminatedVec;
	bool isTerminated;
	std::vector<bool> isTerminatedVec;

	
};


#endif //OBS_NOTTERMINATED_H
