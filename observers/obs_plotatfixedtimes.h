#ifndef OBS_PLOTATFIXEDTIMES_H
#define OBS_PLOTATFIXEDTIMES_H
#include "observer.h"
#include "../configuration.h"

class System;
class Obs_plotAtFixedTimes : public Observer  
{
public:
	//Konstruktor
	Obs_plotAtFixedTimes(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Obs_plotAtFixedTimes(){}
private:
	Obs_plotAtFixedTimes();
	Obs_plotAtFixedTimes(const Obs_plotAtFixedTimes &observer);
	Obs_plotAtFixedTimes &operator=(const Obs_plotAtFixedTimes &observer);
public:
	//Zugriffsmethoden
	bool resumeAble()const{ return true;}
	bool subTimeAble()const{ return true;}

	
	void readParams();
	void init();
	bool observe(double *y_prev, double *y, double *dVmdt, double *temp, double t, double timeStep, bool isResumeAbleStep);
	void finalize(double *y_prev, double *y, double *dVmdt, double *temp, double t){}
	
protected:
	bool forEverySubTimeIndex;
	int SubTimeIndexShift;
	bool saveAll;
	std::vector<double> timeSteps;
	int nextIndex;

};


#endif //OBS_PLOTATFIXEDTIMES_H
