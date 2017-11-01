#ifndef OBSERVER_H
#define OBSERVER_H
#include "../configuration.h"

class System;
class Efield;
class Grid;
class Observer
{
public:
	//Konstruktor
	Observer(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Observer(){}
private:
	Observer();
	Observer(const Observer &observer);
	Observer &operator=(const Observer &observer);
public:
	//Zugriffsmethoden
	virtual bool resumeAble()const{ return false;}
	virtual bool subTimeAble()const{ return false;}
	
	virtual void readParams()=0;
	virtual void init()=0;
	virtual bool observe(double *y_prev, double *y, double *dVmdt, double *temp, double t, double timeStep, bool isResumeAbleStep)=0;
	virtual void finalize(double *y_prev, double *y, double *dVmdt, double *temp, double t)=0;
	
	//Gridcreator
	static Observer* newObserver(const std::string& observerType, System& system, const std::string& configFileName, const std::string& keyPrefix);
	
protected:
	//Configuration
	Configuration cfg;
	//System
	System& system;
	//Efield
	Efield& efield;
	//Gitter
	Grid& grid;
};


#endif //OBSERVER_H
