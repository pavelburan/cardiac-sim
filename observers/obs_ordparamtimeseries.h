#ifndef OBS_ORDPARAMTIMESERIES_H
#define OBS_ORDPARAMTIMESERIES_H
#include "observer.h"
#include "../configuration.h"

class System;
class OrderParameter;
class Obs_ordParamTimeSeries : public Observer  
{
public:
	//Konstruktor
	Obs_ordParamTimeSeries(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Obs_ordParamTimeSeries(){}
private:
	Obs_ordParamTimeSeries();
	Obs_ordParamTimeSeries(const Obs_ordParamTimeSeries &observer);
	Obs_ordParamTimeSeries &operator=(const Obs_ordParamTimeSeries &observer);
public:
	//Zugriffsmethoden
	bool resumeAble()const{ return true;}
	bool subTimeAble()const{ return true;}

	void readParams();
	void init();
	bool observe(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t, double timeStep, bool isResumeAbleStep);
	void finalize(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t);
	
private:
	void loadOrderParameters(const std::string& fileName);
	void saveOrderParameters(const std::string& fileName);
	
protected:
	std::vector<std::string> types;
	std::vector<const OrderParameter*> orderParameters;
	double dt;
	double t_next;
	std::vector<double> timeSteps;
	std::vector< std::vector<double> > data;
};


#endif //OBS_ORDPARAMTIMESERIES_H
