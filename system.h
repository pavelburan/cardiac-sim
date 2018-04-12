#ifndef SYSTEM_H
#define SYSTEM_H
#include "matrvec/matrvec.h"
#include "configuration.h"
#include "grids/grid.h"
#include <string>
#include <vector>
#include <map>

class Efield;
class Grid;
class Observer;
class OrderParameter; 
class System 
{
public:
//Konstruktor 
	System(const std::string& configFileName, const std::string& keyPrefix="");
//Destruktor
	~System();
private:
	System();
	System(const System &area);
	System &operator=(const System &area);

public:
	void readTimeIntegrationParams();
	void readPlotParams();
	void readEfieldParams();
	void readGridParams();
	void readObserverParams();
	void initSystem();
	void initEfield();
	void initGrid();
	void initObservers();
	void calcInitialCondition(double *y0);

	bool obs(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t, double dt);
	void obs_end(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t, double dt);
	void loadState(double *y, const std::string fileName, bool onlyVm=false) const;
	void loadInitialState(double *y0, const std::string fileName) const;
	void saveState(const double *y, const std::string fileName, bool onlyVm=false) const;
	void loadObs1DData(const std::string& fileName);
	void saveObs1DData(const std::string& fileName);
	void saveRandomSeedNumber(const std::string& fileName);
	void loadRandomSeedNumber(const std::string& fileName);
	
//Zugriffsmethoden
	//Zeitintegrationsparameterzugriffsmethoden
	double gett0() const {return t0;}
	double gettb() const {return tb;}
	double gette() const {return te;}
	double getIsTimeCoupledByPulse() const {return isTimeCoupledByPulse;}
	double getdt() const {return dt;}
	std::string getTimeIntScheme() const {return timeIntScheme;}
	std::string getTimeIntSchemeDiffusion() const {return timeIntSchemeDiffusion;}
	//Plotparameterzugriffsmethoden
	double gettp0() const {return tp0;}
	double getdtp() const {return dtp;}
	int getResumeIndex() const {return resumeIndex;}
	//Sonstige Zugriffsmethoden
	Efield& getEfield() const {return *efield;}
	Grid& getGrid() const {return *grid;}
	const OrderParameter* getOrderParameter(std::string index) const {return orderParameters[orderParameterTypeIndices.at(index)];}
	const Observer* getObserver(std::string index) const {return observers[observerTypeIndices.at(index)];}
	bool resume() const {return (resumeIndex >= 0 || cfg.getSubTimeIndex() > 0);}
	int getObsTimeLength(){return obs_time.size();}
	int getn() const {return grid->getn();}
	int getN() const {return N;}
	void setStateStructure(bool newVmVecOrdered, bool newModelVarsVecOrdered);
	int getStateSize() const {return (VmVecOrdered) ? N*(grid->getMaxNumVars()+1) : N*grid->getMaxNumVars();}
	int getStateIndexStride() const {return indexStride;}
	int getStateVarStride() const {return varStride;}
	int getStateFirstVecOrderedVar() const {return firstVecOrderedVar;}
	int getVecIndex(int posIndex, int numVar)const{int l=posIndex%VECLEN; return (N+(posIndex-l)*indexStride+(numVar-firstVecOrderedVar)*varStride+l);}
	void changeStateVector(double* y, int posIndex, double value, int numVar=0)const;
	void changeStateVector(double* y, const std::vector<int> &posIndices, double value, int numVar=0)const;
	void changeStateVector(double* y, const std::vector<int> &posIndices, double* values, int numVar=0)const;
	void changeStateVector(double* y, double value, int numVar=0)const;
	void changeStateVector(double* y, double* values, int numVar=0)const;

protected:
//Configuration
	Configuration cfg;
//Randomparameter
	int randomSeedNumber;
//Zeitintegrationsparameter
	double t0;
	double tb;
	double te;
	bool isTimeCoupledByPulse;
	double dt;
	std::string timeIntScheme;
	std::string timeIntSchemeDiffusion;
	bool VmVecOrdered;
	bool modelVarsVecOrdered;
	int N;
	int indexStride;
	int varStride;
	int firstVecOrderedVar;
//Plotparameter
	double tp0;
	double dtp;
	int resumeIndex;
	bool savePrevious;
	bool saveAll;
//E-Feld
	std::string efieldType;
	Efield *efield;
//Grid
	std::string gridType;
	Grid *grid;
//OrderParameter
	std::vector<std::string> orderParameterTypes;
	std::map<std::string, int> orderParameterTypeIndices;
	std::vector<OrderParameter*> orderParameters;
//Observers
	std::vector<std::string> observerTypes;
	std::map<std::string, int> observerTypeIndices;
	std::vector<Observer*> observers;
	//PLots
	std::vector<double> obs_time;
	std::vector<double> obs_E_data;
	int pcount;
	double tpLast;
	double tpSubTime;
};



#endif //SYSTEM_H

