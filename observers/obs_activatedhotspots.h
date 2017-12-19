#ifndef OBS_ACTIVATEDHOTSPOTS_H
#define OBS_ACTIVATEDHOTSPOTS_H
#include "observer.h"
#include "../configuration.h"
#include <list>

class System;
class Obs_activatedHotSpots : public Observer  
{
public:
	//Konstruktor
	Obs_activatedHotSpots(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Obs_activatedHotSpots(){}
private:
	Obs_activatedHotSpots();
	Obs_activatedHotSpots(const Obs_activatedHotSpots &observer);
	Obs_activatedHotSpots &operator=(const Obs_activatedHotSpots &observer);
	void rekursivRemoveClusterElement(std::vector<int>& posIndicesHet, std::list<int>& posIndicesHetList, std::vector< bool >& isHet, std::vector< std::list<int>::iterator >& hetIterators, std::list<int>::iterator it);
	
public:
	//Zugriffsmethoden
	bool resumeAble()const{ return false;}
	bool subTimeAble()const{ return false;}

	
	void readParams();
	void init();
	bool observe(double *y_prev, double *y, double *dVmdt, double *temp, double t, double timeStep, bool isResumeAbleStep);
	void finalize(double *y_prev, double *y, double *dVmdt, double *temp, double t);
	
protected:
	double maxDistance;
	double minVm;
	double minDt;
	std::vector< std::vector< int > > hetPosIndices;
	std::vector< std::vector< int > > hetBorder1PosIndices;
	std::vector< std::vector< int > > hetBorder2PosIndices;
	std::vector< std::vector< int > > hetBorder3PosIndices;
	std::vector< int > hetActivationStatus;
	std::vector< double > hetActivationTime;

};


#endif //OBS_ACTIVATEDHOTSPOTS
