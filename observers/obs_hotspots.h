#ifndef OBS_HOTSPOTS_H
#define OBS_HOTSPOTS_H
#include "observer.h"
#include "../configuration.h"
#include <list>

class System;
class Obs_hotSpots : public Observer  
{
public:
	//Konstruktor
	Obs_hotSpots(System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Obs_hotSpots(){}
private:
	Obs_hotSpots();
	Obs_hotSpots(const Obs_hotSpots &observer);
	Obs_hotSpots &operator=(const Obs_hotSpots &observer);
	void rekursivRemoveClusterElement(std::vector<int>& posIndicesHet, std::list<int>& posIndicesHetList, std::vector< bool >& isHet, std::vector< std::list<int>::iterator >& hetIterators, std::list<int>::iterator it);
	
public:
	//Zugriffsmethoden
	bool resumeAble()const{ return true;}
	bool subTimeAble()const{ return true;}

	
	void readParams();
	void init();
	bool observe(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t, double timeStep, bool isResumeAbleStep);
	void finalize(double *y_prev, double *y, double *dVmdt_prev, double *dVmdt, double *temp, double t);
	
protected:
	std::vector< std::vector< int > > hetPosIndices;
	std::vector< std::vector< int > > hetBorder1PosIndices;
	std::vector< std::vector< int > > hetBorder2PosIndices;
	std::vector< std::vector< int > > hetBorder3PosIndices;
	std::vector< std::vector< std::vector< int > > > hetBorder2NeighboursBorder1PosIndices;
	std::vector< std::vector< std::vector< int > > > hetBorder2NeighboursBorder3PosIndices;
	//Activated Hot Spots
	double maxDistance;
	double checkDt;
	std::vector< int > hetActivationStatus;
	std::vector< double > hetActivationTime;
	//Excited Hot Spots
	double dt_excitation;
	double t_next_excitation;
	int maxExcitationTimes;
	std::vector< std::vector< double > > hetExcitationTimes;

};


#endif //OBS_HOTSPOTS
