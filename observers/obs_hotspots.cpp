#include "obs_hotspots.h"
#include "../system.h"
#include "../efields/efield.h"
#include "../grids/grid.h"
#include "../models/model.h"
#include "../timeintegration.h"
#include <limits>
#include <algorithm>


Obs_hotSpots::Obs_hotSpots(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),
														isExcitableDt(0.1), state(), a(), b(), Vm_check(0.0),lastVm(), 
														hetBorderWidth(0.0), hetPosIndices(), hetBorder1PosIndices(), hetBorder2PosIndices(), hetBorder3PosIndices(), hetBorder2NeighboursBorder1PosIndices(), hetBorder2NeighboursBorder3PosIndices(), hetBorderLevelsMap(),
														checkHetDt(0.0), minT_excitation(0.0), dt_check_excitation(0.0), t_next_check_excitation(0.0), maxExcitationTimes(0), minFreeHetBorderLevel(0), pointPosindices(), 
														checkActivatedHotSpots(false), hetActivationStatus(), hetActivationTime(), checkExcitedHotSpots(false), hetExcitationStatus(), hetExcitationTime(), hetExcitationTimes(), 
														checkActivationTime(false), activationRatio(0.0), startTimeE(0.0), excitablePosIndices(), excitablePosIndicesList(), activationTime(0.0),
														checkExcitedPoints(false), pointExcitationStatus(), pointExcitationTimes(), checkRefractionPoints(false), pointRefractionStatus(), pointRefractionTimes0(), pointRefractionTimes(){
}

void Obs_hotSpots::rekursivRemoveClusterElement(std::vector<int>& posIndicesHet, std::list<int>& posIndicesHetList, std::vector< bool >& isHet, std::vector< std::list<int>::iterator >& hetIterators, std::list<int>::iterator it){
	int index0 = *it;
	posIndicesHet.push_back(index0);
	isHet[index0] = false;
	it = posIndicesHetList.erase(it);
	
	std::vector<int> neighbours = grid.getPosIndicesSphere(index0, hetBorderWidth);

	for(int i=0;i<neighbours.size();i++){
		if(isHet[neighbours[i]]){
			it = hetIterators[neighbours[i]];
			rekursivRemoveClusterElement(posIndicesHet, posIndicesHetList, isHet, hetIterators, it);
		}
	}
}

bool Obs_hotSpots::isExcitable(double *__restrict__ y, int posIndex, double t){
	const Model &model = grid.getModel(posIndex);
	state[0] = y[posIndex];
	if(state[0] > grid.getVmThresh(posIndex))
		return false;
	else{
		for(int i=1;i<grid.getModel(posIndex).getnVars();i++){
			state[i] = y[system.getVecIndex(posIndex, i)];
		}
		state[0] = grid.getVmThresh(posIndex);
		for(double delta_t=isExcitableDt;delta_t<=5.0*isExcitableDt;delta_t+=isExcitableDt){
			timeInt_nonadaptive(model, state.data(), a.data(), b.data(), t+delta_t, isExcitableDt);
			//std::cerr<<delta_t<<" "<<state[0]<<std::endl;
			if(state[0] > Vm_check)
				return true;
		}
		if(state[0] > grid.getVmThresh(posIndex)){
			double lastVm = state[0];
			for(double delta_t=isExcitableDt;delta_t<=995.0*isExcitableDt;delta_t+=isExcitableDt){
				timeInt_nonadaptive(model, state.data(), a.data(), b.data(), t+delta_t, isExcitableDt);
				//std::cerr<<delta_t<<" "<<state[0]<<std::endl;
				if(state[0] < lastVm)
					return false;
				if(state[0] > Vm_check)
					return true;
				lastVm = state[0];
			}
		}
	}
	return false;
}

void Obs_hotSpots::readParams(){
	//excitability parameters
	cfg.readInto(isExcitableDt, "isExcitableDt");
	cfg.readInto(Vm_check, "Vm_check");
	//Het-border parameters
	cfg.readInto(hetBorderWidth, "hetBorderWidth");
	//Hot Spot parameters
	cfg.readInto(checkHetDt, "checkHetDt");
	//Excitation and refraction paramaters
	cfg.readInto(minT_excitation, "minT_excitation");
	cfg.readInto(dt_check_excitation, "dt_check_excitation");
	cfg.readInto(maxExcitationTimes, "maxExcitationTimes");
	//Observation points Parameters
	cfg.readInto(minFreeHetBorderLevel, "minFreeHetBorderLevel");
	//Activated Hot Spots
	cfg.readInto(checkActivatedHotSpots, "checkActivatedHotSpots");
	//Activation Time
	cfg.readInto(checkActivationTime, "checkActivationTime");
	cfg.readInto(activationRatio, "activationRatio");
	//Excited Hot Spots
	cfg.readInto(checkExcitedHotSpots, "checkExcitedHotSpots");
	//Excited Points
	cfg.readInto(checkExcitedPoints, "checkExcitedPoints");
	//Refraction Time points
	cfg.readInto(checkRefractionPoints, "checkRefractionPoints");
}

void Obs_hotSpots::init(){
	//Excitability
	state.resize(100,0.0);
	a.resize(100,0.0);
	b.resize(100,0.0);
	lastVm.resize(grid.getn(),0.0);	
	
	//Het-border
	std::vector< std::list<int>::iterator > hetIterators(system.getn());
	std::vector< bool > isHet(system.getn());
	
	std::list<int> posIndicesHetList;
	for(int i=0;i<grid.getn();i++){
		if(grid.getIsHet(i)){
			posIndicesHetList.push_back(i);
			hetIterators[i] = --posIndicesHetList.end();
			isHet[i] = true;
		}
		else
			isHet[i] = false;
	}
	int cluster = 0;
	std::list<int>::iterator it;
	while(posIndicesHetList.size() > 0){
		std::vector<int> tempHetPosIndices;
		it = posIndicesHetList.begin();
		rekursivRemoveClusterElement(tempHetPosIndices, posIndicesHetList, isHet, hetIterators, it);
		std::sort(tempHetPosIndices.begin(), tempHetPosIndices.end());
		hetPosIndices.push_back(tempHetPosIndices);
		
		std::vector<int> tempHetBorder1PosIndices;
		for(int j=0;j<tempHetPosIndices.size();j++){
			std::vector<int> neighbours = grid.getPosIndicesSphere(tempHetPosIndices[j], hetBorderWidth);
			for(int k=0;k<neighbours.size();k++)
				if(grid.getIsHet(neighbours[k]) == false)
					tempHetBorder1PosIndices.push_back(neighbours[k]);
		}
		std::sort(tempHetBorder1PosIndices.begin(), tempHetBorder1PosIndices.end());
		std::vector<int>::iterator it1;
		it1 = std::unique(tempHetBorder1PosIndices.begin(), tempHetBorder1PosIndices.end());
		tempHetBorder1PosIndices.resize( std::distance(tempHetBorder1PosIndices.begin(),it1) );
		hetBorder1PosIndices.push_back(tempHetBorder1PosIndices);
		
		std::vector<int> tempHetBorder2PosIndices;
		for(int j=0;j<tempHetBorder1PosIndices.size();j++){
			std::vector<int> neighbours = grid.getPosIndicesSphere(tempHetBorder1PosIndices[j], hetBorderWidth);
			for(int k=0;k<neighbours.size();k++){
				bool isNotHetBorder1 = true;
				for(int l=0;l<tempHetBorder1PosIndices.size();l++){
					if(neighbours[k] == tempHetBorder1PosIndices[l]){
						isNotHetBorder1 = false;
						break;
					}
				}
				if(isNotHetBorder1 && grid.getIsHet(neighbours[k]) == false)
					tempHetBorder2PosIndices.push_back(neighbours[k]);
			}	
		}
		std::sort(tempHetBorder2PosIndices.begin(), tempHetBorder2PosIndices.end());
		std::vector<int>::iterator it2;
		it2 = std::unique(tempHetBorder2PosIndices.begin(), tempHetBorder2PosIndices.end());
		tempHetBorder2PosIndices.resize( std::distance(tempHetBorder2PosIndices.begin(),it2) );
		hetBorder2PosIndices.push_back(tempHetBorder2PosIndices);
		
		std::vector<int> tempHetBorder3PosIndices;
		std::vector< std::vector< int > > tempHetBorder2NeighboursBorder1PosIndices(tempHetBorder2PosIndices.size());
		std::vector< std::vector< int > > tempHetBorder2NeighboursBorder3PosIndices(tempHetBorder2PosIndices.size());
		for(int j=0;j<tempHetBorder2PosIndices.size();j++){
			std::vector<int> neighbours = grid.getPosIndicesSphere(tempHetBorder2PosIndices[j], hetBorderWidth);
			for(int k=0;k<neighbours.size();k++){
				bool isHetBorder1 = false;
				for(int l=0;l<tempHetBorder1PosIndices.size();l++){
					if(neighbours[k] == tempHetBorder1PosIndices[l]){
						isHetBorder1 = true;
						break;
					}
				}
				if(isHetBorder1 == true){
					tempHetBorder2NeighboursBorder1PosIndices[j].push_back(neighbours[k]);
				}
				else{
					bool isHetBorder2 = false;
					for(int l=0;l<tempHetBorder2PosIndices.size();l++){
						if(neighbours[k] == tempHetBorder2PosIndices[l]){
							isHetBorder2 = true;
							break;
						}
					}
					if(isHetBorder2 == false && grid.getIsHet(neighbours[k]) == false){
						tempHetBorder3PosIndices.push_back(neighbours[k]);
					}
					if(grid.getIsHet(neighbours[k]) == false)
						tempHetBorder2NeighboursBorder3PosIndices[j].push_back(neighbours[k]);
				}
			}
			std::sort(tempHetBorder2NeighboursBorder1PosIndices[j].begin(), tempHetBorder2NeighboursBorder1PosIndices[j].end());
			std::sort(tempHetBorder2NeighboursBorder3PosIndices[j].begin(), tempHetBorder2NeighboursBorder3PosIndices[j].end());
		}
		std::sort(tempHetBorder3PosIndices.begin(), tempHetBorder3PosIndices.end());
		std::vector<int>::iterator it3;
		it3 = std::unique(tempHetBorder3PosIndices.begin(), tempHetBorder3PosIndices.end());
		tempHetBorder3PosIndices.resize( std::distance(tempHetBorder3PosIndices.begin(),it3) );
		hetBorder3PosIndices.push_back(tempHetBorder3PosIndices);
		hetBorder2NeighboursBorder1PosIndices.push_back(tempHetBorder2NeighboursBorder1PosIndices);
		hetBorder2NeighboursBorder3PosIndices.push_back(tempHetBorder2NeighboursBorder3PosIndices);
		
		
		cluster++;
	}
	hetBorderLevelsMap.resize(system.getn(),4);
	for(int i=0;i<grid.getn();i++){
		if(grid.getIsHet(i)){
			hetBorderLevelsMap[i] = 0;
		}
	}
	for(int i=0;i<hetBorder1PosIndices.size();i++){
		for(int j=0;j<hetBorder1PosIndices[i].size();j++){
			int posIndex = hetBorder1PosIndices[i][j];
			if(hetBorderLevelsMap[posIndex] > 0)
				hetBorderLevelsMap[posIndex] = 1;
		}
	}
	for(int i=0;i<hetBorder2PosIndices.size();i++){
		for(int j=0;j<hetBorder2PosIndices[i].size();j++){
			int posIndex = hetBorder2PosIndices[i][j];
			if(hetBorderLevelsMap[posIndex] > 1)
				hetBorderLevelsMap[posIndex] = 2;
		}
	}
	for(int i=0;i<hetBorder3PosIndices.size();i++){
		for(int j=0;j<hetBorder3PosIndices[i].size();j++){
			int posIndex = hetBorder3PosIndices[i][j];
			if(hetBorderLevelsMap[posIndex] > 2)
				hetBorderLevelsMap[posIndex] = 3;
		}
	}
	
	
	//Excitation and Refraction
	if(checkExcitedHotSpots || checkExcitedPoints || checkRefractionPoints){
		t_next_check_excitation = system.gettb();
	}
	
	//Observation points
	if(checkExcitedPoints || checkRefractionPoints){
		pointPosindices.reserve(system.getn()); 
		for(int i=0;i<grid.getn();i++){
			if(hetBorderLevelsMap[i] > minFreeHetBorderLevel)
				pointPosindices.push_back(i);
		}
	}
	
	
	//Activated Hot Spots
	if(checkActivatedHotSpots){
		hetActivationStatus.resize(hetBorder1PosIndices.size(), -1);
		hetActivationTime.resize(hetBorder1PosIndices.size(), std::numeric_limits<double>::max());
	}
	
	//Activation Time
	if(checkActivationTime){
		startTimeE = std::numeric_limits<double>::max();
		activationTime = std::numeric_limits<double>::max();
	}
	
	//Excited Hot Spots
	if(checkExcitedHotSpots){	
		hetExcitationStatus.resize(hetBorder1PosIndices.size(), -1);
		hetExcitationTime.resize(hetBorder1PosIndices.size(), std::numeric_limits<double>::max());
		hetExcitationTimes.resize(hetBorder1PosIndices.size());
	}
	
	//Excited Points
	if(checkExcitedPoints)
		pointExcitationStatus.resize(pointPosindices.size(),0);
		pointExcitationTimes.resize(pointPosindices.size());
		
	//Refraction Time
	if(checkRefractionPoints){
		pointRefractionStatus.resize(pointPosindices.size(),0);
		pointRefractionTimes0.resize(pointPosindices.size(),-1.0*std::numeric_limits<double>::max());
		pointRefractionTimes.resize(pointPosindices.size(),-1.0*std::numeric_limits<double>::max());
	}
	
	std::cerr<<"cluster="<<cluster<<std::endl;
}

bool Obs_hotSpots::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt_prev, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	/*
	for(int i=0;i<hetBorderLevelsMap.size();i++){
		if(hetBorderLevelsMap[i] == 1)
			y[i] = -20.0;
		else if(hetBorderLevelsMap[i] == 2)
			y[i] = 0.0;
		else if(hetBorderLevelsMap[i] == 3)
			y[i] = 20.0;
	}*/
	
	/*
	if( t>5.0){
		for(int i=0;i<excitablePosIndices.size();i++){
			int posIndex = excitablePosIndices[i];
			y[posIndex] = -100.0;
		}
	}*/

	//Activated Hot Spots
	if(checkActivatedHotSpots == true && (abs(efield.E(t)) > 0 || abs(efield.E(t-2.0*checkHetDt)))){
		for(int i=0;i<hetActivationStatus.size();i++){
			if(hetActivationStatus[i] == -1){
				for(int j=0;j<hetBorder2PosIndices[i].size();j++){
					int posIndex = hetBorder2PosIndices[i][j];
					if(dVmdt[posIndex] > grid.getdVmdtThresh(posIndex) && y[posIndex] > grid.getVmThresh(posIndex)){
						bool activatedByInside = false;
						for(int k=0;k<hetBorder2NeighboursBorder1PosIndices[i][j].size();k++){
							int posIndex2 = hetBorder2NeighboursBorder1PosIndices[i][j][k];
							if(dVmdt_prev[posIndex2] > grid.getdVmdtThresh(posIndex2) && y_prev[posIndex2] > grid.getVmThresh(posIndex2)){
								activatedByInside = true;
								break;
							}
						}
						
						bool notActivatedByOutside = true;
						for(int k=0;k<hetBorder2NeighboursBorder3PosIndices[i][j].size();k++){
							int posIndex2 = hetBorder2NeighboursBorder3PosIndices[i][j][k];
							if(dVmdt_prev[posIndex2] > grid.getdVmdtThresh(posIndex2) && y_prev[posIndex2] > grid.getVmThresh(posIndex2)){
								notActivatedByOutside = false;
								break;
							}
						}
						if(activatedByInside == true && notActivatedByOutside == true ){
							int count=1;
							for(int k=0;k<hetBorder2PosIndices[i].size();k++){
								int posIndex2 = hetBorder2PosIndices[i][k];
								if(y_prev[posIndex2] > grid.getVmThresh(posIndex2))
									count++;
							}
							if(hetBorder2PosIndices[i].size()/count > 2){
								hetActivationStatus[i] = 0;
								hetActivationTime[i] = t;
								break;
							}
						}
					}
				}
			}
			else if(hetActivationStatus[i] == 0 && (t-hetActivationTime[i]) > checkHetDt){
				/*int count1 = -100;
				for(int j=0;j<hetBorder1PosIndices[i].size();j++){
					int posIndex = hetBorder1PosIndices[i][j];
					if(dVmdt[posIndex] > grid.getdVmdtThresh(posIndex) && y[posIndex] > grid.getVmThresh(posIndex))
						count1++;
				}
				int count2 = -100;
				for(int j=0;j<hetBorder2PosIndices[i].size();j++){
					int posIndex = hetBorder2PosIndices[i][j];
					if(dVmdt[posIndex] > grid.getdVmdtThresh(posIndex) && y[posIndex] > grid.getVmThresh(posIndex))
						count2++;
				}*/
				int count3 = 0;
				for(int j=0;j<hetBorder3PosIndices[i].size();j++){
					int posIndex = hetBorder3PosIndices[i][j];
					if(dVmdt[posIndex] > grid.getdVmdtThresh(posIndex) && y[posIndex] > grid.getVmThresh(posIndex))
						count3++;
				}
				if(count3 > 0)
					hetActivationStatus[i] = 1;
				else{
					hetActivationStatus[i] = -1;
				}
			}
			if(hetActivationStatus[i] == 1){
				for(int j=0;j<hetPosIndices[i].size();j++)
					y[hetPosIndices[i][j]] = 60.0;
			}
			if(hetActivationStatus[i] == 0){
				for(int j=0;j<hetPosIndices[i].size();j++)
					y[hetPosIndices[i][j]] = 0.0;
			}
			if(hetActivationStatus[i] == -1){
				for(int j=0;j<hetPosIndices[i].size();j++)
					y[hetPosIndices[i][j]] = -40.0;
			}
		}
	}
	for(int i=0;i<hetActivationStatus.size();i++){
		if(hetActivationStatus[i] == 1){
			for(int j=0;j<hetPosIndices[i].size();j++)
				y[hetPosIndices[i][j]] = 60.0;
		}
		if(hetActivationStatus[i] == 0){
			for(int j=0;j<hetPosIndices[i].size();j++)
				y[hetPosIndices[i][j]] = 0.0;
		}
		if(hetActivationStatus[i] == -1){
			for(int j=0;j<hetPosIndices[i].size();j++)
				y[hetPosIndices[i][j]] = -40.0;
		}
	}
	
	//Activation Time
	if(checkActivationTime == true){
		if(efield.gettbeg() < t && startTimeE > t){
			startTimeE = t;
			for(int i=0;i<grid.getn();i++){
				if(isExcitable(y, i, t) && grid.getIsHet(i) == false)
					excitablePosIndicesList.push_back(i);
			}
			excitablePosIndices.reserve(excitablePosIndicesList.size());
			std::copy(std::begin(excitablePosIndicesList), std::end(excitablePosIndicesList), std::back_inserter(excitablePosIndices));
		}
		if(excitablePosIndicesList.size() > 0){
			for (std::list<int>::iterator it=excitablePosIndicesList.begin(); it!=excitablePosIndicesList.end(); ++it)
				if(dVmdt[*it] > grid.getdVmdtThresh(*it) && y[*it] > grid.getVmThresh(*it))
					it = excitablePosIndicesList.erase(it);
		}
		if(excitablePosIndices.size() >= 0 && activationTime > t){
			if((1.0 - double(excitablePosIndicesList.size())/excitablePosIndices.size()) > activationRatio)
				activationTime = t - startTimeE;
		}
	}
	
	if(checkExcitedHotSpots || checkExcitedPoints || checkRefractionPoints){
		if(t > t_next_check_excitation - Eps::t()){
			//Excited Hot Spots
			if(checkExcitedHotSpots){
				for(int i=0;i<hetExcitationTimes.size();i++){
					if(hetExcitationStatus[i] == -1){ 
						if(hetExcitationTimes[i].size() < maxExcitationTimes){
							if(hetExcitationTimes[i].size() == 0 || (t-hetExcitationTime[i]) > minT_excitation){
								for(int j=0;j<hetBorder1PosIndices[i].size();j++){
									int posIndex = hetBorder1PosIndices[i][j];
									if(lastVm[posIndex] < Vm_check && y[posIndex] > Vm_check && dVmdt[posIndex]>0.0){
										hetExcitationStatus[i] = 0;
										hetExcitationTime[i] = t;
										break;
									}
								}
							}
						}
					}
					else if(hetExcitationStatus[i] == 0 && (t-hetExcitationTime[i]) > checkHetDt - Eps::t() ){
						int count3 = 0;
						for(int j=0;j<hetBorder3PosIndices[i].size();j++){
							int posIndex = hetBorder3PosIndices[i][j];
							if(lastVm[posIndex] < Vm_check && y[posIndex] > Vm_check && dVmdt[posIndex]>0.0)
								count3++;
						}
						if(count3 > 0)
							hetExcitationTimes[i].push_back(hetExcitationTime[i]);
						hetExcitationStatus[i] = -1;
					}			
				}
			}
			//Excited Points
			if(checkExcitedPoints){
				for(int i=0;i<pointPosindices.size();i++){
					int posIndex = pointPosindices[i];
					if(pointExcitationStatus[i] >= 3)
						;
					else if(pointExcitationStatus[i] == 0){
						//if(y[posIndex] < grid.getVmThresh(posIndex)){
						if(isExcitable(y, posIndex, t)){
							if(pointExcitationTimes[i].size() == 0)
								pointExcitationStatus[i] = 2;
							else
								pointExcitationStatus[i] = 1;					
						}
					}
					else if(pointExcitationStatus[i] == 1){
						if((t-pointExcitationTimes[i].back()) > minT_excitation)
							pointExcitationStatus[i] = 2;
					}
					else if(pointExcitationStatus[i] == 2){
						if(lastVm[posIndex] < Vm_check && y[posIndex] > Vm_check && dVmdt[posIndex]>0.0){
							pointExcitationTimes[i].push_back(t);
							if(pointExcitationTimes[i].size() < maxExcitationTimes)
								pointExcitationStatus[i] = 0;
							else
								pointExcitationStatus[i] = 3;
						}
					}
				}
			}
			//Refraction Points
			if(checkRefractionPoints){
				for(int i=0;i<pointPosindices.size();i++){
					int posIndex = pointPosindices[i];
					if(pointRefractionStatus[i] >= 3)
						;
					else if(pointRefractionStatus[i] == 2){
						/*if(posIndex == 999026){ //LR 500484 LRMod  999026 FK 999592
							if(isExcitable(y, posIndex, t)){
								pointRefractionTimes[i] = t - pointRefractionTimes0[i];
								a[0] = pointRefractionTimes[i];
								pointRefractionStatus[i] = 2;
							}
						}*/
						if(isExcitable(y, posIndex, t)){
							pointRefractionTimes[i] = t - pointRefractionTimes0[i];
							pointRefractionStatus[i] = 3;
						}
					}
					else if(pointRefractionStatus[i] == 1){
						if(lastVm[posIndex] < Vm_check && y[posIndex] > Vm_check && dVmdt[posIndex]>0.0){
							pointRefractionTimes0[i] = t;
							//b[0]++;
							//std::cerr<<posIndex<<std::endl;
							pointRefractionStatus[i] = 2;
						}
					}
					else if(pointRefractionStatus[i] == 0){
						if(y[posIndex] < grid.getVmThresh(posIndex)){
							pointRefractionStatus[i] = 1;
						}
					}
					
					/*
					if(pointRefractionStatus[i] == 1){
						//if(posIndex == 999026){ //LR 500484 LRMod  999026 FK 999592
						//	if(isExcitable(y, posIndex, t)){
						//		pointRefractionTimes[i] = t - pointRefractionTimes0[i];
						//		a[0] = pointRefractionTimes[i];
						//		pointRefractionStatus[i] = 2;
						//	}
						//}
						if(isExcitable(y, posIndex, t)){
							pointRefractionTimes[i] = t - pointRefractionTimes0[i];
							pointRefractionStatus[i] = 2;
						}
					}
					else if(pointRefractionStatus[i] == 0){
						if(lastVm[posIndex] < Vm_check && y[posIndex] > Vm_check && dVmdt[posIndex]>0.0){
							pointRefractionTimes0[i] = t;
							//b[0]++;
							//std::cerr<<posIndex<<std::endl;
							pointRefractionStatus[i] = 1;
						}
					}
					*/
				}
			}
			for(int i=0;i<grid.getn();i++)
				lastVm[i] = y[i];
			t_next_check_excitation += dt_check_excitation;
		}
	}
	//std::cerr<<pointRefractionTimes0.size()<<" test2 "<<a[0]<<" "<<b[0]<<" "<<y[999026]<<std::endl;
	return false;
}

void Obs_hotSpots::finalize(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt_prev, double *__restrict__ dVmdt, double *__restrict__ temp, double t){
	if(cfg.getSubTimeIndex() <= 0){
		if(checkActivatedHotSpots || checkExcitedHotSpots){
			std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("hotSpotPositions.bin");
			std::ofstream file( fileName.c_str(), std::ios::out | std::ios::binary );
			if( file.is_open() ){
				std::vector< int > tempHetPosIndices(hetPosIndices.size(),0.0);
				for(int i=0;i<hetPosIndices.size();i++){
					double x = 0.0;
					double y = 0.0;
					double z = 0.0;
					for(int j=0;j<hetPosIndices[i].size();j++){
						std::vector<double> coordinate = grid.getCoordinate(hetPosIndices[i][j]);
						x += coordinate[0];
						y += coordinate[1];
						z += coordinate[2];
						tempHetPosIndices[i] += hetPosIndices[i][j];
					}
					x /= hetPosIndices[i].size();
					y /= hetPosIndices[i].size();
					z /= hetPosIndices[i].size();
					tempHetPosIndices[i] = grid.getPosIndex(x, y, z);
				}
				file.write(reinterpret_cast<const char*>(tempHetPosIndices.data()),sizeof(int)*tempHetPosIndices.size());
			}
			else{
				std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
				std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
				exit(1);
			}
			file.close();
			cfg.addCleanFile(fileName, 3);
		}
		if(checkExcitedPoints || checkRefractionPoints){
			std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("pointPositions.bin");
			std::ofstream file( fileName.c_str(), std::ios::out | std::ios::binary );
			std::vector< int > tempPointLevels(pointPosindices.size());
			for(int i=0;i<pointPosindices.size();i++)
				tempPointLevels[i] = hetBorderLevelsMap[pointPosindices[i]];
			if( file.is_open() ){
				file.write(reinterpret_cast<const char*>(pointPosindices.data()),sizeof(int)*pointPosindices.size());
				file.write(reinterpret_cast<const char*>(tempPointLevels.data()),sizeof(int)*tempPointLevels.size());
			}
			else{
				std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
				std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
				exit(1);
			}
			file.close();
			cfg.addCleanFile(fileName, 3);
		}
	}
	
	//Activated Hot Spots
	if(checkActivatedHotSpots){
		int anz = std::count(hetActivationStatus.begin(), hetActivationStatus.end(), 1);
		int total = hetActivationStatus.size();
		std::cerr<<"("<<anz<<"/"<<total<<")"<<std::endl;
		std::string fileName = cfg.getPlotFolderSubTimeSavePrefixFileName("anzActivatedHotSpots.txt");
		std::ofstream file( fileName.c_str(), std::ios::app );
		if( file.is_open() ){
			file<<cfg.getSubRepeatIndex()<<" "<<anz<<" "<<total<<"\n";
		}
		else{
			std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
			exit(1);
		}
		file.close();
		cfg.addCleanFile(fileName, 3);
		
		fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("isActivatedHotSpot.bin");
		if(cfg.getSubTimeIndex() <= 0)
			file.open( fileName.c_str(), std::ios::out | std::ios::binary );
		else
			file.open( fileName.c_str(), std::ios::app | std::ios::binary );
		if( file.is_open() ){
			bool *temp = new bool[hetActivationStatus.size()];
			for(int i=0;i<hetActivationStatus.size();i++)
				temp[i] = bool(hetActivationStatus[i] > 0);
			file.write(reinterpret_cast<const char*>(temp),sizeof(bool)*hetActivationStatus.size());
			delete[] temp;
		}
		else{
			std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
			exit(1);
		}
		file.close();
		cfg.addCleanFile(fileName, 3);	
	}
	
	//Activation Time
	if(checkActivationTime == true){
		std::cerr<<(1.0 - double(excitablePosIndicesList.size())/excitablePosIndices.size())<<std::endl;
		std::cerr<<"activationTime"<<activationTime<<" "<<excitablePosIndices.size()<<std::endl;
		std::string fileName = cfg.getPlotFolderSubTimeSavePrefixFileName("activationTime.txt");
		std::ofstream file( fileName.c_str(), std::ios::app );
		if( file.is_open() ){
			file<<cfg.getSubRepeatIndex()<<" "<<activationTime<<" "<<grid.getFraction(excitablePosIndices)<<"\n";
		}
		else{
			std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
			exit(1);
		}
		file.close();
		cfg.addCleanFile(fileName, 3);
	}
	
	if(checkActivatedHotSpots || checkActivationTime){
		std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("activatedHotSpotActivationTime.bin");
		std::ofstream file;
		if(cfg.getSubTimeIndex() <= 0)
			file.open( fileName.c_str(), std::ios::out | std::ios::binary );
		else
			file.open( fileName.c_str(), std::ios::app | std::ios::binary );
		if( file.is_open() ){
			if(checkActivatedHotSpots){
				double temp = double(std::count(hetActivationStatus.begin(), hetActivationStatus.end(), 1));
				file.write(reinterpret_cast<const char*>(&temp),sizeof(double));
			}
			if(checkActivationTime){
				file.write(reinterpret_cast<const char*>(&activationTime),sizeof(double));
				double temp = grid.getFraction(excitablePosIndices);
				file.write(reinterpret_cast<const char*>(&temp),sizeof(double));
			}
		}
		else{
			std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
			exit(1);
		}
		file.close();
		cfg.addCleanFile(fileName, 3);
	}
	
	//Excited Hot Spots
	if(checkExcitedHotSpots){
		std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("hetExcitationTimes.bin");
		std::ofstream file;
		if(cfg.getSubTimeIndex() <= 0)
			file.open( fileName.c_str(), std::ios::out | std::ios::binary );
		else
			file.open( fileName.c_str(), std::ios::app | std::ios::binary );
		if( file.is_open() ){
			for(int i=0;i<hetExcitationTimes.size();i++){
				hetExcitationTimes[i].resize(maxExcitationTimes, std::numeric_limits<double>::min());
				file.write(reinterpret_cast<const char*>(hetExcitationTimes[i].data()),sizeof(double)*hetExcitationTimes[i].size());
			}
		}
		else{
			std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
			exit(1);
		}
		file.close();
		cfg.addCleanFile(fileName, 3);
	}
	
	//Excited Points
	if(checkExcitedPoints){
		std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("pointExcitationTimes.bin");
		std::ofstream file;
		if(cfg.getSubTimeIndex() <= 0)
			file.open( fileName.c_str(), std::ios::out | std::ios::binary );
		else
			file.open( fileName.c_str(), std::ios::app | std::ios::binary );
		if( file.is_open() ){
			for(int i=0;i<pointExcitationTimes.size();i++){
				pointExcitationTimes[i].resize(maxExcitationTimes, -1.0*std::numeric_limits<double>::max());
				file.write(reinterpret_cast<const char*>(pointExcitationTimes[i].data()),sizeof(double)*pointExcitationTimes[i].size());
			}
		}
		else{
			std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
			exit(1);
		}
		file.close();
		cfg.addCleanFile(fileName, 3);
	}
	
	//Refraction Points
	if(checkRefractionPoints){
		std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("pointRefractionTimes.bin");
		std::ofstream file;
		if(cfg.getSubTimeIndex() <= 0)
			file.open( fileName.c_str(), std::ios::out | std::ios::binary );
		else
			file.open( fileName.c_str(), std::ios::app | std::ios::binary );
		if( file.is_open() ){
			for(int i=0;i<pointRefractionTimes.size();i++){
				file.write(reinterpret_cast<const char*>(&pointRefractionTimes0[i]),sizeof(double));
				file.write(reinterpret_cast<const char*>(&pointRefractionTimes[i]),sizeof(double));
			}
		}
		else{
			std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
			exit(1);
		}
		file.close();
		cfg.addCleanFile(fileName, 3);
	}
	
	
}


