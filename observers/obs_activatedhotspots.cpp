#include "obs_activatedhotspots.h"
#include "../system.h"
#include "../efields/efield.h"
#include "../grids/grid.h"
#include <limits>
#include <algorithm>


Obs_activatedHotSpots::Obs_activatedHotSpots(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),
														maxDistance(0.0), minVm(0.0), minDt(0.0), hetPosIndices(), hetBorder1PosIndices(), hetBorder2PosIndices(), hetBorder3PosIndices(), hetActivationStatus(), hetActivationTime(){
}

void Obs_activatedHotSpots::rekursivRemoveClusterElement(std::vector<int>& posIndicesHet, std::list<int>& posIndicesHetList, std::vector< bool >& isHet, std::vector< std::list<int>::iterator >& hetIterators, std::list<int>::iterator it){
	int index0 = *it;
	posIndicesHet.push_back(index0);
	isHet[index0] = false;
	it = posIndicesHetList.erase(it);
	
	std::vector<int> neighbours = grid.getPosIndicesSphere(index0, maxDistance);

	for(int i=0;i<neighbours.size();i++){
		if(isHet[neighbours[i]]){
			it = hetIterators[neighbours[i]];
			rekursivRemoveClusterElement(posIndicesHet, posIndicesHetList, isHet, hetIterators, it);
		}
	}
}

void Obs_activatedHotSpots::readParams(){
	cfg.readInto(maxDistance, "maxDistance");
	cfg.readInto(minVm, "minVm");
	cfg.readInto(minDt, "minDt");
}

void Obs_activatedHotSpots::init(){
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
			std::vector<int> neighbours = grid.getPosIndicesSphere(tempHetPosIndices[j], maxDistance);
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
			std::vector<int> neighbours = grid.getPosIndicesSphere(tempHetBorder1PosIndices[j], maxDistance);
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
		for(int j=0;j<tempHetBorder2PosIndices.size();j++){
			std::vector<int> neighbours = grid.getPosIndicesSphere(tempHetBorder2PosIndices[j], maxDistance);
			for(int k=0;k<neighbours.size();k++){
				bool isNotHetBorder1Border2 = true;
				for(int l=0;l<tempHetBorder1PosIndices.size();l++){
					if(neighbours[k] == tempHetBorder1PosIndices[l]){
						isNotHetBorder1Border2 = false;
						break;
					}
				}
				for(int l=0;l<tempHetBorder2PosIndices.size();l++){
					if(neighbours[k] == tempHetBorder2PosIndices[l]){
						isNotHetBorder1Border2 = false;
						break;
					}
				}
				if(isNotHetBorder1Border2 && grid.getIsHet(neighbours[k]) == false)
					tempHetBorder3PosIndices.push_back(neighbours[k]);
			}	
		}
		std::sort(tempHetBorder3PosIndices.begin(), tempHetBorder3PosIndices.end());
		std::vector<int>::iterator it3;
		it3 = std::unique(tempHetBorder3PosIndices.begin(), tempHetBorder3PosIndices.end());
		tempHetBorder3PosIndices.resize( std::distance(tempHetBorder3PosIndices.begin(),it3) );
		hetBorder3PosIndices.push_back(tempHetBorder3PosIndices);
		
		cluster++;
	}
	hetActivationStatus.resize(hetBorder1PosIndices.size(),-1);
	hetActivationTime.resize(hetBorder1PosIndices.size(),0.0);
	std::cerr<<"cluster="<<cluster<<std::endl;
}

bool Obs_activatedHotSpots::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	if(abs(efield.E(t)) > 0){
		for(int i=0;i<hetActivationStatus.size();i++){
			if(hetActivationStatus[i] == -1){
				for(int j=0;j<hetBorder1PosIndices[i].size();j++){
					int posIndex = hetBorder1PosIndices[i][j];
					if(dVmdt[posIndex] > grid.getdVmdtThresh(posIndex) && y[posIndex] > minVm/*grid.getVmThresh(posIndex)*/){
						bool notActivatedByNeighbour = true;
						for(int k=0;k<hetBorder3PosIndices[i].size();k++){
							int posIndex2 = hetBorder3PosIndices[i][k];
							if(y_prev[posIndex2] > minVm/*grid.getVmThresh(posIndex2)*/){
								notActivatedByNeighbour = false;
								break;
							}
						}
						if(notActivatedByNeighbour){
							hetActivationStatus[i] = 0;
							hetActivationTime[i] = t;
						}
						else
							hetActivationStatus[i] = -1;
						break;
					}
				}
			}
			else if(hetActivationStatus[i] == 0 && (t-hetActivationTime[i]) > minDt){
				int count = 0;
				for(int j=0;j<hetBorder2PosIndices[i].size();j++){
					int posIndex = hetBorder2PosIndices[i][j];
					if(dVmdt[posIndex] > grid.getdVmdtThresh(posIndex) && y[posIndex] > minVm)
						count++;
				}
				if(count > 2)
					hetActivationStatus[i] = 1;
			}
			if(hetActivationStatus[i] == 1){
				for(int j=0;j<hetPosIndices[i].size();j++)
					y[hetPosIndices[i][j]] = 60.0;
			}
			if(hetActivationStatus[i] == 0){
				for(int j=0;j<hetPosIndices[i].size();j++)
					y[hetPosIndices[i][j]] = 0.0;
			}
		}
	}
	return false;
}

void Obs_activatedHotSpots::finalize(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t){
	int anz = std::count(hetActivationStatus.begin(), hetActivationStatus.end(), 1);
	int total = hetActivationStatus.size();
	std::cerr<<"("<<anz<<"/"<<total<<")"<<std::endl;
	std::string fileName = cfg.getPlotFolderSubTimeSavePrefixFileName("activatedHotSpots.txt");
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
}


