#include "obs_activatedhotspots.h"
#include "../system.h"
#include "../efields/efield.h"
#include "../grids/grid.h"
#include <limits>
#include <algorithm>


Obs_activatedHotSpots::Obs_activatedHotSpots(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),
														maxDistance(0.0), minVm(0.0), minDt(0.0), hetPosIndices(), hetBorder1PosIndices(), hetBorder2PosIndices(), hetBorder3PosIndices(), hetBorder2NeighboursBorder1PosIndices(), hetBorder2NeighboursBorder3PosIndices(), hetActivationStatus(), hetActivationTime(){
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
		std::vector< std::vector< int > > tempHetBorder2NeighboursBorder1PosIndices(tempHetBorder2PosIndices.size());
		std::vector< std::vector< int > > tempHetBorder2NeighboursBorder3PosIndices(tempHetBorder2PosIndices.size());
		for(int j=0;j<tempHetBorder2PosIndices.size();j++){
			std::vector<int> neighbours = grid.getPosIndicesSphere(tempHetBorder2PosIndices[j], maxDistance);
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
						tempHetBorder2NeighboursBorder3PosIndices[j].push_back(neighbours[k]);
					}
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
	hetActivationStatus.resize(hetBorder1PosIndices.size(),-1);
	hetActivationTime.resize(hetBorder1PosIndices.size(),0.0);
	std::cerr<<"cluster="<<cluster<<std::endl;
}

bool Obs_activatedHotSpots::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	/*std::vector<bool> isSet(system.getn(),false);
	for(int i=0;i<hetActivationStatus.size();i++){
		for(int j=0;j<hetBorder1PosIndices[i].size();j++){
			int posIndex = hetBorder1PosIndices[i][j];
			if(!isSet[posIndex]){
				y[posIndex] = -20.0;
				isSet[posIndex] = true;
			}
			
		}
	}
	for(int i=0;i<hetActivationStatus.size();i++){
		for(int j=0;j<hetBorder2PosIndices[i].size();j++){
			int posIndex = hetBorder2PosIndices[i][j];
			if(!isSet[posIndex]){
				y[posIndex] = 0.0;
				isSet[posIndex] = true;
			}
		}
	}
	for(int i=0;i<hetActivationStatus.size();i++){
		for(int j=0;j<hetBorder3PosIndices[i].size();j++){
			int posIndex = hetBorder3PosIndices[i][j];
			if(!isSet[posIndex]){
				y[posIndex] = 20.0;
				isSet[posIndex]= true;
			}	
		}
	}
	for(int i=0;i<hetActivationStatus.size();i++){
		for(int j=0;j<hetBorder1PosIndices[i].size();j++){
			for(int k=0;k<hetBorder2PosIndices[i].size();k++){
				for(int l=0;l<hetBorder2PosIndices[i].size();l++){
					int posIndex1 = hetBorder1PosIndices[i][j];
					int posIndex2 = hetBorder2PosIndices[i][k];
					int posIndex3 = hetBorder3PosIndices[i][l];
					if(posIndex1 == posIndex2 || posIndex1 == posIndex3 || posIndex2 == posIndex3)
						std::cerr<<"error"<<std::endl;
				}
			}
		}
	}
	for(int i=0;i<hetActivationStatus.size();i++){
		for(int j=0;j<hetBorder2PosIndices[i].size();j++){
			for(int k=0;k<hetBorder2NeighboursBorder1PosIndices[i][j].size();k++){
				bool found1 = false;
				for(int l=0;l<hetBorder1PosIndices[i].size();l++){
					if(hetBorder1PosIndices[i][l] == hetBorder2NeighboursBorder1PosIndices[i][j][k]){
						found1 = true;
						break;
					}
				}
				if(found1 == false)
				std::cerr<<"error1"<<std::endl;
			}
			for(int k=0;k<hetBorder2NeighboursBorder3PosIndices[i][j].size();k++){
				bool found3 = false;
				for(int l=0;l<hetBorder3PosIndices[i].size();l++){
					if(hetBorder3PosIndices[i][l] == hetBorder2NeighboursBorder3PosIndices[i][j][k]){
						found3 = true;
						break;
					}
				}
				if(found3 == false)
				std::cerr<<"error3"<<std::endl;
			}
		}
	}*/
	//grid.getVmThresh(posIndex)
	if(abs(efield.E(t)) > 0){
		for(int i=0;i<hetActivationStatus.size();i++){
			if(hetActivationStatus[i] == -1){
				for(int j=0;j<hetBorder2PosIndices[i].size();j++){
					int posIndex = hetBorder2PosIndices[i][j];
					if(dVmdt[posIndex] > grid.getdVmdtThresh(posIndex) && y[posIndex] > grid.getVmThresh(posIndex)){
						bool activatedByInside = false;
						for(int k=0;k<hetBorder2NeighboursBorder1PosIndices[i][j].size();k++){
							int posIndex2 = hetBorder2NeighboursBorder1PosIndices[i][j][k];
							if(y_prev[posIndex2] > grid.getVmThresh(posIndex)){
								activatedByInside = true;
								break;
							}
						}
						
						bool notActivatedByOutside = true;
						for(int k=0;k<hetBorder2NeighboursBorder3PosIndices[i][j].size();k++){
							int posIndex2 = hetBorder2NeighboursBorder3PosIndices[i][j][k];
							if(y_prev[posIndex2] > grid.getVmThresh(posIndex)){
								notActivatedByOutside = false;
								break;
							}
						}
						if(activatedByInside == true && notActivatedByOutside == true ){
							//for(int k=0;k<hetBorder1PosIndices[i].size();k++){
							hetActivationStatus[i] = 0;
							hetActivationTime[i] = t;
						}
						else if(activatedByInside == true && notActivatedByOutside == false ){
							hetActivationStatus[i] = -1;
							break;
						}
					}
				}
			}
			else if(hetActivationStatus[i] == 0 && (t-hetActivationTime[i]) > minDt){
				/*int count1 = 10;
				for(int j=0;j<hetBorder1PosIndices[i].size();j++){
					int posIndex = hetBorder1PosIndices[i][j];
					if(y[posIndex] > minVm)
						count1++;
				}
				int count2 = 0;
				for(int j=0;j<hetBorder2PosIndices[i].size();j++){
					int posIndex = hetBorder2PosIndices[i][j];
					if(y[posIndex] > minVm)
						count2++;
				}*/
				int count3 = 0;
				for(int j=0;j<hetBorder3PosIndices[i].size();j++){
					int posIndex = hetBorder3PosIndices[i][j];
					if(dVmdt[posIndex] > grid.getdVmdtThresh(posIndex) && y[posIndex] > grid.getVmThresh(posIndex))
						count3++;
				}
				if(count3 > 1)
					hetActivationStatus[i] = 1;
				else
					hetActivationStatus[i] = -1;
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


