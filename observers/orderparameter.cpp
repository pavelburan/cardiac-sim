#include "orderparameter.h"
#include "../system.h"
#include "../grids/grid.h"
#include <iostream>
#include <list>

OrderParameter* OrderParameter::newOrderParameter(const std::string& orderParameterType, const System& system, const std::string& configFileName, const std::string& keyPrefix){
	if(orderParameterType == "MeanVm")
		return new OP_meanVm(system, configFileName, keyPrefix);
	else if(orderParameterType == "Excitable")
		return new OP_excitable(system, configFileName, keyPrefix);
	else if(orderParameterType == "Excited")
		return new OP_excited(system, configFileName, keyPrefix);
	else if(orderParameterType == "Exciting")
		return new OP_exciting(system, configFileName, keyPrefix);
	else if(orderParameterType == "Cluster")
		return new OP_cluster(system, configFileName, keyPrefix);
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"orderParameterType="<<orderParameterType<<" existiert nicht!"<<std::endl;
		exit(1);
		return NULL;
	}
}

OrderParameter::OrderParameter(const System& system, const std::string& configFileName, const std::string& keyPrefix):system(system),efield(system.getEfield()),grid(system.getGrid()),cfg(configFileName,keyPrefix){
}


OP_meanVm::OP_meanVm(const System& system, const std::string& configFileName, const std::string& keyPrefix):OrderParameter(system,configFileName,keyPrefix){
}

double OP_meanVm::calc(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep)const{
	return grid.getMeanIntegral(y);
}


OP_excitable::OP_excitable(const System& system, const std::string& configFileName, const std::string& keyPrefix):OrderParameter(system,configFileName,keyPrefix){
}

double OP_excitable::calc(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep)const{
	std::vector<bool> isPointOfFraction(grid.getn(), true);
	#pragma omp parallel for //simd
	#pragma ivdep
	for(int i=0;i<grid.getn();i++){
		isPointOfFraction[i] =  y[i] < grid.getVmThresh(i) ? true : false;
	}
	return grid.getFraction(isPointOfFraction);
}


OP_excited::OP_excited(const System& system, const std::string& configFileName, const std::string& keyPrefix):OrderParameter(system,configFileName,keyPrefix){
}

double OP_excited::calc(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep)const{
	std::vector<bool> isPointOfFraction(grid.getn(), true);
	#pragma omp parallel for //simd
	#pragma ivdep
	for(int i=0;i<grid.getn();i++){
		isPointOfFraction[i] =  y[i] > 0.0 ? true : false;
	}
	return grid.getFraction(isPointOfFraction);
}


OP_exciting::OP_exciting(const System& system, const std::string& configFileName, const std::string& keyPrefix):OrderParameter(system,configFileName,keyPrefix){
}

double OP_exciting::calc(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep)const{
	std::vector<bool> isPointOfFraction(grid.getn(), true);
	#pragma omp parallel for //simd
	#pragma ivdep
	for(int i=0;i<grid.getn();i++){
		isPointOfFraction[i] = (dVmdt[i] > grid.getdVmdtThresh(i) && y[i] > grid.getVmThresh(i)) ? true : false;
	}
	return grid.getFraction(isPointOfFraction);
}


OP_cluster::OP_cluster(const System& system, const std::string& configFileName, const std::string& keyPrefix):OrderParameter(system,configFileName,keyPrefix),maxDistance(0.0),itTotalVector(),isExciting(){
}

void OP_cluster::readParams(){
		cfg.readInto(maxDistance, "maxDistance");
}

void OP_cluster::init(){
	itTotalVector.resize(system.getn());
	isExciting.resize(system.getn());
}

void OP_cluster::rekursivRemoveClusterElement(std::list<int>& posIndicesExciting, std::vector<std::list<int>::iterator>& itVector)const{
	std::list<int>::iterator it = itVector.back();
	for(int i=0;i<itVector.size();i++){
		if(itVector[i] == it)
			++itVector[i];
	}
	int index0 = *it;
	isExciting[index0] = false;
	it = posIndicesExciting.erase(it);
	
	std::vector<int> neighbours = grid.getPosIndicesSphere(index0, maxDistance);
	/*std::cerr<<"it0="<<index0<<" size="<<posIndicesExciting.size()<<" vectorSize="<<itVector.size()<<std::endl;
	std::cerr<<"bla "<<neighbours.size()<<std::endl;
	for(int i=0;i<neighbours.size();i++)
		std::cerr<<neighbours[i]<<", ";
	std::cerr<<std::endl;*/
	
	for(int i;i<neighbours.size();i++){
		if(isExciting[neighbours[i]]){
			it = itTotalVector[neighbours[i]];
			itVector.push_back(it);
			rekursivRemoveClusterElement(posIndicesExciting, itVector);
		}
	}
	it = posIndicesExciting.begin();
	while(it != posIndicesExciting.end()){
		//std::cerr<<"it="<<*it<<" it0="<<index0<<" itBeg="<<*posIndicesExciting.begin()<<" itEnd="<<*posIndicesExciting.end()<<" size="<<posIndicesExciting.size()<<" distance="<<grid.getDistance(index0, *it)<<" vectorSize="<<itVector.size()<<std::endl;
		if(grid.getDistance(index0, *it) <= maxDistance){
			itVector.push_back(it);
			rekursivRemoveClusterElement(posIndicesExciting, itVector);
			it = itVector.back();
			itVector.pop_back();
		}
		else
			++it;
	}
}

double OP_cluster::calc(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep)const{
	std::list<int> posIndicesExciting;
	for(int i=0;i<grid.getn();i++){
		if(dVmdt[i] > grid.getdVmdtThresh(i) && y[i] > grid.getVmThresh(i)){
			posIndicesExciting.push_back(i);
			itTotalVector[i] = --posIndicesExciting.end();
			isExciting[i] = true;
		}
		else
			isExciting[i] = false;
	}
	int cluster = 0;
	std::vector<std::list<int>::iterator> itVector;
	itVector.push_back(posIndicesExciting.begin());
	std::list<int>::iterator it = posIndicesExciting.begin();
	while(it != posIndicesExciting.end()){
		rekursivRemoveClusterElement(posIndicesExciting, itVector);
		it = itVector.back();
		cluster++;
	}
	std::cerr<<"cluster="<<cluster<<std::endl;
	return cluster;
}
 
