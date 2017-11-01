#include "obs_ordparamtimeseries.h"
#include "orderparameter.h"
#include "../system.h"
#include "../efields/efield.h"
#include "../grids/grid.h"
#include <limits>

Obs_ordParamTimeSeries::Obs_ordParamTimeSeries(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),
														types(),orderParameters(),dt(0.0),t_next(0.0),timeSteps(),data(){
}

void Obs_ordParamTimeSeries::readParams(){
	cfg.readIntoVector(types, "orderParameters");
	orderParameters.resize(types.size());
	data.resize(types.size());
	for(int i=0;i<orderParameters.size();i++)
		orderParameters[i] = system.getOrderParameter(types[i]);
	
	cfg.readInto(dt, "dt");
	
}

void Obs_ordParamTimeSeries::init(){
	if(system.resume()){
		std::string fileName = cfg.getPlotFolderSubRepeatResumePrefixFileName("ordParamTimeSeries.bin");
		loadOrderParameters(fileName);
		double tb = system.gettb();
		int size = timeSteps.size();
		for(int i=0;i<size;i++){
			if(timeSteps[i] > tb-Eps::t()){
				if(timeSteps[i] < tb+Eps::t()){
					size = i+1;
				}
				else
					size = i;
				break;
			}
		}
		timeSteps.resize(size);
		for(int i=0;i<data.size();i++)
			data[i].resize(size);
		t_next = timeSteps.back() + dt;	
	}
	else
		t_next = system.gettb();
	std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("ordParamTimeSeries.bin");
	cfg.addCleanFile(fileName, 2);
}

bool Obs_ordParamTimeSeries::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	if(t > t_next - Eps::t()){
		timeSteps.push_back(t_next);
		#pragma omp parallel for
		for(int i=0;i<orderParameters.size();i++)
			data[i].push_back(orderParameters[i]->calc(y_prev, y, dVmdt, temp, t, timeStep));
		t_next += dt;
	}
		
	if(isResumeAbleStep){
		std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("ordParamTimeSeries.bin");
		saveOrderParameters(fileName);
	}
	return false;
}

void Obs_ordParamTimeSeries::finalize(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t){
	std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("ordParamTimeSeries.bin");
	saveOrderParameters(fileName);
}

void Obs_ordParamTimeSeries::loadOrderParameters(const std::string& fileName){
	std::ifstream file( fileName.c_str(), std::ios::in | std::ios::binary );
	if( file.is_open() ){
		int loadSize;
		file.seekg(0, file.end);
		loadSize = (file.tellg()/sizeof(double))/(data.size()+1);
		file.seekg(0, file.beg);
		timeSteps.resize(loadSize);
		file.read(reinterpret_cast<char*>(timeSteps.data()),sizeof(double)*loadSize);
		for(int i=0;i<data.size();i++){
			data[i].resize(loadSize);
			file.read(reinterpret_cast<char*>(data[i].data()),sizeof(double)*loadSize);
		}
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be found"<<std::endl;
		exit(1);
	}
	file.close();	
}

void Obs_ordParamTimeSeries::saveOrderParameters(const std::string& fileName){
	std::ofstream file( fileName.c_str(), std::ios::out | std::ios::binary );
	if( file.is_open() ){
		file.write(reinterpret_cast<const char*>(timeSteps.data()),sizeof(double)*timeSteps.size());
		for(int i=0;i<data.size();i++){
			file.write(reinterpret_cast<const char*>(data[i].data()),sizeof(double)*data[i].size());
		}
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
		exit(1);
	}
	file.close();
}
