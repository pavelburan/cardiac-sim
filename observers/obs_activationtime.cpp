#include "obs_activationtime.h"
#include "../system.h"
#include "../efields/efield.h"
#include "../grids/grid.h"
#include <limits>

Obs_activationTime::Obs_activationTime(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),
														ratio(0.0),tau(0.0){
}

void Obs_activationTime::readParams(){
	cfg.readInto(ratio, "ratio");
}

void Obs_activationTime::init(){
	tau = std::numeric_limits<double>::max();
}

bool Obs_activationTime::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	std::vector<bool> isPointOfFraction(grid.getn(), false);
	#pragma omp parallel for
	for(int i=0;i<grid.getn();i++){
		isPointOfFraction[i] = (y[i] >= grid.getVmThresh(i))? true : false;
	}
	if(grid.getFraction(isPointOfFraction) >= ratio){
		tau = t;
		return true;
	}
	return false;
}

void Obs_activationTime::finalize(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t){
	std::string fileName = cfg.getPlotFolderSubTimeSavePrefixFileName("activationTime.bin");
	std::ofstream file( fileName.c_str(), std::ios::app );
	if( file.is_open() ){
		file<<cfg.getSubRepeatIndex()<<" "<<tau<<"\n";
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
		exit(1);
	}
	file.close();
	cfg.addCleanFile(fileName, 3);	
}


