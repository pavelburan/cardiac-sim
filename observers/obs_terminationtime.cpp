#include "obs_terminationtime.h"
#include "../system.h"
#include "../efields/efield.h"
#include "../grids/grid.h"
#include <limits>

Obs_terminationTime::Obs_terminationTime(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),
														finishAfterTermination(false),terminationTime(std::numeric_limits<double>::max()),timeAfterTermination(0.0){
}

void Obs_terminationTime::readParams(){
	cfg.readInto(finishAfterTermination, "finishAfterTermination");
	cfg.readInto(timeAfterTermination, "timeAfterTermination");
}

/*bool Obs_terminationTime::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	double eps = std::numeric_limits< double >::epsilon();
	std::vector<bool> isPointOfFraction(grid.getn(), true);
	#pragma omp parallel for
	#pragma ivdep
	for(int i=0;i<grid.getn();i++){
		isPointOfFraction[i] = (dVmdt[i] > grid.getdVmdtThresh(i) && y[i] > grid.getVmThresh(i))? true : false;
	}
	if(grid.getFraction(isPointOfFraction) < eps && terminationTime > system.gette()){
		terminationTime = t;
		if(finishAfterTermination)
			return true;
	}
	return false;
}*/

bool Obs_terminationTime::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	double eps = std::numeric_limits< double >::epsilon();
	if(terminationTime > system.gette()){
		std::vector<bool> isPointOfFraction(grid.getn(), true);
		#pragma omp parallel for
		#pragma ivdep
		for(int i=0;i<grid.getn();i++){
			isPointOfFraction[i] = (dVmdt[i] > grid.getdVmdtThresh(i) && y[i] > grid.getVmThresh(i))? true : false;
		}
		if(grid.getFraction(isPointOfFraction) < eps)
			terminationTime = t;
	}
	if(finishAfterTermination && (t-terminationTime) > timeAfterTermination - Eps::t())
		return true;
	else
		return false;
}

void Obs_terminationTime::finalize(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t){
	std::string fileName = cfg.getPlotFolderSubTimeSavePrefixFileName(cfg.getSubConfigPrefix()+std::string("terminationTime.txt"));
	std::ofstream file( fileName.c_str(), std::ios::app );
	if( file.is_open() ){
		file<<cfg.getSubRepeatIndex()<<" "<<terminationTime<<"\n";
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
		exit(1);
	}
	file.close();
	cfg.addCleanFile(fileName, 3);
}
