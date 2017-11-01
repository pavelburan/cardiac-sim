#include "obs_notterminated.h"
#include "../system.h"
#include "../efields/efield.h"
#include "../grids/grid.h"
#include <limits>

Obs_notTerminated::Obs_notTerminated(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),
														runAllPulses(false),tSave(0.0),fracNotTerminated(0.0),isTerminated(false){
}

void Obs_notTerminated::readParams(){
	cfg.readInto(runAllPulses, "runAllPulses");
	tSave = system.gette();
	cfg.readInto(tSave, "tSave");
}

void Obs_notTerminated::init(){
	if(system.getIsTimeCoupledByPulse()){
		cfg.printUp("Zeiten hängen vom Puls ab...");
		tSave = tSave + system.gett0();
		cfg.printVar(tSave, "tSave");
		cfg.printDown("ready");
	}
	
	if(cfg.getSubTimeIndex() > 0){
		int subRepeatIndex = -100;
		double temp = 1.0;
		std::string fileName = cfg.getPlotFolderSubTimeSavePrefixFileName("notTerminated.txt",-1);
		std::ifstream file( fileName.c_str() );
		if( file.is_open() ){
			while(file >> subRepeatIndex && file >> temp)
				if(subRepeatIndex == cfg.getSubRepeatIndex())
					break;
		}
		else{
			std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"File "<<fileName<<" could not be found"<<std::endl;
			exit(1);
		}
		file.close();
		if(temp > 0.0-Eps::rel() && temp < 0.0+Eps::rel()){
			fracNotTerminated = 0.0;
			isTerminated = true;
		}	
	}
}

bool Obs_notTerminated::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	if(t > tSave-Eps::t() && t < tSave+Eps::t()){
		std::vector<bool> isPointOfFraction(grid.getn(), true);
		#pragma omp parallel for //simd
		#pragma ivdep
		for(int i=0;i<grid.getn();i++){
			isPointOfFraction[i] = (dVmdt[i] > grid.getdVmdtThresh(i) && y[i] > grid.getVmThresh(i))? true : false;
		}
		fracNotTerminated = grid.getFraction(isPointOfFraction);
		std::cerr<<"fracNotTerminated="<<fracNotTerminated<<std::endl;
		if(fracNotTerminated > 0.0-Eps::rel() && fracNotTerminated < 0.0+Eps::rel())
			isTerminated = true;
	}
	
	if(isTerminated && !runAllPulses){
		cfg.print("Fibrillation is terminated!");
		return true;
	}
	else
		return false;
}

void Obs_notTerminated::finalize(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t){
	std::string fileName = cfg.getPlotFolderSubTimeSavePrefixFileName("notTerminated.txt");
	std::ofstream file( fileName.c_str(), std::ios::app );
	if( file.is_open() ){
		file<<cfg.getSubRepeatIndex()<<" "<<fracNotTerminated<<"\n";
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
		exit(1);
	}
	file.close();
	cfg.addCleanFile(fileName, 3);
}
