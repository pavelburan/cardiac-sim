#include "obs_notterminated.h"
#include "../system.h"
#include "../efields/efield.h"
#include "../grids/grid.h"
#include <limits>

Obs_notTerminated::Obs_notTerminated(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),
														tSave(0.0),tSaveVecRaw(),tSaveVec(),fracNotTerminated(0.0),fracNotTerminatedVec(),isTerminated(false),isTerminatedVec(){
}

void Obs_notTerminated::readParams(){
	//tSave = system.gette();
	tSave = efield.gettend();
	cfg.readInto(tSave, "tSave");
	cfg.readRawInto(tSaveVecRaw, "tSaveVec");
}

void Obs_notTerminated::init(){
	cfg.stringAsVector(tSaveVec, tSaveVecRaw);
	fracNotTerminatedVec.resize(tSaveVec.size(), 0.0);
	isTerminatedVec.resize(tSaveVec.size(), false);
	if(system.getIsTimeCoupledByPulse()){
		cfg.printUp("Zeiten hängen vom Puls ab...");
		tSave = tSave + system.gett0();
		cfg.printVar(tSave, "tSave");
		tSaveVec[0] = tSaveVec[0] + system.gett0(); 
		tSaveVecRaw="[" + std::to_str(tSaveVec[0]);
		for(int i=1;i<tSaveVec.size();i++){
			tSaveVec[i] = tSaveVec[i] + system.gett0(); 
			tSaveVecRaw += "," + std::to_str(tSaveVec[i]);
		}
		tSaveVecRaw += "]";
		cfg.printVar(tSaveVecRaw, "tSaveVec");
		cfg.printDown("ready");
	}
	
}

bool Obs_notTerminated::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt_prev, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
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
	for(int i=0;i<tSaveVec.size();i++){
		if(t > tSaveVec[i]-Eps::t() && t < tSaveVec[i]+Eps::t()){
			std::vector<bool> isPointOfFraction(grid.getn(), true);
			#pragma omp parallel for //simd
			#pragma ivdep
			for(int j=0;j<grid.getn();j++){
				isPointOfFraction[j] = (dVmdt[j] > grid.getdVmdtThresh(j) && y[j] > grid.getVmThresh(j))? true : false;
			}
			fracNotTerminatedVec[i] = grid.getFraction(isPointOfFraction);
			std::cerr<<"fracNotTerminated["<<i<<"]="<<fracNotTerminatedVec[i]<<std::endl;
			if(fracNotTerminated > 0.0-Eps::rel() && fracNotTerminated < 0.0+Eps::rel())
				isTerminatedVec[i] = true;
		}
	}
	
	return false;
}

void Obs_notTerminated::finalize(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt_prev, double *__restrict__ dVmdt, double *__restrict__ temp, double t){
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
	
	if(tSaveVec.size() > 0){
		fileName = cfg.getPlotFolderSubTimeSavePrefixFileName("notTerminatedVec.txt");
		file.open( fileName.c_str(), std::ios::app );
		if( file.is_open() ){
			file<<cfg.getSubRepeatIndex();
			for(int i=0;i<tSaveVec.size();i++){
				file<<" "<<fracNotTerminatedVec[i];
			}
			file<<"\n";
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
