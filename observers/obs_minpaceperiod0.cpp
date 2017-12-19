#include "obs_minpaceperiod0.h"
#include "../system.h"
#include "../efields/efield.h"
#include "../grids/grid.h"
#include <limits>

Obs_minPacePeriod0::Obs_minPacePeriod0(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),
														startT(0.0),dT(0.0),dtCheck(0.0),tCheck(0.0),Ts(),T(0.0){
}

void Obs_minPacePeriod0::readParams(){
	cfg.readInto(startT, "startT");
	cfg.readInto(dT, "dT");
	cfg.readInto(dtCheck, "dtCheck");
}

void Obs_minPacePeriod0::init(){
	T = startT;
	Ts.push_back(T);
	tCheck = system.gett0() + T + dtCheck;
	efield.feedback(system.gett0(), tCheck, T);
}

bool Obs_minPacePeriod0::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	if(t >= tCheck-Eps::t()){
		std::vector<bool> isPointOfFraction(grid.getn(), false);
		#pragma omp parallel for
		#pragma ivdep
		for(int i=0;i<grid.getn();i++){
			isPointOfFraction[i] = (y[i] >= grid.getVmThresh(i))? true : false;
		}
		if(grid.getFraction(isPointOfFraction) >= 0.99){
			T -= dT;
			Ts.push_back(T);
			tCheck = t + T + dtCheck;
			efield.feedback(t, tCheck, T);
			std::vector<int> posIndices(system.getn());
			for(int i=0;i<posIndices.size();i++)
				posIndices[i] = i;
			grid.setRestingState(y, posIndices);
		}
		else
			return true;
	}
	cfg.print("T=" + std::to_str(T));
	return false;
}

void Obs_minPacePeriod0::finalize(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t){
	std::string fileName = cfg.getPlotFolderSubTimeSavePrefixFileName("minPacePeriod0.txt");
	std::ofstream file( fileName.c_str(), std::ios::app );
	if( file.is_open() ){
		T = -1.0;
		if(Ts.size() > 1)
			T = Ts[Ts.size()-2];
		file<<cfg.getSubRepeatIndex()<<" "<<T<<"\n";
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
		exit(1);
	}
	file.close();
	cfg.addCleanFile(fileName, 3);
}


