#include "obs_minpaceperiod.h"
#include "../system.h"
#include "../efields/efield.h"
#include "../grids/grid.h"
#include <limits>

Obs_minPacePeriod::Obs_minPacePeriod(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),
														startT1(0.0),dT1(0.0),startT2(0.0),dT2(0.0),dtCheck(0.0),tPuls(0.0),tCheck(0.0),Ts(),T(0.0){
}

void Obs_minPacePeriod::readParams(){
	cfg.readInto(startT1, "startT1");
	cfg.readInto(dT1, "dT1");
	cfg.readInto(startT2, "startT2");
	cfg.readInto(dT2, "dT2");
	cfg.readInto(tCheck, "dtCheck");
}

void Obs_minPacePeriod::init(){
	T = startT1;
	Ts.push_back(T);
	tPuls = efield.gettbeg();
	tCheck = tPuls + dtCheck;
}

bool Obs_minPacePeriod::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	if(t >= tCheck-Eps::t()){
		std::vector<bool> isPointOfFraction(grid.getn(), false);
		#pragma omp parallel for
		#pragma ivdep
		for(int i=0;i<grid.getn();i++){
			isPointOfFraction[i] = (y[i] >= grid.getVmThresh(i))? true : false;
		}
		if(grid.getFraction(isPointOfFraction) >= 0.99){
			if(T > startT2)
				T -= dT1;
			else{
				T -= dT2;
			}
			tPuls += T;
			Ts.push_back(T);
			tCheck = tPuls + dtCheck;
			efield.feedback(tPuls);

			return false;
		}
		else
			return true;
	}
	cfg.print("T=" + std::to_str(T));
	return false;
}

void Obs_minPacePeriod::finalize(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t){
	std::string fileName = cfg.getPlotFolderSubTimeSavePrefixFileName("minPacePeriod.txt");
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


