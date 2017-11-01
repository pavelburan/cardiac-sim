#include "obs_plotatfixedtimes.h"
#include "../system.h"
#include "../grids/grid.h"
#include "../efields/efield.h"
#include <algorithm>
#include <limits>

Obs_plotAtFixedTimes::Obs_plotAtFixedTimes(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),forEverySubTimeIndex(true),saveAll(false),timeSteps(),nextIndex(0){
}

void Obs_plotAtFixedTimes::readParams(){
	cfg.readInto(forEverySubTimeIndex, "forEverySubTimeIndex");
	cfg.readInto(saveAll, "saveAll");
	cfg.readIntoVector(timeSteps, "timeSteps");
}

void Obs_plotAtFixedTimes::init(){
	if(forEverySubTimeIndex){
		double delta_t = cfg.getSubTimeIndex()*efield.gettPeriod();
		for(int i=0;i<timeSteps.size();i++)
			timeSteps[i] += delta_t;
	}
	/*for(int i=0;i<timeSteps.size();i++)
		std::cerr<<timeSteps[i]<<std::endl;
	}*/
	timeSteps.push_back(std::numeric_limits<double>::infinity());
	std::sort(timeSteps.begin(), timeSteps.end());
	for(int i=0;i<timeSteps.size();i++){
		if(timeSteps[i] > system.gettb())
			nextIndex = i;
			break;
	}	
}

bool Obs_plotAtFixedTimes::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	if(t > timeSteps[nextIndex] - Eps::t()){
		std::string fileName = cfg.getPlotFolderFileName(std::string("(@SR)(@SC)(@ST)y_fixed_") + std::to_str(nextIndex) + std::string(".bin"));
		system.saveState(y, fileName, true);
		std::cerr<<fileName<<std::endl;
		cfg.addCleanFile(fileName, 2);
		nextIndex += 1;
	}
	return false;
}


