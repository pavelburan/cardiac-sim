#include "obs_minperiod.h"
#include "../system.h"
#include "../efields/efield.h"
#include "../grids/grid.h"
#include <limits>

Obs_minPeriod::Obs_minPeriod(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),
														startT1(0.0),dT1(),startT2(0.0),dT2(0.0),pos()/*{.0,.0,.0}*/,posIndex(0),Ts(),tSend(),tReceive(),T(0.0){
	for(int i=0;i<3;i++)
		pos[i] = 0.0;
}

void Obs_minPeriod::readParams(){
	cfg.readInto(startT1, "startT1");
	cfg.readInto(dT1, "dT1");
	cfg.readInto(startT2, "startT2");
	cfg.readInto(dT2, "dT2");
	std::vector<double> posVec;
	cfg.readIntoVector(posVec, "pos");
	for(int i=0;i<posVec.size();i++)	
		pos[i] = posVec[i];
}

void Obs_minPeriod::init(){
	posIndex = grid.getPosIndex(pos[0], pos[1], pos[2]);
	tSend.push_back(system.gettb());
	T = startT1;
	Ts.push_back(T);
}

bool Obs_minPeriod::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt_prev, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	double eps = timeStep/100.0;
	int size = tReceive.size();
	double thresh = grid.getVmThresh(posIndex);
	if(y_prev[posIndex]<thresh && y[posIndex]>=thresh){
		int size = tReceive.size();
		tReceive.push_back(t);
		size += 1;
		if(size > 2){
			if((tReceive[size-1]-tReceive[size-2]) > 1.5*(tReceive[size-2]-tReceive[size-3])){
				return true;
			}
		}
	}
	if(t > tSend.back()+T-eps){
		efield.feedback(t);
		tSend.push_back(t);
		if(T > startT2)
			T -= dT1;
		else{
			T -= dT2;
		}
		Ts.push_back(T);
	}
	if(size > 1)
		cfg.print("minT=" + std::to_str(Ts[size-2]) + " minTsend=" + std::to_str(tSend[size-1]-tSend[size-2]) + " minTreceive=" + std::to_str(tReceive[size-1]-tReceive[size-2]));
	else if(size > 0)
		cfg.print("tSend[0]=" + std::to_str(tSend[size-1]) + " tReceive[0]=" + std::to_str(tReceive[size-1]));
	return false;
}

void Obs_minPeriod::finalize(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt_prev, double *__restrict__ dVmdt, double *__restrict__ temp, double t){
	int size = tReceive.size();
	if(size < 3)
		cfg.print("Not finished");
	else{
		if((tReceive[size-1]-tReceive[size-2]) > 1.5*(tReceive[size-2]-tReceive[size-3])){
			cfg.print("Finished");
			size -= 1;
		}
		else{
			cfg.print("Not finished");
		}
	}
	if(size == 1)
		cfg.print("tSend[0]=" + std::to_str(tSend[size-1]) + " tReceive[0]=" + std::to_str(tReceive[size-1]));
	else if(size > 1)
		cfg.print("minT=" + std::to_str(Ts[size-2]) + " minTreceive=" + std::to_str(tReceive[size-1]-tReceive[size-2]) + " minTsend=" + std::to_str(tSend[size-1]-tSend[size-2]));
}


