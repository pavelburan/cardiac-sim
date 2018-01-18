#include "obs_stimulusrandom.h"
#include "../system.h"
#include "../grids/grid.h"
#include <stdlib.h>

Obs_stimulusRandom::Obs_stimulusRandom(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),t0(0.0),dxi(0.0),deta(0.0),dzeta(0.0),n(0),posIndices(){
}

void Obs_stimulusRandom::readParams(){
	cfg.readInto(t0, "t0");	
	std::vector<double> size(3,0.0);
	cfg.readIntoVector(size, "size");
	dxi = size[0];
	deta = size[1];
	dzeta = size[2];
	cfg.readInto(n, "n");
}

void Obs_stimulusRandom::init(){
	posIndices.resize(n);
	for(int i=0;i<n;i++){
		int posIndex = rand() % grid.getn();
		posIndices[i] = grid.getPosIndicesVolume(posIndex, dxi, deta, dzeta);
	}
}

bool Obs_stimulusRandom::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt_prev, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	if(t > t0-Eps::t() && t < t0+Eps::t()){
		for(int i=0;i<n;i++)
			grid.setExcitedState(y, posIndices[i]);
	}
	return false;
}


