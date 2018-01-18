#include "obs_stimulus.h"
#include "../system.h"
#include "../grids/grid.h"

Obs_stimulus::Obs_stimulus(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),t0(0.0),dt(0.0),te(0.0),xi0(0.0),dxi(0.0),eta0(0.0),deta(0.0),zeta0(0.0),dzeta(0.0),tLast(0.0),posIndices(){
}

void Obs_stimulus::readParams(){
	cfg.readInto(t0, "t0");	
	cfg.readInto(dt, "dt");
	cfg.readInto(te, "te");
	std::vector<double> region(6,0.0);
	cfg.readIntoVector(region, "region");
	xi0 = region[0];
	dxi = region[1];
	eta0 = region[2];
	deta = region[3];
	zeta0 = region[4];
	dzeta = region[5];
}

void Obs_stimulus::init(){
	double sys_tb = system.gettb();
	tLast = (t0>=sys_tb) ? (t0) : (sys_tb - fmod(sys_tb-t0,dt) + dt);
	posIndices = grid.getPosIndicesVolume(xi0, dxi, eta0, deta, zeta0, dzeta);
}

bool Obs_stimulus::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt_prev, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	if(t > tLast-Eps::t() && t < te+Eps::t()){
		grid.setExcitedState(y, posIndices);
		tLast += dt;
	}
	return false;
}


