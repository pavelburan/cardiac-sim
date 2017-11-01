#include "observer.h"
#include "obs_stimulus.h"
#include "obs_stimulusrandom.h"
#include "obs_velocity.h"
#include "obs_minperiod.h"
#include "obs_activationtime.h"
#include "obs_observepoints.h"
#include "obs_notterminated.h"
#include "obs_terminationtime.h"
#include "obs_ordparamtimeseries.h"
#include "obs_plotatfixedtimes.h"
#include "../system.h"
#include <iostream>

Observer* Observer::newObserver(const std::string& observerType, System& system, const std::string& configFileName, const std::string& keyPrefix){
	if(observerType == "Stimulus")
		return new Obs_stimulus(system, configFileName, keyPrefix);
	if(observerType == "StimulusRandom")
		return new Obs_stimulusRandom(system, configFileName, keyPrefix);
	if(observerType == "Velocity")
		return new Obs_velocity(system, configFileName, keyPrefix);
	if(observerType == "MinPeriod")
		return new Obs_minPeriod(system, configFileName, keyPrefix);
	if(observerType == "ActivationTime")
		return new Obs_activationTime(system, configFileName, keyPrefix);
	if(observerType == "ObservationPoints")
		return new Obs_observePoints(system, configFileName, keyPrefix);
	if(observerType == "NotTerminated")
		return new Obs_notTerminated(system, configFileName, keyPrefix);
	if(observerType == "TerminationTime")
		return new Obs_terminationTime(system, configFileName, keyPrefix);
	if(observerType == "OrderParameterTimeSeries")
		return new Obs_ordParamTimeSeries(system, configFileName, keyPrefix);
	if(observerType == "PlotAtFixedTimes")
		return new Obs_plotAtFixedTimes(system, configFileName, keyPrefix);
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"observerType="<<observerType<<" existiert nicht!"<<std::endl;
		exit(1);
		return NULL;
	}
}

Observer::Observer(System& system, const std::string& configFileName, const std::string& keyPrefix):system(system),efield(system.getEfield()),grid(system.getGrid()),cfg(configFileName,keyPrefix){
}
