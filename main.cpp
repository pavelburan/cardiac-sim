#include <iostream>
#include <string.h>
#include <omp.h>
#ifdef PETSC
#include <petscksp.h>
#endif
#include "matrvec/matrvec.h"
#include "configuration.h"
#include "system.h"
#include "timeintegration.h"

int main(int argc, char *argv[]){
	//clustertest
	Eps::setMaxUlps(4);
	Eps::setRelEps(1e-6);
	Configuration::loadConfigParams(argc, argv, "config.py", "  ");
	Configuration cfg(Configuration::getMainConfigFileName());

	//Parameter einlesen//
	System sys(Configuration::getMainConfigFileName());
	sys.readTimeIntegrationParams();
	sys.readPlotParams();
	sys.readEfieldParams();
	sys.readGridParams();
	sys.readObserverParams();
	
	//System initialsieren//
	sys.initSystem();
	
	//Efield initialsieren//
	sys.initEfield();
	
	//Gitter initialisieren//
	sys.initGrid();
	
	//Observer initialisieren//
	sys.initObservers();
	
	//Zeitintegration:
	std::string timeIntScheme = sys.getTimeIntScheme();
	void (*timeInt)(System &sys);
	if(timeIntScheme == "adaptive")
		timeInt = timeIntMono_adaptive;
	else if(timeIntScheme == "nonadaptive")
		timeInt = timeIntMono_nonadaptive;
	
	#ifdef PETSC
	PetscErrorCode ierr;
	PetscInitialize(&argc,&argv,(char*)0,(char*)0);
	#endif
	
	timeInt(sys);
	
	#ifdef PETSC
	ierr = PetscFinalize();
	#endif
	
	std::cout<<"subRepeatPrefix="<<cfg.getSubRepeatIndex()<<" subConfigIndex="<<cfg.getSubConfigIndex()<<" subTimeIndex="<<cfg.getSubTimeIndex()<<" finished :-)"<<std::endl;
	
	return 0;
	
}

