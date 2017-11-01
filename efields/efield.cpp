#include "efield.h"
#include "efield_const.h"
#include "efield_sin.h"
#include "efield_square.h"
#include "../configuration.h"
#include <iostream>
#include <stdlib.h>

Efield* Efield::newEfield(const std::string& efieldType, const System& system, const std::string& configFileName, const std::string& keyPrefix){
	if(efieldType == "const")
		return new Efield_const(system, configFileName, keyPrefix);
	else if(efieldType == "sin")
		return new Efield_sin(system, configFileName, keyPrefix);
	else if(efieldType == "square")
		return new Efield_square(system, configFileName, keyPrefix);
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"efieldType="<<efieldType<<" existiert nicht!"<<std::endl;
		exit(1);
		return NULL;
	}
}

Efield::Efield(const System& system, const std::string& configFileName, const std::string& keyPrefix, double E0, double phi, double theta, double tbeg, double tend, double toff, double Eoff)
													:system(system),cfg(configFileName,keyPrefix),E0(E0),phi(phi),theta(theta),tbeg(tbeg),tend(tend),toff(toff),Eoff(Eoff){
}
	
void Efield::readEfieldParams(){
	cfg.readInto( E0, "E0");
	cfg.readInto( phi, "phi");
	cfg.readInto( theta, "theta");
	cfg.readInto( tbeg, "tbeg");
	cfg.readInto( tend, "tend");
	cfg.readInto( toff, "toff");
	cfg.readInto( Eoff, "Eoff");
}
	
void Efield::feedback(double t){
	tend = tend-tbeg+t; 
	tbeg = t;
}
