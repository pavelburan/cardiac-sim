#include "efield_sin.h"
#include "../fltcomp.h"
#include <math.h>

Efield_sin::Efield_sin(const System& system, const std::string& configFileName, const std::string& keyPrefix, double E0, double phi, double theta, double tbeg, double tend, double toff, double Eoff, double tPeriod)
													:Efield(system,configFileName,keyPrefix,E0,phi,theta,tbeg,tend,toff,Eoff),tPeriod(tPeriod){
}	
void Efield_sin::readEfieldParams(){
	Efield::readEfieldParams();
	cfg.readInto(tPeriod, "tPeriod");
}

double Efield_sin::E(double t)const{
	return (t >= tbeg-Eps::t() && t <= tend+Eps::t()) ? (E0*sin(2*M_PI/tPeriod*(t-tbeg+toff)) + Eoff) : 0.0 ;
}

void Efield_sin::feedback(double newtbeg, double newtend, double newtPeriod){
	tbeg = newtbeg;
	tend = newtend; 
	tPeriod = newtPeriod;
}
