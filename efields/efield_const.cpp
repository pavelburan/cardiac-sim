#include "efield_const.h"

Efield_const::Efield_const(const System& system, const std::string& configFileName, const std::string& keyPrefix, double E0, double phi, double theta, double tbeg, double tend, double toff, double Eoff)
													:Efield(system,configFileName,keyPrefix,E0,phi,theta,tbeg,tend,toff,Eoff){
}	

void Efield_const::readEfieldParams(){
	Efield::readEfieldParams();
}

double Efield_const::E(double t)const{
	return (t >= tbeg-Eps::t() && t <= tend+Eps::t()) ? E0 : 0.0 ;
}

void Efield_const::feedback(double newtbeg, double newtend, double newtPeriod){
	tbeg = newtbeg;
	tend = newtend; 
}
