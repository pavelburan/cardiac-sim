#include "efield_square.h"
#include "../fltcomp.h"
#include <math.h>
Efield_square::Efield_square(const System& system, const std::string& configFileName, const std::string& keyPrefix, double E0, double phi, double theta, double tbeg, double tend, double toff, double Eoff, double a, double b, double c)
													:Efield(system,configFileName,keyPrefix,E0,phi,theta,tbeg,tend,toff,Eoff),a(a),b(b),c(c){
}	
	
void Efield_square::readEfieldParams(){
	Efield::readEfieldParams();
	cfg.readInto(a, "a");
	cfg.readInto(b, "b");
	cfg.readInto(c, "c");
}

double Efield_square::E(double t)const{
	if(t >= tbeg-Eps::t() && t <= tend+Eps::t()){
		double temp = fmod(t-tbeg+toff,a+b+c);
		if(temp <= a+Eps::t())
			return E0 + Eoff;
		else if(temp <= a+b+Eps::t())
			return -1.0*E0 + Eoff;
		else
			return Eoff;
	}
	else
		return 0.0;
}

void Efield_square::feedback(double newtbeg, double newtend, double newtPeriod){
	tbeg = newtbeg;
	tend = newtend;
	/*double scalFac = newtPeriod/(a + b + c);
	a *=scalFac;
	b *=scalFac;
	c *=scalFac;*/
	c = newtPeriod - a - b;
}
