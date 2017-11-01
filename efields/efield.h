#ifndef EFIELD_H
#define EFIELD_H
#include "../configuration.h"
class System;
class Efield
{
public:
	//Konstruktoren
	Efield(const System& system, const std::string& configFileName, const std::string& keyPrefix="", double E0=0.0, double phi=0.0, double theta=0.0, double tbeg=0.0, double tend=0.0, double toff=0.0, double Eoff=0.0);
	//Destruktor
	virtual ~Efield(){}
private:
	Efield();
	Efield(const Efield &efield);
	Efield &operator=(const Efield &efield);
public:
	//Initialisieren
	virtual void readEfieldParams();
	virtual void initEfield(){}
	//Stärke des E-Feldes zurückgeben 
	virtual double E(double t) const = 0;
	//Feedback
	void feedback(double t);
	virtual void feedback(double newtbeg, double newtend, double newtPeriod) = 0;
	
	//Zugriffsmethoden
	const Configuration& getCfg()const{return cfg;}
	const System& getSystem()const{return system;}
	double getE0()const{return E0;}
	double getphi()const{return phi;}
	double gettheta()const{return theta;}
	double gettbeg()const{return tbeg;}
	double gettend()const{return tend;}
	double gettoff()const{return toff;}
	double getEoff()const{return Eoff;}
	virtual double gettPeriod()const = 0;
	void setE0(double newE0){E0 = newE0;}
	void setphi(double newphi){phi = newphi;}
	void settheta(double newtheta){theta = newtheta;}
	void settbeg(double newtbeg){tbeg = newtbeg;}
	void settend(double newtend){tend = newtend;}
	void settoff(double newtoff){Eoff = newtoff;}
	void setEoff(double newEoff){Eoff = newEoff;}
	
	//Efieldcreator
	static Efield* newEfield(const std::string& efieldType, const System& system, const std::string& configFileName, const std::string& keyPrefix);
	
protected:
//Configuration
	Configuration cfg;
//System
	const System& system;
//Efieldparameter
	double E0;
	double phi;
	double theta;
	double tbeg;
	double tend;
	double toff;
	double Eoff;
};
#endif //EFIELD_H

