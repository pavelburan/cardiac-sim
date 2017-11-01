#ifndef EFIELD_CONST_H
#define EFIELD_CONST_H
#include "efield.h"
class System;
class Efield_const : public Efield
{
public:
	//Konstruktoren
	Efield_const(const System& system, const std::string& configFileName, const std::string& keyPrefix="", double E0=0.0, double phi=0.0, double theta=0.0, double tbeg=0.0, double tend=0.0, double toff=0.0, double Eoff=0.0);	
	//Destruktor
	~Efield_const(){}
private:
	Efield_const();
	Efield_const(const Efield_const &efield);
	Efield_const &operator=(const Efield_const &efield);
public:
	//Initialisieren
	void readEfieldParams();
	//Stärke des E-Feldes zurückgeben 
	double E(double t)const;
	//Feedback
	void feedback(double newtbeg, double newtend, double newtPeriod);
	//Zugriffsmethoden
	double gettPeriod()const{return tend-tbeg;}
	
protected:
	
};

#endif //EFIELD_CONST_H

