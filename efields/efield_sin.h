#ifndef EFIELD_SIN_H
#define EFIELD_SIN_H
#include "efield.h"
class System;
class Efield_sin : public Efield
{
public:
	//Konstruktoren
	Efield_sin(const System& system, const std::string& configFileName, const std::string& keyPrefix="", double E0=0.0, double phi=0.0, double theta=0.0, double tbeg=0.0, double tend=0.0, double toff=0.0, double Eoff=0.0, double tPeriod=0.0);
	//Destruktor
	~Efield_sin(){}
private:
	Efield_sin();
	Efield_sin(const Efield_sin &efield);
	Efield_sin &operator=(const Efield_sin &efield);
public:
	//Initialisieren
	void readEfieldParams();
	//Stärke des E-Feldes zurückgeben 
	double E(double t)const;
	//Feedback
	void feedback(double newtbeg, double newtend, double newtPeriod);
	
	//Zugriffsmethoden
	double gettPeriod()const{return tPeriod;}
	void settPeriod(double newtPeriod){tPeriod = newtPeriod;}
protected:
	double tPeriod;
};
#endif //EFIELD_SIN_H

