#ifndef EFIELD_SQUARE_H
#define EFIELD_SQUARE_H
#include "efield.h"
class System;
class Efield_square : public Efield
{
public:
	//Konstruktoren
	Efield_square(const System& system, const std::string& configFileName, const std::string& keyPrefix="", double E0=0.0, double phi=0.0, double theta=0.0, double tbeg=0.0, double tend=0.0, double toff=0.0, double Eoff=0.0, double a=0.0, double b=0.0, double c=1.0);	
	//Destruktor
	~Efield_square(){}
private:
	Efield_square();
	Efield_square(const Efield_square &square);
	Efield_square &operator=(const Efield_square &square);
public:
	//Initialisieren
	void readEfieldParams();
	//Stärke des E-Feldes zurückgeben 
	double E(double t)const;
	//Feedback
	void feedback(double newtbeg, double newtend, double newtPeriod);
	
	//Zugriffsmethoden
	double gettPeriod()const{return a+b+c;}
	double geta()const{return a;}
	double getb()const{return b;}
	double getc()const{return c;}
	void seta(double newa){a = newa;}
	void setb(double newb){b = newb;}
	void setc(double newc){c = newc;}
protected:
	double a;
	double b;
	double c;
};

#endif //EFIELD_SQUARE_H

