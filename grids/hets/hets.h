#ifndef HETS_H
#define HETS_H
#include "../../configuration.h"
#include <vector>
class Grid;
class Hets{
public:
	//Konstruktor
	Hets(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix="", double frac=0.0);
	//Destruktor
	virtual ~Hets(){};
private:
	Hets();
	Hets(const Hets &hets);
	Hets &operator=(const Hets &hets);
public:
	const Configuration& getCfg()const{return cfg;}
	const Grid& getGrid()const{return grid;}
	double getfrac()const{return frac;}

	virtual void readHetsParams(){cfg.readInto(frac, "frac");}
	virtual void initHets(){};
	virtual void setisHet(bool *isHet)const{};
	
	//Hetscreator
	static Hets* newHets(const std::string& hetsType, const Grid& grid, const std::string& configFileName, const std::string& keyPrefix);
protected:
	//Configuration
	Configuration cfg;
	//System
	const Grid& grid;
	double frac;
};
#endif //HETS_H

