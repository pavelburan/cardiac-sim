#include "hets.h"
#include "het_spheres.h"
#include "../../system.h"
#include "../grid.h"
#include <iostream>
#include <stdlib.h>

Hets* Hets::newHets(const std::string& hetsType, const Grid& grid, const std::string& configFileName, const std::string& keyPrefix){
	if(hetsType == "points")
		return new Hets(grid, configFileName, keyPrefix);
	else if(hetsType == "spheres")
		return new Het_spheres(grid, configFileName, keyPrefix);
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"hetsType="<<hetsType<<" existiert nicht!"<<std::endl;
		exit(1);
		return NULL;
	}
}

Hets::Hets(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix, double frac):grid(grid),cfg(configFileName,keyPrefix),frac(frac){
}
