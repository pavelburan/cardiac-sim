#include "model.h"
#include "model_bocf.h"
#include "model_luorudy.h"
#include "model_tnnp2006b.h"
#include "model_mitchschaeff.h"
#include "model_fentonkarma.h"
#include "../configuration.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
//omp//#include <omp.h>

Model* Model::newModel(const std::string& modelType, const Grid& grid, const std::string& configFileName, const std::string& keyPrefix){
	if(modelType == "BOCF")
		return new Model_BOCF(grid, configFileName, keyPrefix);
	else if(modelType == "LuoRudy")
		return new Model_LuoRudy(grid, configFileName, keyPrefix);
	else if(modelType == "TNNP2006b")
		return new Model_TNNP2006b(grid, configFileName, keyPrefix);
	else if(modelType == "MitchSchaeff")
		return new Model_MitchSchaeff(grid, configFileName, keyPrefix);
	else if(modelType == "FentonKarma")
		return new Model_FentonKarma(grid, configFileName, keyPrefix);
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"modelType="<<modelType<<" existiert nicht!"<<std::endl;
		exit(1);
		return NULL;
	}
}

Model::Model(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix, int nVars, int nHH, int nNL, double VmThresh, double dVmdtThresh, double VmMin, double VmMax)
											:grid(grid),cfg(configFileName,keyPrefix),nVars(nVars),nHH(nHH),nNL(nNL),VmThresh(VmThresh),dVmdtThresh(dVmdtThresh),VmMin(VmMin),VmMax(VmMax){
}

void Model::readModelParams(){
	cfg.readInto(VmThresh, "VmThresh");
	cfg.readInto(dVmdtThresh, "dVmdtThresh");
	cfg.readInto(VmMin, "VmMin");
	cfg.readInto(VmMax, "VmMax");
}
