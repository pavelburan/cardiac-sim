#include "grid.h"
#include "grid_fdm2dcartesian.h"
#include "grid_fdm3dcartesian.h"
//#include "bessel.h"
#include "../configuration.h"
#include "../system.h"
#include "../models/model.h"
#include <iostream>
#include <stdlib.h>

Grid* Grid::newGrid(const std::string& gridType, const System& system, const std::string& configFileName, const std::string& keyPrefix){
	if(gridType == "FDM2DCartesian")
		return new Grid_FDM2DCartesian(system, configFileName, keyPrefix);
	else if(gridType == "FDM3DCartesian")
		return new Grid_FDM3DCartesian(system, configFileName, keyPrefix);
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"gridType="<<gridType<<" existiert nicht!"<<std::endl;
		exit(1);
		return NULL;
	}
}

Grid::Grid(const System& system, const std::string& configFileName, const std::string& keyPrefix):system(system),cfg(configFileName,keyPrefix),VmThresh(NULL),dVmdtThresh(NULL){
}

Grid::~Grid(){
	if(VmThresh)
		delete[] VmThresh;	
	if(dVmdtThresh)
		delete[] dVmdtThresh;
}

void Grid::initGrid(){
	VmThresh = new double[getn()];
	dVmdtThresh = new double[getn()];
	for(int i=0;i<getn();i++){
		VmThresh[i] = getModel(i).getVmThresh();
		dVmdtThresh[i] = getModel(i).getdVmdtThresh();
	}
}

/*bool Grid::isExcitable(const double *__restrict__ y, const double *__restrict__ dVmdt, int posIndex)const{
	if(y[posIndex] > VmThresh[posIndex] || dVmdt[posIndex] > dVmdtThresh[posIndex])
		return false;
	else{
		double Vm_check = VmThresh[posIndex]*0.0;
		double dVmdt_check;
		std::cerr<<Vm_check<<std::endl;
		getModel(posIndex).f_Vm(&Vm_check, &y[system.getVecIndex(posIndex,1)], &dVmdt_check, 0.0, system.getStateVarStride());
		std::cerr<<dVmdt_check<<std::endl;
		std::cerr<<y[system.getVecIndex(posIndex,2)]<<" "<<y[system.getVecIndex(posIndex,1)+system.getStateVarStride()]<<std::endl;
		return dVmdt_check > dVmdtThresh[posIndex];
	}
}
*/
