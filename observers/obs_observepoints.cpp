#include "obs_observepoints.h"
#include "../system.h"
#include "../efields/efield.h"
#include "../grids/grid.h"
#include <limits>

Obs_observePoints::Obs_observePoints(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),
														nRandObsPoints(0),fixedObsPoints(""),posIndices(),obsPointsData(){
}

void Obs_observePoints::readParams(){
	cfg.readInto(nRandObsPoints, "nRandObsPoints");
	cfg.readRawInto(fixedObsPoints, "fixedObsPoints");
}

void Obs_observePoints::init(){
	std::vector<double> vec;
	cfg.stringAsVector(vec, fixedObsPoints);
	int dim = grid.getDim();
	int nFixedObsPoints = vec.size() / (dim);
	posIndices.resize(nFixedObsPoints+nRandObsPoints);
	obsPointsData.resize(nFixedObsPoints+nRandObsPoints);
	std::string printString;
	for(int i=0; i<nFixedObsPoints; ++i){
		posIndices[i] = grid.getPosIndex(vec[dim*i], vec[dim*i+1], vec[dim*i+2]);
		printString = "obsPoint" + std::to_str(i) + "=[";
		for(int k=0;k<dim;k++)
			printString += std::to_str(vec[dim*i+k]) + ",";
		printString[printString.size()-1] = ']';
		cfg.print(printString);
	}
	if(system.resume()){
		std::string fileName = cfg.getPlotFolderSubRepeatResumePrefixFileName("observePoints.bin");
		loadObsPoints(fileName);
	}
	else{
		srand (time(NULL));
		for(int i=nFixedObsPoints;i<posIndices.size();i++){
			posIndices[i] = rand()%grid.getn();
		}
	}
	std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("observePoints.bin");
	cfg.addCleanFile(fileName, 1);
}

bool Obs_observePoints::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt_prev, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	#pragma omp parallel for
	for(int i=0;i<posIndices.size();i++)
		obsPointsData[i].push_back(y[posIndices[i]]);
	if(isResumeAbleStep){
		std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("observePoints.bin");
		saveObsPoints(fileName);
	}
	return false;
}

void Obs_observePoints::finalize(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt_prev, double *__restrict__ dVmdt, double *__restrict__ temp, double t){
	std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("observePoints.bin");
	saveObsPoints(fileName);
}

void Obs_observePoints::loadObsPoints(const std::string& fileName){
	std::ifstream file( fileName.c_str(), std::ios::in | std::ios::binary );
	if( file.is_open() ){
		double ind;
		int loadSize;
		int neededSize = system.getObsTimeLength();
		file.seekg(0, file.end);
		loadSize = (file.tellg()/sizeof(double)-posIndices.size())/posIndices.size();
		file.seekg(0, file.beg);
		for(int i=0;i<posIndices.size();i++){
			file.read(reinterpret_cast<char*>(&ind),sizeof(double));
			posIndices[i] = int(ind);
		}
		for(int i=0;i<posIndices.size();i++){
			obsPointsData[i].resize(loadSize);
			file.read(reinterpret_cast<char*>(obsPointsData[i].data()),sizeof(double)*loadSize);
			obsPointsData[i].resize(neededSize);
		}
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be found"<<std::endl;
		exit(1);
	}
	file.close();	
}

void Obs_observePoints::saveObsPoints(const std::string& fileName){
	std::ofstream file( fileName.c_str(), std::ios::out | std::ios::binary );
	if( file.is_open() ){
		double ind;
		for(int i=0;i<posIndices.size();i++){
			ind = double(posIndices[i]);
			file.write(reinterpret_cast<const char*>(&ind),sizeof(double));
		}
		for(int i=0;i<posIndices.size();i++){
			file.write(reinterpret_cast<const char*>(obsPointsData[i].data()),sizeof(double)*obsPointsData[i].size());
		}
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
		exit(1);
	}
	file.close();
}
