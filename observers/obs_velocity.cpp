#include "obs_velocity.h"
#include "../system.h"
#include "../grids/grid.h"
#include <limits>

Obs_velocity::Obs_velocity(System& system, const std::string& configFileName, const std::string& keyPrefix):Observer(system,configFileName,keyPrefix),
														regionBeg()/*{.0,.0,.0,.0,.0,.0}*/,regionEnd()/*{.0,.0,.0,.0,.0,.0}*/,posIndicesBeg(),posIndicesEnd(),tBeg(),tEnd(),nTotal(0),n(0){
	for(int i=0;i<6;i++){
		regionBeg[i] = 0.0;
		regionEnd[i] = 0.0;
	}
}

void Obs_velocity::readParams(){
	std::vector<double> region;
	cfg.readIntoVector(region, "regionBeg");
	for(int i=0;i<region.size();i++)	
		regionBeg[i] = region[i];
	cfg.readIntoVector(region, "regionEnd");
	for(int i=0;i<region.size();i++)	
		regionEnd[i] = region[i];
	
}

void Obs_velocity::init(){
	posIndicesBeg = grid.getPosIndicesVolume(regionBeg[0], regionBeg[1], regionBeg[2], regionBeg[3], regionBeg[4], regionBeg[5]);
	posIndicesEnd = grid.getPosIndicesVolume(regionEnd[0], regionEnd[1], regionEnd[2], regionEnd[3], regionEnd[4], regionEnd[5]);
	tBeg.assign(posIndicesBeg.size(), system.gett0()-1.0);
	tEnd.assign(posIndicesEnd.size(), system.gett0()-1.0);
	nTotal = posIndicesBeg.size() + posIndicesEnd.size();
	n = 0;
}

bool Obs_velocity::observe(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt_prev, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double timeStep, bool isResumeAbleStep){
	double t0 = system.gett0();
	#pragma omp parallel for
	for(int i=0;i<tBeg.size();i++){
		int posIndex = posIndicesBeg[i];
		double thresh = grid.getVmThresh(posIndex);
		if(y_prev[posIndex]<thresh && y[posIndex]>=thresh && tBeg[i]<t0 && t > t0){
			tBeg[i] = t;
			#pragma omp atomic
			++n;
		}
	}
	#pragma omp parallel for
	for(int i=0;i<tEnd.size();i++){
		int posIndex = posIndicesEnd[i];
		double thresh = grid.getVmThresh(posIndex);
		if(y_prev[posIndex]<thresh && y[posIndex]>=thresh && tEnd[i]<t0 && t > t0){
			tEnd[i] = t;
			#pragma omp atomic
			++n;
		}
	}
	if(n==nTotal)
		return true;
	else
		return false;

}

void Obs_velocity::finalize(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt_prev, double *__restrict__ dVmdt, double *__restrict__ temp, double t){
	double t0 = system.gett0();
	double tBegSum = 0.0;
	double tEndSum = 0.0;
	int ntBeg = 0;
	int ntEnd = 0;
	#pragma omp parallel for
	for(int i=0;i<tBeg.size();i++){
		if(tBeg[i] > t0){
			#pragma omp atomic
			tBegSum += tBeg[i];
			#pragma omp atomic
			ntBeg++;
		}
	}
	#pragma omp parallel for
	for(int i=0;i<tEnd.size();i++){
		if(tEnd[i] > t0){
			#pragma omp atomic
			tEndSum += tEnd[i];
			#pragma omp atomic
			ntEnd++;
		}
	}
	
	double distance = 0.0;
	#pragma omp parallel for
	for(int i=0;i<posIndicesEnd.size();i++){
		int posIndexEnd = posIndicesEnd[i];
		double distanceLocal = std::numeric_limits<double>::max();
		for(int j=0;j<posIndicesBeg.size();j++){
			distanceLocal = std::min(grid.getDistance(posIndicesBeg[j], posIndexEnd), distanceLocal);
		}
		#pragma omp atomic
		distance += distanceLocal;
	}
	distance /= posIndicesEnd.size();
		
	double velocity = distance/(tEndSum/ntEnd - tBegSum/ntBeg);
	cfg.print("ntBeg=" + std::to_str(ntBeg) + " tBegSum=" + std::to_str(tBegSum) + " tBeg=" + std::to_str(tBegSum/ntBeg));
	cfg.print("ntEnd=" + std::to_str(ntEnd) + " tEndSum=" + std::to_str(tEndSum) + " tEnd=" + std::to_str(tEndSum/ntEnd));
	cfg.print("v=dx/dt=" + std::to_str(distance) + "/" + std::to_str(tEndSum/ntEnd - tBegSum/ntBeg) + "=" + std::to_str(velocity));
	std::string fileName = cfg.getPlotFolderSubTimeSavePrefixFileName("velocity.txt");
	std::ofstream file( fileName.c_str(), std::ios::app );
	if( file.is_open() ){
		file.precision(16);
		file.setf(std::ios::scientific);
		file<<cfg.getSubRepeatIndex()<<" "<<velocity<<"\n";
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
		exit(1);
	}
	file.close();
	cfg.addCleanFile(fileName, 3);
}


