#include "system.h"
#include "efields/efield.h"
#include "observers/observer.h"
#include "observers/orderparameter.h"
#include <iostream>
#include <fstream>
#include <limits>
#include <stdlib.h>
#include <time.h> 
#include <omp.h>

System::System(const std::string& configFileName, const std::string& keyPrefix):cfg(configFileName,keyPrefix),randomSeedNumber(0),t0(0.0),tb(0.0),te(0.0),isTimeCoupledByPulse(false),dt(0.0),timeIntScheme(),timeIntSchemeDiffusion(),VmVecOrdered(false),modelVarsVecOrdered(false),N(0),indexStride(0),varStride(0),firstVecOrderedVar(0),
														tp0(0.0),dtp(0.0),resumeIndex(-1),savePrevious(false),saveAll(false),efieldType(""),efield(NULL),gridType(""),grid(NULL),orderParameterTypes(),orderParameterTypeIndices(),orderParameters(),observerTypes(),observerTypeIndices(),observers(),
														obs_time(),obs_E_data(),pcount(0),tpLast(0.0),tpSubTime(0.0)
{
	cfg.printUp("Systemparameter werden eingelesen...");
	
	if( resume() ){
		/*std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName(cfg.getSubConfigPrefix()+std::string("randeomSeedNumber.txt"));
		std::ofstream file( fileName.c_str(), std::ios::app );
		if( file.is_open() ){
			file<<cfg.getSubRepeatIndex()<<" "<<terminationTime<<"\n";
		}
		else{
			std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
			exit(1);
		}
		file.close();*/
	}
	else{
		std::string initialRandomSeedFile("NO");
		if(initialRandomSeedFile == "NO"){
			randomSeedNumber = time(NULL) + 86400*(cfg.getSubRepeatIndex() + 1000*cfg.getSubConfigIndex());
		}
		else{
			loadRandomSeedNumber(initialRandomSeedFile);
		}
	}
	cfg.printVar(randomSeedNumber, "randomSeedNumber"); 
	srand(randomSeedNumber);
	std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName(cfg.getSubConfigPrefix()+std::string("randomSeedNumber.txt"));
	cfg.addCleanFile(fileName, 4);

	
	cfg.printDown("ready\n");
}													
	

System::~System(){
	if(grid)
		delete grid;
	if(efield)
		delete efield;
	for(int i=0;i<orderParameters.size();i++)
		delete orderParameters[i];
	for(int i=0;i<observers.size();i++){
		delete observers[i];
	}
}

void System::readTimeIntegrationParams(){
	cfg.printUp("Zeitintegrationsparameter werden eingelesen...");

	cfg.readInto(t0, "t0");
	cfg.readInto(te, "te");
	cfg.readInto(isTimeCoupledByPulse, "isTimeCoupledByPulse");
	cfg.readInto(dt, "dt");
	Eps::setAbsEps_t(Eps::rel()*dt);
	
	cfg.readInto(timeIntScheme, "timeIntScheme");
	if(timeIntScheme != "adaptive" && timeIntScheme != "nonadaptive"){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"timeIntScheme="<<timeIntScheme<<" existiert nicht!"<<std::endl;
		exit(1);
	}
	cfg.readInto(timeIntSchemeDiffusion, "timeIntSchemeDiffusion");
	if(timeIntSchemeDiffusion != "FE" && timeIntSchemeDiffusion != "BE" && timeIntSchemeDiffusion != "CN" && timeIntSchemeDiffusion != "BDF2" && timeIntSchemeDiffusion != "AM2" && timeIntSchemeDiffusion != "MS2"){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"timeIntSchemeDiffusion="<<timeIntSchemeDiffusion<<" existiert nicht!"<<std::endl;
		exit(1);
	}
	cfg.printDown("ready\n");
}

void System::readPlotParams(){
	cfg.printUp("Plotparameter werden eingelesen...");
	
	tpSubTime = 2*te;
	std::vector<int> subIndexShift(2,0);
	subIndexShift[1] += 1;
	cfg.readIntoShifted(subIndexShift, tpSubTime, "t0", "tpSubTime");
	cfg.readInto(tp0, "tp0");
	cfg.readInto(dtp, "dtp");
	cfg.readInto(resumeIndex, "resumeIndex");
	cfg.readInto(savePrevious, "savePrevious");
	cfg.readInto(saveAll, "saveAll");

	cfg.printDown("ready\n");
}

void System::readEfieldParams(){
	if(efield){
		delete efield;
	}
	cfg.printUp("Efieldparameter werden eingelesen...");

	cfg.readInto(efieldType, "efieldType");
	std::string efieldConfigFileName = cfg.getConfigFileName();
	cfg.readInto(efieldConfigFileName, "efieldConfigFileName");
	std::string efieldKeyPrefix = "";
	cfg.readInto(efieldKeyPrefix, "efieldKeyPrefix");

	efield = Efield::newEfield(efieldType, *this, efieldConfigFileName, efieldKeyPrefix);
	efield->readEfieldParams();
	
	cfg.printDown("ready\n");
}

void System::readGridParams(){
	if(grid)
		delete grid;
	cfg.printUp("Gitterparameter werden eingelesen...");
	
	cfg.readInto(gridType, "gridType");
	std::string gridConfigFileName = cfg.getConfigFileName();
	cfg.readInto(gridConfigFileName, "gridConfigFileName");
	std::string gridKeyPrefix = "";
	cfg.readInto(gridKeyPrefix, "gridKeyPrefix");
	
	grid = Grid::newGrid(gridType, *this, gridConfigFileName, gridKeyPrefix);
	grid->readGridParams();
	Eps::setAbsEps_x(Eps::rel()*grid->getminh());
	Eps::setAbsEps_x(Eps::rel());
	cfg.printDown("ready\n");
	
}

void System::readObserverParams(){
	cfg.printUp("Orderparameterparameter und Observerparameter werden eingelesen...");
	
	std::vector<std::string> orderParameterStringVec;
	cfg.readIntoVector(orderParameterStringVec, "orderParameters");
	orderParameters.resize(orderParameterStringVec.size()/3);
	orderParameterTypes.resize(orderParameters.size());
	for(int i=0;i<orderParameters.size();i++){
		orderParameterTypes[i] = orderParameterStringVec[3*i];
		orderParameterTypeIndices[orderParameterTypes[i]] = i;
		orderParameters[i] = OrderParameter::newOrderParameter(orderParameterStringVec[3*i], *this, orderParameterStringVec[3*i+1], orderParameterStringVec[3*i+2]);
		cfg.printUp("orderParameter"+std::to_str(i)+"="+orderParameterTypes[i]+"...");
		orderParameters[i]->readParams();
		cfg.printDown("ready");	
	}

	std::vector<std::string> observerStringVec;
	cfg.readIntoVector(observerStringVec, "observers");
	observers.resize(observerStringVec.size()/3);
	observerTypes.resize(observers.size());
	for(int i=0;i<observers.size();i++){
		observerTypes[i] = observerStringVec[3*i];
		observerTypeIndices[observerTypes[i]] = i;
		observers[i] = Observer::newObserver(observerStringVec[3*i], *this, observerStringVec[3*i+1], observerStringVec[3*i+2]);
		cfg.printUp("Observer"+std::to_str(i)+"="+observerTypes[i]+"...");
		observers[i]->readParams();
		cfg.printDown("ready");	
	}
		
	cfg.printDown("ready\n");
}

void System::initSystem(){
	cfg.printUp("System wird initialsiert...");

	if(isTimeCoupledByPulse){
		cfg.printUp("Zeiten hängen vom Puls ab...");
		double t0_old = t0;
		t0 = t0 + cfg.getSubTimeIndex()*efield->gettPeriod();
		tpSubTime = t0 + efield->gettPeriod();
		te = t0 + ( (te-t0_old > efield->gettPeriod()) ? te-t0_old : efield->gettPeriod() );
		efield->settbeg(efield->gettbeg()+cfg.getSubTimeIndex()*efield->gettPeriod());
		efield->settend(efield->gettbeg()+efield->gettPeriod());
		cfg.printVar(t0, "t0");
		cfg.printVar(te, "te");
		cfg.printVar(tpSubTime, "tpSubTime");
		cfg.printVar(efield->gettbeg(), "E_tbeg");
		cfg.printVar(efield->gettend(), "E_tend");
		cfg.printVar(efield->gettPeriod(), "E_tPeriod");
		cfg.printDown("ready");
	}
	tb = t0;
	if(resumeIndex >= 0)
		tb = tp0 + dtp*resumeIndex;
	
	cfg.printDown("ready\n");
}

void System::initEfield(){
	cfg.printUp("Efield wird initialsiert...");

	efield->initEfield();
	
	cfg.printDown("ready\n");
}

void System::initGrid(){
	cfg.printUp("Gitterparameter werden initialisiert...");
	
	grid->initGrid();

	cfg.printDown("ready\n");
	
	//todo
	setStateStructure(VmVecOrdered, modelVarsVecOrdered);
}

void System::initObservers(){
	cfg.printUp("Orderparameter und Observer werden initialisiert...");
	
	//1D-Plot
	if(resume()){
		std::string fileName = cfg.getPlotFolderSubRepeatResumePrefixFileName("1DPlot.bin");
		loadObs1DData(fileName);
		double tb = gettb();
		for(int i=0;i<obs_time.size();i++){
			if(obs_time[i] >= tb-Eps::t()){
				obs_time.resize(i+1);
				obs_E_data.resize(i+1);
				break;
			}
		}
	}
	std::string fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("1DPlot.bin");
	cfg.addCleanFile(fileName, 2);
	
	for(int i=0;i<orderParameters.size();i++){
		cfg.printUp("Orderparameter"+std::to_str(i)+"="+orderParameterTypes[i]+"...");
		orderParameters[i]->init();
		cfg.printDown("ready");	
	}
	
	for(int i=0;i<observers.size();i++){
		cfg.printUp("Observer"+std::to_str(i)+"="+observerTypes[i]+"...");
		if(resumeIndex >= 0 && !observers[i]->resumeAble()){
			std::cerr<<"Observer"<<i<<"="<<observerTypes[i]<<" is not able to resume!"<<std::endl;
			exit(1);
		}
		if(cfg.getSubTimeIndex() > 0 && !observers[i]->subTimeAble()){
			std::cerr<<"Observer"<<i<<"="<<observerTypes[i]<<" is not able to subTime!"<<std::endl;
			exit(1);
		}
		observers[i]->init();
		cfg.printDown("ready");	
	}
	
	//clean.sh
	cfg.createCleanSH();

	cfg.printDown("ready\n");
}

void System::calcInitialCondition(double *y0){
	cfg.printUp("Anfangsbedingungen werden festgelegt...");
	std::string fileName;
	//tb = t0;
	pcount = 0;
	tpLast = tp0;
	if(cfg.getSubTimeIndex() > 0){
		if(tp0+Eps::t() < tb){
			pcount = int((tb-tp0)/dtp)+1;
			tpLast = tb - fmod(tb-tp0,dtp) + dtp;
		}
		fileName = cfg.getPlotFolderSubRepeatResumePrefixFileName("y.bin");
		loadState(y0,fileName);
	}
	else if(resumeIndex >= 0){
		//tb = tp0 + dtp*resumeIndex;
		pcount = resumeIndex;
		tpLast = tb;
		if(savePrevious){
			bool onlyVm = saveAll ? false : true;
			for(int i=0; i < pcount; i++){
				fileName = cfg.getPlotFolderSubRepeatResumePrefixFileName("y_" + std::to_str(i) + ".bin");
				loadState(y0,fileName,onlyVm);
				fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("y_" + std::to_str(i) + ".bin");
				saveState(y0,fileName,onlyVm);
			}
		}
		fileName = cfg.getPlotFolderSubRepeatResumePrefixFileName("y_" + std::to_str(pcount) + ".bin");
		loadState(y0,fileName);
		fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("y_" + std::to_str(pcount) + ".bin");
		saveState(y0,fileName);
			
		pcount++;
		tpLast += dtp;	
	}
	else{
		std::string initialStateFile("NO");
		cfg.readInto(initialStateFile, "initialStateFile");
		if(initialStateFile == "NO"){
			std::vector<int> posIndices(getn());
			for(int i=0;i<posIndices.size();i++)
				posIndices[i] = i;
			grid->setRestingState(y0, posIndices);
		}
		else{
			fileName = cfg.getInitialFileName(initialStateFile);
			loadInitialState(y0,fileName);	
		}
	}

	std::vector<int> posIndices(getN()-getn(),0);
	for(int i=getn();i<getN();i++)
		posIndices[i-getn()] = i;
	grid->setRestingState(y0, posIndices);
	cfg.printDown("ready\n");
}

bool System::obs(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt_prev, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double dt){
	bool finished = false;
	std::string fileName;
	bool isPlotSaveStep = t >= tpLast-Eps::t() || t >= tpSubTime-Eps::t();
	
	//1DPlot
	obs_time.push_back(t);
	obs_E_data.push_back(efield->E(t));
	
	for(int i=0;i<observers.size();i++){
		finished |= observers[i]->observe(y_prev, y, dVmdt_prev, dVmdt, temp, t, dt, isPlotSaveStep);
	}

	if(isPlotSaveStep){
		//1D-Plot
		fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("1DPlot.bin");
		saveObs1DData(fileName);

		//2D-Plot
		if(t >= tpLast-Eps::t()){ 
			if(!saveAll && pcount > 0 && pcount > resumeIndex){
				fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("y_" + std::to_str(pcount-1) + ".bin");
				loadState(temp, fileName, true);
				saveState(temp, fileName, true);
			}
			fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("y_" + std::to_str(pcount) + ".bin");
			saveState(y, fileName);
			std::cout<<"save("<<pcount<<")"<<std::endl;
			pcount++;
			tpLast += dtp;
		}
		if(t >= tpSubTime-Eps::t()){
			fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("y.bin");
			saveState(y,fileName);
			std::cout<<"save(tpSubTime)"<<std::endl;
			tpSubTime = 2*te;
			cfg.addCleanFile(fileName, 0);
		}
	}
	return finished;
}

void System::obs_end(double *__restrict__ y_prev, double *__restrict__ y, double *__restrict__ dVmdt_prev, double *__restrict__ dVmdt, double *__restrict__ temp, double t, double dt){
	std::string fileName;
	
	cfg.printUp("Observers werden terminiert...");
	
	for(int i=0;i<observers.size();i++){
		cfg.printUp("Observer"+std::to_str(i)+"="+observerTypes[i]+"...");
		observers[i]->finalize(y_prev, y, dVmdt_prev, dVmdt, temp, t);
		cfg.printDown("ready");	
	}
	
	//1D-Plot
	fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("1DPlot.bin");
	saveObs1DData(fileName);
	
	//clean.sh
	cfg.createCleanSH();
	
	cfg.printDown("ready\n");
}

void System::loadState(double *y, const std::string fileName, bool onlyVm) const{
	std::ifstream file(fileName.c_str(), std::ios::in | std::ios::binary );
	if( file.is_open() ){
		int length = 0;
		file.seekg(0, file.end);
		length = file.tellg()/sizeof(double);
		file.seekg(0, file.beg);
		int nVars = grid->getMaxNumVars();
		//std::cerr<<"length="<<length<<" == getn()="<<getn()<<" && onlyVm="<<onlyVm<<" || length="<<length<<" == getn()*nVars="<<getn()*nVars<<std::endl;
		if(length == getn() && onlyVm || length == getn()*nVars){
			file.read(reinterpret_cast<char*>(y),sizeof(double)*getn());
			if(!onlyVm){
				if(VmVecOrdered)
					file.seekg(0, file.beg);
				int N = getN();
				int kStride = getStateIndexStride();
				int mStride = getStateVarStride();
				int m0 = getStateFirstVecOrderedVar();
				for(int m=m0;m<nVars;m++){
					for(int k=0;k<N-VECLEN;k+=VECLEN){
						file.read(reinterpret_cast<char*>(&y[N+k*kStride+(m-m0)*mStride]),sizeof(double)*VECLEN);
					}
					file.read(reinterpret_cast<char*>(&y[N+(N-VECLEN)*kStride+(m-m0)*mStride]),sizeof(double)*(VECLEN-N+getn()));
				}
			}
		}
		else{
			std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"File "<<fileName<<" has not the right size(length="<<length<<" != n="<<getn()<<",nVars="<<nVars<<")"<<std::endl;
			exit(1);
		}
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be found"<<std::endl;
		exit(1);
	}
	file.close();
}

void System::loadInitialState(double *y0, const std::string fileName) const{
	std::ifstream file(fileName.c_str(), std::ios::in | std::ios::binary );
	if( file.is_open() ){
		int initialN = 0;
		file.seekg(0, file.end);
		initialN = file.tellg()/sizeof(double);
		file.seekg(0, file.beg);
		double *temp = new double[initialN];
		file.read(reinterpret_cast<char*>(temp),sizeof(double)*initialN);
		std::cerr<<"loadInitialState test1"<<std::endl;
		grid->interpolateInitialStateOnGrid(y0, temp, initialN);
		std::cerr<<"loadInitialState test2"<<std::endl;
		delete[] temp;
		std::cerr<<"loadInitialState test3"<<std::endl;
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be found"<<std::endl;
		exit(1);
	}
	file.close();
}

void System::saveState(const double *y, const std::string fileName, bool onlyVm)const{
	std::ofstream file( fileName.c_str(), std::ios::out | std::ios::binary );
	if( file.is_open() ){
		int nVars = grid->getMaxNumVars();
		file.write(reinterpret_cast<const char*>(y),sizeof(double)*getn());
		if(!onlyVm){
			int N = getN();
			int kStride = getStateIndexStride();
			int mStride = getStateVarStride();
			int m0 = getStateFirstVecOrderedVar();
			for(int m=1;m<nVars;m++){
				for(int k=0;k<N-VECLEN;k+=VECLEN){
					file.write(reinterpret_cast<const char*>(&y[N+k*kStride+(m-m0)*mStride]),sizeof(double)*VECLEN);
				}
				file.write(reinterpret_cast<const char*>(&y[N+(N-VECLEN)*kStride+(m-m0)*mStride]),sizeof(double)*(VECLEN-N+getn()));
			}
		}
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
		exit(1);
	}
	file.close();
}

void System::loadObs1DData(const std::string& fileName){
	std::ifstream file( fileName.c_str(), std::ios::in | std::ios::binary );
	if( file.is_open() ){
		file.seekg(0, file.end);
		int size = file.tellg()/sizeof(double)/2;
		file.seekg(0, file.beg);
		obs_time.resize(size);
		file.read(reinterpret_cast<char*>(obs_time.data()),sizeof(double)*size);
		obs_E_data.resize(unsigned(size));
		file.read(reinterpret_cast<char*>(obs_E_data.data()),sizeof(double)*size);
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be found"<<std::endl;
		exit(1);
	}
	file.close();
}

void System::saveObs1DData(const std::string& fileName){
	std::ofstream file( fileName.c_str(), std::ios::out | std::ios::binary );
	if( file.is_open() ){
		int size = obs_time.size();
		file.write(reinterpret_cast<const char*>(obs_time.data()),sizeof(double)*size);
		file.write(reinterpret_cast<const char*>(obs_E_data.data()),sizeof(double)*size);
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
		exit(1);
	}
	file.close();
}

void System::saveRandomSeedNumber(const std::string& fileName){
	std::ofstream file( fileName.c_str(), std::ios::out );
	if( file.is_open() ){
		file << randomSeedNumber;
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
		exit(1);
	}
	file.close();
}

void System::loadRandomSeedNumber(const std::string& fileName){
	std::ifstream file( fileName.c_str(), std::ios::in | std::ios::binary );
	if( file.is_open() ){
		file >> randomSeedNumber;
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be found"<<std::endl;
		exit(1);
	}
	file.close();
}

void System::setStateStructure(bool newVmVecOrdered, bool newModelVarsVecOrdered){
	VmVecOrdered =  newVmVecOrdered;
	modelVarsVecOrdered = newModelVarsVecOrdered;
	N = int(ceil(double(grid->getn())/VECLEN))*VECLEN;
	indexStride = (modelVarsVecOrdered) ? ( (VmVecOrdered) ? grid->getMaxNumVars() : (grid->getMaxNumVars()-1) ) : 1;
	varStride = (modelVarsVecOrdered) ? VECLEN : N;
	firstVecOrderedVar = (VmVecOrdered && modelVarsVecOrdered) ? 0 : 1;
}

void System::changeStateVector(double* y, int posIndex, double value, int numVar)const{
	if(numVar == 0)
		y[posIndex] = value;
	if(numVar != 0 || VmVecOrdered)
		y[getVecIndex(posIndex,numVar)] = value;
}

void System::changeStateVector(double* y, const std::vector<int> &posIndices, double value, int numVar)const{
	if(numVar == 0){
		#pragma omp parallel for
		for(int i=0;i<posIndices.size();i++)
			y[posIndices[i]] = value;
	}
	if(numVar != 0 || VmVecOrdered){
		#pragma omp parallel for
		for(int i=0;i<posIndices.size();i++)
			y[getVecIndex(posIndices[i],numVar)] = value;
	}
}

void System::changeStateVector(double* y, const std::vector<int> &posIndices, double* values, int numVar)const{
	if(numVar == 0){
		#pragma omp parallel for
		for(int i=0;i<posIndices.size();i++)
			y[posIndices[i]] = values[i];
	}
	if(numVar != 0 || VmVecOrdered){
		#pragma omp parallel for
		for(int i=0;i<posIndices.size();i++)
			y[getVecIndex(posIndices[i],numVar)] = values[i];
	}
}

void System::changeStateVector(double* y, double value, int numVar)const{
	if(numVar == 0){
		#pragma omp parallel for
		for(int i=0;i<getn();i++)
			y[i] = value;
	}
	if(numVar != 0 || VmVecOrdered){
		#pragma omp parallel for
		for(int i=0;i<getn();i++)
			y[getVecIndex(i,numVar)] = value;
	}
}

void System::changeStateVector(double* y, double* values, int numVar)const{
	if(numVar == 0){
		#pragma omp parallel for
		for(int i=0;i<getn();i++)
			y[i] = values[i];
	}
	if(numVar != 0 || VmVecOrdered){
		#pragma omp parallel for
		for(int i=0;i<getn();i++)
			y[getVecIndex(i,numVar)] = values[i];
	}
}
