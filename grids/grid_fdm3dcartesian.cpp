#include "grid_fdm3dcartesian.h"
//#include "bessel.h"
#include "../configuration.h"
#include "../system.h"
#include "../efields/efield.h"
#include "hets/hets.h"
#include "hets/het_spheres.h"
#include <iostream>
#include <stdlib.h>
#include <limits>

Grid_FDM3DCartesian::Grid_FDM3DCartesian(const System& system, const std::string& configFileName, const std::string& keyPrefix, int nx, int ny, int nz, double lx, double ly, double lz, std::string bcx, std::string bcy, std::string bcz, double lamda[6], double sigma[6], double Dx, double Dy, double Dz)
																	:Grid(system,configFileName,keyPrefix),nx(nx),ny(ny),nz(nz),nxy(nx*ny),n(nx*ny*nz),lx(lx),ly(ly),lz(lz),hx(lx/(nx-1)),hy(ly/(ny-1)),hz(lz/(nz-1)),bcx(bcx),bcy(bcy),bcz(bcz),Dx(Dx),Dy(Dy),Dz(Dz),xx(NULL),yy(NULL),zz(NULL),
																	model(NULL),modelType(""),hetsType(""),hets(NULL),isHet(NULL),nisHet(0),data(NULL),bcVec(NULL),EVec(NULL)
{
	if(bcx == "d"){
		bc_lamda[0] = 0.0;
		bc_lamda[1] = 0.0;
		bc_sigma[0] = 1.0;
		bc_sigma[1] = 1.0;
	}
	else if(bcx == "n"){
		bc_lamda[0] = 1.0;
		bc_lamda[1] = 1.0;
		bc_sigma[0] = 0.0;
		bc_sigma[1] = 0.0;
	}
	else if(bcx == "c" && lamda != NULL && sigma != NULL){
		bc_lamda[0] = lamda[0];
		bc_lamda[1] = lamda[1];
		bc_sigma[0] = sigma[0];
		bc_sigma[1] = sigma[1];
	}
	if(bcy == "d"){
		bc_lamda[2] = 0.0;
		bc_lamda[3] = 0.0;
		bc_sigma[2] = 1.0;
		bc_sigma[3] = 1.0;
	}
	else if(bcy == "n"){
		bc_lamda[2] = 1.0;
		bc_lamda[3] = 1.0;
		bc_sigma[2] = 0.0;
		bc_sigma[3] = 0.0;
	}
	if(bcy == "c" && lamda != NULL && sigma != NULL){
		bc_lamda[2] = lamda[2];
		bc_lamda[3] = lamda[3];
		bc_sigma[2] = sigma[2];
		bc_sigma[3] = sigma[3];
	}
	if(bcz == "d"){
		bc_lamda[4] = 0.0;
		bc_lamda[5] = 0.0;
		bc_sigma[4] = 1.0;
		bc_sigma[5] = 1.0;
	}
	else if(bcz == "n"){
		bc_lamda[4] = 1.0;
		bc_lamda[5] = 1.0;
		bc_sigma[4] = 0.0;
		bc_sigma[5] = 0.0;
	}
	if(bcz == "c" && lamda != NULL && sigma != NULL){
		bc_lamda[4] = lamda[4];
		bc_lamda[5] = lamda[5];
		bc_sigma[4] = sigma[4];
		bc_sigma[5] = sigma[5];
	}
}

Grid_FDM3DCartesian::~Grid_FDM3DCartesian(){
	if(xx)
		delete[] xx;
	if(yy)
		delete[] yy;
	if(zz)
		delete[] yy;
	if(hets)
		delete hets;
	if(isHet)
		delete[] isHet;
	if(data)	
		delete[] data;
	if(bcVec)
		delete[] bcVec;
	if(EVec)
		delete[] EVec;
}

int Grid_FDM3DCartesian::getPosIndex(int i, int j, int k)const{
	i = (i %= nx) >= 0 ? i : (i+nx);
	j = (j %= ny) >= 0 ? j : (j+ny);
	k = (k %= nz) >= 0 ? k : (k+nz);
	return i + nx*j + nxy*k;
}

int Grid_FDM3DCartesian::getPosIndex(double xi, double eta, double zeta)const{
	return int(xi/hx+0.5) + nx*int(eta/hy+0.5) + nxy*int(zeta/hz+0.5);
}

std::vector<int> Grid_FDM3DCartesian::getPosIndicesVolume(int posIndex, double dxi, double deta, double dzeta)const{
	std::vector<int> posIndices;
	int iBeg = posIndex % nx;
	int jBeg = (posIndex / nx) % ny;
	int kBeg = posIndex / nxy;
	int iEnd = iBeg + int(dxi/hx+0.5);
	int jEnd = jBeg + int(deta/hy+0.5);
	int kEnd = kBeg + int(dzeta/hz+0.5);
	
	for(int k=kBeg;k<=kEnd;k++){
		for(int j=jBeg;j<=jEnd;j++){
			for(int i=iBeg;i<=iEnd;i++){
				if(!isHet[getPosIndex(i,j)])
					posIndices.push_back(getPosIndex(i,j,k));
			}
		}
	}

	return posIndices;
}

std::vector<int> Grid_FDM3DCartesian::getPosIndicesSphere(int posIndex, double r)const{
	std::vector<int> posIndices;
	int i0 = posIndex % nx;
	int j0 = (posIndex / nx) % ny;
	int k0 = posIndex / nxy;
	int di = int(r/hx+0.5);
	int dj = int(r/hy+0.5);
	int dk = int(r/hz+0.5);

	int iBeg = i0 - di;
	int iEnd = i0 + di;
	int jBeg = j0 - dj;
	int jEnd = j0 + dj;
	int kBeg = k0 - dk;
	int kEnd = k0 + dk;
	
	for(int k=kBeg;k<=kEnd;k++){
		for(int j=jBeg;j<=jEnd;j++){
			for(int i=iBeg;i<=iEnd;i++){
				if(!isHet[getPosIndex(i,j,k)]){
					if(pow((i-i0)*hx,2.0)+pow((j-j0)*hy,2.0)+pow((k-k0)*hz,2.0)<=pow(r,2.0))
						posIndices.push_back(getPosIndex(i,j,k));
				}
			}
		}
	}

	return posIndices;
}

double &Grid_FDM3DCartesian::getData(int i, int j){//todomaybe
	int k;
	int offset=j-i;
	int offpos=2*n;
	if(offset == -nxy || offset == n-nxy)
		offpos = 0;
	else if(offset == -nx || offset == nxy-nx)
		offpos = 1;
	else if(offset == -1 || offset == nx-1)
		offpos = 2;
	else if(offset == 0)
		offpos = 3;
	else if(offset == 1 || offset == 1-nx)
		offpos = 4;
	else if(offset == nx || offset == nx-nxy)
		offpos = 5;
	else if(offset == nxy || offset == nxy-n)
		offpos = 6;
	
	//#ifdef _ERROR_
	if(offpos==2*n){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Dia doesn't exist!"<<std::endl;
		std::cerr<<offset<<std::endl;
		exit(1);
	}
	//#endif
	return data[offpos*n+i];
}

double Grid_FDM3DCartesian::getDistance(int posIndex1, int posIndex2)const{
	return sqrt(pow(xx[posIndex1]-xx[posIndex2],2) + pow(yy[posIndex1]-yy[posIndex2],2) + pow(zz[posIndex1]-zz[posIndex2],2));
}

double Grid_FDM3DCartesian::getDistanceFromBorder(int posIndex)const{
	int i = posIndex % nx;
	int j = (posIndex / nx) % ny;
	int k = posIndex / nxy;
	double minDistance = std::min(std::min(lx,ly),lz);
	if(bcx != "p"){
		minDistance = std::min(i*hx,minDistance);
		minDistance = std::min((nx-1-i)*hx,minDistance);
	}
	if(bcy != "p"){
		minDistance = std::min(j*hy,minDistance);
		minDistance = std::min((ny-1-j)*hy,minDistance);
	}
	if(bcz != "p"){
		minDistance = std::min(k*hz,minDistance);
		minDistance = std::min((nz-1-k)*hz,minDistance);
	}
	return minDistance;
}

double Grid_FDM3DCartesian::getFraction(const std::vector<bool>& isPointOfFraction)const{
	double sum = 0.0;
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		if(!isHet[i] && isPointOfFraction[i]){
			#pragma omp atomic
			sum += 1.0;
		}
	}
	return sum / (n-nisHet);
}

double Grid_FDM3DCartesian::getFraction(const std::vector<int>& posIndices)const{
	double sum = 0.0;
	#pragma omp parallel for
	for(int i=0;i<posIndices.size();i++){
		int ind = posIndices[i];
		if(!isHet[ind]){
			#pragma omp atomic
			sum += 1.0;
		}
	}
	return sum / (n-nisHet);
}

double Grid_FDM3DCartesian::getMeanIntegral(double *__restrict__ values)const{
	double sum = 0.0;
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		if(!isHet[i]){
			#pragma omp atomic
			sum += values[i];
		}
	}
	return sum / (n-nisHet);
}

double Grid_FDM3DCartesian::getMeanIntegral(double *__restrict__ values, const std::vector<bool>& isPointOfFraction)const{
	double sum = 0.0;
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		if(!isHet[i] && isPointOfFraction[i]){
			#pragma omp atomic
			sum += values[i];
		}
	}
	return sum / (n-nisHet);
}

double Grid_FDM3DCartesian::getMeanIntegral(double *__restrict__ values, const std::vector<int>& posIndices)const{
	int nTotal = 0;
	double sum = 0.0;
	#pragma omp parallel for
	for(int i=0;i<posIndices.size();i++){
		int ind = posIndices[i];
		if(!isHet[ind]){
			#pragma omp atomic
			nTotal += 1;
			#pragma omp atomic
			sum += values[ind];
		}
	}
	return sum / nTotal;
}

void Grid_FDM3DCartesian::setRestingState(double *__restrict__ y, const std::vector<int>& posIndices)const{
	for(int numVar=0;numVar<model->getnVars();numVar++){
		double value = model->getRestingState(numVar);
		system.changeStateVector(y, posIndices, value, numVar);
	}
}

void Grid_FDM3DCartesian::setExcitedState(double *__restrict__ y, const std::vector<int>& posIndices)const{
	for(int numVar=0;numVar<model->getnVars();numVar++){
		double value = model->getExcitedState(numVar);
		system.changeStateVector(y, posIndices, value, numVar);
	}
}

void Grid_FDM3DCartesian::readGridParams(){
	Grid::readGridParams();
	cfg.readInto(nx, "nx");
	cfg.readInto(ny, "ny");
	cfg.readInto(nz, "nz");
	nxy = nx*ny;
	cfg.printVar(nxy, "nxy"); 
	n = nx*ny*nz;
	cfg.printVar(n, "n"); 
	
	cfg.readInto(lx, "lx");
	cfg.readInto(ly, "ly");
	cfg.readInto(lz, "lz");
	hx = lx/(nx-1);
	cfg.printVar(hx, "hx"); 
	hy = ly/(ny-1);
	cfg.printVar(hy, "hy");
	hz = lz/(nz-1);
	cfg.printVar(hz, "hz"); 
	
	cfg.readInto(bcx, "bcx");
	if(bcx == "d"){
		bc_lamda[0] = 0.0;
		bc_lamda[1] = 0.0;
		bc_sigma[0] = 1.0;
		bc_sigma[1] = 1.0;
	}
	else if(bcx == "n"){
		bc_lamda[0] = 1.0;
		bc_lamda[1] = 1.0;
		bc_sigma[0] = 0.0;
		bc_sigma[1] = 0.0;
	}
	else if(bcx == "c"){
		cfg.readInto(bc_lamda[0], "bc_l_lamda");
		cfg.readInto(bc_lamda[1], "bc_r_lamda");
		cfg.readInto(bc_sigma[0], "bc_l_sigma");
		cfg.readInto(bc_sigma[1], "bc_r_sigma");
	}
	cfg.readInto(bcy, "bcy");
	if(bcy == "d"){
		bc_lamda[2] = 0.0;
		bc_lamda[3] = 0.0;
		bc_sigma[2] = 1.0;
		bc_sigma[3] = 1.0;
	}
	else if(bcy == "n"){
		bc_lamda[2] = 1.0;
		bc_lamda[3] = 1.0;
		bc_sigma[2] = 0.0;
		bc_sigma[3] = 0.0;
	}
	else if(bcy == "c"){
		cfg.readInto(bc_lamda[2], "bc_u_lamda");
		cfg.readInto(bc_lamda[3], "bc_d_lamda");
		cfg.readInto(bc_sigma[2], "bc_u_sigma");
		cfg.readInto(bc_sigma[3], "bc_d_sigma");
	}
	cfg.readInto(bcz, "bcz");
	if(bcz == "d"){
		bc_lamda[4] = 0.0;
		bc_lamda[5] = 0.0;
		bc_sigma[4] = 1.0;
		bc_sigma[5] = 1.0;
	}
	else if(bcz == "n"){
		bc_lamda[4] = 1.0;
		bc_lamda[5] = 1.0;
		bc_sigma[4] = 0.0;
		bc_sigma[5] = 0.0;
	}
	else if(bcz == "c"){
		cfg.readInto(bc_lamda[4], "bc_f_lamda");
		cfg.readInto(bc_lamda[5], "bc_b_lamda");
		cfg.readInto(bc_sigma[4], "bc_f_sigma");
		cfg.readInto(bc_sigma[5], "bc_b_sigma");
	}
	
	cfg.readInto(Dx, "Dx");
	cfg.readInto(Dy, "Dy");
	cfg.readInto(Dz, "Dz");
	
	if(model)
		delete model;
	cfg.printUp("Modellparameter werden eingelesen...");
	cfg.readInto(modelType, "modelType");
	std::string modelConfigFileName = cfg.getConfigFileName();
	cfg.readInto(modelConfigFileName, "modelConfigFileName");
	std::string modelKeyPrefix = "";
	cfg.readInto(modelKeyPrefix, "modelKeyPrefix");
	
	model = Model::newModel(modelType, *this, modelConfigFileName, modelKeyPrefix);
	model->readModelParams();
	cfg.printDown("ready");
	
	if(hets)
		delete hets;
	cfg.printUp("Heterogenitätenparameter werden eingelesen...");
	cfg.readInto(hetsType, "hetsType");
	std::string hetsConfigFileName = cfg.getConfigFileName();
	cfg.readInto(hetsConfigFileName, "hetsConfigFileName");
	std::string hetsKeyPrefix = "";
	cfg.readInto(hetsKeyPrefix, "hetsKeyPrefix");
	
	hets = Hets::newHets(hetsType, *this, hetsConfigFileName, hetsKeyPrefix);
	hets->readHetsParams();
	cfg.printDown("ready");
}

void Grid_FDM3DCartesian::initGrid(){
	Grid::initGrid();
	std::string fileName;
	
	//Gitter initialsieren
	if(xx)
		delete[] xx;
	if(yy)
		delete[] yy;
	if(zz)
		delete[] zz;
	if(isHet)
		delete[] isHet;
	xx = new double[n];
	yy = new double[n];
	zz = new double[n];

	isHet = new bool[n];
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			for(int k=0;k<nz;k++){
				xx[getPosIndex(i,j,k)] = i*hx;
				yy[getPosIndex(i,j,k)] = j*hy;
				zz[getPosIndex(i,j,k)] = k*hz;
				isHet[getPosIndex(i,j,k)] = false;
			}
		}
	}
	
	//Heterogenitäten in Gitter einbauen
	bool resume = system.resume();
	cfg.printUp("Heterogenitäten werden initialisiert...");
	hets->initHets();
	cfg.printDown("ready");
	if(resume){
		fileName = cfg.getPlotFolderSubRepeatResumePrefixFileName("Grid.bin");
		loadGrid(fileName);
	}
	else{
		std::string initialGridFile("NO");
		cfg.readInto(initialGridFile, "initialGridFile");
		if(initialGridFile != "NO"){
			fileName = cfg.getInitialFileName(initialGridFile);
			loadGrid(fileName);
		}
		else{
			hets->setisHet(isHet);
			nisHet = 0;
			#pragma omp parallel for
			for(int i=0;i<n;i++){
				if(isHet[i]){
					#pragma omp atomic
					nisHet += 1;
				}
			}
			//Restanteil frac mit Heterogenitäten füllen
			double frac = hets->getfrac();
			if(frac > 0){
				srand(time(NULL));
				int k=0;
				for(int i=0;i<n;i++){
					if(isHet[i])
						k++;
				}
				while(k<int(frac*n)){
					int ind = rand()%n;
					if(!isHet[ind]){
						isHet[ind] = true;
						k++;
					}
				}
			}
		}
	}
	fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("Grid.bin");
	saveGrid(fileName);
	cfg.addCleanFile(fileName, 0);
	
	//Matrizen und Vektoren bestimmen
	if(data)	
		delete[] data;
	if(bcVec)
		delete[] bcVec;
	if(EVec)
		delete[] EVec;
	
	data = new double[7*n];
	bcVec = new double[n];
	EVec = new double[n];
	offsets[0] = -1*nxy;
	offsets[1] = -1*nx;
	offsets[2] = -1;
	offsets[3] = 0;
	offsets[4] = 1;
	offsets[5] = 1*nx;
	offsets[6] = 1*nxy;
	//Alle Punkte von bcVec auf 0 setzen:
	for(int i=0;i<n;i++){
		bcVec[i] = 0.0;
		EVec[i] = 0.0;
	}

	double fakx = Dx/(hx*hx);
	double faky = Dy/(hy*hy);
	double fakz = Dz/(hz*hz);

	double* __restrict__ pdata = data;
	double* __restrict__ pbcVec = bcVec;
	
	//Cauchy RB: lamda*du/dnu + sigma*u = g
	//Für (lamda=0,sigma=1) ergeben sich Dirichlet und für (lamda=1,sigma=0) ergeben sich Neumann RB
	double l_l_div_shxpl = bc_lamda[0]/(bc_sigma[0]*hx+bc_lamda[0]);
	double l_hx_div_shxpl = hx/(bc_sigma[0]*hx+bc_lamda[0]);
	double r_l_div_shxpl = bc_lamda[1]/(bc_sigma[1]*hx+bc_lamda[1]);
	double r_hx_div_shxpl = hx/(bc_sigma[1]*hx+bc_lamda[1]);
	double u_l_div_shypl = bc_lamda[2]/(bc_sigma[2]*hy+bc_lamda[2]);
	double u_hy_div_shypl = hy/(bc_sigma[2]*hy+bc_lamda[2]);
	double d_l_div_shypl = bc_lamda[3]/(bc_sigma[3]*hy+bc_lamda[3]);
	double d_hy_div_shypl = hy/(bc_sigma[3]*hy+bc_lamda[3]);
	double f_l_div_shzpl = bc_lamda[4]/(bc_sigma[4]*hz+bc_lamda[4]);
	double f_hz_div_shzpl = hz/(bc_sigma[4]*hz+bc_lamda[4]);
	double b_l_div_shzpl = bc_lamda[5]/(bc_sigma[5]*hz+bc_lamda[5]);
	double b_hz_div_shzpl = hz/(bc_sigma[5]*hz+bc_lamda[5]);
	
	double *g = new double[n];
	for(int i=0;i<n;i++)
		g[i] = 0.0;
	
	//Diffusionmatrix und bcVec für die Randbedingungen berechnen:
	for(int k=0;k<nz;k++){
		for(int j=0;j<ny;j++){
			for(int i=0;i<nx;i++){
				int index = (k*nxy+j*nx+i);
				double addCenter = 0.0; //Anteil der auf Diagonalelement zusätzlich addiert wird 
				
				//z-Richtung Punkte setzen
				if(k == 0){ //hinterer Rand
					pdata[index+6*n] = fakz;
					if(bcz == "p"){//Periodisch: u(i,j,-1) = u(i,j,nz-1)
						pdata[index] = fakz;
					}
					else{//Cauchy: u(i,j,-1) = b_l_div_shzpl * u(i,j,0) + b_hz_div_shzpl * g(i,j,-1)
						pdata[index] = 0.0;
						pbcVec[index] += fakz * b_hz_div_shzpl * g[index];
						addCenter += fakz * b_l_div_shzpl;
					}
				}
				else if(k+1 == nz){ //vorderer Rand
					pdata[index] = fakz;
					if(bcz== "p"){//Periodisch: u(i,j,nz) = u(i,j,0)
						pdata[index+6*n] = fakz;
					}
					else{ //Cauchy: u(i,j,nz) = f_l_div_shzpl * u(i,j,nz-1) + f_hz_div_shzpl * g(i,j,nz)
						pdata[index+6*n] = 0.0;
						pbcVec[index] += fakz * f_hz_div_shzpl * g[index];
						addCenter += fakz * f_l_div_shzpl;
					}
				}
				else{
					pdata[index] = fakz;
					pdata[index+6*n] = fakz;
				}
				
				//y-Richtung Punkte Setzen
				if(j == 0){ //unterer Rand
					pdata[index+5*n] = faky;
					if(bcy == "p"){//Periodisch: u(i,-1,k) = u(i,ny-1,k)
						pdata[index+n] = faky;
					}
					else{//Cauchy: u(i,-1,k) = d_l_div_shypl * u(i,0,k) + d_hy_div_shypl * g(i,-1,k)
						pdata[index+n] = 0.0;
						pbcVec[index] += faky * d_hy_div_shypl * g[index];
						addCenter += faky * d_l_div_shypl;
					}
				}
				else if(j+1 == ny){ //oberer Rand
					pdata[index+n] = faky;
					if(bcy == "p"){ //Periodisch: u(i,ny,k) = u(i,0,k)
						pdata[index+5*n] = faky;
					}
					else{//Cauchy: u(i,ny,k) = u_l_div_shypl * u(i,ny-1,k) + u_hy_div_shypl * g(i,ny,k)
						pdata[index+5*n] = 0.0;
						pbcVec[index] += faky * u_hy_div_shypl * g[index];
						addCenter += faky * u_l_div_shypl;
					}
				}
				else{
					pdata[index+n] = faky;
					pdata[index+5*n] = faky;
				}
				
				//x-Richtung Punkte setzen
				if(i == 0){ //linker Rand
					pdata[index+4*n] = fakx;
					if(bcx == "p"){//Periodisch: u(-1,j,k) = u(nx-1,j,k)
						pdata[index+2*n] = fakx;
					}
					else{//Cauchy: u(-1,j,k) = l_l_div_shxpl * u(0,j,k) + l_hx_div_shxpl * g(-1,j,k)
						pdata[index+2*n] = 0.0;
						pbcVec[index] += fakx * l_hx_div_shxpl * g[index];
						addCenter += fakx * l_l_div_shxpl;
					}
				}
				else if(i+1 == nx){ // rechter Rand
					pdata[index+2*n] = fakx;
					if(bcx == "p"){//Periodisch: u(nx,j,k) = u(0,j,k)
						pdata[index+4*n] = fakx;
					}
					else{ //Cauchy: u(nx,j,k) = r_div_shxpl * u(nx-1,j,k) + r_hx_div_shxpl * g(nx,j,k)
						pdata[index+4*n] = 0.0;
						pbcVec[index] += fakx * r_hx_div_shxpl * g[index];
						addCenter += fakx * r_l_div_shxpl;
					}
				}
				else{
					pdata[index+2*n] = fakx;
					pdata[index+4*n] = fakx;
				}
				
				//mittleren Punkt setzen
				pdata[index+3*n] = -2.0*(fakx+faky+fakz)+addCenter;
			}
		}
	}
	delete[] g;

	//Heterogenitäten einbauen
	double phi = system.getEfield().getphi();
	double theta = system.getEfield().gettheta();
	for(int I=0;I<n;I++){
		if(isHet[I]){
			int i = I % nx;
			int j = (I / nx) % ny;
			int k = I / nxy;
			int J[6] = {getPosIndex(i,j,k-1),getPosIndex(i,j-1,k),getPosIndex(i-1,j,k),getPosIndex(i+1,j,k),getPosIndex(i,j+1,k),getPosIndex(i,j,k+1)};
			//int J[6] = {I-nxy,I-nx,I-1,I+1,I+nx,I+nxy};
			double addEvec[6] = {-1.0*cos(theta)*hz, -1.0*sin(theta)*sin(phi)*hy, -1.0*sin(theta)*cos(phi)*hx, sin(theta)*cos(phi)*hx, sin(theta)*sin(phi)*hy, cos(theta)*hz};
			for(int m=0;m<6;m++){
				if(!isHet[J[m]])
					EVec[J[m]] += 0.1*getData(J[m],I)*addEvec[m];
				getData(I,I) += getData(I,J[m]);
				getData(J[m],J[m]) += getData(J[m],I);
				getData(I,J[m]) = 0.0;
				getData(J[m],I) = 0.0;
			}
		}
	}

	
	//std::cerr<<"Data:"<<std::endl;
	//pdata=data;
	//for(int i=0;i<n;i++){
	//	for(int j=0;j<7;j++)
	//		std::cerr<<std::scientific<<pdata[i+j*n]<<" ";
	//	std::cerr<<std::endl;
	//}
	
}

void Grid_FDM3DCartesian::interpolateInitialStateOnGrid(double* y0, const double* yInitial, int initialN){
	int initialNx = 0;
	cfg.readInto(initialNx, "initialNx");
	int initialNy = 0;
	cfg.readInto(initialNy, "initialNy");
	int initialNz = 0;
	cfg.readInto(initialNz, "initialNz");
	double initialLx = 0.0;
	cfg.readInto(initialLx, "initialLx");
	double initialLy = 0.0;
	cfg.readInto(initialLy, "initialLy");
	double initialLz = 0.0;
	cfg.readInto(initialLz, "initialLz");
	double initialX0 = 0.0;
	cfg.readInto(initialX0, "initialX0");
	double initialY0 = 0.0;
	cfg.readInto(initialY0, "initialY0");
	double initialZ0 = 0.0;
	cfg.readInto(initialZ0, "initialZ0");
	if(model->getnVars()*initialNx*initialNy*initialNz == initialN){
		if(initialX0 + lx <= initialLx+Eps::x() && initialY0 + ly <= initialLy+Eps::x() && initialZ0 + lz <= initialLz+Eps::x()){
			double initialHx = initialLx/(initialNx-1);
			double initialHy = initialLy/(initialNy-1);
			double initialHz = initialLz/(initialNz-1);
			int i0 = int(initialX0/initialHx+0.5);
			int j0 = int(initialY0/initialHy+0.5);
			int k0 = int(initialZ0/initialHz+0.5);
			
			if(fabs(initialLx/initialNx - lx/nx) <= Eps::x() && fabs(initialLy/initialNy - ly/ny) <= Eps::x() && fabs(initialLz/initialNz - lz/nz) <= Eps::x() || fabs(initialLx/initialNx - hx) <= Eps::x() && fabs(initialLy/initialNy - hy) <= Eps::x() && fabs(initialLz/initialNz - hz) <= Eps::x() || fabs(initialHx - lx/nx) <= Eps::x() && fabs(initialHy - ly/ny) <= Eps::x() && fabs(initialHz - lz/nz) <= Eps::x()){//todo thats cheating 
				initialHx = hx;
				initialHy = hy;
				initialHz = hz;
			}
			
			if(fabs(initialHx-hx) <= Eps::x() && fabs(initialHy-hy) <= Eps::x() && fabs(initialHz-hz) <= Eps::x()){
				int nVars = getMaxNumVars();
				double *temp = new double[getn()];
				for(int m=0;m<nVars;m++){
					for(int k=0;k<nz;k++){
						for(int j=0;j<ny;j++){
							for(int i=0;i<nx;i++){
								temp[k*nxy+j*nx+i] = yInitial[m*initialNx*initialNy*initialNz + (k0+k)*initialNx*initialNz + (j0+j)*initialNx + (i0+i)];
							}
						}
					}
					system.changeStateVector(y0, temp, m);
				}
				delete[] temp;
				/*
				int N = system.getN();
				int nVars = getMaxNumVars();
				int kStride = system.getStateIndexStride();
				int mStride = system.getStateVarStride();
				int m0 = system.getStateFirstVecOrderedVar();
				#pragma omp parallel for
				#pragma ivdep
				for(int k=0;k<N-VECLEN;k+=VECLEN){
					for(int l=0;l<VECLEN;l++){
						int i = (k+l)%nx;
						int j = (k+l)/nx;
						y0[k+l] = yInitial[(j0+j)*initialNx + (i0+i)];
						for(int m=m0;m<nVars;m++){
							y0[N+k*kStride+(m-m0)*mStride+l] = yInitial[m*initialNx*initialNy+(j0+j)*initialNx + (i0+i)];
						}
					}
				}
				int k = (N-VECLEN);
				for(int l=0;l<VECLEN-N+getn();l++){
					int i = (k+l)%nx;
					int j = (k+l)/nx;
					y0[k+l] = yInitial[(j0+j)*initialNx + (i0+i)];
					for(int m=m0;m<nVars;m++){
						y0[N+k*kStride+(m-m0)*mStride+l] = yInitial[m*initialNx*initialNy+(j0+j)*initialNx + (i0+i)];
					}
				}*/
			}
			else{//todo
				std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
				std::cerr<<"initialHx ="<<initialHx <<" != hx="<<hx<<" || initialHy ="<<initialHy <<" != hy="<<hy<<std::endl;
				std::cerr<<"Not supported yet"<<std::endl;
				exit(1);
			}
			std::vector<int> posIndices;
			for(int i=0;i<n;i++)
				if(isHet[i])
					posIndices.push_back(i);
			setRestingState(y0, posIndices);

		}
		else{
			std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"initialX0 + lx ="<<initialX0 + lx <<" > initialLx="<<initialLx<<" || initialY0 + ly="<<initialY0 + ly<<" > initialLy="<<initialLy<<" || initialZ0 + lz="<<initialZ0 + lz<<" > initialLz="<<initialLz<<std::endl;
			exit(1);
		}
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"initialN="<<initialN<<" != model->getnVars()*initialNx*initialNy*initialNz="<<model->getnVars()*initialNx*initialNy*initialNz<<std::endl;
		exit(1);
	}
}

void Grid_FDM3DCartesian::loadGrid(const std::string fileName){
	std::ifstream fileHet( fileName.c_str(), std::ios::in | std::ios::binary );
	if( fileHet.is_open() ){
		fileHet.read(reinterpret_cast<char*>(isHet),sizeof(bool)*n);
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be found"<<std::endl;
		exit(1);
	}
	fileHet.close();
}

void Grid_FDM3DCartesian::saveGrid(const std::string fileName)const{
	std::ofstream fileHet( fileName.c_str(), std::ios::out | std::ios::binary );
	if( fileHet.is_open() ){
		fileHet.write(reinterpret_cast<char*>(isHet),sizeof(bool)*n);
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
		exit(1);
	}
	fileHet.close();
}

DIA<double> Grid_FDM3DCartesian::getCoeffMatrixMonoDIA(double a, double b) const{
	std::string specType;	
	if(bcx == "p" && bcy == "p" && bcz == "p")
		specType = "lex3Dppp";
	else if(bcx == "p" && bcy == "p" && bcz != "p")
		specType = "lex3Dppc";
	else if(bcx == "p" && bcy != "p" && bcz == "p")
		specType = "lex3Dpcp";
	else if(bcx == "p" && bcy != "p" && bcz != "p")
		specType = "lex3Dpcc";
	else if(bcx != "p" && bcy == "p" && bcz == "p")
		specType = "lex3Dcpp";
	else if(bcx != "p" && bcy == "p" && bcz != "p")
		specType = "lex3Dcpc";
	else if(bcx != "p" && bcy != "p" && bcz == "p")
		specType = "lex3Dccp";
	else
		specType = "lex3Dccc";
		
	
	DIA<double> L(n,n,7,data,offsets,specType.c_str());
	L *= b;
	L.addScalToDiagonal(a, 3);
	return L;
}

#ifdef PETSC
void Grid_FDM3DCartesian::getCoeffMatrixMonoPetsc(Mat &L, double a, double b) const{
	MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 7, 0, &L);

	int i_A, j_A;
	double value;
	for(int k=0;k<nz;k++){
		for(int j=0;j<ny;j++){
			for(int i=0;i<nx;i++){
				i_A = i + nx*j + nxy*k;
				
				if(k==0){
					if(bcz == "p"){
						j_A = i_A + nxy*(nz-1);
						value = data[i_A]*b;
						MatSetValues(L,1,&i_A,1,&j_A,&value,INSERT_VALUES);
					}
				}
				else{
					j_A = i_A - nxy;
					value = data[i_A]*b;
					MatSetValues(L,1,&i_A,1,&j_A,&value,INSERT_VALUES);
				}
				
				if(j==0){
					if(bcy == "p"){
						j_A = i_A + nx*(ny-1);
						value = data[i_A+n]*b;
						MatSetValues(L,1,&i_A,1,&j_A,&value,INSERT_VALUES);
					}
				}
				else{
					j_A = i_A - nx;
					value = data[i_A+n]*b;
					MatSetValues(L,1,&i_A,1,&j_A,&value,INSERT_VALUES);
				}
				
				if(i==0){
					if(bcx == "p"){
						j_A = i_A + (nx-1);
						value = data[i_A+2*n]*b;
						MatSetValues(L,1,&i_A,1,&j_A,&value,INSERT_VALUES);
					}
				}
				else{
					j_A = i_A - 1;
					value = data[i_A+2*n]*b;
					MatSetValues(L,1,&i_A,1,&j_A,&value,INSERT_VALUES);
				}
				
				j_A = i_A;
				value = data[i+3*n]*b + a;
				MatSetValues(L,1,&i_A,1,&j_A,&value,INSERT_VALUES);
				
				if(i+1==nx){
					if(bcx == "p"){
						j_A = i_A - (nx-1);
						value = data[i_A+4*n]*b;
						MatSetValues(L,1,&i_A,1,&j_A,&value,INSERT_VALUES);
					}
				}
				else{
					j_A = i_A + 1;
					value = data[i_A+4*n]*b;
					MatSetValues(L,1,&i_A,1,&j_A,&value,INSERT_VALUES);
				}
				
				if(j+1==ny){
					if(bcy == "p"){
						j_A = i_A - nx*(ny-1);
						value = data[i_A+5*n]*b;
						MatSetValues(L,1,&i_A,1,&j_A,&value,INSERT_VALUES);
					}
				}
				else{
					j_A = i_A + nx;
					value = data[i_A+5*n]*b;
					MatSetValues(L,1,&i_A,1,&j_A,&value,INSERT_VALUES);
				}
				
				if(k+1==nz){
					if(bcz == "p"){
						j_A = i_A - nxy*(nz-1);
						value = data[i_A+6*n]*b;
						MatSetValues(L,1,&i_A,1,&j_A,&value,INSERT_VALUES);
					}
				}
				else{
					j_A = i_A + nxy;
					value = data[i_A+6*n]*b;
					MatSetValues(L,1,&i_A,1,&j_A,&value,INSERT_VALUES);
				}
				
				
			}
		}
	}
	
	MatAssemblyBegin(L, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(L, MAT_FINAL_ASSEMBLY);
	MatSetOption(L,MAT_SYMMETRIC,PETSC_TRUE);
}
#endif

void Grid_FDM3DCartesian::getRhsMono(const double *__restrict__ Vm, double *__restrict__ rhsVm, double t, double dt, bool withVm) const{
	double E = system.getEfield().E(t);
	if(withVm){
		#pragma omp parallel for
		#pragma ivdep
		for(int i=0;i<n;i++)
			rhsVm[i] = Vm[i] + dt*(bcVec[i] + E*EVec[i]);
	}
	else{
		#pragma omp parallel for
		#pragma ivdep
		for(int i=0;i<n;i++)
			rhsVm[i] = dt*(bcVec[i] + E*EVec[i]);
	}
}

void Grid_FDM3DCartesian::correctVm(double *__restrict__ u) const{
	double VmMax = model->getVmMax();
	double VmMin = model->getVmMin();
	#pragma omp parallel for
	#pragma ivdep
	for(int i=0;i<n;i++){
		u[i] = (u[i] < VmMax) ? u[i] : VmMax;
		u[i] = (u[i] > VmMin) ? u[i] : VmMin;
	}
}

