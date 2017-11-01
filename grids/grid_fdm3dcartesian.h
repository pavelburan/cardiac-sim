#ifndef GRID_FDM3DCARTESIAN_H
#define GRID_FDM3DCARTESIAN_H
#include "grid.h"
#include "../matrvec/matrvec.h"
#include "../models/model.h"
#include <vector>
#ifdef PETSC
#include <petscksp.h>
#endif
#include <math.h>

class System;
class Model;
class Hets;
class Grid_FDM3DCartesian : public Grid
{
public:
	//Konstruktoren
	Grid_FDM3DCartesian(const System& system, const std::string& configFileName, const std::string& keyPrefix="", int nx=0, int ny=0, int nz=0, double lx=0.0, double ly=0.0, double lz=0.0, std::string bcx="p", std::string bcy="p", std::string bcz="p", double lamda[6]=NULL, double sigma[6]=NULL, double Dx=0.0, double Dy=0.0, double Dz=0.0);
	//Destruktor
	~Grid_FDM3DCartesian();
	
private:
	Grid_FDM3DCartesian();
	Grid_FDM3DCartesian(const Grid_FDM3DCartesian &grid);
	Grid_FDM3DCartesian &operator=(const Grid_FDM3DCartesian &grid);
	int getPosIndex(int i, int j, int k)const;
	double& getData(int i, int j);
public:
	//Zugriffsmethoden:
	int getn() const {return n;}
	int getMaxNumVars() const {return model->getnVars();}
	int getDim()const{return 3;}
	double geth()const{return fmin(fmin(hx,hy),hz);}
	double getVges()const{return lx*ly*lz;}
	double getVSphere(double r)const{return 4.0/3.0*M_PI*r*r*r;};
	int getPosIndex(double xi=0.0, double eta=0.0, double zeta=0.0)const;
	std::vector<int> getPosIndicesVolume(int posIndex, double dxi=0.0, double deta=0.0, double dzeta=0.0)const;
	std::vector<int> getPosIndicesSphere(int posIndex, double r)const;
	double getDistance(int posIndex1, int posIndex2)const;
	double getDistanceFromBorder(int posIndex)const;
	double getFraction(const std::vector<bool>& isPointOfFraction)const;
	double getFraction(const std::vector<int>& posIndices)const;
	double getMeanIntegral(double* values)const;
	double getMeanIntegral(double* values, const std::vector<bool>& isPointOfFraction)const;
	double getMeanIntegral(double* values, const std::vector<int>& posIndices)const;
	
	//Modellzugriff
	Model& getModel(int posIndex) const{return *model;}
	void setRestingState(double* y, const std::vector<int>& posIndices)const;
	void setExcitedState(double* y, const std::vector<int>& posIndices)const;
	
	void readGridParams();
	void initGrid();
	void interpolateInitialStateOnGrid(double* y0, const double* yInitial, int InitialN); 
	void loadGrid(const std::string fileName);
	void saveGrid(const std::string fileName)const;
	
	DIA<double> getCoeffMatrixMonoDIA(double a, double b) const;
	#ifdef PETSC
	void getCoeffMatrixMonoPetsc(Mat &L, double a, double b) const;
	#endif
	void getRhsMono(const double *Vm, double *rhsVm, double t, double dt, bool withVm = false) const;
	void correctVm(double *u) const;

protected:
	int nx;
	int ny;
	int nz;
	int nxy;
	int n;
	double lx;
	double ly;
	double lz;
	double hx;
	double hy;
	double hz;
	std::string bcx;
	std::string bcy;
	std::string bcz;
	//0=links,1=rechts,2=oben,3=unten,4=vorne,5=hinten
	double bc_lamda[6];
	double bc_sigma[6];
	double Dx;
	double Dy;
	double Dz;
	//Gitter
	double *xx;
	double *yy;
	double *zz;
	//Heterogenitäten
	std::string hetsType;
	Hets *hets;
	bool *isHet;
	int nisHet;
	//Eintraege der Matrix
	double *data;     
	int offsets[7];
	double *bcVec;
	double *EVec;
public:
	//Model
	std::string modelType;
	Model* model; 
};

#endif //GRID_FDM3DCARTESIAN_H
