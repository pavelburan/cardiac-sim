#ifndef GRID_FDM2DCARTESIAN_H
#define GRID_FDM2DCARTESIAN_H
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
class Grid_FDM2DCartesian : public Grid
{
public:
	//Konstruktoren
	Grid_FDM2DCartesian(const System& system, const std::string& configFileName, const std::string& keyPrefix="", int nx=0, int ny=0, double lx=0.0, double ly=0.0, std::string bcx="p", std::string bcy="p", double lamda[4]=NULL, double sigma[4]=NULL, double Dx=0.0, double Dy=0.0);
	//Destruktor
	~Grid_FDM2DCartesian();
	
private:
	Grid_FDM2DCartesian();
	Grid_FDM2DCartesian(const Grid_FDM2DCartesian &grid);
	Grid_FDM2DCartesian &operator=(const Grid_FDM2DCartesian &grid);
	int getPosIndex(int i, int j)const;
	double& getData(int i, int j);
public:
	//Zugriffsmethoden:
	int getn() const {return n;}
	int getMaxNumVars() const {return model->getnVars();}
	int getDim()const{return 2;}
	double getminh()const{return fmin(hx,hy);}
	double getmaxh()const{return fmax(hx,hy);}
	double getVges()const{return lx*ly;}
	double getVSphere(double r)const{return M_PI*r*r;};
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
	
	//Hetzugriff
	virtual bool getIsHet(int posIndex) const{ return isHet[posIndex];}
	
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
	int n;
	double lx;
	double ly;
	double hx;
	double hy;
	std::string bcx;
	std::string bcy;
	//0=links,1=rechts,2=oben,3=unten
	double bc_lamda[4];
	double bc_sigma[4];
	double Dx;
	double Dy;
	//Gitter
	double *xx;
	double *yy;
	//Heterogenitäten
	std::string hetsType;
	Hets *hets;
	bool *isHet;
	int nisHet;
	//Eintraege der Matrix
	double *data;     
	int offsets[5];
	double *bcVec;
	double *EVec;
public:
	//Model
	std::string modelType;
	Model* model; 
};

#endif //GRID_FDM2DCARTESIAN_H
