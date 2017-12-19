#ifndef GRID_H
#define GRID_H
#include <vector>
#ifdef PETSC
#include <petscksp.h>
#endif
#include "../matrvec/matrvec.h"
#include "../configuration.h"

class System;
class Model;
class Grid
{
public:
	//Konstruktor
	Grid(const System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~Grid();
private:
	Grid();
	Grid(const Grid &grid);
	Grid &operator=(const Grid &grid);
public:
	//Zugriffsmethoden
	const Configuration& getCfg()const{return cfg;}
	const System& getSystem()const{return system;}
	virtual int getn()const =0;
	virtual int getMaxNumVars()const=0;
	virtual int getDim()const =0;
	virtual double getminh()const=0;
	virtual double getmaxh()const=0;
	virtual double getVges()const =0;
	virtual double getVSphere(double r)const =0;
	virtual int getPosIndex(double xi=0.0, double eta=0.0, double zeta=0.0)const =0;
	virtual std::vector<int> getPosIndicesVolume(int posIndex, double dxi=0.0, double deta=0.0, double dzeta=0.0)const =0;
	virtual std::vector<int> getPosIndicesVolume(double xi0=0.0, double dxi=0.0, double eta0=0.0, double deta=0.0, double zeta0=0.0, double dzeta=0.0)const {return getPosIndicesVolume(getPosIndex(xi0 ,eta0 ,zeta0), dxi, deta, dzeta);}
	virtual std::vector<int> getPosIndicesSphere(int posIndex, double r)const =0;
	virtual std::vector<int> getPosIndicesSphere(double r, double xi0=0.0, double eta0=0.0, double zeta0=0.0)const{return getPosIndicesSphere(getPosIndex(xi0 ,eta0 ,zeta0),r);}
	virtual double getDistance(int posIndex1, int posIndex2)const =0;
	virtual double getDistanceFromBorder(int posIndex)const =0;
	virtual double getFraction(const std::vector<bool>& isPointOfFraction)const =0;
	virtual double getFraction(const std::vector<int>& posIndices)const =0;
	virtual double getMeanIntegral(double* values)const =0;
	virtual double getMeanIntegral(double* values, const std::vector<bool>& isPointOfFraction)const =0;
	virtual double getMeanIntegral(double* values, const std::vector<int>& posIndices)const =0;
	
	//Modellzugriff
	virtual Model& getModel(int posIndex) const =0;
	const double& getVmThresh(int posIndex)const{return VmThresh[posIndex];}
	const double& getdVmdtThresh(int posIndex)const{return dVmdtThresh[posIndex];}
	virtual void setRestingState(double* y, const std::vector<int>& posIndices)const =0;
	virtual void setExcitedState(double* y, const std::vector<int>& posIndices)const =0;
	
	//Hetzugriff
	virtual bool getIsHet(int posIndex) const =0;
	
	virtual void readGridParams(){};
	virtual void initGrid();
	virtual void interpolateInitialStateOnGrid(double* y0, const double* yInitial, int InitialN) =0; 
	virtual void loadGrid(const std::string fileName) =0;
	virtual void saveGrid(const std::string fileName)const =0;
	
	virtual DIA<double> getCoeffMatrixMonoDIA(double a, double b) const =0;
	#ifdef PETSC
	virtual void getCoeffMatrixMonoPetsc(Mat &L, double a, double b) const =0;
	#endif
	virtual void getRhsMono(const double *y, double *rhsVm, double t, double dt, bool withVm = false) const =0;
	virtual void correctVm(double *u) const=0;
	
	//Gridcreator
	static Grid* newGrid(const std::string& gridType, const System& system, const std::string& configFileName, const std::string& keyPrefix);
	
protected:
	//Configuration
	Configuration cfg;
	//System
	const System& system;
	//Thresh
	double* VmThresh;
	double* dVmdtThresh;
	
};
#endif //GRID_H
