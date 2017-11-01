#ifndef ORDERPARAMETER_H
#define ORDERPARAMETER_H
#include "../configuration.h"
#include <list>

class System;
class Efield;
class Grid;
class OrderParameter
{
public:
	//Konstruktor
	OrderParameter(const System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~OrderParameter(){}
private:
	OrderParameter();
	OrderParameter(const OrderParameter &orderParameter);
	OrderParameter &operator=(const OrderParameter &orderParameter);
public:

	virtual void readParams()=0;
	virtual void init()=0;
	virtual double calc(double *y_prev, double *y, double *dVmdt, double *temp, double t, double timeStep)const =0;
	
	//OrderParameterCreator
	static OrderParameter* newOrderParameter(const std::string& orderParameterType, const System& system, const std::string& configFileName, const std::string& keyPrefix);
	
protected:
	//Configuration
	Configuration cfg;
	//System
	const System& system;
	//Efield
	const Efield& efield;
	//Gitter
	const Grid& grid;
};

class OP_meanVm : public OrderParameter
{
public:
	//Konstruktor
	OP_meanVm(const System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~OP_meanVm(){}
private:
	OP_meanVm();
	OP_meanVm(const OrderParameter &orderParameter);
	OP_meanVm &operator=(const OrderParameter &orderParameter);
public:
	//Zugriffsmethoden	
	void readParams(){}
	void init(){};
	double calc(double *y_prev, double *y, double *dVmdt, double *temp, double t, double timeStep)const;
};

class OP_excitable : public OrderParameter
{
public:
	//Konstruktor
	OP_excitable(const System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~OP_excitable(){}
private:
	OP_excitable();
	OP_excitable(const OrderParameter &orderParameter);
	OP_excitable &operator=(const OrderParameter &orderParameter);
public:
	//Zugriffsmethoden	
	void readParams(){}
	void init(){};
	double calc(double *y_prev, double *y, double *dVmdt, double *temp, double t, double timeStep)const;
};

class OP_excited : public OrderParameter
{
public:
	//Konstruktor
	OP_excited(const System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~OP_excited(){}
private:
	OP_excited();
	OP_excited(const OrderParameter &orderParameter);
	OP_excited &operator=(const OrderParameter &orderParameter);
public:
	//Zugriffsmethoden	
	void readParams(){}
	void init(){};
	double calc(double *y_prev, double *y, double *dVmdt, double *temp, double t, double timeStep)const;
};

class OP_exciting : public OrderParameter
{
public:
	//Konstruktor
	OP_exciting(const System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~OP_exciting(){}
private:
	OP_exciting();
	OP_exciting(const OrderParameter &orderParameter);
	OP_exciting &operator=(const OrderParameter &orderParameter);
public:
	//Zugriffsmethoden	
	void readParams(){}
	void init(){};
	double calc(double *y_prev, double *y, double *dVmdt, double *temp, double t, double timeStep)const;
};

class OP_cluster : public OrderParameter
{
public:
	//Konstruktor
	OP_cluster(const System& system, const std::string& configFileName, const std::string& keyPrefix="");
	//Destruktor
	virtual ~OP_cluster(){}
private:
	OP_cluster();
	OP_cluster(const OrderParameter &orderParameter);
	OP_cluster &operator=(const OrderParameter &orderParameter);
	void rekursivRemoveClusterElement(std::list<int>& posIndicesExciting, std::vector<std::list<int>::iterator>& itVector)const;
public:
	//Zugriffsmethoden	
	void readParams();
	void init();
	double calc(double *y_prev, double *y, double *dVmdt, double *temp, double t, double timeStep)const;
protected:
	double maxDistance;
	mutable std::vector< std::list<int>::iterator > itTotalVector;
	mutable std::vector< bool > isExciting;
};

#endif //ORDERPARAMETER_H
