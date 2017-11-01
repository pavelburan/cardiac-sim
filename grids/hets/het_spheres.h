#ifndef HET_SPHERES_H
#define HET_SPHERES_H
#include "hets.h"
#include <vector>
class Grid;
class Het_spheres : public Hets
{
	struct Sphere{
		double r;
		int posIndex;
	};
public:
	//Konstruktor
	Het_spheres(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix="", double rmin=0.0, double rmax=0.0, double dr=0.0, double distance=0.0, double alpha=-2.74, double rho=0.0, double frac=0.0);
	//Destruktor
	~Het_spheres(){};
private:
	Het_spheres();
	Het_spheres(const Het_spheres &hets);
	Het_spheres &operator=(const Het_spheres &hets);
public:
	void loadHets(const std::string fileName);
	void saveHets(const std::string fileName)const;
	void readHetsParams();
	void initHets();
	void setisHet(bool *isHet)const;
	
	double getrmin()const{return rmin;}
	double getrmax()const{return rmax;}
	double getdr()const{return dr;}
	double getdistance()const{return distance;}
	double getalpha()const{return alpha;}
	double getrho()const{return rho;}
	
	const Sphere &operator()(int i) const {return spheres[i];};
	const Sphere &operator[](int i) const {return spheres[i];};
	int size()const{return spheres.size();}
		
protected:
	std::vector<Sphere> spheres;
	std::string fixedSpheres;
	double rmin;
	double rmax;
	double dr;
	double distance;
	double alpha;
	double rho;
};
#endif //HET_SPHERES_H

