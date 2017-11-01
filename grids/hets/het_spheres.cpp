#include "het_spheres.h"
#include "../../system.h"
#include "../grid.h"
#include <iostream>
#include <math.h>
#include <algorithm>

Het_spheres::Het_spheres(const Grid& grid, const std::string& configFileName, const std::string& keyPrefix, double rmin, double rmax, double dr, double distance, double alpha, double rho, double frac)
												:Hets(grid,configFileName,keyPrefix,frac),spheres(),fixedSpheres(),rmin(rmin),rmax(rmax),dr(dr),distance(distance),alpha(alpha),rho(rho){
}

void Het_spheres::loadHets(const std::string fileName){
	std::ifstream file( fileName.c_str(), std::ios::out | std::ios::binary );
	if( file.is_open() ){
		double ind;
		int size = 0;
		file.seekg(0, file.end);
		size = file.tellg()/sizeof(double)/2;
		file.seekg(0, file.beg);
		spheres.resize(size);
		for(int i=0;i<size;i++){
			file.read(reinterpret_cast<char*>(&ind),sizeof(double));
			spheres[i].posIndex = int(ind);
		}
		for(int i=0;i<size;i++)
			file.read(reinterpret_cast<char*>(&(spheres[i].r)),sizeof(double));
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be found"<<std::endl;
		exit(1);
	}
	file.close();	
}

void Het_spheres::saveHets(const std::string fileName)const{
	std::ofstream file( fileName.c_str(), std::ios::out | std::ios::binary );
	if( file.is_open() ){
		double ind;
		for(int i=0;i<spheres.size();i++){
			ind = double(spheres[i].posIndex);
			file.write(reinterpret_cast<const char*>(&ind),sizeof(double));
		}
		for(int i=0;i<spheres.size();i++)
			file.write(reinterpret_cast<const char*>(&(spheres[i].r)),sizeof(double));
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
		exit(1);
	}
	file.close();	
}

void Het_spheres::readHetsParams(){
	Hets::readHetsParams();
	cfg.readRawInto(fixedSpheres, "fixedSpheres");
	cfg.readInto(rmin, "rmin");
	cfg.readInto(rmax, "rmax");
	cfg.readInto(dr, "dr");
	cfg.readInto(distance, "distance");
	cfg.readInto(alpha, "alpha");
	cfg.readInto(rho, "rho");
}

void Het_spheres::initHets(){
	Hets::initHets();
	std::string fileName;
	bool resume = grid.getSystem().resume();
	if(resume){
		fileName = cfg.getPlotFolderSubRepeatResumePrefixFileName("Het_spheres.bin");
		loadHets(fileName);
	}
	else{
		std::string initialHetsFile("NO");
		cfg.readInto(initialHetsFile, "initialHetsFile");
		if(initialHetsFile != "NO"){
			fileName = cfg.getInitialFileName(initialHetsFile);
			loadHets(fileName);
		}
		else{
			std::vector<double> vec;
			Configuration::stringAsVector(vec, fixedSpheres);
			int dim = grid.getDim();
			int nhet = vec.size() / (dim+1);
			spheres.resize(nhet);
			std::string printString;
			for(int i=0; i<nhet; ++i){
				spheres[i].posIndex = grid.getPosIndex(vec[(dim+1)*i], vec[(dim+1)*i+1], vec[(dim+1)*i+2]);
				spheres[i].r = vec[(dim+1)*i+dim];
				printString = "het";
				printString += std::to_str(i) + "=[";
				for(int k=0;k<dim;k++)
					printString += std::to_str(vec[(dim+1)*i+k]) + ",";
				printString += std::to_str(vec[(dim+1)*i+dim]) + "]";
				cfg.print(printString);
			}
			if(rho>0.0){
				double Vges = grid.getVges();
				int nr = int((rmax-rmin)/dr); //number of bins
				double dr = (rmax-rmin)/nr;	  // with of bins
				//Heterogenitätenverteilung berechnen
				double R[65536];
				double Pges = 1.0/(alpha+1)*(pow(rmax,alpha+1)-pow(rmin,alpha+1));
				double ra = rmin;
				double rb, r;
				int nbeg = 0;
				int nend;
				if(nr != 0){
					for(int i=0;i<nr;i++){
						rb = ra + dr;
						r = (ra + rb) / 2;
						double Pdr = 1.0/(alpha+1)*(pow(rb,alpha+1)-pow(ra,alpha+1));
						nend = nbeg + int(65536.0*Pdr/Pges);
						for(int j=nbeg;j<nend;j++)
							R[j] = r;
						ra = rb;
						nbeg = nend;
					}
				}
				else{
					nend = 65536;
					for(int j=0;j<nend;j++)
						R[j] = (rmin+rmax)/2;
				}
				//Heterogenitäten bestimmen
				//srand (time(NULL));
				std::vector<double> hetRadius;
				if(rho < 1.0){//rho is the fraction of tissue, which are part of heterogeneities
					double V = 0.0;
					while(V < Vges*rho){
						r = R[rand()%nend];
						hetRadius.push_back(r);
						V += grid.getVSphere(r);
					}
				}
				else{//rho is the density of heterogeneities
					int Nhet = int(rho * Vges / 100.0);
					for(int i=0;i<Nhet;i++){
						r = R[rand()%nend];
						hetRadius.push_back(r);
					}
				}
				std::sort(hetRadius.begin(), hetRadius.end());
				std::reverse(hetRadius.begin(), hetRadius.end());
				int posIndex;
				Sphere temp;
				for(int i=0;i<int(hetRadius.size());i++){
					r = hetRadius[i];
					bool finished = false;
					int k=0;
					while( !finished && k++ <1000000){
						posIndex = rand() % grid.getn();
						finished = true;
						for(int l=0;l<spheres.size();l++){
							if(fmin(grid.getDistance(posIndex,spheres[l].posIndex)-r-spheres[l].r, grid.getDistanceFromBorder(posIndex)-r) <= distance+Eps::x()){
								finished = false;
								break;
							}
						}
					}
					if(k >= 1000000){
						std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
						std::cerr<<"Heterogeneity "<<i<<" of "<<hetRadius.size()<<" doesn't fit"<<std::endl;
						exit(1);
					}
					temp.r = r;
					temp.posIndex = posIndex;
					spheres.push_back(temp);
				}
			}
		}
	}
	fileName = cfg.getPlotFolderSubRepeatSavePrefixFileName("Het_spheres.bin");
	saveHets(fileName);
	cfg.addCleanFile(fileName, 1);
	cfg.print("Gesamtzahl an Heterogenitäten: "+std::to_str(spheres.size()));
}

void Het_spheres::setisHet(bool *isHet)const{
	for(int k=0;k<spheres.size();k++){
		std::vector<int> posIndicesSphere = grid.getPosIndicesSphere(spheres[k].posIndex, spheres[k].r);
		for(int i=0;i<posIndicesSphere.size();i++){
			isHet[posIndicesSphere[i]] = true;
		}
	}
}
