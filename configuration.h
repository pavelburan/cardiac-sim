#ifndef CONFIGURATION_H
#define CONFIGURATION_H
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "ConfigFile.h"

#define VECLEN 8

namespace std{
	template <typename T>
	string to_str(T value)
	{
		//create an output string stream
		std::ostringstream os ;
		//throw the value into the string stream
		os << value ;
		//convert the string stream into a string and return
		return os.str() ;
	}
	
	template<class T>
	T to_T( const string& s )
	{
		T t;
		//create an input string stream
		std::istringstream is(s);
		//convert the string stream into a value of type T
		is >> t;
		return t;
	}
}

class Eps{
public:
	static int setMaxUlps(int newMaxUlps){mUlps = newMaxUlps;}
	static double setRelEps(double newRelEps){relEps = newRelEps;}
	static double setAbsEps_t(double newAbsEps_t){absEps_t = newAbsEps_t;}
	static double setAbsEps_x(double newAbsEps_x){absEps_x = newAbsEps_x;}
	static double setAbsEps_Vm(double newAbsEps_Vm){absEps_Vm = newAbsEps_Vm;}
	static const int& maxUlps(){return mUlps;}
	static const double& rel(){return relEps;}
	static const double& t(){return absEps_t;}
	static const double& x(){return absEps_x;}
	static const double& Vm(){return absEps_Vm;}
protected:
	static int mUlps;
	static double relEps;
	static double absEps_t;
	static double absEps_x;
	static double absEps_Vm;
};

class Configuration{
public:
	//Konstruktor
	Configuration(const std::string& configFileName, const std::string& keyPrefix=""):configFile(addPrefixesToFileName(configFileName),"=","+=","#",";",",","[","]"),configFileName(configFileName),keyPrefix(keyPrefix){}
	//Destruktor
	~Configuration();
private:
	Configuration();
	Configuration(const Configuration &cfg);
	Configuration &operator=(const Configuration &cfg);
public:
	const std::string& getConfigFileName()const{return configFileName;}
protected:
	ConfigFile configFile;
	std::string configFileName;
	std::string keyPrefix;
	static std::string defaultLevelSpace;
	static std::string actualLevelSpace;
	static std::vector<std::string> levelSpaces;
	static int subRepeatIndex;
	static std::vector<int> subIndex;
	static std::string mainConfigFileName;
	static std::string tempDir;
	static std::string destDir;
	static std::string subRepeatPrefix;
	static std::string subConfigPrefix;
	static std::string subTimePrefix;
	static std::string resumePrefix;
	static std::string savePrefix;
	static std::vector< std::vector<std::string> > cleanFiles;
	
public:
	static void rmSpecialChar(std::string& str, char spezialChar);
	template<class T>
	static int stringAsVector(std::vector<T>& vec, std::string str);
	static std::string addPrefixesToFileName(std::string fileName);
	static int getFileLength(const std::string& fileName, int typeSize=1);

	static void loadConfigParams(int argc, char * const argv[], const std::string& newMainConfigFileName, const std::string& defaultLevelSpace="", const std::vector<std::string>& levelSpaces=std::vector<std::string>(0));
	static void updateLevelSpace();
	static void levelUp();
	static void levelUp(const std::string& levelSpace);
	static void levelDown();
	static void print(const std::string& printString){std::cout<<actualLevelSpace<<printString<<std::endl;}
	static void printUp(const std::string& printString){print(printString); levelUp();}
	static void printUp(const std::string& printString, const std::string& levelSpace){print(printString); levelUp(levelSpace);}
	static void printDown(const std::string& printString){levelDown(); print(printString);}
	
	static void setDefaultLevelSpace(const std::string& newDefaultLevelSpace){defaultLevelSpace = newDefaultLevelSpace;};
	static std::string& getDefaultLevelSpace(){return defaultLevelSpace;};
	static int getSubRepeatIndex(){return subRepeatIndex;}
	static int getSubConfigIndex(){return subIndex[0];}
	static int getSubTimeIndex(){return subIndex[1];}
	static std::string& getMainConfigFileName(){return mainConfigFileName;}
	static std::string& getTempdir(){return tempDir;}
	static std::string& getDestDir(){return destDir;}
	static std::string& getSubRepeatPrefix(){return subRepeatPrefix;}
	static std::string& getSubConfigPrefix(){return subConfigPrefix;}
	static std::string& getSubTimePrefix(){return subTimePrefix;}
	static std::string& getResumePrefix(){return resumePrefix;}
	static std::string& getSavePrefix(){return savePrefix;}
	
	static std::string getDestFolderFileName(const std::string& fileName){ return tempDir + destDir + addPrefixesToFileName(fileName);}
	static std::string getPlotFolderFileName(const std::string& fileName){ return tempDir + destDir + "plotData/" + addPrefixesToFileName(fileName);}
	static std::string getPlotFolderSubRepeatResumePrefixFileName(const std::string& fileName){ return tempDir + destDir + "plotData/" + subRepeatPrefix + resumePrefix + fileName;}
	static std::string getPlotFolderSubRepeatSavePrefixFileName(const std::string& fileName){return tempDir + destDir + "plotData/" + subRepeatPrefix + savePrefix + fileName;}
	static std::string getPlotFolderSubTimeResumePrefixFileName(const std::string& fileName){return tempDir + destDir + "plotData/" + subTimePrefix + resumePrefix + fileName;}
	static std::string getPlotFolderSubTimeResumePrefixFileName(const std::string& fileName, int subTimeIndexShift){return tempDir + destDir + "plotData/" + std::to_str(subIndex[1]+subTimeIndexShift) + "_" + resumePrefix + fileName;}
	static std::string getPlotFolderSubTimeSavePrefixFileName(const std::string& fileName){return tempDir + destDir + "plotData/" + subTimePrefix + savePrefix + fileName;}
	static std::string getPlotFolderSubTimeSavePrefixFileName(const std::string& fileName, int subTimeIndexShift){return tempDir + destDir + "plotData/" + std::to_str(subIndex[1]+subTimeIndexShift) + "_" + savePrefix + fileName;}
	static std::string getInitialFileName(const std::string& fileName){return addPrefixesToFileName(fileName);}
	static void addCleanFile(std::string fileName, int cleanLevel){cleanFiles[cleanLevel].push_back( fileName.substr(fileName.find_last_of("/\\")+1) );}
	static void createCleanSH();
	
	template<class T>
	inline static void printVarStatic(const T& var, const std::string& key){std::cout<<actualLevelSpace<<key<<"="<<var<<std::endl;}
	template<class T>
	inline static void printVarStatic(const T& var, const std::string& key, const bool& found){std::cout<<actualLevelSpace<<key<<"="<<var<<"("<<found<<")"<<std::endl;}
	template<class T>
	inline void printVar(const T& var, std::string key, const std::string& keyName="")const;
	template<class T>
	inline void printVar(const T& var, std::string key, const bool& found, const std::string& keyName="")const;
	template<class T> 
	inline bool readInto(T& var, std::string key, const std::string& keyName="")const;
	inline bool readInto(std::string& var, std::string key, const std::string& keyName="")const;
	inline bool readRawInto(std::string& var, std::string key, const std::string& keyName="")const;
	template<class T>
	inline bool readIntoShifted(std::vector<int> subIndexShift, T& var, std::string key, const std::string& keyName="")const;
	template<class T> 
	bool readIntoVector(std::vector<T>& var, std::string key, const std::string& keyName="")const;
};

template<class T>
int Configuration::stringAsVector(std::vector<T>& vec, std::string str){
	rmSpecialChar(str, '\'');
	rmSpecialChar(str, '[');
	rmSpecialChar(str, '(');
	rmSpecialChar(str, ']');
	rmSpecialChar(str, ')');
	
	int size = 1;
	int pos = -1;
	while((pos = str.find(',',pos+1)) != std::string::npos)
		size++;
	if(vec.size() == 0)	
		vec.resize(size);
	int i = 0;
	int posLast = 0;
	pos = -1;
	while((pos = str.find(',',pos+1)) != std::string::npos && i < vec.size()){
		vec[i++] = std::to_T<T>(str.substr(posLast,pos-posLast));
		posLast = pos+1;
	}
	vec[i] = std::to_T<T>(str.substr(posLast,str.length()-posLast));
	return size;
}

template<class T>
void Configuration::printVar(const T& var, std::string key, const std::string& keyName)const{
	key = keyPrefix + key;
	if(keyName.length() > 0) key = keyName;
	printVarStatic(var, key);
}

template<class T>
void Configuration::printVar(const T& var, std::string key, const bool& found, const std::string& keyName)const{
	key = keyPrefix + key;
	if(keyName.length() > 0) key = keyName;
	printVarStatic(var, key, found);
}
template<class T> 
bool Configuration::readInto(T& var, std::string key, const std::string& keyName)const{
	key = keyPrefix + key;
	bool found = configFile.readElemInto(subIndex, var, key);
	if(keyName.length() > 0) key = keyName;
	printVarStatic(var, key, found);
	return found;
}

bool Configuration::readInto(std::string& var, std::string key, const std::string& keyName)const{
	key = keyPrefix + key;
	bool found = configFile.readElemInto(subIndex, var, key);
	rmSpecialChar(var, '\'');
	if(keyName.length() > 0) key = keyName;
	printVarStatic(var, key, found);	
	return found;
}
bool Configuration::readRawInto(std::string& var, std::string key, const std::string& keyName)const{
	key = keyPrefix + key;
	bool found = configFile.readInto(var, key);
	if(keyName.length() > 0) key = keyName;
	printVarStatic(var, key, found);	
	return found;
}

template<class T>
bool Configuration::readIntoShifted(std::vector<int> subIndexShift, T& var, std::string key, const std::string& keyName)const{
	key = keyPrefix + key;
	bool found = false;
	std::vector<int> shiftedSubIndex(subIndex);
	shiftedSubIndex[0] += subIndexShift[0];
	shiftedSubIndex[1] += subIndexShift[1];
	
	int flag = -1*configFile.elemExist(shiftedSubIndex,key);
	if(flag != 0)
		if(flag < 0 || !((flag&1 && subIndexShift[0]) || (flag&2 && subIndexShift[1]))) 
			found = configFile.readElemInto(shiftedSubIndex, var, key);
	if(keyName.length() > 0) key = keyName;
	printVarStatic(var, key, found);	
	return found;
}

template<class T> 
bool Configuration::readIntoVector(std::vector<T>& var, std::string key, const std::string& keyName)const{
	std::string str;
	bool found = readRawInto(str, key, keyName);
	if(found)
		stringAsVector(var, str);
	return found;
}

#endif //CONFIGURATION_H

