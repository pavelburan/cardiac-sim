#include "configuration.h"
#include <fstream>
#include <stdlib.h>

int Eps::mUlps = 4;
double Eps::relEps = 1e-6;
double Eps::absEps_t = 1e-3;
double Eps::absEps_x = 1e-3;
double Eps::absEps_Vm = 1e-3;

std::string Configuration::defaultLevelSpace = "";
std::string Configuration::actualLevelSpace = "";
std::vector<std::string> Configuration::levelSpaces = std::vector<std::string>(0);
int Configuration::subRepeatIndex = 0;
std::vector<int> Configuration::subIndex = std::vector<int>(2,0);
std::string Configuration::mainConfigFileName = "";
std::string Configuration::tempDir = "";
std::string Configuration::destDir = "";
std::string Configuration::subRepeatPrefix = "";
std::string Configuration::subConfigPrefix = "";
std::string Configuration::subTimePrefix = "";
std::string Configuration::resumePrefix = "";;
std::string Configuration::savePrefix = "";
std::vector< std::vector<std::string> > Configuration::cleanFiles = std::vector< std::vector<std::string> >(5);

Configuration::~Configuration(){}

void Configuration::rmSpecialChar(std::string& str, char spezialChar){
	for(unsigned int i=0; i<str.length(); ++i)
		if(str[i] == spezialChar)
			str.erase(i--,1);
}

std::string Configuration::addPrefixesToFileName(std::string fileName, int subTimeIndexShift){
	int pos = -1;
	while((pos=fileName.find("(@SR)", pos+1)) != string::npos)//SR SC ST
		fileName.replace(pos, 5, subRepeatPrefix);
	pos = -1;
	while((pos=fileName.find("(@SC)", pos+1)) != string::npos)//SR SC ST
		fileName.replace(pos, 5, subConfigPrefix);
	pos = -1;
	std::string subTimePrefixTemp = subTimePrefix;
	if(subTimeIndexShift)
		subTimePrefixTemp = std::to_str(subIndex[1]+subTimeIndexShift) + "_";
	while((pos=fileName.find("(@ST)", pos+1)) != string::npos)//SR SC ST
		fileName.replace(pos, 5, subTimePrefixTemp);
	pos = -1;
	while((pos=fileName.find("(@P)", pos+1)) != string::npos)//SR SC ST
		fileName.replace(pos, 4, subTimePrefix);
	return fileName;
}

int Configuration::getFileLength(const std::string& fileName, int typeSize){
	int length = 0;
	std::ifstream file(fileName.c_str(), std::ios::in | std::ios::binary );
	if( file.is_open() ){
		file.seekg(0, file.end);
		length = file.tellg();
		file.seekg(0, file.beg);
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be found"<<std::endl;
		exit(1);
	}
	file.close();
	return length/typeSize;
}


void Configuration::loadConfigParams(int argc, char * const argv[], const std::string& newMainConfigFileName, const std::string& newDefaultLevelSpace, const std::vector<std::string>& newLevelSpaces){
	bool found;
	mainConfigFileName = newMainConfigFileName;
	defaultLevelSpace = newDefaultLevelSpace;
	levelSpaces = newLevelSpaces;
	updateLevelSpace();
	printUp("Konfigurationsparameter werden eingelesen...");
	
	found = false;
	if(argc > 1){
		mainConfigFileName = argv[1];
		found = true;
	}
	ConfigFile configFile(mainConfigFileName,"=","+=","#",";",",","[","]");

	printVarStatic(mainConfigFileName, "mainConfigFileName", found);
	
	found = false;
	if(argc > 2){
		tempDir = argv[2];
		int length = tempDir.length();
		if(length > 0)
			if(tempDir.at(length-1) != '/' && tempDir.at(length-1) != '\\')
				tempDir.push_back('/');
		found = true;
	}
	printVarStatic(tempDir, "tempDir", found);
	

	found = configFile.readInto(destDir, "destDir");
	rmSpecialChar(destDir, '\'');
	int length = destDir.length();
	if(length > 0)
		if(destDir.at(length-1) != '/' && destDir.at(length-1) != '\\')
			destDir.push_back('/');
	printVarStatic(destDir, "destDir", found);
	
	found = false;
	if(argc > 3){
		subRepeatIndex = atoi(argv[3]);
		if(subRepeatIndex >= 0)
			subRepeatPrefix = argv[3] + std::string("_");
		else
			subRepeatIndex = 0;
		found = true;
	}
	printVarStatic(subRepeatPrefix, "subRepeatPrefix", found);
	
	
	found = false;
	if(argc > 4){
		subIndex[0] = atoi(argv[4]);
		if(subIndex[0] >= 0)
			subConfigPrefix = argv[4] + std::string("_");
		else
			subIndex[0] = 0;
		found = true;
	}
	printVarStatic(subConfigPrefix, "subConfigPrefix", found);
	
	
	found = false;
	if(argc > 5){
		subIndex[1] = atoi(argv[5]);
		if(subIndex[1] >= 0){
			subTimePrefix = argv[5] + std::string("_");
		}
		else
			subIndex[1] = 0;
		found = true;
	}
	printVarStatic(subTimePrefix, "subTimePrefix", found);
	
	
	found = configFile.readElemInto(subIndex, resumePrefix, "resumePrefix");
	rmSpecialChar(resumePrefix, '\'');
	printVarStatic(resumePrefix, "resumePrefix", found);
	
	found = configFile.readElemInto(subIndex, savePrefix, "savePrefix");
	rmSpecialChar(savePrefix, '\'');
	printVarStatic(savePrefix, "savePrefix", found);
	
	printDown("ready\n");
}

void Configuration::updateLevelSpace(){
	actualLevelSpace = "";
	for(int i=0;i<levelSpaces.size();i++)
		actualLevelSpace += levelSpaces[i];
}

void Configuration::levelUp(){
	levelSpaces.push_back(defaultLevelSpace);
	updateLevelSpace();
}

void Configuration::levelUp(const std::string& levelSpace){
	levelSpaces.push_back(levelSpace);
	updateLevelSpace();
}

void Configuration::levelDown(){
	if(levelSpaces.size() > 0)
		levelSpaces.pop_back();
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Reached already printLevel=0"<<std::endl;
		exit(1);
	}
	updateLevelSpace();
}

void Configuration::createCleanSH(){
	std::string fileName = tempDir + destDir + "plotData/" + subRepeatPrefix + subConfigPrefix + subTimePrefix + "clean.sh"; 
	std::ofstream fileClean( fileName.c_str(), std::ios::out );
	if( fileClean.is_open() ){
		fileClean<<"#!/bin/bash \n";
		fileClean<<"cleanLevel=0 \n";
		fileClean<<"if [ $# > 0 ] \n";
		fileClean<<"then \n";
		fileClean<<"cleanLevel=$1 \n";
		fileClean<<"fi \n";
		fileClean<<"if [ $cleanLevel -ge 0 ] \n";
		fileClean<<"then \n";
		for(int i=0;i<cleanFiles[0].size();i++)
			fileClean<<"rm -f "<<cleanFiles[0][i]<<"\n";//0
		fileClean<<"fi \n";
		fileClean<<"if [ $cleanLevel -ge 1 ] \n";
		fileClean<<"then \n";
		for(int i=0;i<cleanFiles[1].size();i++)
			fileClean<<"rm -f "<<cleanFiles[1][i]<<"\n";//1
		fileClean<<"fi \n";
		fileClean<<"if [ $cleanLevel -ge 2 ] \n";
		fileClean<<"then \n";
		for(int i=0;i<cleanFiles[2].size();i++)
			fileClean<<"rm -f "<<cleanFiles[2][i]<<"\n";//2
		fileClean<<"rm -f "<<subRepeatPrefix<<savePrefix<<"y_*.bin"<<"\n"; //2
		fileClean<<"fi \n";
		fileClean<<"if [ $cleanLevel -ge 3 ] \n";
		fileClean<<"then \n";
		for(int i=0;i<cleanFiles[3].size();i++)
			fileClean<<"rm -f "<<cleanFiles[3][i]<<"\n";//3
		fileClean<<"fi \n";
		fileClean<<"if [ $cleanLevel -ge 4 ] \n";
		fileClean<<"then \n";
		for(int i=0;i<cleanFiles[4].size();i++)
			fileClean<<"rm -f "<<cleanFiles[4][i]<<"\n";//4
		fileClean<<"fi \n";
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<fileName<<" could not be created"<<std::endl;
		exit(1);
	}
	fileClean.close();
}

