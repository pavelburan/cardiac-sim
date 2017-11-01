// ConfigFile.cpp

#include "ConfigFile.h"

using std::string;

ConfigFile::ConfigFile( string filename, string delimiter, string addDelimiter,string comment, string terminator,
						string seperator, string LeftBrace, string RightBrace, string sentry )
	: myDelimiter(delimiter), myAddDelimiter(addDelimiter), myComment(comment), myTerminator(terminator), mySeperator(seperator), myLeftBrace(LeftBrace), myRightBrace(RightBrace), mySentry(sentry), myContents()
{
	// Construct a ConfigFile, getting keys and values from given file
	std::ifstream in( filename.c_str() );
	
	if( !in ) throw file_not_found( filename ); 

	in >> (*this);
}


ConfigFile::ConfigFile()
	: myDelimiter( string("=") ), myAddDelimiter( string("+=") ), myComment( string("#") ), myTerminator( string(";") ), mySeperator( string(",") ), myLeftBrace( string("[") ), myRightBrace( string("]") ), mySentry("EndConfigFile"), myContents() 
{
	// Construct a ConfigFile without a file; empty
}

ConfigFile::~ConfigFile(){}

void ConfigFile::remove( const string& key )
{
	// Remove key and its value
	myContents.erase( myContents.find( key ) );
	return;
}


bool ConfigFile::keyExists( const string& key ) const
{
	// Indicate whether key is found
	mapci p = myContents.find( key );
	return ( p != myContents.end() );
}


/* static */
void ConfigFile::trim( string& s )
{
	// Remove leading and trailing whitespace
	static const char whitespace[] = " \n\t\v\r\f";
	s.erase( 0, s.find_first_not_of(whitespace) );
	s.erase( s.find_last_not_of(whitespace) + 1U );
}

/*string ConfigFile::getElemFromArrayString( string str, int index, const string& key ) const
{	
	string elem;
	for(unsigned int i=0; i<str.length(); ++i)
		if(str[i] == '[' || str[i] == ']' || str[i] == '(' || str[i] == ')')
			str.erase(i--,1);
	char tempChar[512]="";
	str.copy(tempChar,512,0);
	char *token;
	token = strtok(tempChar,",");
	elem = token;
	int i = 1;
	token = strtok(0, ",");
	while (token != 0 )
	{
		if(i == index){
			elem = token;
			return elem;
			//break;
		}
		token = strtok(0, ",");
		++i;
	}
	throw index_not_found(index,key);
	return elem;
}*/


int ConfigFile::elemExist( std::vector<int> index, const string& key ) const
{	
	mapci p = myContents.find(key);
	if( p == myContents.end() ) return 0;
	string str(p->second);
	int n = index.size();
	index.push_back(0); 
	std::string delimiterBeg(n,'[');
	std::string delimiterEnd(n,']');
	int pos;
	int posBeg = 0;
	int posEnd = str.length();
	int returnFlag = 0;
	for(int i=n;i>0;i--){
		pos = str.find(delimiterBeg.c_str(),posBeg,i);
		if(pos == std::string::npos){
			if( i!=n ){
				returnFlag |= 1<<i;
			}
			continue;
		}
		posBeg = pos+1;
		for(int j=1;j<=index[i];j++){
			pos = str.find(delimiterBeg.c_str(),posBeg,i);
			if(pos == std::string::npos || pos > posEnd){
				if(j==1){
					returnFlag |= 1<<i;
					break;//return -1;
				}
				else
					return 0;
			}
			posBeg = pos+1;	
		}
		posEnd = str.find(delimiterEnd.c_str(),posBeg,i);
	}
	for(int j=1;j<=index[0];j++){
		pos = str.find(',',posBeg);
		if(pos == std::string::npos || pos >= posEnd){
			if(j==1){
				returnFlag |= 1<<0;
				break;//return -1;
			}
			else
				return 0;
		}
		posBeg = pos+1;
	}
	pos = str.find(',',posBeg);
	if(pos != std::string::npos && pos < posEnd)
			posEnd = pos;
	return (returnFlag > 0) ? (-1*returnFlag) : 1;
}


string ConfigFile::getElemFromArrayString( string str, std::vector<int> index, const string& key ) const
{	
	int n = index.size();
	index.push_back(0); 
	std::string delimiterBeg(n,'[');
	std::string delimiterEnd(n,']');
	int pos;
	int posBeg = 0;
	int posEnd = str.length();
	for(int i=n;i>0;i--){
		pos = str.find(delimiterBeg.c_str(),posBeg,i);
		if(pos == std::string::npos){
			continue;
		}
		posBeg = pos+1;
		for(int j=1;j<=index[i];j++){
			pos = str.find(delimiterBeg.c_str(),posBeg,i);
			if(pos == std::string::npos || pos > posEnd){
				if(j==1)
					break;
				else{
					std::cerr<<"i="<<i<<" j="<<j<<std::endl;
					throw index_not_found(index,key);
				}
			}
			posBeg = pos+1;	
		}
		posEnd = str.find(delimiterEnd.c_str(),posBeg,i);
	}

	for(int j=1;j<=index[0];j++){
		pos = str.find(',',posBeg);
		if(pos == std::string::npos || pos >= posEnd){
			if(j==1)
				break;
			else{
				std::cerr<<"j="<<j<<std::endl;
				throw index_not_found(index,key);
			}
		}
		posBeg = pos+1;
	}
	pos = str.find(',',posBeg);
	if(pos != std::string::npos && pos < posEnd)
			posEnd = pos;
	
	string elem(str,posBeg,posEnd-posBeg);
	return elem;
}


std::ostream& operator<<( std::ostream& os, const ConfigFile& cf )
{
	// Save a ConfigFile to os
	for( ConfigFile::mapci p = cf.myContents.begin();
	     p != cf.myContents.end();
		 ++p )
	{
		os << p->first << " " << cf.myDelimiter << " ";
		os << p->second << std::endl;
	}
	return os;
}


std::istream& operator>>( std::istream& is, ConfigFile& cf )
{
	// Load a ConfigFile from is
	// Read in keys and values, keeping internal whitespace
	typedef string::size_type pos;
	const string& delim  = cf.myDelimiter;  // delimiter
	const string& addDelim   = cf.myAddDelimiter;    // adddelimiter
	const string& comm   = cf.myComment;    // comment
	const string& terminator  = cf.myTerminator;  // terminator
	const string& seperator  = cf.mySeperator;  // seperator
	const string& lBrace  = cf.myLeftBrace;  // leftBrace
	const string& rBrace  = cf.myRightBrace;  // rightBrace
	const string& sentry = cf.mySentry;     // end of file sentry
	const pos delimSkip = delim.length();        // length of delimiter
	const pos addDelimSkip = addDelim.length();        // length of delimiter
	string nextline = "";  // might need to read ahead to see where value ends
	
	while( is || nextline.length() > 0 )
	{
		// Read an entire line at a time
		string line;
		if( nextline.length() > 0 )
		{
			line = nextline;  // we read ahead; use it now
			nextline = "";
		}
		else
		{
			std::getline( is, line );
		}
		
		// Ignore comments
		line = line.substr( 0, line.find(comm) );
		
		// Check for end of file sentry
		if( sentry != "" && line.find(sentry) != string::npos ) return is;
		
		// Parse the line if it contains a delimiter
		pos delimPos = line.find( delim );
		pos addDelimPos = line.find( addDelim );
		if( delimPos != string::npos || addDelimPos != string::npos)
		{
			// Extract the key
			string key;
			if( addDelimPos != string::npos ){
				key = line.substr( 0, addDelimPos );
				line.replace( 0, addDelimPos+addDelimSkip, "" );
			}
			else{
				key = line.substr( 0, delimPos );
				line.replace( 0, delimPos+delimSkip, "" );
			}
			// See if value continues on the next line
			// Stop at blank line, next line with a key, end of stream,
			// or end of file sentry
			bool terminate = false;
			while( !terminate && is )
			{
				std::getline( is, nextline );
				terminate = true;
				
				string nlcopy = nextline;
				ConfigFile::trim(nlcopy);
				if( nlcopy == "" ) continue;
				
				nextline = nextline.substr( 0, nextline.find(comm) );
				if( nextline.find(delim) != string::npos )
					continue;
				if( sentry != "" && nextline.find(sentry) != string::npos )
					continue;
				
				nlcopy = nextline;
				ConfigFile::trim(nlcopy);
				if( nlcopy != "" ) line += "\n";
				line += nextline;
				terminate = false;
			}
			
			// Store key and value
			ConfigFile::trim(key);
			ConfigFile::trim(line);
			pos terminatorPos = line.find( terminator );
			if( terminatorPos == line.length()-terminator.length() && line.length()>=terminator.length())
				line.erase(terminatorPos, terminator.length());
			if(addDelimPos != string::npos){
				if(cf.myContents[key] != "")
					cf.myContents[key] += seperator;
				cf.myContents[key] += line;
			}
			else{
				if(line.find(lBrace) != string::npos && line.find(rBrace) != string::npos && line.length()-lBrace.length()-rBrace.length() == 0 ){
					cf.myContents[key] = "";	// "[]" = ""
				}
				else{
					cf.myContents[key] = line;  // overwrites if key is repeated
				}
			}
		}
	}
	
	return is;
}
