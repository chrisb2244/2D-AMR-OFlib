#ifndef DIAG_List_H
#define DIAG_List_H

#include "dictionary.H"
#include <map>

namespace Foam
{

class DIAG_List
{
	public:
	ClassName("DIAG_List");
	
	std::map<std::string, bool*> myMap;
	
	// Constructor
	DIAG_List(dictionary &diagDict);
	
	void addToList(const std::string lookupName, bool &boolName);
	
	bool lookupValue(const std::string lookupName);
	
	void reportBools();
	
	// The [] operator now returns the value (as a bool) of the given
	// variable.
	const bool& operator[](const std::string name);
	
	private:
	dictionary myDict;
	const bool default_Diag;
	const bool recursive_Search;
	const bool pattern_Matching;
	const bool CONST_FALSE;

};

} // End namespace Foam
#endif
