#include "DIAG_List.H"
#include "IOstreams.H"
#include <map>

Foam::DIAG_List::DIAG_List(dictionary &diagDict):
default_Diag(false),
recursive_Search(false),
pattern_Matching(false),
CONST_FALSE(false)
{
	myDict = diagDict;
}

const bool& Foam::DIAG_List::operator[](const std::string name)
{
	std::map<std::string, bool*>::iterator iter = myMap.find(name);
	if (iter != myMap.end()) {
		return *(iter->second);
	}
	else {
		FatalErrorIn("myReadDiagDict")
			<< "The name passed to the [] operator of the DIAG_List "
			<< "class was invalid."
			<< abort(FatalError);
	}
	return CONST_FALSE;
}
	
void Foam::DIAG_List::reportBools()
{
	Pout<< endl << "******************************" << endl;
	Pout<< "Diagnostic levels are recorded with values: " << endl;
	
	for(std::map<std::string, bool*>::iterator iter = myMap.begin(); iter != myMap.end(); iter++)
	{
		Pout<< iter->first << " = " << *(iter->second) << endl;
	}
	
	Pout<< "******************************" << endl;
}

void Foam::DIAG_List::addToList(const std::string lookupName, bool &boolName)
{
	boolName = lookupValue(lookupName);
	myMap[lookupName] = &boolName;
}

bool Foam::DIAG_List::lookupValue(const std::string lookupName)
{
	return myDict.lookupOrDefault<bool> (lookupName, default_Diag, recursive_Search, pattern_Matching);
}
