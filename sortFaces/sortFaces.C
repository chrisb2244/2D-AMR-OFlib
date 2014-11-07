/*---------------------------------------------------------------------------*\
 libmyDynamicMesh Copyright (C) 2014 Christian Butcher
 chrisb2244@gmail.com

License
	This file is part of a library, libmyDynamicMesh, using and derived 
	from OpenFOAM.
	
	OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    This work is distributed under the same licensing conditions.
    
    You should have received a copy of the GNU General Public License
    along with this library.  If not, see <http://www.gnu.org/licenses/>.
    
\*---------------------------------------------------------------------------*/


#include "sortFaces.H"	// Header file for class
//~ #include <vector>					// Needed for vector of functors
#include "fvCFD.H"					// Completes types for patchFields

namespace Foam
{
	defineTypeNameAndDebug(sortFaces, 0);
}

Foam::sortFaces::sortFaces(DynamicList<Pair<label> >& sourceList)
:
	values_(sourceList)
{
} 


Foam::sortFaces::~sortFaces()
{}

void Foam::sortFaces::printInfo()
{
	Pout<< values_;
}

label Foam::sortFaces::partition(
        const label& left, const label& right, const label& pivotIndex)
{
	label pivotValue = values_[pivotIndex].first();
	swap(values_[pivotIndex], values_[right]);
	label storeIndex = left;
	for (int i = left; i < right; i++)
	{
		if (values_[i].first() <= pivotValue)
		{
			swap(values_[i], values_[storeIndex]);
			storeIndex ++;
		}
	}
	swap(values_[storeIndex], values_[right]);
	return storeIndex;
}

void Foam::sortFaces::quicksort(const label& left, const label& right)
{
	if (left < right)
	{
		label pivotIndex = (left+right)/2;
		label pivotNewIndex = partition(left, right, pivotIndex);
		quicksort(left, pivotNewIndex - 1);
		quicksort(pivotNewIndex + 1, right);
	}
}
	
void Foam::sortFaces::swap(Pair<label> &a, Pair<label> &b)
{
	Pair<label> temp = a;
	a = b;
	b = temp;
}

void Foam::sortFaces::sort()
{
	quicksort(0, values_.size());
}

bool Foam::sortFaces::bCondense()
{
	bool changed = false;
	DynamicList<Pair<label> > originalList(Cvalues_);
	if (originalList.size() == 0)
	{
		Pout<< "No list" << endl;
	}
	Cvalues_.clear();
	if (Cvalues_.size() != 0)
	{
		Pout<< "Didn't clear" << endl;
	}
	for(int i = 0; i < (originalList.size()-1); i++)
	{
		if (originalList[i+1].first() == originalList[i].first())
		{
            originalList[i].second() =
                    (originalList[i].second() + originalList[i+1].second());
			Cvalues_.append(originalList[i]);
			i++;
			changed = true;
			continue;
		}
		Cvalues_.append(originalList[i]);
	}
	if (!changed)
	{
		//~ forAll(originalList, i)
		//~ Cvalues_.append(originalList[i]);
		Cvalues_=originalList;
	}
	return changed;
}

DynamicList<Pair<label> > Foam::sortFaces::condense()
{
	bool again = true;
	Cvalues_ = values_;
	while (again)
	{
		again = bCondense();
	}
	return Cvalues_;
}
