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
    
Description
	This class (regenerateAlphaClass) is used to generate cosine-wave
	based perturbations to an interface between two fluids, for solution
	with solvers such as interDyMFoam or interFoam.
	
	A dictionary named "yPerturbationDict" contains information about
	the shape and extent of the perturbation. An example yPerturbationDict
	would be as follows:
	
	// -------------------------------------- \\

	kStart = 1 	// The lowest wavenumber to perturb
	kEnd = 5	// The highest wavenumber to perturb
	numKs = 5	// The number of wavenumbers between 
				// kStart and kEnd to perturb
	yScalingFactor = 0.15
				// The amplitude to use (in combination
				// with a random (seeded) number)
				// as a prefactor for the perturbation

	\\ -------------------------------------- //

\*---------------------------------------------------------------------------*/



#ifndef sortFaces_H
#define sortFaces_H

#include "DynamicList.H"
#include "UList.H"

namespace Foam
{
// Forward declaration of classes
//~ class Pair<typename>
/*---------------------------------------------------------------------------*\
                    Class sortFaces Declaration
\*---------------------------------------------------------------------------*/

class sortFaces
{
public:
	
	//- Runtime type information
    TypeName("sortFaces");
	
	void printInfo();

	sortFaces(DynamicList<Pair<label> >& sourceList); // Arguments to constructor?
	virtual ~sortFaces();
	void sort();
	DynamicList<Pair<label> > condense();
	
private:
	//- Disallow default bitwise copy construct
	sortFaces(const sortFaces&);

	//- Disallow default bitwise assignment
	void operator=(const sortFaces&);

	List<Pair<label> > values_;
	DynamicList<Pair<label> > Cvalues_;
	bool bCondense();

	
	label partition(const label& left, const label& right, const label& pivotIndex);
	void quicksort(const label& left, const label& right);
	
	void swap(Pair<label> &a, Pair<label> &b);
	

};



} // End of namespace Foam
#endif
