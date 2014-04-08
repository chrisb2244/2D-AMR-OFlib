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



#ifndef regenerateAlphaClass_H
#define regenerateAlphaClass_H

#include "dynamicFvMesh.H"
#include <vector>

namespace Foam
{
// Forward declaration of classes
class polyMesh;

/*---------------------------------------------------------------------------*\
                    Class regenerateAlphaClass Declaration
\*---------------------------------------------------------------------------*/

class regenerateAlphaClass
{
public:
	class cosineFunctor
	{
		public:
			cosineFunctor(double kChosen, double randScaling, double phase);
			double operator()(double x);
		private:
			double k_;
			double scaling_;
			double phase_;
			double PI;
	};
	
	//- Runtime type information
    TypeName("regenerateAlphaClass");

	Foam::volScalarField regenerateAlpha();
	
	void printInfo();

	regenerateAlphaClass(const dynamicFvMesh& mesh, const double seed, const double yMid);
	virtual ~regenerateAlphaClass();
	
private:
	//- Disallow default bitwise copy construct
	regenerateAlphaClass(const regenerateAlphaClass&);

	//- Disallow default bitwise assignment
	void operator=(const regenerateAlphaClass&);

	void calculateProfileHeight();
	void createCosineFunctors
	(
		std::vector<cosineFunctor>& cosineVector,
		const double& k,
		const double& scale,
		const double& phase
	);
	
	// Data values (private)
	const fvMesh& mesh_;
	DynamicList<double> kList_;
	DynamicList<double> scalingList_;
	DynamicList<double> phaseList_;
	std::vector<cosineFunctor> cosineVector;
	double cellDepth_;
	const double randomSeed_;
	const double yMid_;
	double kStart_;
	double kEnd_;
	double kSpacing_;
	int numKs_;
	double yScaling_;
};



} // End of namespace Foam
#endif
