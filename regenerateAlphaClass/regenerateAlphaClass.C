#include "regenerateAlphaClass.H"	// Header file for class
#include <vector>					// Needed for vector of functors
#include <fstream>					// Needed to read the yPerturbationDict
#include "functions.H"				// Needed to read the yPerturbationDict
#include "fvCFD.H"					// Completes types for patchFields

namespace Foam
{
	defineTypeNameAndDebug(regenerateAlphaClass, 0);
}

Foam::regenerateAlphaClass::regenerateAlphaClass(const dynamicFvMesh& mesh, const double seed, const double yMid) 
:
	mesh_(mesh),
	randomSeed_(seed),
	yMid_(yMid)
{
	// Read these values from a dictionary
	std::string yScalingFactorLine, kStartLine, kEndLine, numKsLine, line;
	std::ifstream yPert;
	yPert.open("system/yPerturbationDict");

	while (yPert.good())
	{
		getline(yPert, line);
		lineFinder(yScalingFactorLine,"yScalingFactor",line);
		lineFinder(kStartLine,"kStart",line);
		lineFinder(kEndLine,"kEnd",line);
		lineFinder(numKsLine,"numKs",line);
	}

	extractNumber(kStart_, "kStart =", kStartLine);
	extractNumber(kEnd_, "kEnd =", kEndLine);
	extractNumber(yScaling_, "yScalingFactor =", yScalingFactorLine);
	extractNumber(numKs_, "numKs =", numKsLine);
	
	if (debug)
	{
		Pout<< "The values of kStart, kEnd, yScalingFactor are " << kStart_
			<< ", " << kEnd_ << ", " << yScaling_ << ".\n"
			<< numKs_ << " k values will be used to form the perturbation." << endl;
	}
	
	// Need to get this still.
	cellDepth_ = 1.0;

	// Calculated values
	kSpacing_ = (kEnd_ - kStart_)/numKs_;

	
	// Function calls
	srand(randomSeed_);
	if (1) printInfo();
	calculateProfileHeight();
} 

Foam::regenerateAlphaClass::~regenerateAlphaClass()
{}

Foam::volScalarField Foam::regenerateAlphaClass::regenerateAlpha()
{
	if (debug) Pout<< "Calling regenerateAlphaClass::regenerateAlpha()" << endl;
	// Going to set values of the data points in a constructor.
	// The constructor also sets up a vector of cosineFunctors to get heights for given x.
	// This needs to work out which cells should be given which values,
	// and write an alpha1 file.
	volScalarField aBound
    (
        IOobject
        (
            "alpha.water.boundaryFields",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_
    );
    
	wordList wantedTypes;
	forAll(aBound.boundaryField(), i)
	{
		wantedTypes.append(aBound.boundaryField().types()[i]);
	}
	
	if (debug) Pout<< "wantedTypes = " << wantedTypes << endl;
	
	volScalarField alpha1
	(
		IOobject
		(
			"alpha.water",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE,
			true // registry
		),
		mesh_,
		dimensionedScalar
		(
			"alpha.water",
			dimensionSet(0, 0, 0, 0, 0, 0, 0),
			0.0
		),
		wantedTypes
	);
	if (debug) Pout<< "Field created " << endl;

	double yHeight = yMid_; // This value is passed in at construction from the calling program
	forAll(mesh_.cellCentres(), i)
	{
		double cellHeight = sqrt(mesh_.cellVolumes()[i] / cellDepth_);
		yHeight = yMid_;
		Foam::vector position(mesh_.cellCentres()[i]);
		for (unsigned int k=0; k < cosineVector.size(); k++)
		{
			yHeight += cosineVector[k](position[0]);
		}
		double deltaY = yHeight - position[1];
		if (deltaY < -cellHeight)
		{
			alpha1.internalField()[i] = 0.0;
		}
		else if (deltaY > cellHeight)
		{
			alpha1.internalField()[i] = 1.0;
		}
		else
		{
			double value_Alpha = (0.5+(0.5 * (deltaY/cellHeight)));
			alpha1.internalField()[i] = value_Alpha;
		}
	}
	
	return alpha1;
}

void Foam::regenerateAlphaClass::printInfo()
{
	Pout<< nl;
	Pout<< "Lowest perturbed wavenumber (kStart) = " << kStart_ << endl;
	Pout<< "Highest perturbed wavenumber (kEnd) = " << kEnd_ << endl;
	Pout<< "The spacing between perturbed wavenumbers (kSpacing) = " << (kEnd_ - kStart_)/numKs_ << " or " << kSpacing_ << endl;
	Pout<< "The scaling factor applied to the perturbation (yScaling) = " << yScaling_ << endl;
	Pout<< endl;
}

class cosineFunctor
{
	private:
		double k_;
		double scaling_;
		double phase_;
		double PI;
};

Foam::regenerateAlphaClass::cosineFunctor::cosineFunctor(double kChosen, double randScaling, double phase)
:
	k_(kChosen),
	scaling_(randScaling),
	phase_(phase),
	PI(3.14159265359)
	{}
	
double Foam::regenerateAlphaClass::cosineFunctor::operator()(double x)
{
	return (scaling_ * cos((2*PI*k_*x) + phase_));
}
	
void Foam::regenerateAlphaClass::calculateProfileHeight()
{
	const double PI = 3.14159265359;
	for (int i=0; i<numKs_; i++)
	{
		kList_.append(kStart_ + (kSpacing_ * i));
		scalingList_.append((yScaling_*rand())/RAND_MAX);
		//~ phaseList_.append(0);
		phaseList_.append((2*PI*rand())/RAND_MAX);
	}
	forAll (kList_, k)
	{
		createCosineFunctors(cosineVector, kList_[k], scalingList_[k], phaseList_[k]);
	}
}

void Foam::regenerateAlphaClass::createCosineFunctors
(
	std::vector<cosineFunctor>& cosineVector,
	const double& k,
	const double& scale,
	const double& phase
)
{
	cosineFunctor cF(k, scale, phase);
	cosineVector.push_back(cF);
}
