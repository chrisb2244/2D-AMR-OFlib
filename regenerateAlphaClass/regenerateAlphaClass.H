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
