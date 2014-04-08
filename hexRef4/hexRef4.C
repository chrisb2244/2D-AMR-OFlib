/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

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

#include "hexRef4.H"

#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyAddFace.H"
#include "polyAddPoint.H"
#include "polyAddCell.H"
#include "polyModifyFace.H"
#include "syncTools.H"
#include "faceSet.H"
#include "cellSet.H"
#include "pointSet.H"
#include "OFstream.H"
#include "Time.H"
#include "FaceCellWave.H"
#include "mapDistributePolyMesh.H"
#include "refinementData.H"
#include "refinementDistanceData.H"
#include "degenerateMatcher.H"

#include "primitiveMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hexRef4, 0);

    //- Reduction class. If x and y are not equal assign value.
    template<int value>
    class ifEqEqOp
    {
        public:
        void operator()(label& x, const label y) const
        {
            x = (x == y) ? x : value;
        }
    };
}
bool cellLevel2Debug = false;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::vector Foam::hexRef4::calcSingleFaceNormal(const pointField& p, const label& faceI) const
{
	const labelList& f = mesh_.faces()[faceI];
	point fCentre = p[f[0]];
	
	vector sumN = vector::zero;
	scalar sumA = 0.0;
	
	label nPoints = f.size();

	for (label pi = 1; pi < nPoints; pi++)
	{
		fCentre += p[f[pi]];
	}
	fCentre /= nPoints;
	
	for (label pi = 0; pi < nPoints; pi++)
	{
		const point& nextPoint = p[f[(pi + 1) % nPoints]];
		vector n = (nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);
		scalar a = mag(n);
		sumN += n;
		sumA += a;
	}
	// This is to deal with zero-area faces. Mark very small faces
	// to be detected in e.g., processorPolyPatch.
	if (sumA < ROOTVSMALL)
	{
		return vector::zero;
	}
	else
	{
		return 0.5*sumN;
	}
}

double Foam::hexRef4::getCellDepth()
{
	// Assumes all the cells are equally deep.
	// Might not be true, but likely much larger problems if it isn't!
	DynamicList<scalar> depths(0);
	label point0 = mesh_.cellPoints()[0][0];
	point pt0 = mesh_.points()[point0];
	forAll(mesh_.cellPoints()[0], pt)
	{
		label pointI = mesh_.cellPoints()[0][pt];
		point ptI = mesh_.points()[pointI];
		vector vec(ptI - pt0);
		if (mag(vec[2]) > 1E-14)
		{
			if (depths.size() == 0)
			{
				depths.append(mag(vec[2]));
			}
			else if (depths.size() == 1)
			{
				if (depths[0] != mag(vec[2]))
				{
					FatalErrorIn("getCellDepth()")
						<< "The value of z computed between point "
						<< pointI
						<< " and point 0 of cell 0 was not the same as for an earlier point "
						<< "with non-zero difference in z."
						<< abort(FatalError);
				}
			}
		}
	}
	if (debug) Pout<< "Cell depth = " << depths[0] << endl;
	return depths[0];
}

double Foam::hexRef4::getAspectRatio()
{
	// Assumes all cells have the same aspect ratio (although not size)
	// Could also be wrong, but shouldn't be (and if so, probably problems elsewhere)
	double area = mesh_.cellVolumes()[0] / cellDepth_;
	if (debug) Pout<< "Cell surface area (for cell0) is " << area << endl;
	label point0 = mesh_.cellPoints()[0][0];
	point pt0 = mesh_.points()[point0];
	forAll(mesh_.cellPoints()[0], pt)
	{
		label pointI = mesh_.cellPoints()[0][pt];
		point ptI = mesh_.points()[pointI];
		vector vec(ptI - pt0);
		if (mag(vec[0]) > 1E-14 && mag(vec[1]) > 1E-14)
		{
			if (debug)
			{
				Pout<< "Cell aspect ratio of X:Y = " << mag(vec[0]) << ":" << mag(vec[1]) << endl;
				Pout<< "Setting aspect value to " << (mag(vec[1])/mag(vec[0])) << endl;
			}
			return (mag(vec[1])/mag(vec[0]));
		}
	}
	FatalErrorIn("getAspectRatio()")
		<< "Couldn't find a point with both X and Y larger than point0 of cell0 "
		<< "by at least 1E-10"
		<< abort(FatalError);
	return -1;
}


int Foam::hexRef4::calcRelevantDirs(const Foam::vector& dir, const bool Fatal=true)
{
	DynamicList<label> relevDirs;
	scalar relevantDir = 12345;
	for (int i = 0; i < 3; i++)
	{
		if (mag(dir[i]) > 1E-10)
		{
			relevDirs.append(i+1); // 1 for x, 2 for y, 3 for z.
		}
	}
	if (relevDirs.size() != 1)
	{
		if (Fatal)
		{
			FatalErrorIn("hexRef4::calcRelevantDirs(..)")
				<< "Error finding the relevant direction. "
				<< "The number of non-zero(ish) directions != 1. "
				<< nl
				<< "mag(dir) = " << vector(mag(dir[0]), mag(dir[1]), mag(dir[2])) 
				<< nl 
				<< "dir = " << dir
				<< abort(FatalError);
		}
		else
		{
			return 100;
		}
	}
	else
	{
		relevantDir = ((dir[relevDirs[0]-1] > 0) ? relevDirs[0] : (-1*relevDirs[0]));
	}
	
	if (relevantDir == 12345)
	{
		FatalErrorIn("hexRef4::calcRelevantDirs(..)")
			<< "relevantDir was not set during the function. This error should not have been reached"
			<< abort(FatalError);
	}
	return relevantDir;
}

void Foam::hexRef4::setIndicators(const label& own, const label& nei, label& ownInd, label& neiInd, const bool refinementBool[])
{
	if (refinementBool[0] && refinementBool[1])
	{
		if (debug) Pout<< "Both owner and neighbour passed to setIndicators are refined" << endl;
		neiInd = cellLevel_[nei] > 0 ? history_.myParentCell(nei) : nei;
		ownInd = cellLevel_[own] > 0 ? history_.myParentCell(own) : own;
	}
	else if (refinementBool[0] && !refinementBool[1])
	{
		if (debug) Pout<< "Owner passed to setIndicators is refined" << endl;
		neiInd = nei;
		ownInd = cellLevel_[own] > 0 ? history_.myParentCell(own) : own;
	}
	else if (!refinementBool[0] && refinementBool[1])
	{
		if (debug) Pout<< "Neighbour passed to setIndicators is refined" << endl;
		neiInd = cellLevel_[nei] > 0 ? history_.myParentCell(nei) : nei;
		ownInd = own;
	}
	else
	{
		FatalErrorIn("hexRef4::setIndicators(..)")
			<< "It seems that neither cell is refined "
			<< "(or at least, that neither bool is true)"
			<< nl
			<< abort(FatalError);
	}
}

void Foam::hexRef4::setBoundaryCellInfo(const label& faceI, const labelList& parentAddedCells, label& own, const int a)
{
	if (parentAddedCells.size() == 0) // Not sure why cells that are not refined are being called here.
	{
		Pout<< "Warning: setBoundaryCellInfo called with parentAddedCells.size() == 0, "
			<< "ie, the parent cell is not refined in this timestep." << endl;
		return;
	}
	if (parentAddedCells.size() != 4)
	{
		FatalErrorIn("hexRef4::setBoundaryCellInfo(..)")
			<< "The cellAddedCells list was not of size 4. "
			<< parentAddedCells
			<< abort(FatalError);
	}

	vector faceNormal = calcSingleFaceNormal(mesh_.points(), faceI);
	int relevantDir = calcRelevantDirs(faceNormal);


	if (debug) Pout << "setBoundaryCellInfo for a = " << a << ", case " << relevantDir << endl;
	if (a==1)
	{
		
		switch (relevantDir)
		{
			case -1: // left edge
			own = parentAddedCells[2];
			break;
			
			case +1: // right edge
			own = parentAddedCells[3];
			break;
			
			case -2: // bottom edge
			own = parentAddedCells[1];
			break;
			
			case +2: // top edge
			own = parentAddedCells[3];
			break;
			
			default:
			FatalErrorIn("hexRef4::setRefinement(..)")
				<< "Could not match the relevantDir found with an expected case (external face)"
				<< nl
				<< "faceI : " << faceI
				<< ", relevantDir : " << relevantDir
				<< ", faceNormal = " << faceNormal
				<< abort(FatalError);
			break;
		}
	}
	else if (a==0)
	{
		switch (relevantDir)
		{
			case -1: // left edge
			own = parentAddedCells[0];
			break;
			
			case +1: // right edge
			own = parentAddedCells[1];
			break;
			
			case -2: // bottom edge
			own = parentAddedCells[0];
			break;
			
			case +2: // top edge
			own = parentAddedCells[2];
			break;
			
			default:
			FatalErrorIn("hexRef4::setRefinement(..)")
				<< "Could not match the relevantDir found with an expected case (external face)"
				<< nl
				<< "faceI : " << faceI
				<< ", relevantDir : " << relevantDir
				<< ", faceNormal = " << faceNormal
				<< abort(FatalError);
			break;
		}
	}
	else
	{
		FatalErrorIn("setBoundaryCellInfo(..)")
			<< "The value of a passed was not equal to 0 or 1. (a = " << a
			<< ")" << abort(FatalError);
	}
}

Foam::vector Foam::hexRef4::setInternalDir
(
	const labelListList& cellAddedCells,
	const bool refinementBool[],
	const label& own,
	const label& nei,
	const label& ownInd,
	const label& neiInd,
	const polyTopoChange& meshMod
)
{
	vector dir(0,0,0);
	if (refinementBool[2] || refinementBool[3])
	{
		// Need to set vector dir more carefully.
		if (refinementBool[2]) // ownMoreRefined
		{
			Pout<< "ownMoreRefined" << endl;
			label indexToNei = 0;
			label parentCell = history_.myParentCell(nei);
			forAll(cellAddedCells[parentCell], i)
			{
				if (cellAddedCells[parentCell][i] == nei)
				{
					indexToNei = i;
					break;
				}
			}
			if (indexToNei == 0)
			{
				FatalErrorIn("setInternalFaceInfoA0(..)")
					<< "Found the neighbour cell to be its own parent when claiming ownMoreRefined"
					<< abort(FatalError);
			}
			point neiCentre(mesh_.cellCentres()[parentCell]);
			double quartX = (mesh_.cellVolumes()[parentCell] / cellDepth_);
			double quartY = (sqrt(quartX * aspectX_to_Y_))/4.0;
			quartX = (sqrt(quartX/aspectX_to_Y_))/4.0;
			neiCentre[0] = (indexToNei == 2 ? neiCentre[0] - quartX : neiCentre[0] + quartX);
			neiCentre[1] = (indexToNei == 1 ? neiCentre[1] - quartY : neiCentre[1] + quartY);
			if (debug)
			{
				Pout<< "parentCell = " << parentCell << endl;
				Pout<< "ownInd = " << ownInd << " and cellCentres()[ownInd] = " << mesh_.cellCentres()[ownInd] << endl;
				Pout<< "Original centre : " << mesh_.cellCentres()[parentCell] << endl;
				Pout<< "Using neiCentre = " << neiCentre << endl;
			}

			dir = neiCentre - mesh_.cellCentres()[ownInd];
		}
		else // neiMoreRefined
		{
			Pout<< "neiMoreRefined" << endl;
			label indexToOwn = 0;
			label parentCell = history_.myParentCell(own);
			forAll(cellAddedCells[parentCell], i)
			{
				if (cellAddedCells[parentCell][i] == own)
				{
					indexToOwn = i;
					if (debug) Pout<< "i = " << i << endl;
					break;
				}
			}
			if (indexToOwn == 0)
			{
				FatalErrorIn("setInternalFaceInfoA0(..)")
					<< "Found the owner cell to be its own parent when claiming neiMoreRefined"
					<< abort(FatalError);
			}
			if (debug) Pout<< "cellVolume for parentCell " << parentCell << " : " << mesh_.cellVolumes()[parentCell] << endl;
			point ownCentre(mesh_.cellCentres()[parentCell]);
			double quartX = (mesh_.cellVolumes()[parentCell] / cellDepth_);
			double quartY = (sqrt(quartX * aspectX_to_Y_))/4.0;
			quartX = (sqrt(quartX/aspectX_to_Y_))/4.0;
			ownCentre[0] = (indexToOwn == 2 ? ownCentre[0] - quartX : ownCentre[0] + quartX);
			ownCentre[1] = (indexToOwn == 1 ? ownCentre[1] - quartY : ownCentre[1] + quartY);
			if (debug)
			{
				Pout<< "parentCell = " << parentCell << endl;
				Pout<< "neiInd = " << neiInd << " and cellCentres()[neiInd] = " << mesh_.cellCentres()[neiInd] << endl;
				Pout<< "Original centre : " << mesh_.cellCentres()[parentCell] << endl;
				Pout<< "Using ownCentre = " << ownCentre << endl;
			}

			dir = mesh_.cellCentres()[neiInd] - ownCentre;
		}
	}
	else
	{
		dir = (mesh_.cellCentres()[neiInd] - mesh_.cellCentres()[ownInd]);
	}
	return dir;
}

void Foam::hexRef4::setInternalFaceInfoA1
(
	const labelListList& cellAddedCells, 
	const bool refinementBool[],
	label& own, 
	label& nei,
	label& ownInd,
	label& neiInd,
	const polyTopoChange& meshMod
)
{
	if (debug) Pout<< "Entered setInternalFaceInfoA1" << endl;
	vector dir(setInternalDir(cellAddedCells, refinementBool, own, nei, ownInd, neiInd, meshMod));
	int relevantDir = calcRelevantDirs(dir);
	if (debug) Pout<< "A1: Case " << relevantDir << endl;
	
	switch (relevantDir)
	{
		case -1:
		if (refinementBool[0]) // 0, +1
		{
			own = cellAddedCells[ownInd][2];
		}
		if (refinementBool[1]) // 0, -1
		{
			nei = cellAddedCells[neiInd][3];
		}
		break;
		
		case +1:
		if (refinementBool[0]) // 0, +1
		{
			own = cellAddedCells[ownInd][3];
		}
		if (refinementBool[1]) // 0, -1
		{
			nei = cellAddedCells[neiInd][2];
		}
		break;
		
		case -2:
		if (refinementBool[0] && refinementBool[1]) // 0
		{
			nei = cellAddedCells[neiInd][3];
			own = cellAddedCells[ownInd][1];
		}
		else if (!refinementBool[0] && refinementBool[1]) // -1
		{
			nei = cellAddedCells[neiInd][3];
		}
		else // +1
		{
			own = cellAddedCells[ownInd][1];
		}
		break;
		
		case +2:
		if (refinementBool[0] && refinementBool[1]) // 0
		{
			own = cellAddedCells[ownInd][3];
			nei = cellAddedCells[neiInd][1];
		}
		else if (!refinementBool[0] && refinementBool[1]) // -1
		{
			own = ownInd;
			nei = cellAddedCells[neiInd][1];
		}
		else // +1
		{
			own = cellAddedCells[ownInd][3];
		}
		break;

		default:
		FatalErrorIn("hexRef4::setInternalFaceInfoA1(..)")
			<< "Could not match the relevantDir found with an expected case"
			<< nl
			<< "ownInd : " << ownInd
			<< ", neiInd : " << neiInd
			<< ", relevantDir : " << relevantDir
			<< ", dir = " << dir
			<< abort(FatalError);
		break;
	}
}

void Foam::hexRef4::setInternalFaceInfoA0
(
	const labelListList& cellAddedCells, 
	const bool refinementBool[],
	label& own, 
	label& nei,
	label& ownInd,
	label& neiInd,
	const polyTopoChange& meshMod
)
{
	if (debug) Pout<< "Entered setInternalFaceInfoA0" << endl;
	vector dir(setInternalDir(cellAddedCells, refinementBool, own, nei, ownInd, neiInd, meshMod));
	int relevantDir = calcRelevantDirs(dir);
	
	if (debug) Pout<< "A0: Case " << relevantDir << endl;
	switch (relevantDir)
	{
		case -1:
		if (refinementBool[0]) // 0, +1
		{
			own = cellAddedCells[ownInd][0];
		}
		if (refinementBool[1]) // 0, -1
		{
			nei = cellAddedCells[neiInd][1];
		}
		break;
		
		case +1:
		if (refinementBool[0]) // 0, +1
		{
			own = cellAddedCells[ownInd][1];
		}
		if (refinementBool[1]) // 0, -1
		{
			nei = cellAddedCells[neiInd][0];
		}
		break;
		
		case -2:
		if (refinementBool[0] && refinementBool[1]) // 0
		{
			nei = cellAddedCells[neiInd][2];
			own = cellAddedCells[ownInd][0];
		}
		else if (!refinementBool[0] && refinementBool[1]) // -1
		{
			nei = cellAddedCells[neiInd][2];
		}
		else // +1
		{
			own = cellAddedCells[ownInd][0];
		}
		break;
		
		case +2:
		if (refinementBool[0] && refinementBool[1]) // 0
		{
			own = cellAddedCells[ownInd][2];
			nei = cellAddedCells[neiInd][0];
		}
		else if (!refinementBool[0] && refinementBool[1]) // -1
		{
			own = ownInd;
			nei = cellAddedCells[neiInd][0];
		}
		else // +1
		{
			own = cellAddedCells[ownInd][2];
		}
		break;

		default:
		FatalErrorIn("hexRef4::setInternalFaceInfoA0(..)")
			<< "Could not match the relevantDir found with an expected case"
			<< nl
			<< "ownInd : " << ownInd
			<< ", neiInd : " << neiInd
			<< ", relevantDir : " << relevantDir
			<< ", dir = " << dir
			<< abort(FatalError);
		break;
	}
}

void Foam::hexRef4::calcFaceNormalVector(const pointField& p, vectorField& fAreas) const
{
	const faceList& fs = mesh_.faces();

    forAll(fs, facei)
    {
        const labelList& f = fs[facei];
        label nPoints = f.size();

        // If the face is a triangle, do a direct calculation for efficiency
        // and to avoid round-off error-related problems
        if (nPoints == 3)
        {
            fAreas[facei] = 0.5*((p[f[1]] - p[f[0]])^(p[f[2]] - p[f[0]]));
        }
        else
        {
            vector sumN = vector::zero;
            scalar sumA = 0.0;

            point fCentre = p[f[0]];
            for (label pi = 1; pi < nPoints; pi++)
            {
                fCentre += p[f[pi]];
            }

            fCentre /= nPoints;

            for (label pi = 0; pi < nPoints; pi++)
            {
                const point& nextPoint = p[f[(pi + 1) % nPoints]];

                vector n = (nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);
                scalar a = mag(n);

                sumN += n;
                sumA += a;
            }

            // This is to deal with zero-area faces. Mark very small faces
            // to be detected in e.g., processorPolyPatch.
            (sumA < ROOTVSMALL) ? fAreas[facei] = vector::zero : fAreas[facei] = 0.5*sumN;
            //~ {
                //~ fAreas[facei] = vector::zero;
            //~ }
            //~ else
            //~ {
                //~ fAreas[facei] = 0.5*sumN;
            //~ }
        }
    }
}

void Foam::hexRef4::mySplitSideFaces
(
	const labelListList& cellAddedCells,
	const labelList& faceMidPoint,
	const labelList& edgeMidPoint,
	const label& faceI,
	const label& oldOwn,
	const label& oldNei,
	polyTopoChange& meshMod
)
{
	// Values of oldOwn and oldNei are those given by getFaceNeighbours for faceI
	DynamicList<label> fEdgesStorage;
	
	const face& f = mesh_.faces()[faceI];
	const labelList& fEdges = mesh_.faceEdges
	(
		faceI,
		fEdgesStorage
	);
	
	DynamicList<label> newFaceVertsA[2];
	{ // set two loops of vertices for the side faces.
		label fp = 0;
		newFaceVertsA[0].append(f[fp]);
		if (edgeMidPoint[fEdges[fp]] >= 0)
		{
			newFaceVertsA[0].append(edgeMidPoint[fEdges[fp]]);
			newFaceVertsA[1].append(edgeMidPoint[fEdges[fp]]);
			fp = f.fcIndex(fp);
			newFaceVertsA[1].append(f[fp]);
			fp = f.fcIndex(fp);
			newFaceVertsA[1].append(f[fp]);
			if (edgeMidPoint[fEdges[fp]] >= 0)
			{
				newFaceVertsA[0].append(edgeMidPoint[fEdges[fp]]);
				newFaceVertsA[1].append(edgeMidPoint[fEdges[fp]]);
				fp = f.fcIndex(fp);
				newFaceVertsA[0].append(f[fp]);
			}
			else
			{
				FatalErrorIn("mySplitSideFaces(..)")
					<< "A mid point was not found opposite an already found mid point"
					<< abort(FatalError);
			}
		}
		else
		{
			fp = f.fcIndex(fp);
			newFaceVertsA[0].append(f[fp]);
			if (edgeMidPoint[fEdges[fp]] >= 0)
			{
				newFaceVertsA[0].append(edgeMidPoint[fEdges[fp]]);
				newFaceVertsA[1].append(edgeMidPoint[fEdges[fp]]);
				fp = f.fcIndex(fp);
				newFaceVertsA[1].append(f[fp]);
				fp = f.fcIndex(fp);
				newFaceVertsA[1].append(f[fp]);
				if (edgeMidPoint[fEdges[fp]] >= 0)
				{
					newFaceVertsA[0].append(edgeMidPoint[fEdges[fp]]);
					newFaceVertsA[1].append(edgeMidPoint[fEdges[fp]]);
				}
				else
				{
					FatalErrorIn("mySplitSideFaces(..)")
					<< "A mid point was not found opposite an already found mid point"
					<< abort(FatalError);
				}
			}
			else
			{
				FatalErrorIn("mySplitSideFaces(..)")
					<< "A mid point was not found on either of two consecutive edges"
					<< abort(FatalError);
			}
		}
	}
	
	if ((newFaceVertsA[0].size() != 4) 
	 || (newFaceVertsA[1].size() != 4))
	{
		FatalErrorIn("hexRef4::setRefinement(..)")
			<< "newFaceVertsA does not have 4 vertices for one of its loops"
			<< nl
			<< "[0] : " << newFaceVertsA[0] << endl
			<< "[1] : " << newFaceVertsA[1] << endl
			<< abort(FatalError);
	}

	vector centre0(0,0,0), centre1(0,0,0);
	forAll(newFaceVertsA[0], i)
	{
		centre0 += meshMod.points()[newFaceVertsA[0][i]];
	}
	forAll(newFaceVertsA[1], i)
	{
		centre1 += meshMod.points()[newFaceVertsA[1][i]];
	}
	centre0 /= 4;
	centre1 /= 4;
	
	vector dir(centre1 - centre0);
	if (dir[0] < -1E-10 || dir[1] < -1E-10 || dir[2] < -1E-10)
	{
		if (debug) Pout<< "Switching vectors newFaceVertsA - original dir = " << dir << endl;
		DynamicList<label> tempVerts = newFaceVertsA[0];
		newFaceVertsA[0] = newFaceVertsA[1];
		newFaceVertsA[1] = tempVerts;
		vector tempCentre = centre0;
		centre0 = centre1;
		centre1 = tempCentre;
	}
	
	label own = oldOwn;
	label nei = oldNei;
	if (!mesh_.isInternalFace(faceI))
	{
		nei = -1;
		
		if (debug) Pout<< "Boundary cell - Initial own/nei are " << oldOwn << ", " << oldNei << endl;
		const labelList& pACells = (cellLevel_[oldOwn] > 0) ? cellAddedCells[history_.myParentCell(oldOwn)]:cellAddedCells[oldOwn];
		//~ const labelList& pACells = ((cellLevel_[oldOwn] - oldCellLevel_[oldOwn]) != 0) ? cellAddedCells[history_.myParentCell(oldOwn)]:cellAddedCells[oldOwn];
		if (debug) Pout<< "pACells = " << pACells << endl;
		setBoundaryCellInfo(faceI, pACells, own, 1);
		face newFace;
		newFace.transfer(newFaceVertsA[1]);
		if (debug) Pout<< "Adding boundary face (a=1) with final own/nei = " << own << ", " << nei << endl;
		addFace(meshMod, faceI, newFace, own, nei);
		
		setBoundaryCellInfo(faceI, pACells, own, 0);
		face sameFace;
		sameFace.transfer(newFaceVertsA[0]);
		if (debug) Pout<< "Modding boundary face (a=0) with final own/nei = " << own << ", " << nei << endl;
		modFace(meshMod, faceI, sameFace, own, nei);
	}
	else // internal cell
	{
		if (debug) Pout<< "Internal cell - Initial own/nei are " << oldOwn << ", " << oldNei << endl;
		int cellLevelDiff = cellLevel_[oldOwn] - cellLevel_[oldNei]; // care with own/nei here
		label ownInd, neiInd;
		bool refinementBool[4];
		bool& ownRefined = refinementBool[0];
		bool& neiRefined = refinementBool[1];
		bool& ownMoreRefined = refinementBool[2];
		bool& neiMoreRefined = refinementBool[3];
		ownRefined = false;
		neiRefined = false;
		ownMoreRefined = false;
		neiMoreRefined = false;
		switch (cellLevelDiff)
		{
			case 0:
			ownRefined = true;
			//~ refinementBool[0]=1;
			neiRefined = true;
			//~ refinementBool[1]=1;
			break;
			
			case -1:
			neiRefined = true;
			//~ refinementBool[1]=1;
			// new section
			forAll(cellAddedCells, cellI)
			{
				forAll(cellAddedCells[cellI], i)
				{
					if ((cellAddedCells[cellI][i] == oldOwn) && (cellLevel_[oldOwn] > 0) && (i != 0))
					{
						neiMoreRefined = true;
						break;
					}
				}
			}
			break;
			
			case +1:
			ownRefined = true;
			//~ refinementBool[0]=1;
			// new section
			forAll(cellAddedCells, cellI)
			{
				forAll(cellAddedCells[cellI], i)
				{
					if ((cellAddedCells[cellI][i] == oldNei) && (cellLevel_[oldNei] > 0) && (i != 0)) // I think this only works if nei is being refined this timestep.
					{
						ownMoreRefined = true;
						break;
					}
				}
			}
			break;
			
			default:
			FatalErrorIn("hexRef4::setRefinement(..)")
				<< "Change in cell level greater than 1, "
				<< "and not a border cell. Fails to obey "
				<< "the 2-1 refinement rule"
				<< abort(FatalError);
			break;
		}
		if (refinementBool[2] && refinementBool[3])
		{
			FatalErrorIn("hexRef4::setRefinement(..)")
				<< "Claims that both owner and neighbour are "
				<< "more refined than the other. "
				<< "refinementBool[4] = "
				<< refinementBool[0]
				<< refinementBool[1]
				<< refinementBool[2]
				<< refinementBool[3]
				<< abort(FatalError);
		}
		
		bool requiresSwitch = false;		
		
		setIndicators(own, nei, ownInd, neiInd, refinementBool);
		
		// A1
		//~ setInternalFaceInfoA1(faceI, cellAddedCells, ownRefined, neiRefined, own, nei, requiresSwitch, ownInd, neiInd);
		setInternalFaceInfoA1(cellAddedCells, refinementBool, own, nei, ownInd, neiInd, meshMod);
		if (own > nei)
		{
			label temp = own;
			own = nei;
			nei = temp;
			temp = ownInd;
			ownInd = neiInd;
			neiInd = temp;
			//~ if (requiresSwitch)
			//~ {
				//~ requiresSwitch = false;
			//~ }
			//~ else
			//~ {
				if (debug) Pout<< "A1. Setting requiresSwitch = true, because own > nei" << endl;
				requiresSwitch = true;
			//~ }
		}
		if (requiresSwitch)
		{
			if (debug) Pout<< "Switching due to requiresSwitch = true" << endl;
			label temp = newFaceVertsA[1][0];
			newFaceVertsA[1][0] = newFaceVertsA[1][2];
			newFaceVertsA[1][2] = temp;
		}
		face newFace;
		newFace.transfer(newFaceVertsA[1]);
		if (debug)
		{
			if (refinementBool[2])
			{
			}
			else if (refinementBool[3])
			{
			}
			else
			{		
				checkInternalOrientation
				(
					meshMod,
					own,
					faceI,
					mesh_.cellCentres()[ownInd],
					mesh_.cellCentres()[neiInd],
					newFace
				);
			}
		}
		if (debug) Pout<< "Adding face (a=1) with final own/nei = " << own << ", " << nei << endl;
		addFace(meshMod, faceI, newFace, own, nei);
		
		// A0
		own = oldOwn; // reset own and nei
		nei = oldNei;
		requiresSwitch = false; // reset requiresSwitch
		
		setIndicators(own, nei, ownInd, neiInd, refinementBool);
		
		//~ setInternalFaceInfoA0(faceI, cellAddedCells, ownRefined, neiRefined, own, nei, requiresSwitch, ownInd, neiInd);
		setInternalFaceInfoA0(cellAddedCells, refinementBool, own, nei, ownInd, neiInd, meshMod);
		if (own > nei)
		{
			label temp = own;
			own = nei;
			nei = temp;
			temp = ownInd;
			ownInd = neiInd;
			neiInd = temp;
			//~ if (requiresSwitch)
			//~ {
				//~ requiresSwitch = false;
			//~ }
			//~ else
			//~ {
				if (debug) Pout<< "A0. Setting requiresSwitch = true, because own > nei" << endl;
				requiresSwitch = true;
			//~ }
		}
		if (requiresSwitch)
		{
			if (debug) Pout<< "2. Switching due to requiresSwitch = true" << endl;
			label temp = newFaceVertsA[0][0];
			newFaceVertsA[0][0] = newFaceVertsA[0][2];
			newFaceVertsA[0][2] = temp;
		}
		face sameFace;
		sameFace.transfer(newFaceVertsA[0]);
		if (debug)
		{
			if (refinementBool[2])
			{
			}
			else if (refinementBool[3])
			{
			}
			else
			{	
				checkInternalOrientation
				(
					meshMod,
					own,
					faceI,
					mesh_.cellCentres()[ownInd],
					mesh_.cellCentres()[neiInd],
					sameFace
				);
			}
		}
		if (debug) Pout<< "Modding face (a=0) with final own/nei = " << own << ", " << nei << endl;
		modFace(meshMod, faceI, sameFace, own, nei);
	}
}

void Foam::hexRef4::myCreateInternalFaces
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& faceMidPoint,
    const labelList& edgeMidPoint,
    const label& cellI,
    polyTopoChange& meshMod,
    const Foam::vector& normalDir
)
{
	label nFacesAdded = 0;
	DynamicList<label> storage;
	
	if (debug) Pout<< nl << "cellI = " << cellI << endl;
	
	const cell& cFaces = mesh_.cells()[cellI];

	label facesProcessed = 0;
	DynamicList<label> cF(2);
	DynamicList<label> ncF(4);
	forAll(cFaces, i)
	{
		label faceI = cFaces[i];
		Foam::vector faceNormal = calcSingleFaceNormal(mesh_.points(), faceI);
		faceNormal /= mag(faceNormal);
		
		if (mag(faceNormal & normalDir) > 0.5)
		{
			cF.append(faceI);
			facesProcessed++;
		}
		else
		{
			ncF.append(faceI);
			facesProcessed++;
		}
		
		if (facesProcessed > 6)
		{
			if (cF.size() != 2)
			{
				FatalErrorIn("hexRef4::myCreateInternalFaces(..)")
					<< "There were not two faces chosen as front and back. "
					<< "cF.size() = " 
					<< cF.size()
					<< abort(FatalError);
			}
			Pout<< "Warning: More than 6 faces processed. Check for errors due to this!" << endl;
		}
	}

	label edgeMidPoint1 = -1;
	label edgeMidPoint2 = -1;
	const label faceMidPoint1 = faceMidPoint[cF[0]];
	const label faceMidPoint2 = faceMidPoint[cF[1]];
	const labelList& edgesOnFace1 = mesh_.faceEdges(cF[0], storage);
	const labelList& edgesOnFace2 = mesh_.faceEdges(cF[1], storage);
	
	if (debug)
	{
		Pout<< "edgesOnFace1 = " << edgesOnFace1 << endl;
		Pout<< "edgesOnFace2 = " << edgesOnFace2 << endl;
	}
	
	DynamicList<label> newFaceVerts(4);
	DynamicList<label> facesToCombine;
	
	forAll(ncF, i)
	{
		label faceI = ncF[i];
		vector directionToFaceFromCentre(mesh_.faceCentres()[faceI] - mesh_.cellCentres()[cellI]);
		// The false here causes a return 100, instead of FatalError, if the relevantDir cannot be found
		label relevantDir = calcRelevantDirs(directionToFaceFromCentre, false);
		if (relevantDir == 100)
		{
			// face is offset
			facesToCombine.append(faceI);
			if (debug) Pout<< "Added faceI = " << faceI << " to a list of faces to combine. List: " << facesToCombine << endl;
			continue;
		}
			
		label eP = -1;
		const labelList& fEdges = mesh_.faceEdges(faceI, storage);
		
		if (debug)
		{
			Pout<< "The points forming the target face are " << mesh_.faces()[faceI] << endl;
			Pout<< "The edges of the target face are " << fEdges << endl;
		}
		
		forAll(edgesOnFace1, eI)
		{
			eP = findEdge(fEdges, edgesOnFace1[eI]);
			if (eP != -1)
			{
				edgeMidPoint1 = edgeMidPoint[fEdges[eP]];
				if (debug) Pout<< "edgeMidPoint1 = " << edgeMidPoint1 << endl;
				break;
			}
		}

		forAll(edgesOnFace2, eI)
		{
			eP = findEdge(fEdges, edgesOnFace2[eI]);
			if (eP != -1)
			{
				edgeMidPoint2 = edgeMidPoint[fEdges[eP]];
				if (debug) Pout<< "edgeMidPoint2 = " << edgeMidPoint2 << endl;
				break;
			}
		}
	
		newFaceVerts.append(edgeMidPoint1);
		newFaceVerts.append(edgeMidPoint2);
		newFaceVerts.append(faceMidPoint2);
		newFaceVerts.append(faceMidPoint1);
		
		myCreateInternalFacesSubSection
		(
			newFaceVerts,
			relevantDir,
			cellI,
			faceI,
			cellAddedCells,
			meshMod,
			nFacesAdded
		);
	}
	
	// Add face for the split sides
	if (facesToCombine.size()%2 != 0)
	{
		FatalErrorIn("myCreateInternalFaces(..)")
			<< "The list of faces to combine does not contain an even number of faces. "
			<< "List: " << facesToCombine << abort(FatalError);
	}
	label numPairedFaces = facesToCombine.size() / 2;
	
	// Better handling now
	DynamicList<Pair<label> > pairedFaces;
	if (numPairedFaces == 0)
	{
		goto MY_CREATE_INTERNAL_FACES_END;
	}
	
	forAll(facesToCombine, i)
	{
		label refDir = calcRelevantDirs(calcSingleFaceNormal(mesh_.points(), facesToCombine[i]));
		forAll(pairedFaces, k)
		{
			if (facesToCombine[i] == pairedFaces[k].second())
			{
				if (debug) Pout<< "facesToCombine[i] == pairedFaces[k].second() : " << facesToCombine[i] << endl;
				break;
			}
		}
		for(int j = i+1; j<facesToCombine.size(); j++)
		{
			if
			(
				calcRelevantDirs(calcSingleFaceNormal(mesh_.points(), facesToCombine[j])) == refDir
			 &&	(
					mesh_.faceCentres()[facesToCombine[i]][0] == mesh_.faceCentres()[facesToCombine[j]][0]
				 || mesh_.faceCentres()[facesToCombine[i]][1] == mesh_.faceCentres()[facesToCombine[j]][1]
				)
			)
			{
				pairedFaces.append(Pair<label>(facesToCombine[i], facesToCombine[j]));
			}
		}
	}	
	if (pairedFaces.size() != numPairedFaces)
	{
		FatalErrorIn("myCreateInternalFaces()")
			<< "The list of paired faces does not have the same number of entries " 
			<< "as the value of numPairedFaces. "
			<< nl
			<< "numPairedFaces = " << numPairedFaces << nl
			<< "pairedFaces = " << pairedFaces << nl
			<< abort(FatalError);
	}
	
	if (debug) Pout<< "Adding faces between paired faces: numPairs = " << numPairedFaces << endl;
	forAll(pairedFaces, i)
	{
		DynamicList<label> sharedPoints(2);
		const labelList& pF1 = mesh_.faces()[pairedFaces[i].first()];
		const labelList& pF2 = mesh_.faces()[pairedFaces[i].second()];

		forAll(pF1, pt)
		{
			
			forAll(pF2, pt2)
			{
				if (pF1[pt] == pF2[pt2])
				{
					sharedPoints.append(pF1[pt]);
				}
			}
		}
		if (sharedPoints.size() != 2)
		{
			FatalErrorIn("myCreateInternalFaces(..)")
				<< "The number of shared points found was not equal to 2. "
				<< "Shared points: "
				<< sharedPoints
				<< abort(FatalError);
		}
		newFaceVerts.append(sharedPoints[0]);
		newFaceVerts.append(sharedPoints[1]);
		newFaceVerts.append(faceMidPoint2);
		newFaceVerts.append(faceMidPoint1);

		point newCentre(0,0,0);
		
		label faceI = pairedFaces[i].first();
		forAll(mesh_.faces()[faceI], pt)
		{
			label point = mesh_.faces()[faceI][pt];
			newCentre += mesh_.points()[point];
			Pout<< "Added point " << mesh_.points()[point] << endl;
		}
		faceI = pairedFaces[i].second();
		forAll(mesh_.faces()[faceI], pt)
		{
			label point = mesh_.faces()[faceI][pt];
			newCentre += mesh_.points()[point];
			Pout<< "Added point " << mesh_.points()[point] << endl;
		}
		
		//~ forAll(facesToCombine, i)
		//~ forAll(pairedFaces[i], k)
		//~ {
			//~ label faceI = facesToCombine[i];
			//~ forAll(mesh_.faces()[faceI], pt)
			//~ {
				//~ label point = mesh_.faces()[faceI][pt];
				//~ newCentre += mesh_.points()[point];
				//~ Pout<< "Added point " << mesh_.points()[point] << endl;
			//~ }
		//~ }
		// Fix this ^^
		
		newCentre /= 8;
		if (debug)
		{
			Pout<< "newFaceVerts = " << newFaceVerts
				<< " and newCentre = " << newCentre << endl;
		}
		vector directionToFaceFromCentre(newCentre - mesh_.cellCentres()[cellI]);
		label relevantDir = calcRelevantDirs(directionToFaceFromCentre);
		
		myCreateInternalFacesSubSection
		(
			newFaceVerts,
			relevantDir,
			cellI,
			pairedFaces[i].first(),
			cellAddedCells,
			meshMod,
			nFacesAdded
		);

		if (debug) Pout<< "directionToFaceFromCentre = " << directionToFaceFromCentre << endl;
	}
		
	MY_CREATE_INTERNAL_FACES_END:
	
	if (nFacesAdded != 4)
	{
		FatalErrorIn("myCreateInternalFaces(..)")
			<< "The number of faces added by myCreateInternalFaces "
			<< "was not equal to 4. Instead, " << nFacesAdded
			<< " were added."
			<< abort(FatalError);
	}
}

void Foam::hexRef4::myCreateInternalFacesSubSection
(
	DynamicList<label>& newFaceVerts,
	const label& relevantDir,
	const label& cellI,
	const label& faceI,
	const labelListList& cellAddedCells,
	polyTopoChange& meshMod,
	label& nFacesAdded
)
{
	if (debug)
	{
		Pout<< "The newFaceVerts (for added internal faces) are " << newFaceVerts << endl;
		Pout<< "relevantDir = " << relevantDir << endl;
	}
	bool requiresSwitch = false;
	label own = -50, nei = -50;
	switch (relevantDir)
	{
		case +1:
		own = cellAddedCells[cellI][1];
		nei = cellAddedCells[cellI][3];
		break;
		
		case -1:
		own = cellAddedCells[cellI][0];
		nei = cellAddedCells[cellI][2];
		requiresSwitch = true;
		break;
		
		case +2:
		own = cellAddedCells[cellI][2];
		nei = cellAddedCells[cellI][3];
		requiresSwitch = true;
		break;
		
		case -2:
		own = cellAddedCells[cellI][0];
		nei = cellAddedCells[cellI][1];
		break;
		
		default:
		FatalErrorIn("myCreateInternalFaces(..)")
			<< "Could not match the relevantDir : " << relevantDir
			<< " to an expected direction (x, y)"
			<< abort(FatalError);
		break;
	}

	if (requiresSwitch)
	{
		label temp = newFaceVerts[0];
		newFaceVerts[0] = newFaceVerts[2];
		newFaceVerts[2] = temp;
	}
	
	face newFace;
	newFace.transfer(newFaceVerts);
	
	if (debug)
	{
		face compactFace(identity(newFace.size()));
		pointField compactPoints(meshMod.points(), newFace);
		Pout<< "The coordinates of the vertices are: " << compactPoints << endl;
	}
	label faceLabel = -1;
	faceLabel = meshMod.setAction
	(
		polyAddFace
		(
			newFace,                    // face
			own,                        // owner
			nei,                        // neighbour
			-1,                         // master point
			-1,                         // master edge
			faceI,                  	// master face for addition
			false,                      // flux flip
			-1,                         // patch for face
			-1,                         // zone for face
			false                       // face zone flip
		)
	);
	if (own == -50 || nei == -50)
	{
		FatalErrorIn("hexRef4::myCreateInternalFaces(..)") 
			<< "own and nei were not changed by the switch options"
			<< abort(FatalError);
	}
	nFacesAdded++;
	if (debug) Pout<< "A face was added between new cells own/nei = " << meshMod.faceOwner()[faceLabel] << ", " << meshMod.faceNeighbour()[faceLabel] << endl;
}

Foam::label Foam::hexRef4::findEdge
(
	const labelList& targetList,
	const label targetLabel
)
{
	forAll(targetList, i)
	{
		if(targetList[i] == targetLabel)
		{
			return i;
		}
	}
	return -1;
}

void Foam::hexRef4::reorder
(
    const labelList& map,
    const label len,
    const label null,
    labelList& elems
)
{
    labelList newElems(len, null);

    forAll(elems, i)
    {
        label newI = map[i];

        if (newI >= len)
        {
            FatalErrorIn("hexRef4::reorder(..)") << abort(FatalError);
        }

        if (newI >= 0)
        {
            newElems[newI] = elems[i];
        }
    }

    elems.transfer(newElems);
}


void Foam::hexRef4::getFaceInfo
(
    const label faceI,
    label& patchID,
    label& zoneID,
    label& zoneFlip
) const
{
    patchID = -1;

    if (!mesh_.isInternalFace(faceI))
    {
        patchID = mesh_.boundaryMesh().whichPatch(faceI);
    }

    zoneID = mesh_.faceZones().whichZone(faceI);

    zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
    }
}


// Adds a face on top of existing faceI.
Foam::label Foam::hexRef4::addFace
(
    polyTopoChange& meshMod,
    const label faceI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(faceI, patchID, zoneID, zoneFlip);

    label newFaceI = -1;

    if ((nei == -1) || (own < nei))
    {
        // Ordering ok.
        newFaceI = meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                faceI,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
    else
    {
        // Reverse owner/neighbour
        newFaceI = meshMod.setAction
        (
            polyAddFace
            (
                newFace.reverseFace(),      // face
                nei,                        // owner
                own,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                faceI,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
    return newFaceI;
}


// Adds an internal face from an edge. Assumes orientation correct.
// Problem is that the face is between four new vertices. So what do we provide
// as master? The only existing mesh item we have is the edge we have split.
// Have to be careful in only using it if it has internal faces since otherwise
// polyMeshMorph will complain (because it cannot generate a sensible mapping
// for the face)
Foam::label Foam::hexRef4::addInternalFace
(
    polyTopoChange& meshMod,
    const label meshFaceI,
    const label meshPointI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    if (mesh_.isInternalFace(meshFaceI))
    {
        return meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                meshFaceI,                  // master face for addition
                false,                      // flux flip
                -1,                         // patch for face
                -1,                         // zone for face
                false                       // face zone flip
            )
        );
    }
    else
    {
        // Two choices:
        // - append (i.e. create out of nothing - will not be mapped)
        //   problem: field does not get mapped.
        // - inflate from point.
        //   problem: does interpolative mapping which constructs full
        //   volPointInterpolation!

        // For now create out of nothing

        return meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                -1,                         // master face for addition
                false,                      // flux flip
                -1,                         // patch for face
                -1,                         // zone for face
                false                       // face zone flip
            )
        );


        ////- Inflate-from-point:
        //// Check if point has any internal faces we can use.
        //label masterPointI = -1;
        //
        //const labelList& pFaces = mesh_.pointFaces()[meshPointI];
        //
        //forAll(pFaces, i)
        //{
        //    if (mesh_.isInternalFace(pFaces[i]))
        //    {
        //        // meshPoint uses internal faces so ok to inflate from it
        //        masterPointI = meshPointI;
        //
        //        break;
        //    }
        //}
        //
        //return meshMod.setAction
        //(
        //    polyAddFace
        //    (
        //        newFace,                    // face
        //        own,                        // owner
        //        nei,                        // neighbour
        //        masterPointI,               // master point
        //        -1,                         // master edge
        //        -1,                         // master face for addition
        //        false,                      // flux flip
        //        -1,                         // patch for face
        //        -1,                         // zone for face
        //        false                       // face zone flip
        //    )
        //);
    }
}


// Modifies existing faceI for either new owner/neighbour or new face points.
void Foam::hexRef4::modFace
(
    polyTopoChange& meshMod,
    const label faceI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(faceI, patchID, zoneID, zoneFlip);

    if
    (
        (own != mesh_.faceOwner()[faceI])
     || (mesh_.isInternalFace(faceI) && (nei != mesh_.faceNeighbour()[faceI]))
     || (newFace != mesh_.faces()[faceI])
    )
    {
        if ((nei == -1) || (own < nei))
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace,            // modified face
                    faceI,              // label of face being modified
                    own,                // owner
                    nei,                // neighbour
                    false,              // face flip
                    patchID,            // patch for face
                    false,              // remove from zone
                    zoneID,             // zone for face
                    zoneFlip            // face flip in zone
                )
            );
        }
        else
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace.reverseFace(),  // modified face
                    faceI,                  // label of face being modified
                    nei,                    // owner
                    own,                    // neighbour
                    false,                  // face flip
                    patchID,                // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }
    }
}


// Bit complex way to determine the unrefined edge length.
Foam::scalar Foam::hexRef4::getLevel0EdgeLength() const
{
    if (cellLevel_.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "hexRef4::getLevel0EdgeLength() const"
        )   << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size()
            << endl
            << "This might be because of a restart with inconsistent cellLevel."
            << abort(FatalError);
    }

    // Determine minimum edge length per refinement level
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar GREAT2 = sqr(GREAT);

    label nLevels = gMax(cellLevel_)+1;

    scalarField typEdgeLenSqr(nLevels, GREAT2);


    // 1. Look only at edges surrounded by cellLevel cells only.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Per edge the cellLevel of connected cells. -1 if not set,
        // labelMax if different levels, otherwise levels of connected cells.
        labelList edgeLevel(mesh_.nEdges(), -1);

        forAll(cellLevel_, cellI)
        {
            const label cLevel = cellLevel_[cellI];

            const labelList& cEdges = mesh_.cellEdges(cellI);

            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];

                if (edgeLevel[edgeI] == -1)
                {
                    edgeLevel[edgeI] = cLevel;
                }
                else if (edgeLevel[edgeI] == labelMax)
                {
                    // Already marked as on different cellLevels
                }
                else if (edgeLevel[edgeI] != cLevel)
                {
                    edgeLevel[edgeI] = labelMax;
                }
            }
        }

        // Make sure that edges with different levels on different processors
        // are also marked. Do the same test (edgeLevel != cLevel) on coupled
        // edges.
        syncTools::syncEdgeList
        (
            mesh_,
            edgeLevel,
            ifEqEqOp<labelMax>(),
            labelMin
        );

        // Now use the edgeLevel with a valid value to determine the
        // length per level.
        forAll(edgeLevel, edgeI)
        {
            const label eLevel = edgeLevel[edgeI];

            if (eLevel >= 0 && eLevel < labelMax)
            {
                const edge& e = mesh_.edges()[edgeI];

                scalar edgeLenSqr = magSqr(e.vec(mesh_.points()));

                typEdgeLenSqr[eLevel] = min(typEdgeLenSqr[eLevel], edgeLenSqr);
            }
        }
    }

    // Get the minimum per level over all processors. Note minimum so if
    // cells are not cubic we use the smallest edge side.
    Pstream::listCombineGather(typEdgeLenSqr, minEqOp<scalar>());
    Pstream::listCombineScatter(typEdgeLenSqr);

    if (debug)
    {
        Pout<< "hexRef4::getLevel0EdgeLength() :"
            << " After phase1: Edgelengths (squared) per refinementlevel:"
            << typEdgeLenSqr << endl;
    }


    // 2. For any levels where we haven't determined a valid length yet
    //    use any surrounding cell level. Here we use the max so we don't
    //    pick up levels between celllevel and higher celllevel (will have
    //    edges sized according to highest celllevel)
    //    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField maxEdgeLenSqr(nLevels, -GREAT2);

    forAll(cellLevel_, cellI)
    {
        const label cLevel = cellLevel_[cellI];

        const labelList& cEdges = mesh_.cellEdges(cellI);

        forAll(cEdges, i)
        {
            const edge& e = mesh_.edges()[cEdges[i]];

            scalar edgeLenSqr = magSqr(e.vec(mesh_.points()));

            maxEdgeLenSqr[cLevel] = max(maxEdgeLenSqr[cLevel], edgeLenSqr);
        }
    }

    Pstream::listCombineGather(maxEdgeLenSqr, maxEqOp<scalar>());
    Pstream::listCombineScatter(maxEdgeLenSqr);

    if (debug)
    {
        Pout<< "hexRef4::getLevel0EdgeLength() :"
            << " Crappy Edgelengths (squared) per refinementlevel:"
            << maxEdgeLenSqr << endl;
    }


    // 3. Combine the two sets of lengths
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(typEdgeLenSqr, levelI)
    {
        if (typEdgeLenSqr[levelI] == GREAT2 && maxEdgeLenSqr[levelI] >= 0)
        {
            typEdgeLenSqr[levelI] = maxEdgeLenSqr[levelI];
        }
    }

    if (debug)
    {
        Pout<< "hexRef4::getLevel0EdgeLength() :"
            << " Final Edgelengths (squared) per refinementlevel:"
            << typEdgeLenSqr << endl;
    }

    // Find lowest level present
    scalar level0Size = -1;

    forAll(typEdgeLenSqr, levelI)
    {
        scalar lenSqr = typEdgeLenSqr[levelI];

        if (lenSqr < GREAT2)
        {
            level0Size = Foam::sqrt(lenSqr)*(1<<levelI);

            if (debug)
            {
                Pout<< "hexRef4::getLevel0EdgeLength() :"
                    << " For level:" << levelI
                    << " have edgeLen:" << Foam::sqrt(lenSqr)
                    << " with equivalent level0 len:" << level0Size
                    << endl;
            }
            break;
        }
    }

    if (level0Size == -1)
    {
        FatalErrorIn("hexRef4::getLevel0EdgeLength()")
            << "Problem : typEdgeLenSqr:" << typEdgeLenSqr << abort(FatalError);
    }

    return level0Size;
}


// Check whether pointI is an anchor on cellI.
// If it is not check whether any other point on the face is an anchor cell.
Foam::label Foam::hexRef4::getAnchorCell
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label cellI,
    const label faceI,
    const label pointI
) const
{
    if (cellAnchorPoints[cellI].size())
    {
        label index = (findIndex(cellAnchorPoints[cellI], pointI)%4) ;

        if (index != -1)
        {
            return (cellAddedCells[cellI][index]);
        }


        // pointI is not an anchor cell.
        // Maybe we are already a refined face so check all the face
        // vertices.
        const face& f = mesh_.faces()[faceI];

        forAll(f, fp)
        {
            label index = findIndex(cellAnchorPoints[cellI], f[fp]);

            if (index != -1)
            {
                return cellAddedCells[cellI][index];
            }
        }

        // Problem.
        dumpCell(cellI);
        Perr<< "cell:" << cellI << " anchorPoints:" << cellAnchorPoints[cellI]
            << endl;

        FatalErrorIn("hexRef4::getAnchorCell(..)")
            << "Could not find point " << pointI
            << " in the anchorPoints for cell " << cellI << endl
            << "Does your original mesh obey the 2:1 constraint and"
            << " did you use consistentRefinement to make your cells to refine"
            << " obey this constraint as well?"
            << abort(FatalError);

        return -1;
    }
    else
    {
        return cellI;
    }
}


// Get new owner and neighbour
void Foam::hexRef4::getFaceNeighbours
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label faceI,
    const label pointI,

    label& own,
    label& nei
) const
{
    // Is owner split?
    own = getAnchorCell
    (
        cellAnchorPoints,
        cellAddedCells,
        mesh_.faceOwner()[faceI],
        faceI,
        pointI
    );

    if (mesh_.isInternalFace(faceI))
    {
        nei = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            mesh_.faceNeighbour()[faceI],
            faceI,
            pointI
        );
    }
    else
    {
        nei = -1;
    }
}


// Get point with the lowest pointLevel
Foam::label Foam::hexRef4::findMinLevel(const labelList& f) const
{
    label minLevel = labelMax;
    label minFp = -1;

    forAll(f, fp)
    {
        label level = pointLevel_[f[fp]];

        if (level < minLevel)
        {
            minLevel = level;
            minFp = fp;
        }
    }

    return minFp;
}


// Get point with the highest pointLevel
Foam::label Foam::hexRef4::findMaxLevel(const labelList& f) const
{
    label maxLevel = labelMin;
    label maxFp = -1;

    forAll(f, fp)
    {
        label level = pointLevel_[f[fp]];

        if (level > maxLevel)
        {
            maxLevel = level;
            maxFp = fp;
        }
    }

    return maxFp;
}


Foam::label Foam::hexRef4::countAnchors
(
    const labelList& f,
    const label anchorLevel
) const
{
    label nAnchors = 0;

    forAll(f, fp)
    {
        if (pointLevel_[f[fp]] <= anchorLevel)
        {
            nAnchors++;
        }
    }
    return nAnchors;
}


void Foam::hexRef4::dumpCell(const label cellI) const
{
    OFstream str(mesh_.time().path()/"cell_" + Foam::name(cellI) + ".obj");
    Pout<< "hexRef4 : Dumping cell as obj to " << str.name() << endl;

    const cell& cFaces = mesh_.cells()[cellI];

    Map<label> pointToObjVert;
    label objVertI = 0;

    forAll(cFaces, i)
    {
        const face& f = mesh_.faces()[cFaces[i]];

        forAll(f, fp)
        {
            if (pointToObjVert.insert(f[fp], objVertI))
            {
                meshTools::writeOBJ(str, mesh_.points()[f[fp]]);
                objVertI++;
            }
        }
    }

    forAll(cFaces, i)
    {
        const face& f = mesh_.faces()[cFaces[i]];

        forAll(f, fp)
        {
            label pointI = f[fp];
            label nexPointI = f[f.fcIndex(fp)];

            str << "l " << pointToObjVert[pointI]+1
                << ' ' << pointToObjVert[nexPointI]+1 << nl;
        }
    }
}


// Find point with certain pointLevel. Skip any higher levels.
Foam::label Foam::hexRef4::findLevel
(
    const label faceI,
    const face& f,
    const label startFp,
    const bool searchForward,
    const label wantedLevel
) const
{
    label fp = startFp;

    forAll(f, i)
    {
        label pointI = f[fp];

        if (pointLevel_[pointI] < wantedLevel)
        {
            dumpCell(mesh_.faceOwner()[faceI]);
            if (mesh_.isInternalFace(faceI))
            {
                dumpCell(mesh_.faceNeighbour()[faceI]);
            }

            FatalErrorIn("hexRef4::findLevel(..)")
                << "face:" << f
                << " level:" << UIndirectList<label>(pointLevel_, f)()
                << " startFp:" << startFp
                << " wantedLevel:" << wantedLevel
                << abort(FatalError);
        }
        else if (pointLevel_[pointI] == wantedLevel)
        {
            return fp;
        }

        if (searchForward)
        {
            fp = f.fcIndex(fp);
        }
        else
        {
            fp = f.rcIndex(fp);
        }
    }

    dumpCell(mesh_.faceOwner()[faceI]);
    if (mesh_.isInternalFace(faceI))
    {
        dumpCell(mesh_.faceNeighbour()[faceI]);
    }

    FatalErrorIn("hexRef4::findLevel(..)")
        << "face:" << f
        << " level:" << UIndirectList<label>(pointLevel_, f)()
        << " startFp:" << startFp
        << " wantedLevel:" << wantedLevel
        << abort(FatalError);

    return -1;
}


// Gets cell level such that the face has four points <= level.
Foam::label Foam::hexRef4::getAnchorLevel(const label faceI) const
{
    const face& f = mesh_.faces()[faceI];

    if (f.size() <= 4)
    {
        return pointLevel_[f[findMaxLevel(f)]];
    }
    else
    {
        label ownLevel = cellLevel_[mesh_.faceOwner()[faceI]];

        if (countAnchors(f, ownLevel) == 4)
        {
            return ownLevel;
        }
        else if (countAnchors(f, ownLevel+1) == 4)
        {
            return ownLevel+1;
        }
        else
        {
            return -1;
        }
    }
}


void Foam::hexRef4::checkInternalOrientation
(
    polyTopoChange& meshMod,
    const label cellI,
    const label faceI,
    const point& ownPt,
    const point& neiPt,
    const face& newFace
)
{
    face compactFace(identity(newFace.size()));
    pointField compactPoints(meshMod.points(), newFace);

    vector n(compactFace.normal(compactPoints));

    vector dir(neiPt - ownPt);

    if ((dir & n) < 0)
    {
        FatalErrorIn("checkInternalOrientation(..)")
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace << endl
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " neiPt:" << neiPt
			<< " dir:" << dir
			<< " n:" << n
            << abort(FatalError);
    }

    vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    scalar s = (fcToOwn&n) / (dir&n);

    if (s < 0.1 || s > 0.9)
    {
        FatalErrorIn("checkInternalOrientation(..)")
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace << endl
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " neiPt:" << neiPt
            << " s:" << s
            << abort(FatalError);
    }
}


void Foam::hexRef4::checkBoundaryOrientation
(
    polyTopoChange& meshMod,
    const label cellI,
    const label faceI,
    const point& ownPt,
    const point& boundaryPt,
    const face& newFace
)
{
    face compactFace(identity(newFace.size()));
    pointField compactPoints(meshMod.points(), newFace);

    vector n(compactFace.normal(compactPoints));

    vector dir(boundaryPt - ownPt);

    if ((dir & n) < 0)
    {
        FatalErrorIn("checkBoundaryOrientation(..)")
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << abort(FatalError);
    }

    vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    scalar s = (fcToOwn&dir) / magSqr(dir);

    if (s < 0.7 || s > 1.3)
    {
        WarningIn("checkBoundaryOrientation(..)")
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << " s:" << s
            << endl;
            //<< abort(FatalError);
    }
}


// If p0 and p1 are existing vertices check if edge is split and insert
// splitPoint.
void Foam::hexRef4::insertEdgeSplit
(
    const labelList& edgeMidPoint,
    const label p0,
    const label p1,
    DynamicList<label>& verts
) const
{
    if (p0 < mesh_.nPoints() && p1 < mesh_.nPoints())
    {
        label edgeI = meshTools::findEdge(mesh_, p0, p1);

        if (edgeI != -1 && edgeMidPoint[edgeI] != -1)
        {
            verts.append(edgeMidPoint[edgeI]);
        }
    }
}


// Internal faces are one per edge between anchor points. So one per midPoint
// between the anchor points. Here we store the information on the midPoint
// and if we have enough information:
// - two anchors
// - two face mid points
// we add the face. Note that this routine can get called anywhere from
// two times (two unrefined faces) to four times (two refined faces) so
// the first call that adds the information creates the face.
Foam::label Foam::hexRef4::storeMidPointInfo
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& edgeMidPoint,
    const label cellI,
    const label faceI,
    const bool faceOrder,
    const label edgeMidPointI,
    const label anchorPointI,
    const label faceMidPointI,

    Map<edge>& midPointToAnchors,
    Map<edge>& midPointToFaceMids,
    polyTopoChange& meshMod
) const
{
    // See if need to store anchors.

    bool changed = false;
    bool haveTwoAnchors = false;

    Map<edge>::iterator edgeMidFnd = midPointToAnchors.find(edgeMidPointI);

    if (edgeMidFnd == midPointToAnchors.end())
    {
        midPointToAnchors.insert(edgeMidPointI, edge(anchorPointI, -1));
    }
    else
    {
        edge& e = edgeMidFnd();

        if (anchorPointI != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = anchorPointI;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoAnchors = true;
        }
    }

    bool haveTwoFaceMids = false;

    Map<edge>::iterator faceMidFnd = midPointToFaceMids.find(edgeMidPointI);

    if (faceMidFnd == midPointToFaceMids.end())
    {
        midPointToFaceMids.insert(edgeMidPointI, edge(faceMidPointI, -1));
    }
    else
    {
        edge& e = faceMidFnd();

        if (faceMidPointI != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = faceMidPointI;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoFaceMids = true;
        }
    }

    // Check if this call of storeMidPointInfo is the one that completed all
    // the necessary information.

    if (changed && haveTwoAnchors && haveTwoFaceMids)
    {
        const edge& anchors = midPointToAnchors[edgeMidPointI];
        const edge& faceMids = midPointToFaceMids[edgeMidPointI];

        label otherFaceMidPointI = faceMids.otherVertex(faceMidPointI);

        // Create face consistent with anchorI being the owner.
        // Note that the edges between the edge mid point and the face mids
        // might be marked for splitting. Note that these edge splits cannot
        // be between cellMid and face mids.


		// This section might still be needed for 2D?
        DynamicList<label> newFaceVerts(4);
        if (faceOrder == (mesh_.faceOwner()[faceI] == cellI))
        {
            newFaceVerts.append(faceMidPointI);

            // Check & insert edge split if any
            insertEdgeSplit
            (
                edgeMidPoint,
                faceMidPointI,  // edge between faceMid and
                edgeMidPointI,  // edgeMid
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointI,
                otherFaceMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(otherFaceMidPointI);
            newFaceVerts.append(cellMidPoint[cellI]);
        }
        else
        {
            newFaceVerts.append(otherFaceMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                otherFaceMidPointI,
                edgeMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointI,
                faceMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(faceMidPointI);
            newFaceVerts.append(cellMidPoint[cellI]);
        }

        face newFace;
        newFace.transfer(newFaceVerts);

        label anchorCell0 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            cellI,
            faceI,
            anchorPointI
        );
        label anchorCell1 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            cellI,
            faceI,
            anchors.otherVertex(anchorPointI)
        );


        label own, nei;
        point ownPt, neiPt;

        if (anchorCell0 < anchorCell1)
        {
            own = anchorCell0;
            nei = anchorCell1;

            ownPt = mesh_.points()[anchorPointI];
            neiPt = mesh_.points()[anchors.otherVertex(anchorPointI)];

        }
        else
        {
            own = anchorCell1;
            nei = anchorCell0;
            newFace.flip();

            ownPt = mesh_.points()[anchors.otherVertex(anchorPointI)];
            neiPt = mesh_.points()[anchorPointI];
        }

        if (debug)
        {
            point ownPt, neiPt;

            if (anchorCell0 < anchorCell1)
            {
                ownPt = mesh_.points()[anchorPointI];
                neiPt = mesh_.points()[anchors.otherVertex(anchorPointI)];
            }
            else
            {
                ownPt = mesh_.points()[anchors.otherVertex(anchorPointI)];
                neiPt = mesh_.points()[anchorPointI];
            }

            checkInternalOrientation
            (
                meshMod,
                cellI,
                faceI,
                ownPt,
                neiPt,
                newFace
            );
        }

        return addInternalFace
        (
            meshMod,
            faceI,
            anchorPointI,
            newFace,
            own,
            nei
        );
    }
    else
    {
        return -1;
    }
}

void Foam::hexRef4::walkFaceToMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label faceI,
    const label startFp,
    DynamicList<label>& faceVerts
) const
{
    const face& f = mesh_.faces()[faceI];
    const labelList& fEdges = mesh_.faceEdges(faceI);

    label fp = startFp;

    // Starting from fp store all (1 or 2) vertices until where the face
    // gets split

    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] >= 0)
        {
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (pointLevel_[f[fp]] <= cLevel)
        {
            // Next anchor. Have already append split point on edge in code
            // above.
            return;
        }
        else if (pointLevel_[f[fp]] == cLevel+1)
        {
            // Mid level
            faceVerts.append(f[fp]);
            return;
        }
        else if (pointLevel_[f[fp]] == cLevel+2)
		{
            // Store and continue to cLevel+1.
            faceVerts.append(f[fp]);
        }
    }
}


// Same as walkFaceToMid but now walk back.
void Foam::hexRef4::walkFaceFromMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label faceI,
    const label startFp,
    DynamicList<label>& faceVerts
) const
{
    const face& f = mesh_.faces()[faceI];
    const labelList& fEdges = mesh_.faceEdges(faceI);

    label fp = f.rcIndex(startFp);

    while (true)
    {
        if (pointLevel_[f[fp]] <= cLevel)
        {
            // anchor.
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel+1)
        {
            // Mid level
            faceVerts.append(f[fp]);
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel+2)
        {
            // Continue to cLevel+1.
        }
        fp = f.rcIndex(fp);
    }

    // Store
    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] >= 0)
        {
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (fp == startFp)
        {
            break;
        }
        faceVerts.append(f[fp]);
    }
}


// Updates refineCell (cells marked for refinement) so across all faces
// there will be 2:1 consistency after refinement.
Foam::label Foam::hexRef4::faceConsistentRefinement
(
    const bool maxSet,
    PackedBoolList& refineCell
) const
{
    label nChanged = 0;

    // Internal faces.
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        label nei = mesh_.faceNeighbour()[faceI];
        label neiLevel = cellLevel_[nei] + refineCell.get(nei);

        if (ownLevel > (neiLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(nei);
            }
            else
            {
                refineCell.unset(own);
            }
            nChanged++;
        }
        else if (neiLevel > (ownLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(own);
            }
            else
            {
                refineCell.unset(nei);
            }
            nChanged++;
        }
    }


    // Coupled faces. Swap owner level to get neighbouring cell level.
    // (only boundary faces of neiLevel used)
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

        neiLevel[i] = cellLevel_[own] + refineCell.get(own);
    }

    // Swap to neighbour
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);

    // Now we have neighbour value see which cells need refinement
    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        if (ownLevel > (neiLevel[i]+1))
        {
            if (!maxSet)
            {
                refineCell.unset(own);
                nChanged++;
            }
        }
        else if (neiLevel[i] > (ownLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(own);
                nChanged++;
            }
        }
    }

    return nChanged;
}


// Debug: check if wanted refinement is compatible with 2:1
void Foam::hexRef4::checkWantedRefinementLevels
(
    const labelList& cellsToRefine
) const
{
    PackedBoolList refineCell(mesh_.nCells());
    forAll(cellsToRefine, i)
    {
        refineCell.set(cellsToRefine[i]);
    }

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        label nei = mesh_.faceNeighbour()[faceI];
        label neiLevel = cellLevel_[nei] + refineCell.get(nei);

        if (mag(ownLevel-neiLevel) > 1)
        {
            dumpCell(own);
            dumpCell(nei);
            FatalErrorIn
            (
                "hexRef4::checkWantedRefinementLevels(const labelList&)"
            )   << "cell:" << own
                << " current level:" << cellLevel_[own]
                << " level after refinement:" << ownLevel
                << nl
                << "neighbour cell:" << nei
                << " current level:" << cellLevel_[nei]
                << " level after refinement:" << neiLevel
                << nl
                << "which does not satisfy 2:1 constraints anymore."
                << abort(FatalError);
        }
    }

    // Coupled faces. Swap owner level to get neighbouring cell level.
    // (only boundary faces of neiLevel used)
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

        neiLevel[i] = cellLevel_[own] + refineCell.get(own);
    }

    // Swap to neighbour
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);

    // Now we have neighbour value see which cells need refinement
    forAll(neiLevel, i)
    {
        label faceI = i + mesh_.nInternalFaces();

        label own = mesh_.faceOwner()[faceI];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        if (mag(ownLevel - neiLevel[i]) > 1)
        {
            label patchI = mesh_.boundaryMesh().whichPatch(faceI);

            dumpCell(own);
            FatalErrorIn
            (
                "hexRef4::checkWantedRefinementLevels(const labelList&)"
            )   << "Celllevel does not satisfy 2:1 constraint."
                << " On coupled face "
                << faceI
                << " on patch " << patchI << " "
                << mesh_.boundaryMesh()[patchI].name()
                << " owner cell " << own
                << " current level:" << cellLevel_[own]
                << " level after refinement:" << ownLevel
                << nl
                << " (coupled) neighbour cell will get refinement "
                << neiLevel[i]
                << abort(FatalError);
        }
    }
}


// Set instance for mesh files
void Foam::hexRef4::setInstance(const fileName& inst)
{
    if (debug)
    {
        Pout<< "hexRef4::setInstance(const fileName& inst) : "
            << "Resetting file instance to " << inst << endl;
    }

    cellLevel_.instance() = inst;
    pointLevel_.instance() = inst;
    level0Edge_.instance() = inst;
    history_.instance() = inst;
}


void Foam::hexRef4::collectLevelPoints
(
    const labelList& f,
    const label level,
    DynamicList<label>& points
) const
{
    forAll(f, fp)
    {
        if (pointLevel_[f[fp]] <= level)
        {
            points.append(f[fp]);
        }
    }
}


void Foam::hexRef4::collectLevelPoints
(
    const labelList& meshPoints,
    const labelList& f,
    const label level,
    DynamicList<label>& points
) const
{
    forAll(f, fp)
    {
        label pointI = meshPoints[f[fp]];
        if (pointLevel_[pointI] <= level)
        {
            points.append(pointI);
        }
    }
}


// Return true if we've found 6 quads. faces guaranteed to be outwards pointing.
bool Foam::hexRef4::matchHexShape
(
    const label cellI,
    const label cellLevel,
    DynamicList<face>& quads
) const
{
    const cell& cFaces = mesh_.cells()[cellI];

    // Work arrays
    DynamicList<label> verts(4);
    quads.clear();


    // 1. pick up any faces with four cellLevel points

    forAll(cFaces, i)
    {
        label faceI = cFaces[i];
        const face& f = mesh_.faces()[faceI];

        verts.clear();
        collectLevelPoints(f, cellLevel, verts);
        if (verts.size() == 4)
        {
            if (mesh_.faceOwner()[faceI] != cellI)
            {
                reverse(verts);
            }
            quads.append(face(0));
            labelList& quadVerts = quads.last();
            quadVerts.transfer(verts);
        }
    }


    if (quads.size() < 6)
    {
        Map<labelList> pointFaces(2*cFaces.size());

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];
            const face& f = mesh_.faces()[faceI];

            // Pick up any faces with only one level point.
            // See if there are four of these where the commont point
            // is a level+1 point. This common point is then the mid of
            // a split face.

            verts.clear();
            collectLevelPoints(f, cellLevel, verts);
            if (verts.size() == 1)
            {
                // Add to pointFaces for any level+1 point (this might be
                // a midpoint of a split face)
                forAll(f, fp)
                {
                    label pointI = f[fp];
                    if (pointLevel_[pointI] == cellLevel+1)
                    {
                        Map<labelList>::iterator iter =
                            pointFaces.find(pointI);
                        if (iter != pointFaces.end())
                        {
                            labelList& pFaces = iter();
                            if (findIndex(pFaces, faceI) == -1)
                            {
                                pFaces.append(faceI);
                            }
                        }
                        else
                        {
                            pointFaces.insert
                            (
                                pointI,
                                labelList(1, faceI)
                            );
                        }
                    }
                }
            }
        }

        // 2. Check if we've collected any midPoints.
        forAllConstIter(Map<labelList>, pointFaces, iter)
        {
            const labelList& pFaces = iter();

            if (pFaces.size() == 4)
            {
                // Collect and orient.
                faceList fourFaces(pFaces.size());
                forAll(pFaces, pFaceI)
                {
                    label faceI = pFaces[pFaceI];
                    const face& f = mesh_.faces()[faceI];
                    if (mesh_.faceOwner()[faceI] == cellI)
                    {
                        fourFaces[pFaceI] = f;
                    }
                    else
                    {
                        fourFaces[pFaceI] = f.reverseFace();
                    }
                }

                primitivePatch bigFace
                (
                    SubList<face>(fourFaces, fourFaces.size()),
                    mesh_.points()
                );
                const labelListList& edgeLoops = bigFace.edgeLoops();

                if (edgeLoops.size() == 1)
                {
                    // Collect the 4 cellLevel points
                    verts.clear();
                    collectLevelPoints
                    (
                        bigFace.meshPoints(),
                        bigFace.edgeLoops()[0],
                        cellLevel,
                        verts
                    );

                    if (verts.size() == 4)
                    {
                        quads.append(face(0));
                        labelList& quadVerts = quads.last();
                        quadVerts.transfer(verts);
                    }
                }
            }
        }
    }

    return (quads.size() == 6);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh, read refinement data
Foam::hexRef4::hexRef4(const polyMesh& mesh, const bool readHistory)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh_.nCells(), 0)
    ),
    oldCellLevel_(cellLevel_),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh_.nPoints(), 0)
    ),
    oldPointLevel_(pointLevel_),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar("level0Edge", dimLength, getLevel0EdgeLength())
    ),
    history_
    (
        IOobject
        (
            "refinementTree",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (readHistory ? mesh_.nCells() : 0)  // All cells visible if not be read
    ),
    cellDepth_(getCellDepth()),
    aspectX_to_Y_(getAspectRatio()),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if (readHistory)
    {
        // Make sure we don't use the master-only reading. Bit of a hack for
        // now.
        regIOobject::fileCheckTypes oldType =
            regIOobject::fileModificationChecking;
        regIOobject::fileModificationChecking = regIOobject::timeStamp;
        history_.readOpt() = IOobject::READ_IF_PRESENT;
        if (history_.headerOk())
        {
            history_.read();
        }
        regIOobject::fileModificationChecking = oldType;
    }

    if (history_.active() && history_.visibleCells().size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "hexRef4::hexRef4(const polyMesh&)"
        )   << "History enabled but number of visible cells in it "
            << history_.visibleCells().size()
            << " is not equal to the number of cells in the mesh "
            << mesh_.nCells()
            << abort(FatalError);
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorIn
        (
            "hexRef4::hexRef4(const polyMesh&)"
        )   << "Restarted from inconsistent cellLevel or pointLevel files."
            << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }


    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));


    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}


// Construct from components
Foam::hexRef4::hexRef4
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const refinementTree& history,
    const scalar level0Edge
)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        cellLevel
    ),
    oldCellLevel_(cellLevel_),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointLevel
    ),
    oldPointLevel_(pointLevel_),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar
        (
            "level0Edge",
            dimLength,
            (level0Edge >= 0 ? level0Edge : getLevel0EdgeLength())
        )
    ),
    history_
    (
        IOobject
        (
            "refinementTree",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        history
    ),
     cellDepth_(getCellDepth()),
    aspectX_to_Y_(getAspectRatio()),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if (history_.active() && history_.visibleCells().size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "hexRef4::hexRef4(const polyMesh&, const labelList&"
            ", const labelList&, const refinementTree&)"
        )   << "History enabled but number of visible cells in it "
            << history_.visibleCells().size()
            << " is not equal to the number of cells in the mesh "
            << mesh_.nCells() << abort(FatalError);
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorIn
        (
            "hexRef4::hexRef4(const polyMesh&, const labelList&"
            ", const labelList&, const refinementTree&)"
        )   << "Incorrect cellLevel or pointLevel size." << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }

    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));


    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}


// Construct from components
Foam::hexRef4::hexRef4
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const scalar level0Edge
)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        cellLevel
    ),
    oldCellLevel_(cellLevel_),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointLevel
    ),
    oldPointLevel_(pointLevel_),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar
        (
            "level0Edge",
            dimLength,
            (level0Edge >= 0 ? level0Edge : getLevel0EdgeLength())
        )
    ),
    history_
    (
        IOobject
        (
            "refinementTree",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        List<refinementTree::splitCell4>(0),
        labelList(0)
    ),
     cellDepth_(getCellDepth()),
    aspectX_to_Y_(getAspectRatio()),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorIn
        (
            "hexRef4::hexRef4(const polyMesh&, const labelList&"
            ", const labelList&)"
        )   << "Incorrect cellLevel or pointLevel size." << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }

    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));

    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::hexRef4::consistentRefinement
(
    const labelList& cellsToRefine,
    const bool maxSet
) const
{
    // Loop, modifying cellsToRefine, until no more changes to due to 2:1
    // conflicts.
    // maxSet = false : unselect cells to refine
    // maxSet = true  : select cells to refine

    // Go to straight boolList.
    PackedBoolList refineCell(mesh_.nCells());
    forAll(cellsToRefine, i)
    {
        refineCell.set(cellsToRefine[i]);
    }

    while (true)
    {
        label nChanged = faceConsistentRefinement(maxSet, refineCell);

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef4::consistentRefinement : Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }
    }


    // Convert back to labelList.
    label nRefined = 0;

    forAll(refineCell, cellI)
    {
        if (refineCell.get(cellI))
        {
            nRefined++;
        }
    }

    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(refineCell, cellI)
    {
        if (refineCell.get(cellI))
        {
            newCellsToRefine[nRefined++] = cellI;
        }
    }

    if (debug)
    {
        checkWantedRefinementLevels(newCellsToRefine);
    }

    return newCellsToRefine;
}


// Given a list of cells to refine determine additional cells to refine
// such that the overall refinement:
// - satisfies maxFaceDiff (e.g. 2:1) across neighbouring faces
// - satisfies maxPointDiff (e.g. 4:1) across selected point connected
//   cells. This is used to ensure that e.g. cells on the surface are not
//   point connected to cells which are 8 times smaller.
/*
Foam::labelList Foam::hexRef4::consistentSlowRefinement
(
    const label maxFaceDiff,
    const labelList& cellsToRefine,
    const labelList& facesToCheck,
    const label maxPointDiff,
    const labelList& pointsToCheck
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();


    if (maxFaceDiff <= 0)
    {
        FatalErrorIn
        (
            "hexRef4::consistentSlowRefinement"
            "(const label, const labelList&, const labelList&"
            ", const label, const labelList&)"
        )   << "Illegal maxFaceDiff " << maxFaceDiff << nl
            << "Value should be >= 1" << exit(FatalError);
    }


    // Bit tricky. Say we want a distance of three cells between two
    // consecutive refinement levels. This is done by using FaceCellWave to
    // transport out the new refinement level. It gets decremented by one
    // every cell it crosses so if we initialize it to maxFaceDiff
    // we will get a field everywhere that tells us whether an unselected cell
    // needs refining as well.


    // Initial information about (distance to) cellLevel on all cells
    List<refinementData> allCellInfo(mesh_.nCells());

    // Initial information about (distance to) cellLevel on all faces
    List<refinementData> allFaceInfo(mesh_.nFaces());

    forAll(allCellInfo, cellI)
    {
        // maxFaceDiff since refinementData counts both
        // faces and cells.
        allCellInfo[cellI] = refinementData
        (
            maxFaceDiff*(cellLevel_[cellI]+1),// when cell is to be refined
            maxFaceDiff*cellLevel_[cellI]     // current level
        );
    }

    // Cells to be refined will have cellLevel+1
    forAll(cellsToRefine, i)
    {
        label cellI = cellsToRefine[i];

        allCellInfo[cellI].count() = allCellInfo[cellI].refinementCount();
    }


    // Labels of seed faces
    DynamicList<label> seedFaces(mesh_.nFaces()/100);
    // refinementLevel data on seed faces
    DynamicList<refinementData> seedFacesInfo(mesh_.nFaces()/100);

    // Dummy additional info for FaceCellWave
    int dummyTrackData = 0;


    // Additional buffer layer thickness by changing initial count. Usually
    // this happens on boundary faces. Bit tricky. Use allFaceInfo to mark
    // off thus marked faces so they're skipped in the next loop.
    forAll(facesToCheck, i)
    {
        label faceI = facesToCheck[i];

        if (allFaceInfo[faceI].valid(dummyTrackData))
        {
            // Can only occur if face has already gone through loop below.
            FatalErrorIn
            (
                "hexRef4::consistentSlowRefinement"
                "(const label, const labelList&, const labelList&"
                ", const label, const labelList&)"
            )   << "Argument facesToCheck seems to have duplicate entries!"
                << endl
                << "face:" << faceI << " occurs at positions "
                << findIndices(facesToCheck, faceI)
                << abort(FatalError);
        }


        const refinementData& ownData = allCellInfo[faceOwner[faceI]];

        if (mesh_.isInternalFace(faceI))
        {
            // Seed face if neighbouring cell (after possible refinement)
            // will be refined one more than the current owner or neighbour.

            const refinementData& neiData = allCellInfo[faceNeighbour[faceI]];

            label faceCount;
            label faceRefineCount;
            if (neiData.count() > ownData.count())
            {
                faceCount = neiData.count() + maxFaceDiff;
                faceRefineCount = faceCount + maxFaceDiff;
            }
            else
            {
                faceCount = ownData.count() + maxFaceDiff;
                faceRefineCount = faceCount + maxFaceDiff;
            }

            seedFaces.append(faceI);
            seedFacesInfo.append
            (
                refinementData
                (
                    faceRefineCount,
                    faceCount
                )
            );
            allFaceInfo[faceI] = seedFacesInfo.last();
        }
        else
        {
            label faceCount = ownData.count() + maxFaceDiff;
            label faceRefineCount = faceCount + maxFaceDiff;

            seedFaces.append(faceI);
            seedFacesInfo.append
            (
                refinementData
                (
                    faceRefineCount,
                    faceCount
                )
            );
            allFaceInfo[faceI] = seedFacesInfo.last();
        }
    }


    // Just seed with all faces inbetween different refinement levels for now
    // (alternatively only seed faces on cellsToRefine but that gives problems
    //  if no cells to refine)
    forAll(faceNeighbour, faceI)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[faceI].valid(dummyTrackData))
        {
            label own = faceOwner[faceI];
            label nei = faceNeighbour[faceI];

            // Seed face with transported data from highest cell.

            if (allCellInfo[own].count() > allCellInfo[nei].count())
            {
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    own,
                    allCellInfo[own],
                    FaceCellWave<refinementData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFaces.append(faceI);
                seedFacesInfo.append(allFaceInfo[faceI]);
            }
            else if (allCellInfo[own].count() < allCellInfo[nei].count())
            {
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    nei,
                    allCellInfo[nei],
                    FaceCellWave<refinementData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFaces.append(faceI);
                seedFacesInfo.append(allFaceInfo[faceI]);
            }
        }
    }

    // Seed all boundary faces with owner value. This is to make sure that
    // they are visited (probably only important for coupled faces since
    // these need to be visited from both sides)
    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[faceI].valid(dummyTrackData))
        {
            label own = faceOwner[faceI];

            // Seed face with transported data from owner.
            refinementData faceData;
            faceData.updateFace
            (
                mesh_,
                faceI,
                own,
                allCellInfo[own],
                FaceCellWave<refinementData, int>::propagationTol(),
                dummyTrackData
            );
            seedFaces.append(faceI);
            seedFacesInfo.append(faceData);
        }
    }


    // face-cell-face transport engine
    FaceCellWave<refinementData, int> levelCalc
    (
        mesh_,
        allFaceInfo,
        allCellInfo,
        dummyTrackData
    );

    while (true)
    {
        if (debug)
        {
            Pout<< "hexRef4::consistentSlowRefinement : Seeded "
                << seedFaces.size() << " faces between cells with different"
                << " refinement level." << endl;
        }

        // Set seed faces
        levelCalc.setFaceInfo(seedFaces.shrink(), seedFacesInfo.shrink());
        seedFaces.clear();
        seedFacesInfo.clear();

        // Iterate until no change. Now 2:1 face difference should be satisfied
        levelCalc.iterate(mesh_.globalData().nTotalFaces()+1);


        // Now check point-connected cells (face-connected cells already ok):
        // - get per point max of connected cells
        // - sync across coupled points
        // - check cells against above point max

        if (maxPointDiff == -1)
        {
            // No need to do any point checking.
            break;
        }

        // Determine per point the max cell level. (done as count, not
        // as cell level purely for ease)
        labelList maxPointCount(mesh_.nPoints(), 0);

        forAll(maxPointCount, pointI)
        {
            label& pLevel = maxPointCount[pointI];

            const labelList& pCells = mesh_.pointCells(pointI);

            forAll(pCells, i)
            {
                pLevel = max(pLevel, allCellInfo[pCells[i]].count());
            }
        }

        // Sync maxPointCount to neighbour
        syncTools::syncPointList
        (
            mesh_,
            maxPointCount,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Update allFaceInfo from maxPointCount for all points to check
        // (usually on boundary faces)

        // Per face the new refinement data
        Map<refinementData> changedFacesInfo(pointsToCheck.size());

        forAll(pointsToCheck, i)
        {
            label pointI = pointsToCheck[i];

            // Loop over all cells using the point and check whether their
            // refinement level is much less than the maximum.

            const labelList& pCells = mesh_.pointCells(pointI);

            forAll(pCells, pCellI)
            {
                label cellI = pCells[pCellI];

                refinementData& cellInfo = allCellInfo[cellI];

                if
                (
                   !cellInfo.isRefined()
                 && (
                        maxPointCount[pointI]
                      > cellInfo.count() + maxFaceDiff*maxPointDiff
                    )
                )
                {
                    // Mark cell for refinement
                    cellInfo.count() = cellInfo.refinementCount();

                    // Insert faces of cell as seed faces.
                    const cell& cFaces = mesh_.cells()[cellI];

                    forAll(cFaces, cFaceI)
                    {
                        label faceI = cFaces[cFaceI];

                        refinementData faceData;
                        faceData.updateFace
                        (
                            mesh_,
                            faceI,
                            cellI,
                            cellInfo,
                            FaceCellWave<refinementData, int>::propagationTol(),
                            dummyTrackData
                        );

                        if (faceData.count() > allFaceInfo[faceI].count())
                        {
                            changedFacesInfo.insert(faceI, faceData);
                        }
                    }
                }
            }
        }

        label nChanged = changedFacesInfo.size();
        reduce(nChanged, sumOp<label>());

        if (nChanged == 0)
        {
            break;
        }


        // Transfer into seedFaces, seedFacesInfo
        seedFaces.setCapacity(changedFacesInfo.size());
        seedFacesInfo.setCapacity(changedFacesInfo.size());

        forAllConstIter(Map<refinementData>, changedFacesInfo, iter)
        {
            seedFaces.append(iter.key());
            seedFacesInfo.append(iter());
        }
    }


    if (debug)
    {
        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label own = mesh_.faceOwner()[faceI];
            label ownLevel =
                cellLevel_[own]
              + (allCellInfo[own].isRefined() ? 1 : 0);

            label nei = mesh_.faceNeighbour()[faceI];
            label neiLevel =
                cellLevel_[nei]
              + (allCellInfo[nei].isRefined() ? 1 : 0);

            if (mag(ownLevel-neiLevel) > 1)
            {
                dumpCell(own);
                dumpCell(nei);
                FatalErrorIn
                (
                    "hexRef4::consistentSlowRefinement"
                    "(const label, const labelList&, const labelList&"
                    ", const label, const labelList&)"
                )   << "cell:" << own
                    << " current level:" << cellLevel_[own]
                    << " current refData:" << allCellInfo[own]
                    << " level after refinement:" << ownLevel
                    << nl
                    << "neighbour cell:" << nei
                    << " current level:" << cellLevel_[nei]
                    << " current refData:" << allCellInfo[nei]
                    << " level after refinement:" << neiLevel
                    << nl
                    << "which does not satisfy 2:1 constraints anymore." << nl
                    << "face:" << faceI << " faceRefData:" << allFaceInfo[faceI]
                    << abort(FatalError);
            }
        }


        // Coupled faces. Swap owner level to get neighbouring cell level.
        // (only boundary faces of neiLevel used)

        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList neiCount(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList neiRefCount(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
            neiLevel[i] = cellLevel_[own];
            neiCount[i] = allCellInfo[own].count();
            neiRefCount[i] = allCellInfo[own].refinementCount();
        }

        // Swap to neighbour
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);
        syncTools::swapBoundaryFaceList(mesh_, neiCount);
        syncTools::swapBoundaryFaceList(mesh_, neiRefCount);

        // Now we have neighbour value see which cells need refinement
        forAll(neiLevel, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            label own = mesh_.faceOwner()[faceI];
            label ownLevel =
                cellLevel_[own]
              + (allCellInfo[own].isRefined() ? 1 : 0);

            label nbrLevel =
                neiLevel[i]
              + ((neiCount[i] >= neiRefCount[i]) ? 1 : 0);

            if (mag(ownLevel - nbrLevel) > 1)
            {
                dumpCell(own);
                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                FatalErrorIn
                (
                    "hexRef4::consistentSlowRefinement"
                    "(const label, const labelList&, const labelList&"
                    ", const label, const labelList&)"
                )   << "Celllevel does not satisfy 2:1 constraint."
                    << " On coupled face "
                    << faceI
                    << " refData:" << allFaceInfo[faceI]
                    << " on patch " << patchI << " "
                    << mesh_.boundaryMesh()[patchI].name() << nl
                    << "owner cell " << own
                    << " current level:" << cellLevel_[own]
                    << " current count:" << allCellInfo[own].count()
                    << " current refCount:"
                    << allCellInfo[own].refinementCount()
                    << " level after refinement:" << ownLevel
                    << nl
                    << "(coupled) neighbour cell"
                    << " has current level:" << neiLevel[i]
                    << " current count:" << neiCount[i]
                    << " current refCount:" << neiRefCount[i]
                    << " level after refinement:" << nbrLevel
                    << abort(FatalError);
            }
        }
    }

    // Convert back to labelList of cells to refine.

    label nRefined = 0;

    forAll(allCellInfo, cellI)
    {
        if (allCellInfo[cellI].isRefined())
        {
            nRefined++;
        }
    }

    // Updated list of cells to refine
    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(allCellInfo, cellI)
    {
        if (allCellInfo[cellI].isRefined())
        {
            newCellsToRefine[nRefined++] = cellI;
        }
    }

    if (debug)
    {
        Pout<< "hexRef4::consistentSlowRefinement : From "
            << cellsToRefine.size() << " to " << newCellsToRefine.size()
            << " cells to refine." << endl;
    }

    return newCellsToRefine;
}


Foam::labelList Foam::hexRef4::consistentSlowRefinement2
(
    const label maxFaceDiff,
    const labelList& cellsToRefine,
    const labelList& facesToCheck
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    if (maxFaceDiff <= 0)
    {
        FatalErrorIn
        (
            "hexRef4::consistentSlowRefinement2"
            "(const label, const labelList&, const labelList&)"
        )   << "Illegal maxFaceDiff " << maxFaceDiff << nl
            << "Value should be >= 1" << exit(FatalError);
    }

    const scalar level0Size = 2*maxFaceDiff*level0EdgeLength();


    // Bit tricky. Say we want a distance of three cells between two
    // consecutive refinement levels. This is done by using FaceCellWave to
    // transport out the 'refinement shell'. Anything inside the refinement
    // shell (given by a distance) gets marked for refinement.

    // Initial information about (distance to) cellLevel on all cells
    List<refinementDistanceData> allCellInfo(mesh_.nCells());

    // Initial information about (distance to) cellLevel on all faces
    List<refinementDistanceData> allFaceInfo(mesh_.nFaces());

    // Dummy additional info for FaceCellWave
    int dummyTrackData = 0;


    // Mark cells with wanted refinement level
    forAll(cellsToRefine, i)
    {
        label cellI = cellsToRefine[i];

        allCellInfo[cellI] = refinementDistanceData
        (
            level0Size,
            mesh_.cellCentres()[cellI],
            cellLevel_[cellI]+1             // wanted refinement
        );
    }
    // Mark all others with existing refinement level
    forAll(allCellInfo, cellI)
    {
        if (!allCellInfo[cellI].valid(dummyTrackData))
        {
            allCellInfo[cellI] = refinementDistanceData
            (
                level0Size,
                mesh_.cellCentres()[cellI],
                cellLevel_[cellI]           // wanted refinement
            );
        }
    }


    // Labels of seed faces
    DynamicList<label> seedFaces(mesh_.nFaces()/100);
    // refinementLevel data on seed faces
    DynamicList<refinementDistanceData> seedFacesInfo(mesh_.nFaces()/100);

    const pointField& cc = mesh_.cellCentres();

    forAll(facesToCheck, i)
    {
        label faceI = facesToCheck[i];

        if (allFaceInfo[faceI].valid(dummyTrackData))
        {
            // Can only occur if face has already gone through loop below.
            FatalErrorIn
            (
                "hexRef4::consistentSlowRefinement2"
                "(const label, const labelList&, const labelList&)"
            )   << "Argument facesToCheck seems to have duplicate entries!"
                << endl
                << "face:" << faceI << " occurs at positions "
                << findIndices(facesToCheck, faceI)
                << abort(FatalError);
        }

        label own = faceOwner[faceI];

        label ownLevel =
        (
            allCellInfo[own].valid(dummyTrackData)
          ? allCellInfo[own].originLevel()
          : cellLevel_[own]
        );

        if (!mesh_.isInternalFace(faceI))
        {
            // Do as if boundary face would have neighbour with one higher
            // refinement level.
            const point& fc = mesh_.faceCentres()[faceI];

            refinementDistanceData neiData
            (
                level0Size,
                2*fc - cc[own],    // est'd cell centre
                ownLevel+1
            );

            allFaceInfo[faceI].updateFace
            (
                mesh_,
                faceI,
                own,        // not used (should be nei)
                neiData,
                FaceCellWave<refinementDistanceData, int>::propagationTol(),
                dummyTrackData
            );
        }
        else
        {
            label nei = faceNeighbour[faceI];

            label neiLevel =
            (
                allCellInfo[nei].valid(dummyTrackData)
              ? allCellInfo[nei].originLevel()
              : cellLevel_[nei]
            );

            if (ownLevel == neiLevel)
            {
                // Fake as if nei>own or own>nei (whichever one 'wins')
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel+1),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel+1),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
            }
            else
            {
                // Difference in level anyway.
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
            }
        }
        seedFaces.append(faceI);
        seedFacesInfo.append(allFaceInfo[faceI]);
    }


    // Create some initial seeds to start walking from. This is only if there
    // are no facesToCheck.
    // Just seed with all faces inbetween different refinement levels for now
    forAll(faceNeighbour, faceI)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[faceI].valid(dummyTrackData))
        {
            label own = faceOwner[faceI];

            label ownLevel =
            (
                allCellInfo[own].valid(dummyTrackData)
              ? allCellInfo[own].originLevel()
              : cellLevel_[own]
            );

            label nei = faceNeighbour[faceI];

            label neiLevel =
            (
                allCellInfo[nei].valid(dummyTrackData)
              ? allCellInfo[nei].originLevel()
              : cellLevel_[nei]
            );

            if (ownLevel > neiLevel)
            {
                // Set face to owner data. (since face not yet would be copy)
                seedFaces.append(faceI);
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFacesInfo.append(allFaceInfo[faceI]);
            }
            else if (neiLevel > ownLevel)
            {
                seedFaces.append(faceI);
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFacesInfo.append(allFaceInfo[faceI]);
            }
        }
    }

    seedFaces.shrink();
    seedFacesInfo.shrink();

    // face-cell-face transport engine
    FaceCellWave<refinementDistanceData, int> levelCalc
    (
        mesh_,
        seedFaces,
        seedFacesInfo,
        allFaceInfo,
        allCellInfo,
        mesh_.globalData().nTotalCells()+1,
        dummyTrackData
    );


    //if (debug)
    //{
    //    // Dump wanted level
    //    volScalarField wantedLevel
    //    (
    //        IOobject
    //        (
    //            "wantedLevel",
    //            fMesh.time().timeName(),
    //            fMesh,
    //            IOobject::NO_READ,
    //            IOobject::AUTO_WRITE,
    //            false
    //        ),
    //        fMesh,
    //        dimensionedScalar("zero", dimless, 0)
    //    );
    //
    //    forAll(wantedLevel, cellI)
    //    {
    //        wantedLevel[cellI] = allCellInfo[cellI].wantedLevel(cc[cellI]);
    //    }
    //
    //    Pout<< "Writing " << wantedLevel.objectPath() << endl;
    //    wantedLevel.write();
    //}


    // Convert back to labelList of cells to refine.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // 1. Force original refinement cells to be picked up by setting the
    // originLevel of input cells to be a very large level (but within range
    // of 1<< shift inside refinementDistanceData::wantedLevel)
    forAll(cellsToRefine, i)
    {
        label cellI = cellsToRefine[i];

        //allCellInfo[cellI].originLevel() = sizeof(label)*8-2; ****PROBLEM?****
        allCellInfo[cellI].originLevel() = sizeof(label)*4-2;
        allCellInfo[cellI].origin() = cc[cellI];
    }

    // 2. Extend to 2:1. I don't understand yet why this is not done
    // 2. Extend to 2:1. For non-cube cells the scalar distance does not work
    // so make sure it at least provides 2:1.
    PackedBoolList refineCell(mesh_.nCells());
    forAll(allCellInfo, cellI)
    {
        label wanted = allCellInfo[cellI].wantedLevel(cc[cellI]);

        if (wanted > cellLevel_[cellI]+1)
        {
            refineCell.set(cellI);
        }
    }
    faceConsistentRefinement(true, refineCell);

    while (true)
    {
        label nChanged = faceConsistentRefinement(true, refineCell);

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef4::consistentSlowRefinement2 : Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }
    }

    // 3. Convert back to labelList.
    label nRefined = 0;

    forAll(refineCell, cellI)
    {
        if (refineCell.get(cellI))
        {
            nRefined++;
        }
    }

    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(refineCell, cellI)
    {
        if (refineCell.get(cellI))
        {
            newCellsToRefine[nRefined++] = cellI;
        }
    }

    if (debug)
    {
        Pout<< "hexRef4::consistentSlowRefinement2 : From "
            << cellsToRefine.size() << " to " << newCellsToRefine.size()
            << " cells to refine." << endl;

        // Check that newCellsToRefine obeys at least 2:1.

        {
            cellSet cellsIn(mesh_, "cellsToRefineIn", cellsToRefine);
            Pout<< "hexRef4::consistentSlowRefinement2 : writing "
                << cellsIn.size() << " to cellSet "
                << cellsIn.objectPath() << endl;
            cellsIn.write();
        }
        {
            cellSet cellsOut(mesh_, "cellsToRefineOut", newCellsToRefine);
            Pout<< "hexRef4::consistentSlowRefinement2 : writing "
                << cellsOut.size() << " to cellSet "
                << cellsOut.objectPath() << endl;
            cellsOut.write();
        }

        // Extend to 2:1
        PackedBoolList refineCell(mesh_.nCells());
        forAll(newCellsToRefine, i)
        {
            refineCell.set(newCellsToRefine[i]);
        }
        const PackedBoolList savedRefineCell(refineCell);

        label nChanged = faceConsistentRefinement(true, refineCell);

        {
            cellSet cellsOut2
            (
                mesh_, "cellsToRefineOut2", newCellsToRefine.size()
            );
            forAll(refineCell, cellI)
            {
                if (refineCell.get(cellI))
                {
                    cellsOut2.insert(cellI);
                }
            }
            Pout<< "hexRef4::consistentSlowRefinement2 : writing "
                << cellsOut2.size() << " to cellSet "
                << cellsOut2.objectPath() << endl;
            cellsOut2.write();
        }

        if (nChanged > 0)
        {
            forAll(refineCell, cellI)
            {
                if (refineCell.get(cellI) && !savedRefineCell.get(cellI))
                {
                    dumpCell(cellI);
                    FatalErrorIn
                    (
                        "hexRef4::consistentSlowRefinement2"
                        "(const label, const labelList&, const labelList&)"
                    )   << "Cell:" << cellI << " cc:"
                        << mesh_.cellCentres()[cellI]
                        << " was not marked for refinement but does not obey"
                        << " 2:1 constraints."
                        << abort(FatalError);
                }
            }
        }
    }

    return newCellsToRefine;
}
*/

// Top level driver to insert topo changes to do all refinement.
Foam::labelListList Foam::hexRef4::setRefinement
(
    const labelList& cellsToRefine,
    polyTopoChange& meshMod,
    const vector& normalDir
)
{
	bool title_debug = debug;
	if (debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Checking initial mesh just to make sure" << endl;

        checkMesh();
        // Cannot call checkRefinementlevels since hanging points might
        // get triggered by the mesher after subsetting.
        //checkRefinementLevels(-1, labelList(0));
    }

    // Clear any saved point/cell data.
    savedPointLevel_.clear();
    savedCellLevel_.clear();

    // New point/cell level. Copy of pointLevel for existing points.
    DynamicList<label> newCellLevel(cellLevel_.size());
    forAll(cellLevel_, cellI)
    {
        newCellLevel.append(cellLevel_[cellI]);
    }
    DynamicList<label> newPointLevel(pointLevel_.size());
    forAll(pointLevel_, pointI)
    {
        newPointLevel.append(pointLevel_[pointI]);
    }
    
    if (title_debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Allocating " << cellsToRefine.size() << " cell midpoints."
            << endl;
    }

    // Mid point per refined cell.
    // -1 : not refined
    // >=0: label of mid point.
    labelList cellMidPoint(mesh_.nCells(), -1);

    forAll(cellsToRefine, i)
    {
        label cellI = cellsToRefine[i];
        cellMidPoint[cellI] = 12345;
    }

    if (debug)
    {
        cellSet splitCells(mesh_, "splitCells", cellsToRefine.size());

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
                splitCells.insert(cellI);
            }
        }

        Pout<< "hexRef4::setRefinement : Dumping " << splitCells.size()
            << " cells to split to cellSet " << splitCells.objectPath()
            << endl;

        splitCells.write();
    }

    // Split edges
    // ~~~~~~~~~~~
    
    if (title_debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Allocating edge midpoints."
            << endl;
    }

    // Unrefined edges are ones between cellLevel or lower points.
    // If any cell using this edge gets split then the edge needs to be split.

    // -1  : no need to split edge
    // >=0 : label of introduced mid point
    labelList edgeMidPoint(mesh_.nEdges(), -1);

    forAll(cellMidPoint, cellI)
    {
        if (cellMidPoint[cellI] >= 0)
        {
            const labelList& cEdges = mesh_.cellEdges(cellI);

            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];
				
				// Test to remove edges parallel to refinement from edge-splitting
				if (mag(normalDir & meshTools::normEdgeVec(mesh_, edgeI)) > 0.1)
				{
					continue; // go to the next edge.
				}
				else
				{
					const edge& e = mesh_.edges()[edgeI];

					if
					(
						pointLevel_[e[0]] <= cellLevel_[cellI]
					 && pointLevel_[e[1]] <= cellLevel_[cellI]
					)
					{
						edgeMidPoint[edgeI] = 12345;    // mark need for splitting
					}
			    }
            }
        }
    }

    // Synchronize edgeMidPoint across coupled patches. Take max so that
    // any split takes precedence.
    syncTools::syncEdgeList
    (
        mesh_,
        edgeMidPoint,
        maxEqOp<label>(),
        labelMin
    );


    // Introduce edge points
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        // Phase 1: calculate midpoints and sync.
        // This needs doing for if people do not write binary and we slowly
        // get differences.

        pointField edgeMids(mesh_.nEdges(), point(-GREAT, -GREAT, -GREAT));

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                // Edge marked to be split.
                edgeMids[edgeI] = mesh_.edges()[edgeI].centre(mesh_.points());
                if (0)
                {
					Pout<< "edgeI = " << edgeI << " was split in the middle, at "
					<< edgeMids[edgeI] << endl;
				}
            }
        }
        syncTools::syncEdgePositions
        (
            mesh_,
            edgeMids,
            maxEqOp<vector>(),
            point(-GREAT, -GREAT, -GREAT)
        );

        // Phase 2: introduce points at the synced locations.
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                // Edge marked to be split. Replace edgeMidPoint with actual
                // point label.

                const edge& e = mesh_.edges()[edgeI];

                edgeMidPoint[edgeI] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        edgeMids[edgeI],            // point
                        e[0],                       // master point
                        -1,                         // zone for point
                        true                        // supports a cell
                    )
                );

                newPointLevel(edgeMidPoint[edgeI]) =
                    max(pointLevel_[e[0]], pointLevel_[e[1]]) + 1;
            }
        }
    }

    if (debug)
    {
        OFstream str(mesh_.time().path()/"edgeMidPoint.obj");

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                const edge& e = mesh_.edges()[edgeI];

                meshTools::writeOBJ(str, e.centre(mesh_.points()));
            }
        }

        Pout<< "hexRef4::setRefinement :"
            << " Dumping edge centres to split to file " << str.name() << endl;
    }


    // Calculate face level
    // ~~~~~~~~~~~~~~~~~~~~
    // (after splitting)

    if (title_debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Allocating face midpoints."
            << endl;
    }

    // Face anchor level. There are guaranteed 4 points with level
    // <= anchorLevel. These are the corner points.
    labelList faceAnchorLevel(mesh_.nFaces());

    for (label faceI = 0; faceI < mesh_.nFaces(); faceI++)
    {
        faceAnchorLevel[faceI] = getAnchorLevel(faceI);
    }

    // -1  : no need to split face
    // >=0 : label of introduced mid point
    labelList faceMidPoint(mesh_.nFaces(), -1);
    PackedBoolList frontAndBackFaces(mesh_.nFaces());

	{ // scope faceAreas, freeing the memory after frontAndBackFaces are set.
		Foam::vectorField faceAreas(mesh_.nFaces());
		calcFaceNormalVector(mesh_.points(), faceAreas);
		
		forAll(faceAreas, faceI)
		{
			Foam::vector faceNormal = (faceAreas[faceI] / mag(faceAreas[faceI]));
			if (mag(faceNormal & normalDir) > 0.5)
			{
				// Face is perpendicular to the direction not to be refined.
				// It is a member of the frontAndBack patch.
				frontAndBackFaces.set(faceI);
			}
		}
	}
	
    // Internal faces: look at cells on both sides. Uniquely determined since
    // face itself guaranteed to be same level as most refined neighbour.
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
		if (faceAnchorLevel[faceI] >= 0)
		{
			label own = mesh_.faceOwner()[faceI];	// the label of the cell which owns faceI
			label ownLevel = cellLevel_[own];		// the level of the cell which owns faceI
			label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0); // if the cell is being refined, newOwnLevel = ownLevel + 1

			label nei = mesh_.faceNeighbour()[faceI];
			label neiLevel = cellLevel_[nei];
			label newNeiLevel = neiLevel + (cellMidPoint[nei] >= 0 ? 1 : 0); // same with the neighbouring cell

			if
			(
				(newOwnLevel > faceAnchorLevel[faceI]
			 || newNeiLevel > faceAnchorLevel[faceI])
			 && frontAndBackFaces.get(faceI) // Take only faces perpendicular to the normalDirection
			)
			{
				faceMidPoint[faceI] = 12345;    // mark to be split
			}
		}
    }

    // Coupled patches handled like internal faces except now all information
    // from neighbour comes from across processor.
    // Boundary faces are more complicated since the boundary face can
    // be more refined than its owner (or neighbour for coupled patches)
    // (does not happen if refining/unrefining only, but does e.g. when
    //  refinining and subsetting)

    {
        labelList newNeiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(newNeiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
            label ownLevel = cellLevel_[own];
            label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

            newNeiLevel[i] = newOwnLevel;
        }

        // Swap.
        syncTools::swapBoundaryFaceList(mesh_, newNeiLevel);

        // So now we have information on the neighbour.

        forAll(newNeiLevel, i)
        {
            label faceI = i+mesh_.nInternalFaces();

			if (faceAnchorLevel[faceI] >= 0)
			{
				label own = mesh_.faceOwner()[faceI];
				label ownLevel = cellLevel_[own];
				label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

				if
				(
					(newOwnLevel > faceAnchorLevel[faceI]
				|| newNeiLevel[i] > faceAnchorLevel[faceI])
				&& frontAndBackFaces.get(faceI) // Take only faces perpendicular to the normalDirection
				)
				{
					faceMidPoint[faceI] = 12345;    // mark to be split
				}
			}
        }
    }

    // Synchronize faceMidPoint across coupled patches. (logical or)
    syncTools::syncFaceList
    (
        mesh_,
        faceMidPoint,
        maxEqOp<label>()
    );



    // Introduce face points
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        // Phase 1: determine mid points and sync. See comment for edgeMids
        // above
        pointField bFaceMids
        (
            mesh_.nFaces()-mesh_.nInternalFaces(),
            point(-GREAT, -GREAT, -GREAT)
        );

        forAll(bFaceMids, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            if (faceMidPoint[faceI] >= 0)
            {
                bFaceMids[i] = mesh_.faceCentres()[faceI];
            }
        }
        syncTools::syncBoundaryFacePositions
        (
            mesh_,
            bFaceMids,
            maxEqOp<vector>()
        );

        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] >= 0)
            {
                // Face marked to be split. Replace faceMidPoint with actual
                // point label.

                const face& f = mesh_.faces()[faceI];

                faceMidPoint[faceI] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        (
                            faceI < mesh_.nInternalFaces()
                          ? mesh_.faceCentres()[faceI]
                          : bFaceMids[faceI-mesh_.nInternalFaces()]
                        ),                          // point
                        f[0],                       // master point
                        -1,                         // zone for point
                        true                        // supports a cell
                    )
                );

                // Determine the level of the corner points and midpoint will
                // be one higher.
                newPointLevel(faceMidPoint[faceI]) = faceAnchorLevel[faceI]+1;
            }
        }
    }

    if (debug)
    {
        faceSet splitFaces(mesh_, "splitFaces", cellsToRefine.size());

        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] >= 0)
            {
                splitFaces.insert(faceI);
            }
        }

        Pout<< "hexRef4::setRefinement : Dumping " << splitFaces.size()
            << " faces to split to faceSet " << splitFaces.objectPath() << endl;

        splitFaces.write();
    }


    // Information complete
    // ~~~~~~~~~~~~~~~~~~~~
    // At this point we have all the information we need. We should no
    // longer reference the cellsToRefine to refine. All the information is:
    // - cellMidPoint >= 0 : cell needs to be split
    // - faceMidPoint >= 0 : face needs to be split
    // - edgeMidPoint >= 0 : edge needs to be split

    // Get the corner/anchor points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (title_debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Finding cell anchorPoints (8 per cell)"
            << endl;
    }

    // There will always be 8 points on the hex that have were introduced
    // with the hex and will have the same or lower refinement level.

    // Per cell the 8 corner points.
    labelListList cellAnchorPoints(mesh_.nCells());

    {
        labelList nAnchorPoints(mesh_.nCells(), 0);

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
                cellAnchorPoints[cellI].setSize(8);
            }
        }

        forAll(pointLevel_, pointI)
        {
            const labelList& pCells = mesh_.pointCells(pointI);

            forAll(pCells, pCellI)
            {
                label cellI = pCells[pCellI];

                if
                (
                    cellMidPoint[cellI] >= 0
                 && pointLevel_[pointI] <= cellLevel_[cellI]
                )
                {
                    if (nAnchorPoints[cellI] == 8)
                    {
                        dumpCell(cellI);
                        FatalErrorIn
                        (
                            "hexRef4::setRefinement(const labelList&"
                            ", polyTopoChange&)"
                        )   << "cell " << cellI
                            << " of level " << cellLevel_[cellI]
                            << " uses more than 8 points of equal or"
                            << " lower level" << nl
                            << "Points so far:" << cellAnchorPoints[cellI]
                            << abort(FatalError);
                    }

                    cellAnchorPoints[cellI][nAnchorPoints[cellI]++] = pointI;
                }
            }
        }

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
                if (nAnchorPoints[cellI] != 8)
                {
                    dumpCell(cellI);

                    const labelList& cPoints = mesh_.cellPoints(cellI);

                    FatalErrorIn
                    (
                        "hexRef4::setRefinement(const labelList&"
                        ", polyTopoChange&)"
                    )   << "cell " << cellI
                        << " of level " << cellLevel_[cellI]
                        << " does not seem to have 8 points of equal or"
                        << " lower level" << endl
                        << "cellPoints:" << cPoints << endl
                        << "pointLevels:"
                        << UIndirectList<label>(pointLevel_, cPoints)() << endl
                        << abort(FatalError);
                }
            }
        }
    }


    // Add the cells
    // ~~~~~~~~~~~~~

    if (title_debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Adding cells (1/2 per anchorPoint)"
            << endl;
    }
    // Per cell the 7 added cells (+ original cell)
    // Becomes 3 added cells for 2D
    labelListList cellAddedCells(mesh_.nCells());

    forAll(cellAnchorPoints, cellI)
    {
        const labelList& cAnchors = cellAnchorPoints[cellI];
		
        if (cAnchors.size() == 8)
        {
            labelList& cAdded = cellAddedCells[cellI];
            // cAdded is a reference to cellAddedCells[cellI]. i.e. replace everywhere and no change
            cAdded.setSize(4);

            // Original cell at 0
            cAdded[0] = cellI;
	
            // Update cell level
            newCellLevel[cellI] = cellLevel_[cellI]+1;

            for (label i = 1; i < 4; i++) // Changed from 8 - only adding 3 new cells
            {
                cAdded[i] = meshMod.setAction
                (
                    polyAddCell								// cellI is the only non-negative value, so the (new) cell's master is set by cellID.
                    (
                        -1,                                 // master point
                        -1,                                 // master edge
                        -1,                                 // master face
                        cellI,                              // master cell
                        mesh_.cellZones().whichZone(cellI)  // zone for cell
                    )
                );
                newCellLevel(cAdded[i]) = cellLevel_[cellI]+1;
            }
            
            if (debug)
            {
				Pout<< "cellAddedCells[" << cellI
					<< "] = " << cellAddedCells[cellI] << endl;
			}
        }
    }
    
    // This was previously much lower, but it seems like only the labels
    // for the created cells need to exist - information about location etc
    // is not stored by history_.
    
    
    pointLevel_.transfer(newPointLevel);
    cellLevel_.transfer(newCellLevel);
    
    
    // Update the live split cells tree.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // New unrefinement structure
    if (history_.active())
    {
        if (title_debug) // should be if(debug)
        {
            Pout<< "hexRef4::setRefinement :"
                << " Updating refinement history to " << cellLevel_.size()
                << " cells" << endl;
        }

        // Extend refinement history for new cells
        history_.resize(cellLevel_.size());

        forAll(cellAddedCells, cellI)
        {
            const labelList& addedCells = cellAddedCells[cellI];

            if (addedCells.size())
            {
                // Cell was split.
                history_.storeSplit(cellI, addedCells);
            }
        }
    }
    
    // Faces
    // ~~~~~
    // 1. existing faces that get split (into four always)
    // 2. existing faces that do not get split but only edges get split
    // 3. existing faces that do not get split but get new owner/neighbour
    // 4. new internal faces inside split cells.

    if (title_debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Marking faces to be handled"
            << endl;
    }

    // Get all affected faces.
    PackedBoolList affectedFace(mesh_.nFaces());

    {
        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
				if (0)
				{
					Pout<< nl << "cellI = " << cellI << endl;
				}
                const cell& cFaces = mesh_.cells()[cellI];

                forAll(cFaces, i) // marks every face of cellI for refinement.
                {
                    affectedFace.set(cFaces[i]);
                    if (0)
                    {
						Pout<< "cFaces[" << i << "] = " << cFaces[i] << endl;
					}
                }
            }
        }

        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] >= 0)
            {
                affectedFace.set(faceI);
                if (0)
                {
					Pout<< "faceI = " << faceI << ", set for (4-split) refinement" << endl;
				}
            }
        }

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
				if (0)
				{
					Pout << nl << "edgeI = " << edgeI << endl;
				}
                const labelList& eFaces = mesh_.edgeFaces(edgeI);

                forAll(eFaces, i)
                {
                    affectedFace.set(eFaces[i]);
                    if (0)
                    {
						Pout<< "eFaces[i] = " << eFaces[i] << ", set for refinement" << endl;
					}
                }
            }
        }
    }


    // 1. Faces that get split
    // ~~~~~~~~~~~~~~~~~~~~~~~

    if (title_debug)
    {
        Pout<< "hexRef4::setRefinement : Splitting faces" << endl;
    }

    forAll(faceMidPoint, faceI)
	{
        if (faceMidPoint[faceI] >= 0 && affectedFace.get(faceI))
        {
			if (!frontAndBackFaces.get(faceI))
			{
				FatalErrorIn("hexRef4::setRefinement(..)")
					<< "A non-front/back face is being split 4 ways. faceI = " 
					<< faceI
					<< abort(FatalError);
			}
            // Face needs to be split and hasn't yet been done in some way
            // (affectedFace - is impossible since this is first change but
            //  just for completeness)

            const face& f = mesh_.faces()[faceI];

            // Has original faceI been used (three faces added, original gets
            // modified)
            bool modifiedFace = false;

            label anchorLevel = faceAnchorLevel[faceI];

            face newFace(4); // Creates four 'face's? with value -1?

            forAll(f, fp)
            // for all fp in mesh_.faces()[faceI]
            {
                label pointI = f[fp];

                //~ if (pointLevel_[pointI] <= anchorLevel)
                if (pointLevel_[pointI] <= anchorLevel)
                {
                    // point is anchor. Start collecting face.
					
                    DynamicList<label> faceVerts(4);

                    faceVerts.append(pointI);

                    // Walk forward to mid point.
                    // - if next is +2 midpoint is +1
                    // - if next is +1 it is midpoint
                    // - if next is +0 there has to be edgeMidPoint

                    walkFaceToMid
                    (
                        edgeMidPoint,
                        anchorLevel,
                        faceI,
                        fp,
                        faceVerts
                    );

                    faceVerts.append(faceMidPoint[faceI]);

                    walkFaceFromMid
                    (
                        edgeMidPoint,
                        anchorLevel,
                        faceI,
                        fp,
                        faceVerts
                    );
                    
                    if (debug)
                    {
						Pout<< nl << "Info from part 1. faceI = " << faceI << ", faceVerts = "
						<< faceVerts << endl;
						Pout<< "Point locations are: ";
						for (int a=0; a<faceVerts.size(); a++)
						{
							Pout<< meshMod.points()[faceVerts[a]] << ", ";
						}
						Pout<< endl;
					}
					
					point newFaceCentre(0,0,0);
					for (int a = 0; a < faceVerts.size(); a++)
					{
						//~ if (pointLevel_...
						newFaceCentre += meshMod.points()[faceVerts[a]];
					}
					newFaceCentre /= 4;

                    newFace.transfer(faceVerts);

                    // Get new owner/neighbour
                    label own, nei;
                    getFaceNeighbours
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        faceI,
                        pointI,          // Anchor point

                        own,
                        nei
                    );

                    if (cellLevel_[own] > 0)
                    {
						vector newCentreFromOldCentre(newFaceCentre - mesh_.faceCentres()[faceI]);
						if (newCentreFromOldCentre[0] < 0 && newCentreFromOldCentre[1] < 0)
						{
							own = cellAddedCells[history_.myParentCell(own)][0];
						}
						else if (newCentreFromOldCentre[0] < 0 && newCentreFromOldCentre[1] > 0)
						{
							own = cellAddedCells[history_.myParentCell(own)][2];
						}
						else if (newCentreFromOldCentre[0] > 0 && newCentreFromOldCentre[1] < 0)
						{
							own = cellAddedCells[history_.myParentCell(own)][1];
						}
						else if (newCentreFromOldCentre[0] > 0 && newCentreFromOldCentre[1] > 0)
						{
							own = cellAddedCells[history_.myParentCell(own)][3];
						}
						else
						{
							// How would you possibly get here?
							FatalErrorIn("hexRef4::setRefinement(..)")
								<< "vector newCentreFromOldCentre had components that were in some way magical "
								<< "(probably equal to zero):"
								<< newCentreFromOldCentre
								<< abort(FatalError);
						}
					}

                    if (debug)
                    {
						if (mesh_.isInternalFace(faceI))
                        {
                            label oldOwn = mesh_.faceOwner()[faceI];
                            label oldNei = mesh_.faceNeighbour()[faceI];

                            checkInternalOrientation
                            (
                                meshMod,
                                oldOwn,
                                faceI,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.cellCentres()[oldNei],
                                newFace
                            );
                        }
                        else
                        {
                            label oldOwn = mesh_.faceOwner()[faceI];

                            checkBoundaryOrientation
                            (
                                meshMod,
                                oldOwn,
                                faceI,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.faceCentres()[faceI],
                                newFace
                            );
                        }
                    }

                    if (!modifiedFace)
                    {
                        modifiedFace = true;
						if (debug)
						{
							Pout<< "1. Modding face with own, nei : " << own << nei << endl;
						}
                        modFace(meshMod, faceI, newFace, own, nei);
                    }
                    else
                    {
						if (debug)
						{
							Pout<< "2. Adding face with own, nei : " << own << nei << endl;
						}
                        addFace(meshMod, faceI, newFace, own, nei);
                    }
                }
            }

            // Mark face as having been handled
            affectedFace.unset(faceI);
            if (0)
            {
				Pout<< "Face number " << faceI << " has been handled by "
					<< "the first part of the face-splitting section." << endl;
			}
        }
    }



    // 2. faces that do not get split but use edges that get split
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (title_debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Adding edge splits to unsplit faces"
            << endl;
    }
    

    DynamicList<label> eFacesStorage;
    DynamicList<label> fEdgesStorage;

    forAll(edgeMidPoint, edgeI)
    {
        if (edgeMidPoint[edgeI] >= 0)
        {
            // Split edge. Check that face not already handled above.

            const labelList& eFaces = mesh_.edgeFaces(edgeI, eFacesStorage);

            forAll(eFaces, i)
            {
                label faceI = eFaces[i];

                if (frontAndBackFaces.get(faceI) && affectedFace.get(faceI)) // Add points to unrefined bordering refined (front and back faces)
				{
					if (faceMidPoint[faceI] >= 0) // Fatal error
					{
						FatalErrorIn("hexRef4::setRefinement(..)")
			                << "The section of code designed to handle unrefined front "
			                << "and back faces which border refined areas was triggered "
			                << "for a face marked for refinement. "
			                << "faceI : " << faceI
			                << abort(FatalError);
					}
					if (debug)
					{
						Pout<< "faceI = " << faceI << " is a front or back face." << endl;
					}
					// need to modFace with new (5th) point.
					const face& f = mesh_.faces()[faceI];
					const labelList& fEdges = mesh_.faceEdges
					(
						faceI,
						fEdgesStorage
					);
					
					DynamicList<label> newFaceVerts(f.size());
					
					label fp = 0;
					newFaceVerts.clear();
					
					while (true)
					{
						newFaceVerts.append(f[fp]);
						label edgeI = fEdges[fp];
						if (edgeMidPoint[edgeI] >= 0)
						{
							newFaceVerts.append(edgeMidPoint[edgeI]);
						}
						fp = f.fcIndex(fp);
						if (fp == 0)
						{
							break;
						}
					}
					//~ 
					//~ forAll(f, fp)
                    //~ {
                        //~ newFaceVerts.append(f[fp]);
//~ 
                        //~ label edgeI = fEdges[fp];
//~ 
                        //~ if (edgeMidPoint[edgeI] >= 0)
                        //~ {
                            //~ newFaceVerts.append(edgeMidPoint[edgeI]);
                        //~ }
                    //~ }
					
					face newFace;
					if (debug) Pout<< "newFaceVerts = " << newFaceVerts << endl;
					newFace.transfer(newFaceVerts);

					// The point with the lowest level should be an anchor
                    // point of the neighbouring cells.
                    label anchorFp = findMinLevel(f);

                    label own, nei;
                    getFaceNeighbours
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        faceI,
                        f[anchorFp],          // Anchor point

                        own,
                        nei
                    );

                    if (debug)
                    {
					label oldOwn = meshMod.faceOwner()[faceI];
					label oldNei = meshMod.faceNeighbour()[faceI];	
                        if (mesh_.isInternalFace(faceI))
                        {
                            checkInternalOrientation
                            (
                                meshMod,
                                oldOwn,
                                faceI,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.cellCentres()[oldNei],
                                newFace
                            );
                        }
                        else
                        {
                            checkBoundaryOrientation
                            (
                                meshMod,
                                oldOwn,
                                faceI,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.faceCentres()[faceI],
                                newFace
                            );
                        }
                    }

                    meshMod.setAction
		            (
		                polyModifyFace
		                (
		                    newFace,            // modified face
		                    faceI,              // label of face being modified
		                    own,                // owner
		                    nei,                // neighbour
		                    false,              // face flip
		                    mesh_.boundaryMesh().whichPatch(faceI),            // patch for face
		                    false,              // remove from zone
		                    -1,             // zone for face
		                    false            // face flip in zone
		                )
		            );
                    affectedFace.unset(faceI);
				}
				
				if (faceMidPoint[faceI] < 0 && affectedFace.get(faceI)) // Split side faces
                {
                    // Unsplit face. Add edge splits to face.
                    // Side faces fulfill this criteria.
                    if (frontAndBackFaces.get(faceI))
                    {
						FatalErrorIn("hexRef4::setRefinement(..)")
			                << "The section of code designed to split side faces was triggered for a front or back face"
			                << nl
			                << "faceI : " << faceI
			                << abort(FatalError);
					}

                    const face& f = mesh_.faces()[faceI];

					label oldOwn, oldNei;
					getFaceNeighbours
					(
						cellAnchorPoints,
						cellAddedCells,
						faceI,
						f[findMinLevel(f)],
						oldOwn,
						oldNei
					);

					mySplitSideFaces(cellAddedCells, faceMidPoint, edgeMidPoint, faceI, oldOwn, oldNei, meshMod);

					affectedFace.unset(faceI);
                }
            }
        }
    }
    
    // 3. faces that do not get split but whose owner/neighbour change
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // I think this section is never called - have set to if(1) at the bottom to check for that
    // 
    if (title_debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Changing owner/neighbour for otherwise unaffected faces"
            << endl;
    }
    forAll(affectedFace, faceI)
    {
        if (affectedFace.get(faceI))
        {
            const face& f = mesh_.faces()[faceI];

            // The point with the lowest level should be an anchor
            // point of the neighbouring cells.
            label anchorFp = findMinLevel(f);

            label own, nei;
            getFaceNeighbours
            (
                cellAnchorPoints,
                cellAddedCells,
                faceI,
                f[anchorFp],          // Anchor point
                own,
                nei
            );

            modFace(meshMod, faceI, f, own, nei);

            // Mark face as having been handled
            if (1)
            {
				Pout<< "faceI = " << faceI << " was handled by part 3 of the face splitting section" << endl;
			}
            affectedFace.unset(faceI);
        }
    }

    // 4. new internal faces inside split cells.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (title_debug)
    {
        Pout<< "hexRef4::setRefinement :"
            << " Create new internal faces for split cells"
            << endl;
    }

    forAll(cellMidPoint, cellI)
    {
        if (cellMidPoint[cellI] >= 0)
        {
            myCreateInternalFaces
            (
                cellAnchorPoints,
                cellAddedCells,
                cellMidPoint,
                faceMidPoint,
                edgeMidPoint,
                cellI,
                meshMod,
                normalDir
            );
        }
    }

    // Extend pointLevels and cellLevels for the new cells. Could also be done
    // in updateMesh but saves passing cellAddedCells out of this routine.

    // Check
    if (debug)
    {
        label minPointI = labelMax;
        label maxPointI = labelMin;

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
                minPointI = min(minPointI, cellMidPoint[cellI]);
                maxPointI = max(maxPointI, cellMidPoint[cellI]);
            }
        }
        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] >= 0)
            {
                minPointI = min(minPointI, faceMidPoint[faceI]);
                maxPointI = max(maxPointI, faceMidPoint[faceI]);
            }
        }
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                minPointI = min(minPointI, edgeMidPoint[edgeI]);
                maxPointI = max(maxPointI, edgeMidPoint[edgeI]);
            }
        }

        if (minPointI != labelMax && minPointI != mesh_.nPoints())
        {
            FatalErrorIn("hexRef4::setRefinement(..)")
                << "Added point labels not consecutive to existing mesh points."
                << nl
                << "mesh_.nPoints():" << mesh_.nPoints()
                << " minPointI:" << minPointI
                << " maxPointI:" << maxPointI
                << abort(FatalError);
        }
    }

	//~ if (transferLater)
	//~ {
		//~ pointLevel_.transfer(newPointLevel);
		//~ cellLevel_.transfer(newCellLevel);
	//~ }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Compact cellAddedCells.

    labelListList refinedCells(cellsToRefine.size());

    forAll(cellsToRefine, i)
    {
        label cellI = cellsToRefine[i];

        refinedCells[i].transfer(cellAddedCells[cellI]);
    }
    
    if (debug) Pout<< "Final size of meshMod.points(): " << meshMod.points().size() << endl;

    return refinedCells;
}

void Foam::hexRef4::storeData
(
    const labelList& pointsToStore,
    const labelList& facesToStore,
    const labelList& cellsToStore
)
{
    savedPointLevel_.resize(2*pointsToStore.size());
    forAll(pointsToStore, i)
    {
        label pointI = pointsToStore[i];
        savedPointLevel_.insert(pointI, pointLevel_[pointI]);
    }

    savedCellLevel_.resize(2*cellsToStore.size());
    forAll(cellsToStore, i)
    {
        label cellI = cellsToStore[i];
        savedCellLevel_.insert(cellI, cellLevel_[cellI]);
    }
}


// Gets called after the mesh change. setRefinement will already have made
// sure the pointLevel_ and cellLevel_ are the size of the new mesh so we
// only need to account for reordering.
void Foam::hexRef4::updateMesh(const mapPolyMesh& map)
{
    Map<label> dummyMap(0);

    updateMesh(map, dummyMap, dummyMap, dummyMap);
}


// Gets called after the mesh change. setRefinement will already have made
// sure the pointLevel_ and cellLevel_ are the size of the new mesh so we
// only need to account for reordering.
void Foam::hexRef4::updateMesh
(
    const mapPolyMesh& map,
    const Map<label>& pointsToRestore,
    const Map<label>& facesToRestore,
    const Map<label>& cellsToRestore
)
{
    // Update celllevel
    if (debug)
    {
        Pout<< "hexRef4::updateMesh :"
            << " Updating various lists"
            << endl;
    }

    {
        const labelList& reverseCellMap = map.reverseCellMap();

        if (debug)
        {
            Pout<< "hexRef4::updateMesh :"
                << " reverseCellMap:" << map.reverseCellMap().size()
                << " cellMap:" << map.cellMap().size()
                << " nCells:" << mesh_.nCells()
                << " nOldCells:" << map.nOldCells()
                << " cellLevel_:" << cellLevel_.size()
                << " reversePointMap:" << map.reversePointMap().size()
                << " pointMap:" << map.pointMap().size()
                << " nPoints:" << mesh_.nPoints()
                << " nOldPoints:" << map.nOldPoints()
                << " pointLevel_:" << pointLevel_.size()
                << endl;
        }

        if (reverseCellMap.size() == cellLevel_.size())
        {
            // Assume it is after hexRef4 that this routine is called.
            // Just account for reordering. We cannot use cellMap since
            // then cells created from cells would get cellLevel_ of
            // cell they were created from.
            reorder(reverseCellMap, mesh_.nCells(), -1, cellLevel_);
        }
        else
        {
            // Map data
            const labelList& cellMap = map.cellMap();

            labelList newCellLevel(cellMap.size());
            forAll(cellMap, newCellI)
            {
                label oldCellI = cellMap[newCellI];

                if (oldCellI == -1)
                {
                    newCellLevel[newCellI] = -1;
                }
                else
                {
                    newCellLevel[newCellI] = cellLevel_[oldCellI];
                }
            }
            cellLevel_.transfer(newCellLevel);
        }

        // See if any cells to restore. This will be for some new cells
        // the corresponding old cell.
        forAllConstIter(Map<label>, cellsToRestore, iter)
        {
            label newCellI = iter.key();
            label storedCellI = iter();

            Map<label>::iterator fnd = savedCellLevel_.find(storedCellI);

            if (fnd == savedCellLevel_.end())
            {
                FatalErrorIn("hexRef4::updateMesh(const mapPolyMesh&)")
                    << "Problem : trying to restore old value for new cell "
                    << newCellI << " but cannot find old cell " << storedCellI
                    << " in map of stored values " << savedCellLevel_
                    << abort(FatalError);
            }
            cellLevel_[newCellI] = fnd();
        }

        //if (findIndex(cellLevel_, -1) != -1)
        //{
        //    WarningIn("hexRef4::updateMesh(const mapPolyMesh&)")
        //        << "Problem : "
        //        << "cellLevel_ contains illegal value -1 after mapping
        //        << " at cell " << findIndex(cellLevel_, -1) << endl
        //        << "This means that another program has inflated cells"
        //        << " (created cells out-of-nothing) and hence we don't know"
        //        << " their cell level. Continuing with illegal value."
        //        << abort(FatalError);
        //}
    }


    // Update pointlevel
    {
        const labelList& reversePointMap = map.reversePointMap();

        if (reversePointMap.size() == pointLevel_.size())
        {
            // Assume it is after hexRef4 that this routine is called.
            reorder(reversePointMap, mesh_.nPoints(), -1,  pointLevel_);
        }
        else
        {
            // Map data
            const labelList& pointMap = map.pointMap();

            labelList newPointLevel(pointMap.size());

            forAll(pointMap, newPointI)
            {
                label oldPointI = pointMap[newPointI];

                if (oldPointI == -1)
                {
                    //FatalErrorIn("hexRef4::updateMesh(const mapPolyMesh&)")
                    //    << "Problem : point " << newPointI
                    //    << " at " << mesh_.points()[newPointI]
                    //    << " does not originate from another point"
                    //    << " (i.e. is inflated)." << nl
                    //    << "Hence we cannot determine the new pointLevel"
                    //    << " for it." << abort(FatalError);
                    newPointLevel[newPointI] = -1;
                }
                else
                {
                    newPointLevel[newPointI] = pointLevel_[oldPointI];
                }
            }
            pointLevel_.transfer(newPointLevel);
        }

        // See if any points to restore. This will be for some new points
        // the corresponding old point (the one from the call to storeData)
        forAllConstIter(Map<label>, pointsToRestore, iter)
        {
            label newPointI = iter.key();
            label storedPointI = iter();

            Map<label>::iterator fnd = savedPointLevel_.find(storedPointI);

            if (fnd == savedPointLevel_.end())
            {
                FatalErrorIn("hexRef4::updateMesh(const mapPolyMesh&)")
                    << "Problem : trying to restore old value for new point "
                    << newPointI << " but cannot find old point "
                    << storedPointI
                    << " in map of stored values " << savedPointLevel_
                    << abort(FatalError);
            }
            pointLevel_[newPointI] = fnd();
        }

        //if (findIndex(pointLevel_, -1) != -1)
        //{
        //    WarningIn("hexRef4::updateMesh(const mapPolyMesh&)")
        //        << "Problem : "
        //        << "pointLevel_ contains illegal value -1 after mapping"
        //        << " at point" << findIndex(pointLevel_, -1) << endl
        //        << "This means that another program has inflated points"
        //        << " (created points out-of-nothing) and hence we don't know"
        //        << " their point level. Continuing with illegal value."
        //        //<< abort(FatalError);
        //}
    }

    // Update refinement tree
    if (history_.active())
    {
        history_.updateMesh(map);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Update face removal engine
    faceRemover_.updateMesh(map);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


// Gets called after mesh subsetting. Maps are from new back to old.
void Foam::hexRef4::subset
(
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap
)
{
    // Update celllevel
    if (debug)
    {
        Pout<< "hexRef4::subset :"
            << " Updating various lists"
            << endl;
    }

    if (history_.active())
    {
        WarningIn
        (
            "hexRef4::subset(const labelList&, const labelList&"
            ", const labelList&)"
        )   << "Subsetting will not work in combination with unrefinement."
            << nl
            << "Proceed at your own risk." << endl;
    }


    // Update celllevel
    {
        labelList newCellLevel(cellMap.size());

        forAll(cellMap, newCellI)
        {
            newCellLevel[newCellI] = cellLevel_[cellMap[newCellI]];
        }

        cellLevel_.transfer(newCellLevel);

        if (findIndex(cellLevel_, -1) != -1)
        {
            FatalErrorIn("hexRef4::subset(..)")
                << "Problem : "
                << "cellLevel_ contains illegal value -1 after mapping:"
                << cellLevel_
                << abort(FatalError);
        }
    }

    // Update pointlevel
    {
        labelList newPointLevel(pointMap.size());

        forAll(pointMap, newPointI)
        {
            newPointLevel[newPointI] = pointLevel_[pointMap[newPointI]];
        }

        pointLevel_.transfer(newPointLevel);

        if (findIndex(pointLevel_, -1) != -1)
        {
            FatalErrorIn("hexRef4::subset(..)")
                << "Problem : "
                << "pointLevel_ contains illegal value -1 after mapping:"
                << pointLevel_
                << abort(FatalError);
        }
    }

    // Update refinement tree
    if (history_.active())
    {
        history_.subset(pointMap, faceMap, cellMap);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Nothing needs doing to faceRemover.
    //faceRemover_.subset(pointMap, faceMap, cellMap);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


// Gets called after the mesh distribution
void Foam::hexRef4::distribute(const mapDistributePolyMesh& map)
{
    if (debug)
    {
        Pout<< "hexRef4::distribute :"
            << " Updating various lists"
            << endl;
    }

    // Update celllevel
    map.distributeCellData(cellLevel_);
    // Update pointlevel
    map.distributePointData(pointLevel_);

    // Update refinement tree
    if (history_.active())
    {
        history_.distribute(map);
    }

    // Update face removal engine
    faceRemover_.distribute(map);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


void Foam::hexRef4::checkMesh() const
{
    const scalar smallDim = 1e-6 * mesh_.bounds().mag();

    if (debug)
    {
        Pout<< "hexRef4::checkMesh : Using matching tolerance smallDim:"
            << smallDim << endl;
    }

    // Check owner on coupled faces.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // There should be only one coupled face between two cells. Why? Since
    // otherwise mesh redistribution might cause multiple faces between two
    // cells
    {
        labelList nei(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(nei, i)
        {
            nei[i] = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, nei);

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.coupled())
            {
                // Check how many faces between owner and neighbour. Should
                // be only one.
                HashTable<label, labelPair, labelPair::Hash<> >
                    cellToFace(2*pp.size());

                label faceI = pp.start();

                forAll(pp, i)
                {
                    label own = mesh_.faceOwner()[faceI];
                    label bFaceI = faceI-mesh_.nInternalFaces();

                    if (!cellToFace.insert(labelPair(own, nei[bFaceI]), faceI))
                    {
                        dumpCell(own);
                        FatalErrorIn("hexRef4::checkMesh()")
                            << "Faces do not seem to be correct across coupled"
                            << " boundaries" << endl
                            << "Coupled face " << faceI
                            << " between owner " << own
                            << " on patch " << pp.name()
                            << " and coupled neighbour " << nei[bFaceI]
                            << " has two faces connected to it:"
                            << faceI << " and "
                            << cellToFace[labelPair(own, nei[bFaceI])]
                            << abort(FatalError);
                    }

                    faceI++;
                }
            }
        }
    }

    // Check face areas.
    // ~~~~~~~~~~~~~~~~~

    {
        scalarField neiFaceAreas(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(neiFaceAreas, i)
        {
            neiFaceAreas[i] = mag(mesh_.faceAreas()[i+mesh_.nInternalFaces()]);
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, neiFaceAreas);

        forAll(neiFaceAreas, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            const scalar magArea = mag(mesh_.faceAreas()[faceI]);

            if (mag(magArea - neiFaceAreas[i]) > smallDim)
            {
                const face& f = mesh_.faces()[faceI];
                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                dumpCell(mesh_.faceOwner()[faceI]);

                FatalErrorIn("hexRef4::checkMesh()")
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << faceI
                    << " on patch " << patchI
                    << " " << mesh_.boundaryMesh()[patchI].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has face area:" << magArea
                    << " (coupled) neighbour face area differs:"
                    << neiFaceAreas[i]
                    << " to within tolerance " << smallDim
                    << abort(FatalError);
            }
        }
    }


    // Check number of points on faces.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        labelList nVerts(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(nVerts, i)
        {
            nVerts[i] = mesh_.faces()[i+mesh_.nInternalFaces()].size();
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, nVerts);

        forAll(nVerts, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            const face& f = mesh_.faces()[faceI];

            if (f.size() != nVerts[i])
            {
                dumpCell(mesh_.faceOwner()[faceI]);

                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                FatalErrorIn("hexRef4::checkMesh()")
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << faceI
                    << " on patch " << patchI
                    << " " << mesh_.boundaryMesh()[patchI].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has size:" << f.size()
                    << " (coupled) neighbour face has size:"
                    << nVerts[i]
                    << abort(FatalError);
            }
        }
    }


    // Check points of face
    // ~~~~~~~~~~~~~~~~~~~~
    {
        // Anchor points.
        pointField anchorPoints(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(anchorPoints, i)
        {
            label faceI = i+mesh_.nInternalFaces();
            const point& fc = mesh_.faceCentres()[faceI];
            const face& f = mesh_.faces()[faceI];
            const vector anchorVec(mesh_.points()[f[0]] - fc);

            anchorPoints[i] = anchorVec;
        }

        // Replace data on coupled patches with their neighbour ones. Apply
        // rotation transformation (but not separation since is relative vector
        // to point on same face.
        syncTools::swapBoundaryFaceList(mesh_, anchorPoints);

        forAll(anchorPoints, i)
        {
            label faceI = i+mesh_.nInternalFaces();
            const point& fc = mesh_.faceCentres()[faceI];
            const face& f = mesh_.faces()[faceI];
            const vector anchorVec(mesh_.points()[f[0]] - fc);

            if (mag(anchorVec - anchorPoints[i]) > smallDim)
            {
                dumpCell(mesh_.faceOwner()[faceI]);

                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                FatalErrorIn("hexRef4::checkMesh()")
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << faceI
                    << " on patch " << patchI
                    << " " << mesh_.boundaryMesh()[patchI].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has anchor vector:" << anchorVec
                    << " (coupled) neighbour face anchor vector differs:"
                    << anchorPoints[i]
                    << " to within tolerance " << smallDim
                    << abort(FatalError);
            }
        }
    }

    if (debug)
    {
        Pout<< "hexRef4::checkMesh : Returning" << endl;
    }
}


void Foam::hexRef4::checkRefinementLevels
(
    const label maxPointDiff,
    const labelList& pointsToCheck
) const
{
    if (debug)
    {
        Pout<< "hexRef4::checkRefinementLevels :"
            << " Checking 2:1 refinement level" << endl;
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorIn("hexRef4::checkRefinementLevels(const label)")
            << "cellLevel size should be number of cells"
            << " and pointLevel size should be number of points."<< nl
            << "cellLevel:" << cellLevel_.size()
            << " mesh.nCells():" << mesh_.nCells() << nl
            << "pointLevel:" << pointLevel_.size()
            << " mesh.nPoints():" << mesh_.nPoints()
            << abort(FatalError);
    }


    // Check 2:1 consistency.
    // ~~~~~~~~~~~~~~~~~~~~~~

    {
        // Internal faces.
        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label own = mesh_.faceOwner()[faceI];
            label nei = mesh_.faceNeighbour()[faceI];

            if (mag(cellLevel_[own] - cellLevel_[nei]) > 1)
            {
                dumpCell(own);
                dumpCell(nei);

                FatalErrorIn
                (
                    "hexRef4::checkRefinementLevels(const label)"
                )   << "Celllevel does not satisfy 2:1 constraint." << nl
                    << "On face " << faceI << " owner cell " << own
                    << " has refinement " << cellLevel_[own]
                    << " neighbour cell " << nei << " has refinement "
                    << cellLevel_[nei]
                    << abort(FatalError);
            }
        }

        // Coupled faces. Get neighbouring value
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

            neiLevel[i] = cellLevel_[own];
        }

        // No separation
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);

        forAll(neiLevel, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            label own = mesh_.faceOwner()[faceI];

            if (mag(cellLevel_[own] - neiLevel[i]) > 1)
            {
                dumpCell(own);

                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                FatalErrorIn
                (
                    "hexRef4::checkRefinementLevels(const label)"
                )   << "Celllevel does not satisfy 2:1 constraint."
                    << " On coupled face " << faceI
                    << " on patch " << patchI << " "
                    << mesh_.boundaryMesh()[patchI].name()
                    << " owner cell " << own << " has refinement "
                    << cellLevel_[own]
                    << " (coupled) neighbour cell has refinement "
                    << neiLevel[i]
                    << abort(FatalError);
            }
        }
    }


    // Check pointLevel is synchronized
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        labelList syncPointLevel(pointLevel_);

        // Get min level
        syncTools::syncPointList
        (
            mesh_,
            syncPointLevel,
            minEqOp<label>(),
            labelMax
        );


        forAll(syncPointLevel, pointI)
        {
            if (pointLevel_[pointI] != syncPointLevel[pointI])
            {
                FatalErrorIn
                (
                    "hexRef4::checkRefinementLevels(const label)"
                )   << "PointLevel is not consistent across coupled patches."
                    << endl
                    << "point:" << pointI << " coord:" << mesh_.points()[pointI]
                    << " has level " << pointLevel_[pointI]
                    << " whereas the coupled point has level "
                    << syncPointLevel[pointI]
                    << abort(FatalError);
            }
        }
    }


    // Check 2:1 across points (instead of faces)
    if (maxPointDiff != -1)
    {
        // Determine per point the max cell level.
        labelList maxPointLevel(mesh_.nPoints(), 0);

        forAll(maxPointLevel, pointI)
        {
            const labelList& pCells = mesh_.pointCells(pointI);

            label& pLevel = maxPointLevel[pointI];

            forAll(pCells, i)
            {
                pLevel = max(pLevel, cellLevel_[pCells[i]]);
            }
        }

        // Sync maxPointLevel to neighbour
        syncTools::syncPointList
        (
            mesh_,
            maxPointLevel,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Check 2:1 across boundary points
        forAll(pointsToCheck, i)
        {
            label pointI = pointsToCheck[i];

            const labelList& pCells = mesh_.pointCells(pointI);

            forAll(pCells, i)
            {
                label cellI = pCells[i];

                if
                (
                    mag(cellLevel_[cellI]-maxPointLevel[pointI])
                  > maxPointDiff
                )
                {
                    dumpCell(cellI);

                    FatalErrorIn
                    (
                        "hexRef4::checkRefinementLevels(const label)"
                    )   << "Too big a difference between"
                        << " point-connected cells." << nl
                        << "cell:" << cellI
                        << " cellLevel:" << cellLevel_[cellI]
                        << " uses point:" << pointI
                        << " coord:" << mesh_.points()[pointI]
                        << " which is also used by a cell with level:"
                        << maxPointLevel[pointI]
                        << abort(FatalError);
                }
            }
        }
    }


    //- Gives problems after first splitting off inside mesher.
    //// Hanging points
    //{
    //    // Any patches with points having only two edges.
    //
    //    boolList isHangingPoint(mesh_.nPoints(), false);
    //
    //    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    //
    //    forAll(patches, patchI)
    //    {
    //        const polyPatch& pp = patches[patchI];
    //
    //        const labelList& meshPoints = pp.meshPoints();
    //
    //        forAll(meshPoints, i)
    //        {
    //            label pointI = meshPoints[i];
    //
    //            const labelList& pEdges = mesh_.pointEdges()[pointI];
    //
    //            if (pEdges.size() == 2)
    //            {
    //                isHangingPoint[pointI] = true;
    //            }
    //        }
    //    }
    //
    //    syncTools::syncPointList
    //    (
    //        mesh_,
    //        isHangingPoint,
    //        andEqOp<bool>(),        // only if all decide it is hanging point
    //        true,                   // null
    //        false                   // no separation
    //    );
    //
    //    //OFstream str(mesh_.time().path()/"hangingPoints.obj");
    //
    //    label nHanging = 0;
    //
    //    forAll(isHangingPoint, pointI)
    //    {
    //        if (isHangingPoint[pointI])
    //        {
    //            nHanging++;
    //
    //            Pout<< "Hanging boundary point " << pointI
    //                << " at " << mesh_.points()[pointI]
    //                << endl;
    //            //meshTools::writeOBJ(str, mesh_.points()[pointI]);
    //        }
    //    }
    //
    //    if (returnReduce(nHanging, sumOp<label>()) > 0)
    //    {
    //        FatalErrorIn
    //        (
    //            "hexRef4::checkRefinementLevels(const label)"
    //        )   << "Detected a point used by two edges only (hanging point)"
    //            << nl << "This is not allowed"
    //            << abort(FatalError);
    //    }
    //}
}


const Foam::cellShapeList& Foam::hexRef4::cellShapes() const
{
    if (cellShapesPtr_.empty())
    {
        if (debug)
        {
            Pout<< "hexRef4::cellShapes() : calculating splitHex cellShapes."
                << " cellLevel:" << cellLevel_.size()
                << " pointLevel:" << pointLevel_.size()
                << endl;
        }

        const cellShapeList& meshShapes = mesh_.cellShapes();
        cellShapesPtr_.reset(new cellShapeList(meshShapes));

        label nSplitHex = 0;
        label nUnrecognised = 0;

        forAll(cellLevel_, cellI)
        {
            if (meshShapes[cellI].model().index() == 0)
            {
                label level = cellLevel_[cellI];

                // Return true if we've found 6 quads
                DynamicList<face> quads;
                bool haveQuads = matchHexShape
                (
                    cellI,
                    level,
                    quads
                );

                if (haveQuads)
                {
                    faceList faces(quads.xfer());
                    cellShapesPtr_()[cellI] = degenerateMatcher::match(faces);
                    nSplitHex++;
                }
                else
                {
                    nUnrecognised++;
                }
            }
        }
        if (debug)
        {
            Pout<< "hexRef4::cellShapes() :"
                << " nCells:" << mesh_.nCells() << " of which" << nl
                << "    primitive:" << (mesh_.nCells()-nSplitHex-nUnrecognised)
                << nl
                << "    split-hex:" << nSplitHex << nl
                << "    poly     :" << nUnrecognised << nl
                << endl;
        }
    }
    return cellShapesPtr_();
}



//
// Unrefinement
// ~~~~~~~~~~~~
//


Foam::labelList Foam::hexRef4::getSplitPoints() const
{
    if (debug)
    {
        checkRefinementLevels(-1, labelList(0));
    }

    if (debug)
    {
        Pout<< "hexRef4::getSplitPoints :"
            << " Calculating unrefineable points" << endl;
    }


    if (!history_.active())
    {
        FatalErrorIn("hexRef4::getSplitPoints()")
            << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    // Master cell
    // -1 undetermined
    // -2 certainly not split point
    // >= label of master cell
    labelList splitMaster(mesh_.nPoints(), -1);
    labelList splitMasterLevel(mesh_.nPoints(), 0);

    // Unmark all with not 8 cells
    //const labelListList& pointCells = mesh_.pointCells();

    for (label pointI = 0; pointI < mesh_.nPoints(); pointI++)
    {
        const labelList& pCells = mesh_.pointCells(pointI);

        if (pCells.size() != 8)
        {
            splitMaster[pointI] = -2;
        }
    }

    // Unmark all with different master cells
    const labelList& visibleCells = history_.visibleCells();

    forAll(visibleCells, cellI)
    {
        const labelList& cPoints = mesh_.cellPoints(cellI);

        if (visibleCells[cellI] != -1 && history_.parentIndex(cellI) >= 0)
        {
            label parentIndex = history_.parentIndex(cellI);

            // Check same master.
            forAll(cPoints, i)
            {
                label pointI = cPoints[i];

                label masterCellI = splitMaster[pointI];

                if (masterCellI == -1)
                {
                    // First time visit of point. Store parent cell and
                    // level of the parent cell (with respect to cellI). This
                    // is additional guarantee that we're referring to the
                    // same master at the same refinement level.

                    splitMaster[pointI] = parentIndex;
                    splitMasterLevel[pointI] = cellLevel_[cellI] - 1;
                }
                else if (masterCellI == -2)
                {
                    // Already decided that point is not splitPoint
                }
                else if
                (
                    (masterCellI != parentIndex)
                 || (splitMasterLevel[pointI] != cellLevel_[cellI] - 1)
                )
                {
                    // Different masters so point is on two refinement
                    // patterns
                    splitMaster[pointI] = -2;
                }
            }
        }
        else
        {
            // Either not visible or is unrefined cell
            forAll(cPoints, i)
            {
                label pointI = cPoints[i];

                splitMaster[pointI] = -2;
            }
        }
    }

    // Unmark boundary faces
    for
    (
        label faceI = mesh_.nInternalFaces();
        faceI < mesh_.nFaces();
        faceI++
    )
    {
        const face& f = mesh_.faces()[faceI];

        forAll(f, fp)
        {
            splitMaster[f[fp]] = -2;
        }
    }


    // Collect into labelList

    label nSplitPoints = 0;

    forAll(splitMaster, pointI)
    {
        if (splitMaster[pointI] >= 0)
        {
            nSplitPoints++;
        }
    }

    labelList splitPoints(nSplitPoints);
    nSplitPoints = 0;

    forAll(splitMaster, pointI)
    {
        if (splitMaster[pointI] >= 0)
        {
            splitPoints[nSplitPoints++] = pointI;
        }
    }

    return splitPoints;
}


//void Foam::hexRef4::markIndex
//(
//    const label maxLevel,
//    const label level,
//    const label index,
//    const label markValue,
//    labelList& indexValues
//) const
//{
//    if (level < maxLevel && indexValues[index] == -1)
//    {
//        // Mark
//        indexValues[index] = markValue;
//
//        // Mark parent
//        const splitCell4& split = history_.splitCells()[index];
//
//        if (split.parent_ >= 0)
//        {
//            markIndex
//            (
//              maxLevel, level+1, split.parent_, markValue, indexValues);
//            )
//        }
//    }
//}
//
//
//// Get all cells which (down to level) originate from the same cell.
//// level=0 returns cell only, level=1 returns the 8 cells this cell
//// originates from, level=2 returns 64 cells etc.
//// If the cell does not originate from refinement returns just itself.
//void Foam::hexRef4::markCellClusters
//(
//    const label maxLevel,
//    labelList& cluster
//) const
//{
//    cluster.setSize(mesh_.nCells());
//    cluster = -1;
//
//    const DynamicList<splitCell4>& splitCells = history_.splitCells();
//
//    // Mark all splitCells down to level maxLevel with a cell originating from
//    // it.
//
//    labelList indexLevel(splitCells.size(), -1);
//
//    forAll(visibleCells, cellI)
//    {
//        label index = visibleCells[cellI];
//
//        if (index >= 0)
//        {
//            markIndex(maxLevel, 0, index, cellI, indexLevel);
//        }
//    }
//
//    // Mark cells with splitCell
//}


Foam::labelList Foam::hexRef4::consistentUnrefinement
(
    const labelList& pointsToUnrefine,
    const bool maxSet
) const
{
    if (debug)
    {
        Pout<< "hexRef4::consistentUnrefinement :"
            << " Determining 2:1 consistent unrefinement" << endl;
    }

    if (maxSet)
    {
        FatalErrorIn
        (
            "hexRef4::consistentUnrefinement(const labelList&, const bool"
        )   << "maxSet not implemented yet."
            << abort(FatalError);
    }

    // Loop, modifying pointsToUnrefine, until no more changes to due to 2:1
    // conflicts.
    // maxSet = false : unselect points to refine
    // maxSet = true: select points to refine

    // Maintain boolList for pointsToUnrefine and cellsToUnrefine
    PackedBoolList unrefinePoint(mesh_.nPoints());

    forAll(pointsToUnrefine, i)
    {
        label pointI = pointsToUnrefine[i];

        unrefinePoint.set(pointI);
    }


    while (true)
    {
        // Construct cells to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        PackedBoolList unrefineCell(mesh_.nCells());

        forAll(unrefinePoint, pointI)
        {
            if (unrefinePoint.get(pointI))
            {
                const labelList& pCells = mesh_.pointCells(pointI);

                forAll(pCells, j)
                {
                    unrefineCell.set(pCells[j]);
                }
            }
        }


        label nChanged = 0;


        // Check 2:1 consistency taking refinement into account
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Internal faces.
        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label own = mesh_.faceOwner()[faceI];
            label ownLevel = cellLevel_[own] - unrefineCell.get(own);

            label nei = mesh_.faceNeighbour()[faceI];
            label neiLevel = cellLevel_[nei] - unrefineCell.get(nei);

            if (ownLevel < (neiLevel-1))
            {
                // Since was 2:1 this can only occur if own is marked for
                // unrefinement.

                if (maxSet)
                {
                    unrefineCell.set(nei);
                }
                else
                {
                    // could also combine with unset:
                    // if (!unrefineCell.unset(own))
                    // {
                    //     FatalErrorIn("hexRef4::consistentUnrefinement(..)")
                    //         << "problem cell already unset"
                    //         << abort(FatalError);
                    // }
                    if (unrefineCell.get(own) == 0)
                    {
                        FatalErrorIn("hexRef4::consistentUnrefinement(..)")
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(own);
                }
                nChanged++;
            }
            else if (neiLevel < (ownLevel-1))
            {
                if (maxSet)
                {
                    unrefineCell.set(own);
                }
                else
                {
                    if (unrefineCell.get(nei) == 0)
                    {
                        FatalErrorIn("hexRef4::consistentUnrefinement(..)")
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(nei);
                }
                nChanged++;
            }
        }


        // Coupled faces. Swap owner level to get neighbouring cell level.
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

            neiLevel[i] = cellLevel_[own] - unrefineCell.get(own);
        }

        // Swap to neighbour
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);

        forAll(neiLevel, i)
        {
            label faceI = i+mesh_.nInternalFaces();
            label own = mesh_.faceOwner()[faceI];
            label ownLevel = cellLevel_[own] - unrefineCell.get(own);

            if (ownLevel < (neiLevel[i]-1))
            {
                if (!maxSet)
                {
                    if (unrefineCell.get(own) == 0)
                    {
                        FatalErrorIn("hexRef4::consistentUnrefinement(..)")
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(own);
                    nChanged++;
                }
            }
            else if (neiLevel[i] < (ownLevel-1))
            {
                if (maxSet)
                {
                    if (unrefineCell.get(own) == 1)
                    {
                        FatalErrorIn("hexRef4::consistentUnrefinement(..)")
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.set(own);
                    nChanged++;
                }
            }
        }

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef4::consistentUnrefinement :"
                << " Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }


        // Convert cellsToUnrefine back into points to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Knock out any point whose cell neighbour cannot be unrefined.
        forAll(unrefinePoint, pointI)
        {
            if (unrefinePoint.get(pointI))
            {
                const labelList& pCells = mesh_.pointCells(pointI);

                forAll(pCells, j)
                {
                    if (!unrefineCell.get(pCells[j]))
                    {
                        unrefinePoint.unset(pointI);
                        break;
                    }
                }
            }
        }
    }


    // Convert back to labelList.
    label nSet = 0;

    forAll(unrefinePoint, pointI)
    {
        if (unrefinePoint.get(pointI))
        {
            nSet++;
        }
    }

    labelList newPointsToUnrefine(nSet);
    nSet = 0;

    forAll(unrefinePoint, pointI)
    {
        if (unrefinePoint.get(pointI))
        {
            newPointsToUnrefine[nSet++] = pointI;
        }
    }

    return newPointsToUnrefine;
}


void Foam::hexRef4::setUnrefinement
(
    const labelList& splitPointLabels,
    polyTopoChange& meshMod
)
{
    if (!history_.active())
    {
        FatalErrorIn
        (
            "hexRef4::setUnrefinement(const labelList&, polyTopoChange&)"
        )   << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "hexRef4::setUnrefinement :"
            << " Checking initial mesh just to make sure" << endl;

        checkMesh();

        forAll(cellLevel_, cellI)
        {
            if (cellLevel_[cellI] < 0)
            {
                FatalErrorIn
                (
                    "hexRef4::setUnrefinement"
                    "("
                        "const labelList&, "
                        "polyTopoChange&"
                    ")"
                )   << "Illegal cell level " << cellLevel_[cellI]
                    << " for cell " << cellI
                    << abort(FatalError);
            }
        }


        // Write to sets.
        pointSet pSet(mesh_, "splitPoints", splitPointLabels);
        pSet.write();

        cellSet cSet(mesh_, "splitPointCells", splitPointLabels.size());

        forAll(splitPointLabels, i)
        {
            const labelList& pCells = mesh_.pointCells(splitPointLabels[i]);

            forAll(pCells, j)
            {
                cSet.insert(pCells[j]);
            }
        }
        cSet.write();

        Pout<< "hexRef4::setRefinement : Dumping " << pSet.size()
            << " points and "
            << cSet.size() << " cells for unrefinement to" << nl
            << "    pointSet " << pSet.objectPath() << nl
            << "    cellSet " << cSet.objectPath()
            << endl;
    }


    labelList cellRegion;
    labelList cellRegionMaster;
    labelList facesToRemove;

    {
		// There are 12 faces per split point in 3D for internal cells (fewer if on a border)
		// For 2D, there are 8 per splitPoint (always on a border)
        //~ labelHashSet splitFaces(12*splitPointLabels.size());
        labelHashSet splitFaces(8*splitPointLabels.size());

        forAll(splitPointLabels, i)
        {
            const labelList& pFaces = mesh_.pointFaces()[splitPointLabels[i]];

            forAll(pFaces, j)
            {
                splitFaces.insert(pFaces[j]);
            }
        }

        // Check with faceRemover what faces will get removed. Note that this
        // can be more (but never less) than splitFaces provided.
        faceRemover_.compatibleRemoves
        (
            splitFaces.toc(),   // pierced faces
            cellRegion,         // per cell -1 or region it is merged into
            cellRegionMaster,   // per region the master cell
            facesToRemove       // new faces to be removed.
        );

        if (facesToRemove.size() != splitFaces.size())
        {
            FatalErrorIn
            (
                "hexRef4::setUnrefinement(const labelList&, polyTopoChange&)"
            )   << "Ininitial set of split points to unrefine does not"
                << " seem to be consistent or not mid points of refined cells"
                << abort(FatalError);
        }
    }

    // Redo the region master so it is consistent with our master.
    // This will guarantee that the new cell (for which faceRemover uses
    // the region master) is already compatible with our refinement structure.

    forAll(splitPointLabels, i)
    {
        label pointI = splitPointLabels[i];

        // Get original cell label

        const labelList& pCells = mesh_.pointCells(pointI);

        // Check
        // It seems like cells on the unrefined side of a refinement line would fail this?
        // However, the same would presumably be true in 3D, so I assume there is a reason
        // that they don't (they're not listed?)
        if (pCells.size() != 4)
        {
            FatalErrorIn
            (
                "hexRef4::setUnrefinement(const labelList&, polyTopoChange&)"
            )   << "splitPoint " << pointI
                << " should have 4 cells using it. It has " << pCells
                << abort(FatalError);
        }


        // Check that the lowest numbered pCells is the master of the region
        // (should be guaranteed by directRemoveFaces)
        //if (debug)
        {
            label masterCellI = min(pCells);

            forAll(pCells, j)
            {
                label cellI = pCells[j];

                label region = cellRegion[cellI];

                if (region == -1)
                {
                    FatalErrorIn("hexRef4::setUnrefinement(..)")
                        << "Ininitial set of split points to unrefine does not"
                        << " seem to be consistent or not mid points"
                        << " of refined cells" << nl
                        << "cell:" << cellI << " on splitPoint " << pointI
                        << " has no region to be merged into"
                        << abort(FatalError);
                }

                if (masterCellI != cellRegionMaster[region])
                {
                    FatalErrorIn("hexRef4::setUnrefinement(..)")
                        << "cell:" << cellI << " on splitPoint:" << pointI
                        << " in region " << region
                        << " has master:" << cellRegionMaster[region]
                        << " which is not the lowest numbered cell"
                        << " among the pointCells:" << pCells
                        << abort(FatalError);
                }
            }
        }
    }

    // Insert all commands to combine cells. Never fails so don't have to
    // test for success.
    faceRemover_.setRefinement
    (
        facesToRemove,
        cellRegion,
        cellRegionMaster,
        meshMod
    );

    // Remove the 8 cells that originated from merging around the split point
    // and adapt cell levels (not that pointLevels stay the same since points
    // either get removed or stay at the same position.
    
    // Although the comment here specifically mentions 8 cells, I think the code
    // should work for 4 cells as well.
    forAll(splitPointLabels, i)
    {
        label pointI = splitPointLabels[i];

        const labelList& pCells = mesh_.pointCells(pointI);

        label masterCellI = min(pCells);

        forAll(pCells, j)
        {
            cellLevel_[pCells[j]]--;
        }
		// combineCells also appears to work independently of number of children
        history_.combineCells(masterCellI, pCells);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // history_.updateMesh will take care of truncating.
}


// Write refinement to polyMesh directory.
bool Foam::hexRef4::write() const
{
    bool writeOk =
        cellLevel_.write()
     && pointLevel_.write()
     && level0Edge_.write();

    if (history_.active())
    {
        writeOk = writeOk && history_.write();
    }

    return writeOk;
}


// ************************************************************************* //
