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

Class
    Foam::hexRef4

Description
    Refinement of (split) hexes using polyTopoChange.

SourceFiles
    hexRef4.C

\*---------------------------------------------------------------------------*/

#ifndef hexRef4_H
#define hexRef4_H

#include "labelIOList.H"
#include "face.H"
#include "HashSet.H"
#include "DynamicList.H"
#include "primitivePatch.H"
#include "removeFaces.H"
#include "refinementTree.H"
#include "PackedBoolList.H"
#include "uniformDimensionedFields.H"
#include "cellShapeList.H"

#include "fvMesh.H"
#include "primitiveMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyPatch;
class polyTopoChange;
class mapPolyMesh;
class mapDistributePolyMesh;

/*---------------------------------------------------------------------------*\
                           Class hexRef4 Declaration
\*---------------------------------------------------------------------------*/

class hexRef4
{
    // Private data

        //- Reference to underlying mesh.
        const polyMesh& mesh_;

        //- Per cell the refinement level
        labelIOList cellLevel_;
        labelIOList oldCellLevel_;

        //- Per point the refinement level
        labelIOList pointLevel_;
        labelIOList oldPointLevel_;

        //- Typical edge length between unrefined points
        uniformDimensionedScalarField level0Edge_;

        //- Refinement history
        refinementTree history_;
        
        // cell Depth and aspect ratio
        const vector cell0Lengths_;
        const double aspectX_to_Y_; // Multiply x per cell by this to get y per cell
        

        //- Face remover engine
        removeFaces faceRemover_;

        //- Level of saved points
        Map<label> savedPointLevel_;

        //- Level of saved cells
        Map<label> savedCellLevel_;

        //- cell shapes when seen as split hexes
        mutable autoPtr<cellShapeList> cellShapesPtr_;
        
    // Private Member Functions

		Foam::vector calcSingleFaceNormal(const pointField&, const label&) const;
		int calcRelevantDirs(const Foam::vector&, const bool Fatal);
		int calcRelevantDirs(const Foam::vector&, const label& owner, const label& neighbour, const bool Fatal);

		//- Set up a vectorField holding the vector normals of each face.
		void calcFaceNormalVector(const pointField&, vectorField&) const;
		
		void setIndicators(const label& own, const label& nei, label& ownInd, label& neiInd, const bool refinementBool[]);
		
		double getCellLength(const label& dimension);
		double getAspectRatio();
		
		void setBoundaryCellInfo
		(
			const label& faceI,
			const labelList& parentAddedCells, 
			label& own,
			const int a
		);
		
		//~ void setInternalFaceInfo
		//~ (
			//~ const labelListList& cellAddedCells,
			//~ const label& set,
			//~ label& own,
			//~ label& nei,
			//~ const polyTopoChange& meshMod,
			//~ bool& isAFace
		//~ );
		
		void setInternalFaceInfo
		(
			const labelListList& cellAddedCells,
			const label& set,
			label& own,
			label& nei,
			const polyTopoChange& meshMod
		);
		
		void myGetFaceVertsLists
		(
			DynamicList<label> (&newFaceVertsA)[2],
			const labelList& edgeMidPoint,
			const label& faceI,
			polyTopoChange& meshMod
		);
		
		void mySwitchFaceVerts
		(
			labelList& faceVerts,
			std::string reason = "Switching vertices for an unspecified reason"
		);
		
		void mySplitSideFaces
		(
			const labelListList& cellAddedCells,
			const labelList& faceMidPoint,
			const labelList& edgeMidPoint,
			const label& faceI,
			const label& oldOwn,
			const label& oldNei,
			polyTopoChange& meshMod
		);
		
		void myCreateInternalFaces
		(
			const labelListList& cellAnchorPoints,
			const labelListList& cellAddedCells,
			const labelList& faceMidPoint,
			const labelList& edgeMidPoint,
			const label& cellI,

			polyTopoChange& meshMod,
			const Foam::vector& normalDir
		);
		
		void myCreateInternalFacesSubSection
		(
			DynamicList<label>& newFaceVerts,
			const label& relevantDir,
			const label& cellI,
			const label& faceI,
			const labelListList& cellAddedCells,
			polyTopoChange& meshMod,
			label& nFacesAdded
		);
		
		void myCombineFaces
		(
			const labelList& facesToCombine,
			DynamicList<Pair<label> >& pairedFaces
		);
		
		Foam::label findEdge
		(
			const labelList& targetList,
			const label targetLabel
		);
    
        //- Reorder according to map.
        static void reorder
        (
            const labelList& map,
            const label len,
            const label null,
            labelList& elems
        );

        //- Get patch and zone info
        void getFaceInfo
        (
            const label faceI,
            label& patchID,
            label& zoneID,
            label& zoneFlip
        ) const;

        //- Adds a face on top of existing faceI. Reverses if necessary.
        label addFace
        (
            polyTopoChange& meshMod,
            const label faceI,
            const face& newFace,
            const label own,
            const label nei
        ) const;

        //- Adds internal face from point. No checks on reversal.
        label addInternalFace
        (
            polyTopoChange& meshMod,
            const label meshFaceI,
            const label meshPointI,
            const face& newFace,
            const label own,
            const label nei
        ) const;

        //- Modifies existing faceI for either new owner/neighbour or new face
        //  points. Reverses if necessary.
        void modFace
        (
            polyTopoChange& meshMod,
            const label faceI,
            const face& newFace,
            const label own,
            const label nei
        ) const;

        scalar getLevel0EdgeLength() const;

        //- Get cell added to point of cellI (if any)
        label getAnchorCell
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const label cellI,
            const label faceI,
            const label pointI
        ) const;

        //- Get new owner and neighbour (in unspecified order) of pointI
        //  on faceI.
        void getFaceNeighbours
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const label faceI,
            const label pointI,

            label& own,
            label& nei
        ) const;


        //- Get index of minimum pointlevel.
        label findMinLevel(const labelList& f) const;
        //- Get maximum pointlevel.
        label findMaxLevel(const labelList& f) const;
        //- Count number of vertices <= anchorLevel
        label countAnchors(const labelList&, const label) const;
        //- Debugging: dump cell as .obj file
        void dumpCell(const label cellI) const;
        //- Find index of point with wantedLevel, starting from fp.
        label findLevel
        (
            const label faceI,
            const face& f,
            const label startFp,
            const bool searchForward,
            const label wantedLevel
        ) const;

        ////- Print levels of list of points.
        //void printLevels(Ostream&, const labelList&) const;

        //- debug:check orientation of added internal face
        static void checkInternalOrientation
        (
            polyTopoChange& meshMod,
            const label cellI,
            const label faceI,
            const point& ownPt,
            const point& neiPt,
            const face& newFace
        );

        //- debug:check orientation of new boundary face
        static void checkBoundaryOrientation
        (
            polyTopoChange& meshMod,
            const label cellI,
            const label faceI,
            const point& ownPt,
            const point& boundaryPt,
            const face& newFace
        );

        //- If p0 and p1 are existing vertices check if edge is split and insert
        // splitPoint.
        void insertEdgeSplit
        (
            const labelList& edgeMidPoint,
            const label p0,
            const label p1,
            DynamicList<label>& verts
        ) const;

        //- Store in maps correspondence from midpoint to anchors and faces.
        label storeMidPointInfo
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const labelList& cellMidPoint,
            const labelList& edgeMidPoint,
            const label cellI,
            const label faceI,
            const bool faceOrder,
            const label midPointI,
            const label anchorPointI,
            const label faceMidPointI,

            Map<edge>& midPointToAnchors,
            Map<edge>& midPointToFaceMids,
            polyTopoChange& meshMod
        ) const;

        //- Create all internal faces from an unsplit face.
        void createInternalFromSplitFace
        (
            const labelListList& cellAnchorPoints,
            const labelListList& cellAddedCells,
            const labelList& cellMidPoint,
            const labelList& faceMidPoint,
            const labelList& edgeMidPoint,
            const label cellI,
            const label faceI,

            Map<edge>& midPointToAnchors,
            Map<edge>& midPointToFaceMids,
            polyTopoChange& meshMod,
            label& nFacesAdded
        ) const;

        //- Store vertices from startFp upto face split point.
        //  Used when splitting face into 4.
        void walkFaceToMid
        (
            const labelList& edgeMidPoint,
            const label cLevel,
            const label faceI,
            const label startFp,
            DynamicList<label>& faceVerts
        ) const;

        //- Same as walkFaceToMid but now walk back.
        void walkFaceFromMid
        (
            const labelList& edgeMidPoint,
            const label cLevel,
            const label faceI,
            const label startFp,
            DynamicList<label>& faceVerts
        ) const;

        //- Updates refineCell so consistent 2:1 refinement. Returns local
        //  number of cells changed.
        label faceConsistentRefinement
        (
            const bool maxSet,
            PackedBoolList& refineCell
        ) const;

        //- Check wanted refinement for 2:1 consistency
        void checkWantedRefinementLevels(const labelList&) const;


        // Cellshape recognition

            //- Collect all points on face of certain level
            void collectLevelPoints
            (
                const labelList& f,
                const label level,
                DynamicList<label>& points
            ) const;

            //- Collect all points on face (in local numbering) of certain level
            void collectLevelPoints
            (
                const labelList& meshPoints,
                const labelList& f,
                const label level,
                DynamicList<label>& points
            ) const;

            //- Collect all faces with four corner points and return true if
            //  hex was matched (6 faces of each four corner points)
            bool matchHexShape
            (
                const label cellI,
                const label cellLevel,
                DynamicList<face>& quads
            ) const;


        //- Disallow default bitwise copy construct
        hexRef4(const hexRef4&);

        //- Disallow default bitwise assignment
        void operator=(const hexRef4&);


public:

    //- Runtime type information
    ClassName("hexRef4");

    static bool extended_debug;

    // Constructors

        //- Construct from mesh, read_if_present refinement data
        //  (from write below). If readHistory is true does read_if_present
        //  of refinement history. If false clears all history
        hexRef4(const polyMesh& mesh, const bool readHistory = true);

        //- Construct from mesh and un/refinement data and optional size of
        //  starting cells
        hexRef4
        (
            const polyMesh& mesh,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const refinementTree& history,
            const scalar level0Edge = -GREAT
        );

        //- Construct from mesh and refinement data and optional size of
        //  starting cells
        hexRef4
        (
            const polyMesh& mesh,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const scalar level0Edge = -GREAT
        );


    // Member Functions

        // Access

            const polyMesh& mesh() const
            {
                return mesh_;
            }

            const labelIOList& cellLevel() const
            {
                return cellLevel_;
            }

            const labelIOList& pointLevel() const
            {
                return pointLevel_;
            }

            const refinementTree& history() const
            {
                return history_;
            }

            //- Typical edge length between unrefined points
            scalar level0EdgeLength() const
            {
                return level0Edge_.value();
            }

        // Refinement

            //- Gets level such that the face has four points <= level.
            label getAnchorLevel(const label faceI) const;
			
            //- Given valid mesh and current cell level and proposed
            //  cells to refine calculate any clashes (due to 2:1) and return
            //  ok list of cells to refine.
            //  Either adds cells to refine to set (maxSet = true) or
            //  removes cells to refine (maxSet = false)
            labelList consistentRefinement
            (
                const labelList& cellsToRefine,
                const bool maxSet
            ) const;
            
			/*
            //- Like consistentRefinement but slower:
            //
            //  - specify number of cells between consecutive refinement levels
            //    (consistentRefinement equivalent to 1)
            //  - specify max level difference between point-connected cells.
            //    (-1 to disable) Note that with normal 2:1 limitation
            //    (maxFaceDiff=1) there can be 8:1 size difference across point
            //    connected cells so maxPointDiff allows you to make that less.
            //  cellsToRefine : cells we're thinking about refining. It will
            //                  extend this set. All refinement levels will be
            //                  at least maxFaceDiff layers thick.
            //  facesToCheck  : additional faces where to implement the
            //                  maxFaceDiff thickness (usually only boundary
            //                  faces)
            labelList consistentSlowRefinement
            (
                const label maxFaceDiff,
                const labelList& cellsToRefine,
                const labelList& facesToCheck,
                const label maxPointDiff,
                const labelList& pointsToCheck
            ) const;

            //- Like consistentSlowRefinement but uses different meshWave
            //  (proper distance instead of toplogical count). No point checks
            //  yet.
            labelList consistentSlowRefinement2
            (
                const label maxFaceDiff,
                const labelList& cellsToRefine,
                const labelList& facesToCheck
            ) const;*/

            //- Insert refinement. All selected cells will be split into 8.
            //  Returns per element in cells the 8 cells they were split into.
            //  Guarantees that the 0th element is the original cell label.
            //  Mapping:
            //  -split cells: 7 new ones get added from original
            //  -split faces: original gets modified; 3 new ones get added
            //               from original
            //  -added internal faces: added from original cell face(if
            //   that was internal) or created out-of-nothing (so will not
            //   get mapped!). Note: could make this inflate from point but
            //   that will allocate interpolation.
            //  -points added to split edge: added from edge start()
            //  -midpoints added: added from cellPoints[0].
            labelListList setRefinement
            (
                const labelList& cells,
                polyTopoChange&,
                const vector& //,polyMesh&
            );

            //- Update local numbering for changed mesh.
            void updateMesh(const mapPolyMesh&);


            // Restoring : is where other processes delete and reinsert data.
            //             These callbacks allow this to restore the cellLevel
            //             and pointLevel for reintroduced points.
            //             Is not related to undoing my refinement

                //- Signal points/face/cells for which to store data
                void storeData
                (
                    const labelList& pointsToStore,
                    const labelList& facesToStore,
                    const labelList& cellsToStore
                );

                //- Update local numbering + undo
                //  Data to restore given as new pointlabel + stored pointlabel
                //  (i.e. what was in pointsToStore)
                void updateMesh
                (
                    const mapPolyMesh&,
                    const Map<label>& pointsToRestore,
                    const Map<label>& facesToRestore,
                    const Map<label>& cellsToRestore
                );


            //- Update local numbering for subsetted mesh.
            //  Gets new-to-old maps. Not compatible with unrefinement.
            void subset
            (
                const labelList& pointMap,
                const labelList& faceMap,
                const labelList& cellMap
            );

            //- Update local numbering for mesh redistribution
            void distribute(const mapDistributePolyMesh&);

            //- Debug: Check coupled mesh for correctness
            void checkMesh() const;

            //- Debug: Check 2:1 consistency across faces.
            //  maxPointDiff==-1 : only check 2:1 across faces
            //  maxPointDiff!=-1 : check point-connected cells.
            void checkRefinementLevels
            (
                const label maxPointDiff,
                const labelList& pointsToCheck
            ) const;

            //- Utility: get hexes as cell shapes
            const cellShapeList& cellShapes() const;


        // Unrefinement (undoing refinement, not arbitrary coarsening)

            //- Return the points at the centre of top-level split cells
            //  that can be unsplit.
            labelList getSplitPoints() const;

            //- Given proposed
            //  splitPoints to unrefine according to calculate any clashes
            //  (due to 2:1) and return ok list of points to unrefine.
            //  Either adds points to refine to set (maxSet = true) or
            //  removes points to refine (maxSet = false)
            labelList consistentUnrefinement
            (
                const labelList& pointsToUnrefine,
                const bool maxSet
            ) const;

            //- Remove some refinement. Needs to be supplied output of
            //  consistentUnrefinement. Only call if undoable set.
            //  All 8 pointCells of a split point will be combined into
            //  the lowest numbered cell of those 8.
            void setUnrefinement
            (
                const labelList& splitPointLabels,
                polyTopoChange&
            );

        // Write

            // Set instance for mesh files
            void setInstance(const fileName& inst);

            //- Force writing refinement+history to polyMesh directory.
            bool write() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
