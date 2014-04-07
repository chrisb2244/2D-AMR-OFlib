/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::refinementHistory

Description
    All refinement history. Used in unrefinement.

    - visibleCells: valid for the current mesh and contains per cell -1
    (cell unrefined) or an index into splitCells_.
    - splitCells: for every split contains the parent (also index into
      splitCells) and optionally a subsplit as 8 indices into splitCells.
      Note that the numbers in splitCells are not cell labels, they are purely
      indices into splitCells.

    E.g. 2 cells, cell 1 (second cell!) gets refined so end up with 9 cells:
    \verbatim
        // splitCells
        9
        (
        -1 (1 2 3 4 5 6 7 8)
        0 0()
        0 0()
        0 0()
        0 0()
        0 0()
        0 0()
        0 0()
        0 0()
        )

        // visibleCells
        9(-1 1 2 3 4 5 6 7 8)
    \endverbatim


    So cell0 (visibleCells=-1) is unrefined.
    Cells 1-8 have all valid splitCells entries which are:
      - parent:0
      - subsplits:0()

    The parent 0 refers back to the splitcell entries.


SourceFiles
    refinementHistory.C

\*---------------------------------------------------------------------------*/

#ifndef refinementHistory_H
#define refinementHistory_H

#include "DynamicList.H"
#include "labelList.H"
#include "FixedList.H"
#include "SLList.H"
#include "autoPtr.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class mapPolyMesh;
class mapDistributePolyMesh;

/*---------------------------------------------------------------------------*\
                           Class refinementHistory Declaration
\*---------------------------------------------------------------------------*/

class refinementTree
:
    public regIOobject
{
public:

    class splitCell4
    {
    public:

        // Index to original splitCell this cell was refined off from
        // -1: top level cell
        // -2: free splitCell (so should also be in freeSplitCells_)
        label parent_;

        //- cells this cell was refined into
        autoPtr<FixedList<label, 4> > addedCellsPtr_;

        //- Construct null (parent = -1)
        splitCell4();

        //- Construct from parent
        splitCell4(const label parent);

        //- Construct from Istream
        splitCell4(Istream& is);

        //- Construct as deep copy
        splitCell4(const splitCell4&);

        //- Copy operator since autoPtr otherwise 'steals' storage.
        void operator=(const splitCell4& s)
        {
            // Check for assignment to self
            if (this == &s)
            {
                FatalErrorIn("splitCell4::operator=(const Foam::splitCell4&)")
                    << "Attempted assignment to self"
                    << abort(FatalError);
            }

            parent_ = s.parent_;

            addedCellsPtr_.reset
            (
                s.addedCellsPtr_.valid()
              ? new FixedList<label, 4>(s.addedCellsPtr_())
              : NULL
            );
        }

        bool operator==(const splitCell4& s) const
        {
            if (addedCellsPtr_.valid() != s.addedCellsPtr_.valid())
            {
                return false;
            }
            else if (parent_ != s.parent_)
            {
                return false;
            }
            else if (addedCellsPtr_.valid())
            {
                return addedCellsPtr_() == s.addedCellsPtr_();
            }
            else
            {
                return true;
            }
        }

        bool operator!=(const splitCell4& s) const
        {
            return !operator==(s);
        }

        friend Istream& operator>>(Istream&, splitCell4&);
        friend Ostream& operator<<(Ostream&, const splitCell4&);
    }; // end of splitCell4 class


private:

    // Private data

        //- Storage for splitCells
        DynamicList<splitCell4> splitCells_;

        //- Unused indices in splitCells
        DynamicList<label> freeSplitCells_;

        //- Currently visible cells. Indices into splitCells.
        labelList visibleCells_;


    // Private Member Functions

        //- Debug write
        static void writeEntry
        (
            const List<splitCell4>&,
            const splitCell4&
        );
        //- Debug write
        static void writeDebug
        (
            const labelList&,
            const List<splitCell4>&
        );

        //- Check consistency of structure, i.e. indices into splitCells_.
        void checkIndices() const;

        //- Allocate a splitCell. Return index in splitCells_.
        label allocateSplitCell(const label parent, const label i);

        //- Free a splitCell.
        void freeSplitCell(const label index);

        //- Mark entry in splitCells. Recursively mark its parent and subs.
        void markSplit
        (
            const label,
            labelList& oldToNew,
            DynamicList<splitCell4>&
        ) const;

        void countProc
        (
            const label index,
            const label newProcNo,
            labelList& splitCellProc,
            labelList& splitCellNum
        ) const;

public:

    // Declare name of the class and its debug switch
    TypeName("refinementTree");


    // Constructors

        //- Construct (read) given an IOobject
        refinementTree(const IOobject&);

        //- Construct (read) or construct null
        // Seems like this gets cut in the refinementTree in pdf
        refinementTree
        (
            const IOobject&,
            const List<splitCell4>& splitCells,
            const labelList& visibleCells
        );

        //- Construct (read) or construct from initial number of cells
        //  (all visible)
        refinementTree(const IOobject&, const label nCells);

        //- Construct as copy
        refinementTree(const IOobject&, const refinementTree&);

        //- Construct from Istream
        refinementTree(const IOobject&, Istream&);


    // Member Functions


        //- Per cell in the current mesh (i.e. visible) either -1 (unrefined)
        //  or an index into splitCells.
        const labelList& visibleCells() const
        {
            return visibleCells_;
        }

        //- Storage for splitCell8s.
        const DynamicList<splitCell4>& splitCells() const
        {
            return splitCells_;
        }

        //- Cache of unused indices in splitCells
        const DynamicList<label>& freeSplitCells() const
        {
            return freeSplitCells_;
        }

        //- Is there unrefinement history. Note that this will fall over if
        //  there are 0 cells in the mesh. But this gives problems with
        //  lots of other programs anyway.
        bool active() const
        {
            return visibleCells_.size() > 0;
        }

        //- Get parent of cell
        label parentIndex(const label cellI) const
        {
            label index = visibleCells_[cellI];

            if (index < 0)
            {
                FatalErrorIn("refinementHistory::parentIndex(const label)")
                    << "Cell " << cellI << " is not refined and visible"
                    << abort(FatalError);
            }
            return splitCells_[index].parent_;
        }
        
        label myParentCell(const label cellI) const;
        DynamicList<Foam::label> parentList(const label cellI) const;
        Foam::label findInVis(const label splitIndex) const;

        //- Store splitting of cell into 4
        void storeSplit
        (
            const label cellI,
            const labelList& addedCells
        );

        //- Store combining 4 cells into master
        void combineCells
        (
            const label masterCellI,
            const labelList& combinedCells
        );


        //- Update numbering for mesh changes
        void updateMesh(const mapPolyMesh&);

        //- Update numbering for subsetting
        void subset
        (
            const labelList& pointMap,
            const labelList& faceMap,
            const labelList& cellMap
        );

        //- Update local numbering for mesh redistribution.
        //  Can only distribute clusters sent across in one go; cannot
        //  handle parts recombined in multiple passes.
        void distribute(const mapDistributePolyMesh&);


        //- Compact splitCells_. Removes all freeSplitCells_ elements.
        void compact();

        //- Extend/shrink storage. additional visibleCells_ elements get
        //  set to -1.
        void resize(const label nCells);

        //- Debug write
        void writeDebug() const;


        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, refinementTree&);
        friend Ostream& operator<<(Ostream&, const refinementTree&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
