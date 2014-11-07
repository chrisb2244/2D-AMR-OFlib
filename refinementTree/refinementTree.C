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

#include "DynamicList.H"
#include "refinementTree.H"
#include "ListOps.H"
#include "mapPolyMesh.H"
#include "mapDistributePolyMesh.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(refinementTree, 0);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::refinementTree::writeEntry
(
    const List<splitCell4>& splitCells,
    const splitCell4& split
)
{
    // Write me:
    if (split.addedCellsPtr_.valid())
    {
        Pout<< "parent:" << split.parent_
            << " subCells:" << split.addedCellsPtr_()
            << endl;
    }
    else
    {
        Pout<< "parent:" << split.parent_
            << " no subcells"
            << endl;
    }

    if (split.parent_ >= 0)
    {
        Pout<< "parent data:" << endl;
        // Write my parent
        string oldPrefix = Pout.prefix();
        Pout.prefix() = "  " + oldPrefix;
        writeEntry(splitCells, splitCells[split.parent_]);
        Pout.prefix() = oldPrefix;
    }
}


void Foam::refinementTree::writeDebug
(
    const labelList& visibleCells,
    const List<splitCell4>& splitCells
)
{
    string oldPrefix = Pout.prefix();
    Pout.prefix() = "";

    forAll(visibleCells, cellI)
    {
        label index = visibleCells[cellI];

        if (index >= 0)
        {
            Pout<< "Cell from refinement:" << cellI << " index:" << index
                << endl;

            string oldPrefix = Pout.prefix();
            Pout.prefix() = "  " + oldPrefix;
            writeEntry(splitCells, splitCells[index]);
            Pout.prefix() = oldPrefix;
        }
        else
        {
            Pout<< "Unrefined cell:" << cellI << " index:" << index << endl;
        }
    }
    Pout.prefix() = oldPrefix;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
Foam::refinementTree::splitCell4::splitCell4()
:
    parent_(-1),
    addedCellsPtr_(NULL)
{}


//- Construct as child element of parent
Foam::refinementTree::splitCell4::splitCell4(const label parent)
:
    parent_(parent),
    addedCellsPtr_(NULL)
{}


//- Construct from Istream
Foam::refinementTree::splitCell4::splitCell4(Istream& is)
{
    is >> *this;
}


//- Construct as (deep) copy.
Foam::refinementTree::splitCell4::splitCell4(const splitCell4& sc)
:
    parent_(sc.parent_),
    addedCellsPtr_
    (
        sc.addedCellsPtr_.valid()
      ? new FixedList<label, 4>(sc.addedCellsPtr_())
      : NULL
    )
{}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, refinementTree::splitCell4& sc)
{
    labelList addedCells;

    is >> sc.parent_ >> addedCells;

    if (addedCells.size())
    {
        sc.addedCellsPtr_.reset(new FixedList<label, 4>(addedCells));
    }
    else
    {
        sc.addedCellsPtr_.reset(NULL);
    }

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const refinementTree::splitCell4& sc
)
{
    // Output as labelList so we can have 0 sized lists. Alternative is to
    // output as fixedlist with e.g. -1 elements and check for this upon
    // reading. However would cause much more data to be transferred.

    if (sc.addedCellsPtr_.valid())
    {
        return os
            << sc.parent_
            << token::SPACE
            << labelList(sc.addedCellsPtr_());
    }
    else
    {
        return os << sc.parent_ << token::SPACE << labelList(0);
    }
}


void Foam::refinementTree::checkIndices() const
{
    // Check indices.
    forAll(visibleCells_, i)
    {
        if (visibleCells_[i] < 0 && visibleCells_[i] >= splitCells_.size())
        {
            FatalErrorIn("refinementTree::checkIndices() const")
                << "Illegal entry " << visibleCells_[i]
                << " in visibleCells at location" << i << nl
                << "It points outside the range of splitCells : 0.."
                << splitCells_.size()-1
                << abort(FatalError);
        }
    }
}


Foam::label Foam::refinementTree::allocateSplitCell
(
    const label parent,
    const label i
)
{
    label index = -1;

    if (freeSplitCells_.size())
    {
		//Pout<< "freeSplitCells_.size() = " << freeSplitCells_.size() << endl;
        index = freeSplitCells_.remove();

        splitCells_[index] = splitCell4(parent);
    }
    else
    {
		//Pout<< "splitCells_.size() = " << splitCells_.size() << endl;
        index = splitCells_.size();

        splitCells_.append(splitCell4(parent));
    }


    // Update the parent field
    if (parent >= 0)
    {
        splitCell4& parentSplit = splitCells_[parent];

        if (parentSplit.addedCellsPtr_.empty())
        {
            // Allocate storage on parent for the 4 subcells.
            parentSplit.addedCellsPtr_.reset(new FixedList<label, 4>(-1));
        }


        // Store me on my parent
        FixedList<label, 4>& parentSplits = parentSplit.addedCellsPtr_();

        parentSplits[i] = index;
    }

    return index;
}


void Foam::refinementTree::freeSplitCell(const label index)
{
	Pout<< "unlikely?" << endl;
    splitCell4& split = splitCells_[index];

    // Make sure parent does not point to me anymore.
    if (split.parent_ >= 0)
    {
        autoPtr<FixedList<label, 4> >& subCellsPtr =
            splitCells_[split.parent_].addedCellsPtr_;

        if (subCellsPtr.valid())
        {
            FixedList<label, 4>& subCells = subCellsPtr();

            label myPos = findIndex(subCells, index);

            if (myPos == -1)
            {
                FatalErrorIn("refinementTree::freeSplitCell")
                    << "Problem: cannot find myself in"
                    << " parents' children" << abort(FatalError);
            }
            else
            {
                subCells[myPos] = -1;
            }
        }
    }

    // Mark splitCell as free
    split.parent_ = -2;

    // Add to cache of free splitCells
    freeSplitCells_.append(index);
}


// Mark entry in splitCells. Recursively mark its parent and subs.
void Foam::refinementTree::markSplit
(
    const label index,
    labelList& oldToNew,
    DynamicList<splitCell4>& newSplitCells
) const
{
    if (oldToNew[index] == -1)
    {
        // Not yet compacted.

        const splitCell4& split = splitCells_[index];

        oldToNew[index] = newSplitCells.size();
        //Pout<< "split = " << split << endl;
        newSplitCells.append(split);

        if (split.parent_ >= 0)
        {
            markSplit(split.parent_, oldToNew, newSplitCells);
        }
        if (split.addedCellsPtr_.valid())
        {
            const FixedList<label, 4>& splits = split.addedCellsPtr_();

            forAll(splits, i)
            {
                if (splits[i] >= 0)
                {
                    markSplit(splits[i], oldToNew, newSplitCells);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementTree::refinementTree(const IOobject& io)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn
        (
            "refinementTree::refinementTree(const IOobject&)"
        )   << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
            << " does not support automatic rereading."
            << endl;
    }

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }

    if (debug)
    {
        Pout<< "refinementTree::refinementTree :"
            << " constructed history from IOobject :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


//- Read or construct
Foam::refinementTree::refinementTree
(
    const IOobject& io,
    const List<splitCell4>& splitCells,
    const labelList& visibleCells
)
:
    regIOobject(io),
    splitCells_(splitCells),
    freeSplitCells_(0),
    visibleCells_(visibleCells)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn
        (
            "refinementTree::refinementTree"
            "(const IOobject&, const List<splitCell4>&, const labelList&)"
        )   << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
            << " does not support automatic rereading."
            << endl;
    }

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }

    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementTree::refinementTree :"
            << " constructed history from IOobject or components :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// Construct from initial number of cells (all visible)
Foam::refinementTree::refinementTree
(
    const IOobject& io,
    const label nCells
)
:
    regIOobject(io),
    freeSplitCells_(0)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn
        (
            "refinementTree::refinementTree"
            "(const IOobject&, const label)"
        )   << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
            << " does not support automatic rereading."
            << endl;
    }

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
    else
    {
		if (debug)
		{
			Pout<< "Constructed from nCells, not IOobject." << endl;
		}
        visibleCells_.setSize(nCells);
        splitCells_.setCapacity(nCells);

        for (label cellI = 0; cellI < nCells; cellI++)
        {
            visibleCells_[cellI] = cellI;
            splitCells_.append(splitCell4());
        }
    }

    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementTree::refinementTree :"
            << " constructed history from IOobject or initial size :"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// Construct as copy
Foam::refinementTree::refinementTree
(
    const IOobject& io,
    const refinementTree& rh
)
:
    regIOobject(io),
    splitCells_(rh.splitCells()),
    freeSplitCells_(rh.freeSplitCells()),
    visibleCells_(rh.visibleCells())
{
    if (debug)
    {
        Pout<< "refinementTree::refinementTree : constructed initial"
            << " history." << endl;
    }
}


// Construct from Istream
Foam::refinementTree::refinementTree(const IOobject& io, Istream& is)
:
    regIOobject(io),
    splitCells_(is),
    freeSplitCells_(0),
    visibleCells_(is)
{
    // Check indices.
    checkIndices();

    if (debug)
    {
        Pout<< "refinementTree::refinementTree :"
            << " constructed history from Istream"
            << " splitCells:" << splitCells_.size()
            << " visibleCells:" << visibleCells_.size()
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::refinementTree::resize(const label size)
{
    label oldSize = visibleCells_.size();

    if (debug)
    {
        Pout<< "refinementTree::resize from " << oldSize << " to " << size
            << " cells" << endl;
    }

    visibleCells_.setSize(size);

    // Set additional elements to -1.
    for (label i = oldSize; i < visibleCells_.size(); i++)
    {
        visibleCells_[i] = -1;
    }
}


void Foam::refinementTree::updateMesh(const mapPolyMesh& map)
{
	Pout<< "updateMesh called" << endl;
    if (active())
    {
        const labelList& reverseCellMap = map.reverseCellMap();

        // Note that only the live cells need to be renumbered.

        labelList newVisibleCells(map.cellMap().size(), -1);

        forAll(visibleCells_, cellI)
        {
            if (visibleCells_[cellI] != -1)
            {
                label index = visibleCells_[cellI];

                // Check not already set
                if (splitCells_[index].addedCellsPtr_.valid())
                {
                    FatalErrorIn
                    (
                        "refinementTree::updateMesh(const mapPolyMesh&)"
                    )   << "Problem" << abort(FatalError);
                }

                label newCellI = reverseCellMap[cellI];

                if (newCellI >= 0)
                {
                    newVisibleCells[newCellI] = index;
                }
            }
        }

        if (debug)
        {
            Pout<< "refinementTree::updateMesh : from "
                << visibleCells_.size()
                << " to " << newVisibleCells.size()
                << " cells" << endl;
        }

        visibleCells_.transfer(newVisibleCells);
    }
}


// Update numbering for subsetting
void Foam::refinementTree::subset
(
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap
)
{
    if (active())
    {
        labelList newVisibleCells(cellMap.size(), -1);

        forAll(newVisibleCells, cellI)
        {
            label oldCellI = cellMap[cellI];

            label index = visibleCells_[oldCellI];

            // Check that cell is live (so its parent has no refinement)
            if (index >= 0 && splitCells_[index].addedCellsPtr_.valid())
            {
                FatalErrorIn
                (
                    "refinementTree::subset"
                    "(const labelList&, const labelList&, const labelList&)"
                )   << "Problem" << abort(FatalError);
            }

            newVisibleCells[cellI] = index;
        }

        if (debug)
        {
            Pout<< "refinementTree::updateMesh : from "
                << visibleCells_.size()
                << " to " << newVisibleCells.size()
                << " cells" << endl;
        }

        visibleCells_.transfer(newVisibleCells);
    }
}


void Foam::refinementTree::countProc
(
    const label index,
    const label newProcNo,
    labelList& splitCellProc,
    labelList& splitCellNum
) const
{
    if (splitCellProc[index] != newProcNo)
    {
        // Different destination processor from other cells using this
        // parent. Reset count.
        splitCellProc[index] = newProcNo;
        splitCellNum[index] = 1;
    }
    else
    {
        splitCellNum[index]++;

        // Increment parent if whole splitCell moves to same processor
        if (splitCellNum[index] == 4)
        {
            Pout<< "Moving " << splitCellNum[index]
                << " cells originating from cell " << index
                << " from processor " << Pstream::myProcNo()
                << " to processor " << splitCellProc[index]
                << endl;

            label parent = splitCells_[index].parent_;

            if (parent >= 0)
            {
                string oldPrefix = Pout.prefix();
                Pout.prefix() = "  " + oldPrefix;

                countProc(parent, newProcNo, splitCellProc, splitCellNum);

                Pout.prefix() = oldPrefix;
            }
        }
    }
}


void Foam::refinementTree::distribute(const mapDistributePolyMesh& map)
{
    if (!active())
    {
        FatalErrorIn
        (
            "refinementTree::distribute(const mapDistributePolyMesh&)"
        )   << "Calling distribute on inactive history" << abort(FatalError);
    }


    if (!Pstream::parRun())
    {
        return;
    }

    // Remove unreferenced history.
    Pout<< "Call in refinementTree.C to compact()" << endl;
    compact();

    //Pout<< nl << "--BEFORE:" << endl;
    //writeDebug();
    //Pout<< "---------" << nl << endl;


    // Distribution is only partially functional.
    // If all 4 cells resulting from a single parent are sent across in one
    // go it will also send across that part of the refinement history.
    // If however e.g. first 1 and then the other 7 are sent across the
    // history will not be reconstructed.

    // Determine clusters. This is per every entry in splitCells_ (that is
    // a parent of some refinement) a label giving the processor it goes to
    // if all its children are going to the same processor.

    // Per visible cell the processor it goes to.
    labelList destination(visibleCells_.size());

    const labelListList& subCellMap = map.cellMap().subMap();

    forAll(subCellMap, procI)
    {
        const labelList& newToOld = subCellMap[procI];

        forAll(newToOld, i)
        {
            label oldCellI = newToOld[i];

            destination[oldCellI] = procI;
        }
    }

//Pout<< "refinementTree::distribute :"
//    << " destination:" << destination << endl;

    // Per splitCell entry the processor it moves to
    labelList splitCellProc(splitCells_.size(), -1);
    // Per splitCell entry the number of live cells that move to that processor
    labelList splitCellNum(splitCells_.size(), 0);

    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            countProc
            (
                splitCells_[index].parent_,
                destination[cellI],
                splitCellProc,
                splitCellNum
            );
        }
    }

    //Pout<< "refinementTree::distribute :"
    //    << " splitCellProc:" << splitCellProc << endl;
    //
    //Pout<< "refinementTree::distribute :"
    //    << " splitCellNum:" << splitCellNum << endl;


    // Create subsetted refinement tree consisting of all parents that
    // move in their whole to other processor.
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        //Pout<< "-- Subetting for processor " << procI << endl;

        // From uncompacted to compacted splitCells.
        labelList oldToNew(splitCells_.size(), -1);

        // Compacted splitCells. Similar to subset routine below.
        Pout<< "distribute, defo unlikely" << endl;
        DynamicList<splitCell4> newSplitCells(splitCells_.size());

        // Loop over all entries. Note: could recurse like countProc so only
        // visit used entries but is probably not worth it.

        forAll(splitCells_, index)
        {
//            Pout<< "oldCell:" << index
//                << " proc:" << splitCellProc[index]
//                << " nCells:" << splitCellNum[index]
//                << endl;

// -----------------------------------------------------------------------------
            if (splitCellProc[index] == procI && splitCellNum[index] == 4)
            {
                // Entry moves in its whole to procI
                oldToNew[index] = newSplitCells.size();
                newSplitCells.append(splitCells_[index]);

                //Pout<< "Added oldCell " << index
                //    << " info " << newSplitCells.last()
                //    << " at position " << newSplitCells.size()-1
                //    << endl;
            }
        }

        // Add live cells that are subsetted.
        forAll(visibleCells_, cellI)
        {
            label index = visibleCells_[cellI];

            if (index >= 0 && destination[cellI] == procI)
            {
                label parent = splitCells_[index].parent_;

                //Pout<< "Adding refined cell " << cellI
                //    << " since moves to "
                //    << procI << " old parent:" << parent << endl;

                // Create new splitCell with parent
                oldToNew[index] = newSplitCells.size();
                newSplitCells.append(splitCell4(parent));
            }
        }

        //forAll(oldToNew, index)
        //{
        //    Pout<< "old:" << index << " new:" << oldToNew[index]
        //        << endl;
        //}

        newSplitCells.shrink();

        // Renumber contents of newSplitCells
        forAll(newSplitCells, index)
        {
            splitCell4& split = newSplitCells[index];

            if (split.parent_ >= 0)
            {
                split.parent_ = oldToNew[split.parent_];
            }
            if (split.addedCellsPtr_.valid())
            {
                FixedList<label, 4>& splits = split.addedCellsPtr_();

                forAll(splits, i)
                {
                    if (splits[i] >= 0)
                    {
                        splits[i] = oldToNew[splits[i]];
                    }
                }
            }
        }


        const labelList& subMap = subCellMap[procI];

        // New visible cells.
        labelList newVisibleCells(subMap.size(), -1);

        forAll(subMap, newCellI)
        {
            label oldCellI = subMap[newCellI];

            label oldIndex = visibleCells_[oldCellI];

            if (oldIndex >= 0)
            {
                newVisibleCells[newCellI] = oldToNew[oldIndex];
            }
        }

        //Pout<< nl << "--Subset for domain:" << procI << endl;
        //writeDebug(newVisibleCells, newSplitCells);
        //Pout<< "---------" << nl << endl;


        // Send to neighbours
        OPstream toNbr(Pstream::blocking, procI);
        toNbr << newSplitCells << newVisibleCells;
    }


    // Receive from neighbours and merge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Remove all entries. Leave storage intact.
    splitCells_.clear();

    visibleCells_.setSize(map.mesh().nCells());
    visibleCells_ = -1;

    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        IPstream fromNbr(Pstream::blocking, procI);
        List<splitCell4> newSplitCells(fromNbr);
        labelList newVisibleCells(fromNbr);

        //Pout<< nl << "--Received from domain:" << procI << endl;
        //writeDebug(newVisibleCells, newSplitCells);
        //Pout<< "---------" << nl << endl;


        // newSplitCells contain indices only into newSplitCells so
        // renumbering can be done here.
        label offset = splitCells_.size();

        //Pout<< "**Renumbering data from proc " << procI << " with offset "
        //    << offset << endl;

        forAll(newSplitCells, index)
        {
            splitCell4& split = newSplitCells[index];

            if (split.parent_ >= 0)
            {
                split.parent_ += offset;
            }
            if (split.addedCellsPtr_.valid())
            {
                FixedList<label, 4>& splits = split.addedCellsPtr_();

                forAll(splits, i)
                {
                    if (splits[i] >= 0)
                    {
                        splits[i] += offset;
                    }
                }
            }

            splitCells_.append(split);
        }


        // Combine visibleCell.
        const labelList& constructMap = map.cellMap().constructMap()[procI];

        forAll(newVisibleCells, i)
        {
            visibleCells_[constructMap[i]] = newVisibleCells[i] + offset;
        }
    }
    splitCells_.shrink();

    //Pout<< nl << "--AFTER:" << endl;
    //writeDebug();
    //Pout<< "---------" << nl << endl;
}


void Foam::refinementTree::compact()
{
	Pout<< "compact called" << endl;
    if (debug)
    {
        Pout<< "refinementTree::compact() Entering with:"
            << " freeSplitCells_:" << freeSplitCells_.size()
            << " splitCells_:" << splitCells_.size()
            << " visibleCells_:" << visibleCells_.size()
            << endl;
            
        writeDebug(visibleCells_, splitCells_);

        // Check all free splitCells are marked as such
        forAll(freeSplitCells_, i)
        {
            label index = freeSplitCells_[i];

            if (splitCells_[index].parent_ != -2)
            {
                FatalErrorIn("refinementTree::compact()")
                    << "Problem index:" << index
                    << abort(FatalError);
            }
        }

        // Check none of the visible cells are marked as free
        forAll(visibleCells_, cellI)
        {
            if
            (
                visibleCells_[cellI] >= 0
             && splitCells_[visibleCells_[cellI]].parent_ == -2
            )
            {
                FatalErrorIn("refinementTree::compact()")
                    << "Problem : visible cell:" << cellI
                    << " is marked as being free." << abort(FatalError);
            }
        }
    }
	
	
    DynamicList<splitCell4> newSplitCells(splitCells_.size());

    // From uncompacted to compacted splitCells.
    labelList oldToNew(splitCells_.size(), -1);

    // Mark all used splitCell entries. These are either indexed by visibleCells
    // or indexed from other splitCell entries.

    // Mark from visibleCells
    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            // Make sure we only mark visible indices if they either have a
            // parent or subsplits.
            if
            (
                splitCells_[index].parent_ != -1
                    // catches cells that DO have a parent, ie they are refined
             || splitCells_[index].addedCellsPtr_.valid()
                    // catches cells that have children, ie they are refined
            )
            {
                markSplit(index, oldToNew, newSplitCells);
            }
        }
    }

    // Mark from splitCells
    forAll(splitCells_, index)
    {
        if (splitCells_[index].parent_ == -2)
        {
            // freed cell.
        }
        else if
        (
            splitCells_[index].parent_ == -1
                // is a cellLevel 0 cell, it has no parent
         && splitCells_[index].addedCellsPtr_.empty()
                // cell has no children, it is unrefined
        )
        {
            // recombined cell. No need to keep since no parent and no subsplits
            // Note that gets marked if reachable from other index!
        }
        else
        {
            // Is used element.
            markSplit(index, oldToNew, newSplitCells);
        }
    }
    
    
    //Pout<< "Complete list of newSplitCells, followed by oldToNew: " << nl
    //<< newSplitCells << nl << nl << nl << oldToNew << nl << nl<< endl;


    // Now oldToNew is fully complete and compacted elements are in
    // newSplitCells.
    // Renumber contents of newSplitCells and visibleCells.
    forAll(newSplitCells, index)
    {

        splitCell4& split = newSplitCells[index];

        if (split.parent_ >= 0)
        {
            split.parent_ = oldToNew[split.parent_];
        }
        if (split.addedCellsPtr_.valid())
        {
            FixedList<label, 4>& splits = split.addedCellsPtr_();

            forAll(splits, i)
            {
                if (splits[i] >= 0)
                {
                    splits[i] = oldToNew[splits[i]];
                }
            }
        }
    }


    if (debug)
    {
        Pout<< "refinementTree::compact : compacted splitCells from "
            << splitCells_.size() << " to " << newSplitCells.size() << endl;
    }
    
    //Pout<< "Final list of splitCells (newSplitCells): " << nl
    //<< newSplitCells << nl << nl << nl << endl;

    splitCells_.transfer(newSplitCells);
    freeSplitCells_.clearStorage();


    if (debug)
    {
        Pout<< "refinementTree::compact() NOW:"
            << " freeSplitCells_:" << freeSplitCells_.size()
            << " splitCells_:" << splitCells_.size()
            << " newSplitCells:" << newSplitCells.size()
            << " visibleCells_:" << visibleCells_.size()
            << endl;
    }


    // Adapt indices in visibleCells_
    forAll(visibleCells_, cellI)
    {
        label index = visibleCells_[cellI];

        if (index >= 0)
        {
            // Note that oldToNew can be -1 so it resets newVisibleCells.
            visibleCells_[cellI] = oldToNew[index];
        }
        else
        {
            // Keep -1 value.
        }
    }
}
Foam::label Foam::refinementTree::myParentCell(const label cellI) const
{
	if (debug) Pout<< "cellI = " << cellI << endl;
	label visIndex = visibleCells_[cellI];
    if (visIndex < 0) // visibleCells_ gives -1 for the root,
                      // base level parent cells (initial mesh cells)
	{
		if (visIndex != -1)
		{
			// Unexpected...
			Pout<< "Warning: visIndex = " << visIndex << endl;
		}
		if (debug) Pout<< "visIndex = -1" << endl;
		return cellI;
	}
	else
	{
		label splitIndex = splitCells_[visIndex].parent_;
        if (debug) {
            Pout<< "splitCells_[splitIndex] = splitCells_["
                << splitIndex << "] = "
                << splitCells_[splitIndex]
                << endl;
        }
        FixedList<label, 4> splitList =
                splitCells_[splitIndex].addedCellsPtr_();
		forAll (visibleCells_, i)
		{
			if (visibleCells_[i] == splitList[0])
			{
				if (debug) Pout<< "parentCell (i) = " << i << endl;
				return i;
			}
		}
		DynamicList<label> parentListI(parentList(cellI));
		if (parentListI.size() > 0) return parentListI[1];
		FatalErrorIn("myParentCell(..)")
			<< "Couldn't find splitList[0] = " << splitList[0]
			<< " in visibleCells_"
			<< nl << nl
			//~ << visibleCells_
			//~ << nl << nl
			<< "parentList = " << parentList(cellI)
			<< abort(FatalError);
	}
	return -1;
}

Foam::label Foam::refinementTree::findInVis(const label splitIndex) const
{
	forAll(visibleCells_, i)
	{
		if (visibleCells_[i] == splitIndex)
		{
			return i;
		}
	}
	FatalErrorIn("findInVis(..)")
		<< "Could not find the index " << splitIndex
		<< " in visibleCells_"
		<< abort(FatalError);
	return -1;
}


// Returns a list starting with the cell given and moving up through its parents.
// The same label can appear repeatedly if the cell was refined from itself 
// (ie it is the lower left corner of a refined cell)
Foam::DynamicList<Foam::label> Foam::refinementTree::parentList(
        const label cellI) const
{
	DynamicList<Foam::label> parentList;
	label X = visibleCells_[cellI];
	parentList.append(cellI);
	if (X < 0)
	{
		return parentList;
	}
	else
	{
		label Y = splitCells_[X].parent_;
		while (splitCells_[Y].parent_ != -1)
		{
            FixedList<label, 4> splitList = splitCells_[Y].addedCellsPtr_();
            // This is my (A B C D)
			parentList.append(findInVis(splitList[0]));
			Y = splitCells_[Y].parent_;
		}
		FixedList<label, 4> splitList = splitCells_[Y].addedCellsPtr_();
		if (splitCells_[splitList[0]].parent_ == Y)
		{
			// search for more (another?) parent
			if (!splitCells_[splitList[0]].addedCellsPtr_.valid())
			{
				parentList.append(findInVis(splitList[0]));
				return parentList;
			}
			splitList = splitCells_[splitList[0]].addedCellsPtr_();
			parentList.append(findInVis(splitList[0]));
		}
		else
		{
			parentList.append(findInVis(splitList[0]));
		}
		return parentList;
	}
	FatalErrorIn("parentList(..)")
        << "Reached the end of the function without returning "
        << "a list of parent cells"
		<< abort(FatalError);
	return parentList;
}

	
	

void Foam::refinementTree::writeDebug() const
{
    writeDebug(visibleCells_, splitCells_);
}


void Foam::refinementTree::storeSplit
(
    const label cellI,
    const labelList& addedCells
)
{
    label parentIndex = -1;

    if (visibleCells_[cellI] != -1)
    {
        // Was already live. The current live cell becomes the
        // parent of the cells split off from it.

        parentIndex = visibleCells_[cellI];

        // It is no longer live (note that actually cellI gets alive
        // again below since is addedCells[0])
        visibleCells_[cellI] = -1;
    }
    else
    {
		Pout<< "Don't expect this to be called from storeSplit." << endl;
        Pout<< "The value of visibleCells_[cellI] for cellI = " << cellI
            << " is ";
		Pout<< visibleCells_[cellI] << endl;
        // Create 0th level. -1 parent to denote this.
        parentIndex = allocateSplitCell(-1, -1);
    }

    // Create live entries for added cells that point to the
    // cell they were created from (parentIndex)
    forAll(addedCells, i)
    {
        label addedCellI = addedCells[i];

        // Create entries for the split off cells. All of them
        // are visible.
        //Pout<< "allocateSplitCell called with args parentIndex: "
        //    << parentIndex << ", i: "<< i << endl;
        visibleCells_[addedCellI] = allocateSplitCell(parentIndex, i);
    }
}


void Foam::refinementTree::combineCells
(
    const label masterCellI,
    const labelList& combinedCells
)
{
    // Save the parent structure
    label parentIndex = splitCells_[visibleCells_[masterCellI]].parent_;

    // Remove the information for the combined cells
    forAll(combinedCells, i)
    {
        label cellI = combinedCells[i];

        freeSplitCell(visibleCells_[cellI]);
        visibleCells_[cellI] = -1;
    }

    splitCell4& parentSplit = splitCells_[parentIndex];
    parentSplit.addedCellsPtr_.reset(NULL);
    visibleCells_[masterCellI] = parentIndex;
}


bool Foam::refinementTree::readData(Istream& is)
{
    is >> *this;
    return !is.bad();
}


bool Foam::refinementTree::writeData(Ostream& os) const
{
    os << *this;

    return os.good();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, refinementTree& rh)
{
    rh.freeSplitCells_.clearStorage();

    is >> rh.splitCells_ >> rh.visibleCells_;

    // Check indices.
    rh.checkIndices();

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const refinementTree& rh)
{
	//const_cast<refinementTree&>(rh).compact();
    

    return os   << "// splitCells" << nl
                << rh.splitCells_ << nl
                << "// visibleCells" << nl
                << rh.visibleCells_;
}


// ************************************************************************* //
