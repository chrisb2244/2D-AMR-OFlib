2D-AMR-OFlib
============

Note: This library may require OpenFOAM version 2.3.0 (or greater?).
Errors are reported with OF-2.1.1.


Acknowledgements : 

	Vincent Rivola (vinz @ cfd-online forums) - for sending a 2D damBreak with obstacle case to test.
	James Carow - for presenting the same error as VR and alerting me to possible version dependency (Currently unresolved)


Library set based on dynamicRefineFvMesh and the libdynamicfvmesh.so for OpenFOAM to allow 2D AMR

The dynamicMeshDict in \<case\>/constant should be changed to read:

	dynamicFvMesh dynamicRefineFvMeshHexRef4;
	dynamicRefineFvMeshHexRef4Coeffs
	{
		...
		// Standard options here
		...
	}
	toggle = 0; // or 1, see toggle_Explanation file

Contains:

	dynamicRefineFvMeshHexRef4:
	Class based on the dynamicRefineFvMesh. Only a few changes - 
	some hardcoded numbers changed from 7 to 3 or 8 to 4, added
	new functions refineAtZero() and updateAtZero() to try and
	leverage regenerateAlphaClass in the solver. These functions are
	also why dynamicFvMesh.H is included locally rather than from
	the original OF sources.

	hexRef4:
	Main class. Implements mesh cutting, and is based on hexRef8.
	Several new functions, mostly named with a my~ prefix, eg
	mySplitSideFaces(). The bulk of the meaningful changes are within
	the setRefinement() function.
	Initially, the function is similar to the hexRef8 version, and
	creates lists of the cells to be refined, splitting faces and 
	edges as needed. Unlike for hexRef8, not all faces of cells to
	be refined are split - there is a test for orientation of the face
	using calcFaceNormalVector() and comparing it with the direction
	to not refine. Also changed, the cellMidPoints' points no longer
	support a cell.
	Only 3 new cells are added in the next section, down from 7.
	They are stored in cellAddedCells, just like in hexRef8. After
	points have been added, the cellLevel and pointLevels are updated.
	This is much earlier than the updates to history_ in hexRef8. The
	values of cellLevel and pointLevel are needed to properly match up
	face owners and neighbours in the current implementation, but a
	copy of the original (pre-refinement) cell levels and point levels
	are held in the variables oldCellLevel_ and oldPointLevel_. These
	presumably are a waste of space if not used, and so I will probably
	remove them once the class is tidier.
	The splitting of cells is still done in several sections, as in
	hexRef8. Initially, the chosen faces are split (those with
	faceMidPoint >= 0). The code for this section is largely the same
	as for hexRef8.
	The second section holds most of the changes, and is split into two
	main parts. In the first part, additional points are added to front
	and back (the sides being cut) faces which are not themselves refined.
	In the second part, the 'side faces', that is, those which are cut in
	one direction, but not along their length, are handled. The handling
	of these faces is done by mySplitSideFaces(..).
	The third section is (hopefully) never triggered, but is also
	unchanged from hexRef8, apart from a diagnostic message at the bottom.
	The fourth section adds the internal faces between new cells. This is
	done by myCreateInternalFaces(..). myCreateInternalFaces() seems to be
	pretty solid, and so I will probably not look at it more until I have
	fixed other stuff - especially mySplitSideFaces, which needs some work
	for interfaces that become complicated.

	refinementTree:
	Class based on refinementHistory. Few changes, mostly the addition of
	some functions. myParentCell(..) is important in numerous places
	within hexRef4 during the splitting of cells (in determining owner 
	and neighbour). A parentList(..) function is also written, but that is
	mostly just to help with diagnostics when things don't work.

	regenerateAlphaClass:
	This class is new, and not a part of the 2D AMR. Really, it should
	probably be in a separate library, but the possibility of linking
	it more closely might be something worth looking at later.
	The class takes some input parameters (yMid should not be one of them,
	but it's currently simpler!) and then recalculates a cosine-based 
	perturbation to alpha. This is then remapped and returned to the
	solver. The aim is to allow for refinement and remapping of alpha
	before the start of the time loop, so that the mesh can be usefully
	refined at the start of solving, and since the seed for 'random'ness
	is passed by the solver, the alpha perturbation will have the same
	shape each time, so long as it is given the same seed. (ie. don't 
	generate a new seed each time you initialise an object of 
	regenerateAlphaClass type within the solver).

	hexRef4/DIAG_List:
	Used to store pointers to each of the many DIAG_XYZ booleans being
	held in hexRef4 now, which allow easy customization of the output
	to log files of diagnostic/info messages. These can be (currently) 
	changed using the format 
	*((hexRef4object).diagList.myMap["nameOfVar"]) = true/false;
	This is likely to be changed to a function, to make it simpler to 
	access and allow some automated checking of the validity of the 
	name passed.

	sortFaces:
	Used in hexRef4 to sort and print lists of faces to be refined, along
	with a rudimentary description of why faces are set to be refined.
	Will be removed when everything works well (at least for dependencies)
	Might leave the files in just so that they can be used if wanted.

ToDo:

	Test and improve for the damBreak 2D case.
	Determine reasoning for errors in pairedFaces for part 3 face
	splitting.
	Implement unrefinement.
	Check through flux corrections for updateAtZero/refineAtZero. Not all 
	fields 	exist at that point, so relatively few need correction, but at
	least some might (read - alpha).
	
