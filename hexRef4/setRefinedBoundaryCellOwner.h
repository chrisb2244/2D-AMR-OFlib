vector faceNormal = calcSingleFaceNormal(mesh_.points(), faceI);
int relevantDir = calcRelevantDirs(faceNormal);
if (debug) Pout<< "faceNormal = " << faceNormal << " and relevantDir = " << relevantDir << endl;

switch (relevantDir)
{
	case -1: // left edge
	own = cellAddedCells[parentOwnCell][2];
	break;
	
	case +1: // right edge
	own = cellAddedCells[parentOwnCell][3];
	break;
	
	case -2: // bottom edge
	own = cellAddedCells[parentOwnCell][1];
	break;
	
	case +2: // top edge
	own = cellAddedCells[parentOwnCell][3];
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
