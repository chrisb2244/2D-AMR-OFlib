#ifndef VersionCompat_DynRef_H
#define VersionCompat_DynRef_H

#if OF_MAJOR == 2
	#if OF_MINOR < 3
		// 'Old' functions
		Pout<< "OF_MINOR is less than 3: " << OF_MINOR << endl;
		bool Foam::dynamicRefineFvMeshHexRef4::topoChanging(const bool c)
		{
			return Foam::polyMesh::changing(c);
		}
		bool Foam::dynamicRefineFvMeshHexRef4::topoChanging()
		{
			return Foam::polyMesh::changing();
		}
	#endif
	
	#if OF_MINOR >= 3
		// 'New' style.
		bool Foam::dynamicRefineFvMeshHexRef4::topoChanging(const bool c)
		{
			return Foam::polyMesh::topoChanging(c);
		}
		bool Foam::dynamicRefineFvMeshHexRef4::topoChanging()
		{
			return Foam::polyMesh::topoChanging();
		}
	#endif
#endif

#if OF_MAJOR != 2
	// Do something here
#endif

#endif
