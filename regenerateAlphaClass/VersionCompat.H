#if OF_MAJOR == 2
	#if OF_MINOR < 3
		// Need to make lots of changes, old multiphase interface
		Pout<< "OF_MINOR is less than 3: " << OF_MINOR << nl
			<< "Making necessary changes to variable names etc" << endl;
		#define BOUNDARY_FIELDS_NAME "alpha1.boundaryFields"
		#define FIELD_NAME "alpha1"
	#endif
	
	#if OF_MINOR >= 3
		// 'New' style.
		#define BOUNDARY_FIELDS_NAME "alpha.water.boundaryFields"
		#define FIELD_NAME "alpha.water"
	#endif
#endif
