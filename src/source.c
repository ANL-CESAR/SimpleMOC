#include"SimpleMOC_header.h"

Source * initialize_sources( Input I, size_t * nbytes )
{
	// Allocate space
	Source * sources = (Source *) malloc( I.n_source_regions_per_node * sizeof(Source) );
	*nbytes += I.n_source_regions_per_node * sizeof(Source);

	// determine number of cross section and coarse axial regions
	long n_xs_regions = I.n_source_regions_per_node / 8;
	
	// Allocate scattering matrix matrix ptrs
	double *** s_matrices = (double ***) malloc( n_xs_regions * sizeof(double**) );
	*nbytes += n_xs_regions * sizeof(double **);

	// Allocate space for ALL scattering matrix ptrs
	double ** s_matrix_ptrs = (double **) malloc( n_xs_regions * I.n_egroups * sizeof(double *));
	*nbytes += n_xs_regions * sizeof(double **);

	// Allocate space for ALL scattering data
	if(I.mype==0) printf("Scattering data requires %zu MB of data...\n", n_xs_regions * I.n_egroups * I.n_egroups * sizeof(double) / 1024 / 1024);
	double * s_matrix_data = (double *) malloc( n_xs_regions * I.n_egroups * I.n_egroups * 
			sizeof(double));
	*nbytes += n_xs_regions * I.n_egroups * I.n_egroups * sizeof(double);

	// Stitch allocation ptrs together
	for( long i = 0; i < n_xs_regions; i++ )
		s_matrices[i] = &s_matrix_ptrs[i * I.n_egroups];

	for( long i = 0; i < n_xs_regions; i++ )
		for( long j = 0; j < I.n_egroups; j++ )
			s_matrices[i][j] = &s_matrix_data[i * I.n_egroups * I.n_egroups + j * I.n_egroups];

	// Iniitalize Scattering Matrix Values
	for( long i = 0; i < n_xs_regions; i++ )
		for( long j = 0; j < I.n_egroups; j++ )
			for( long k = 0; k < I.n_egroups; k++ )
				s_matrices[i][j][k] = urand();

	/*
	 * Create data scrtucture for storing XS data (and chi) as follows:
	 * An array is created which stores in contigious memory as
	 * [ ..., Total_XS, nu*SigmaF, SigmaA, Chi, ...]
	 */
	
	if(I.mype==0) printf("Beginning XS Allocation...\n");

	// Allocate space for XS ptrs (by region)
	double *** XS = (double ***) malloc( n_xs_regions * sizeof(double **) );
	*nbytes += n_xs_regions * sizeof(double **);

	// Allocate space for each XS type of interest (total, nu*SigmaF, and chi)
	double ** XS_arrays = (double **) malloc (n_xs_regions * I.n_egroups * sizeof(double *) );
	*nbytes += n_xs_regions * I.n_egroups * sizeof(double *);

	// Allocate space for total XS data
	double * XS_data = (double *) malloc( n_xs_regions * I.n_egroups * 4 * sizeof(double) );
	*nbytes += n_xs_regions * I.n_egroups * 4 * sizeof(double);

	// stitch allocation ptrs together for XS data
	for( long i = 0; i < n_xs_regions; i++)
		XS[i] = &XS_arrays[i * I.n_egroups];

	for( long i = 0; i < n_xs_regions; i++)
		for(long j = 0; j < I.n_egroups; j++)
			XS[i][j] = &XS_data[i * I.n_egroups * 4 + j * 4];

	// Initialize XS data
	for( long i = 0; i < n_xs_regions; i++)
		for( int j = 0; j < I.n_egroups; j++)
			for( int k = 0; k < 4; k++)
				XS[i][j][k] = urand();

	/////////////////////////////////////////////////////////////////////////////////

	if(I.mype==0) printf("Beginning Source Parameter Allocation...\n");

	// Allocate space for source parameters (quadratic axially)
	double *** fineSource = (double ***) malloc( I.n_source_regions_per_node
		   	* sizeof(double **) );
	*nbytes += I.n_source_regions_per_node * sizeof(double **);

	// Allocate space for array pointers to parameters
	double ** fineSourcePtrs = (double **) malloc ( I.n_source_regions_per_node
		   	* I.fai * sizeof(double *) );
	*nbytes += I.n_source_regions_per_node * I.fai * sizeof(double *);

	// Allocate space for parameter data
	double * fineSourceData = (double *) malloc( I.n_source_regions_per_node
		   	* I.fai * I.n_egroups * sizeof(double) );
	*nbytes += I.n_source_regions_per_node * I.fai * I.n_egroups * sizeof(double);

	// stitch allocation ptrs together for source parameter data
	for( long i = 0; i < I.n_source_regions_per_node; i++)
		fineSource[i] = &fineSourcePtrs[i * I.fai];

	for( long i = 0; i < I.n_source_regions_per_node; i++)
		for(long j = 0; j < I.fai; j++)
			fineSource[i][j] = &fineSourceData[i * I.fai * I.n_egroups
			   	+ j * I.n_egroups];

	// Initialize source parameters
	for( long i = 0; i < I.n_source_regions_per_node; i++)
		for( int j = 0; j < I.fai; j++)
			for( int k = 0; k < I.n_egroups; k++)
				fineSource[i][j][k] = urand();

	////////////////////////////////////////////////////////////////////////////////////

	if(I.mype==0) printf("Beginning Fine Flux Allocation...\n");

	// Allocate space for source parameters (quadratic axially)
	double *** fineFlux = (double ***) malloc( I.n_source_regions_per_node
		   	* sizeof(double **) );
	*nbytes += I.n_source_regions_per_node * sizeof(double **);

	// Allocate space for array pointers to parameters
	double ** fineFluxPtrs = (double **) malloc ( I.n_source_regions_per_node
		   	* I.fai * sizeof(double *) );
	*nbytes += I.n_source_regions_per_node * I.fai * sizeof(double *);

	// Allocate space for parameter data
	double * fineFluxData = (double *) malloc( I.n_source_regions_per_node
		   	* I.fai * I.n_egroups * sizeof(double) );
	*nbytes += I.n_source_regions_per_node * I.fai * I.n_egroups * sizeof(double);

	// stitch allocation ptrs together for source parameter data
	for( long i = 0; i < I.n_source_regions_per_node; i++)
		fineFlux[i] = &fineFluxPtrs[i * I.fai];

	for( long i = 0; i < I.n_source_regions_per_node; i++)
		for(int j = 0; j < I.fai; j++)
			fineFlux[i][j] = &fineFluxData[i * I.n_egroups * I.fai + j
			   	* I.n_egroups];

	// Initialize source parameters
	for( long i = 0; i < I.n_source_regions_per_node; i++)
		for( int j = 0; j < I.fai; j++)
			for( int k = 0; k < I.n_egroups; k++)
				fineFlux[i][j][k] = 0;

	////////////////////////////////////////////////////////////////////////////////////

	// Assign to source regions
	for( long i = 0; i < I.n_source_regions_per_node; i++ )
	{
		long idx;
		if( i == 0 )
			idx = 0;
		else
			idx = rand() % n_xs_regions;

		sources[i].scattering_matrix = s_matrices[idx];
		sources[i].XS = XS[idx];
		sources[i].fine_flux = fineFlux[i];
		sources[i].fine_source = fineSource[i]; 

		// initialize FSR volume
		sources[i].vol = urand();
	}

	// free memory of temporary variables
	free( s_matrices );
	free( XS );
	free( fineFlux );
	free( fineSource);

	return sources;
}

void free_sources( Input I, Source * sources )
{
	// Free XS's
	free( sources[0].XS[0] );
	free( sources[0].XS );
	// Free Flux's
	free( sources[0].fine_flux[0] );
	free( sources[0].fine_flux );
	// Free scattering matrices
	free( sources[0].scattering_matrix[0] );
	free( sources[0].scattering_matrix );
	// Free source values
	free( sources[0].fine_source[0] );
	free( sources[0].fine_source );
}
