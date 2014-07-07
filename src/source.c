#include"SimpleMOC_header.h"

Source * initialize_sources( Input I )
{
	// Allocate space
	Source * sources = (Source *) malloc( I.n_source_regions_per_node * sizeof(Source) );

	// determine number of cross section regions
	long n_xs_regions = I.n_source_regions_per_node / 8;
	
	// Allocate scattering matrix matrix ptrs
	double *** s_matrices = (double ***) malloc( n_xs_regions * sizeof(double**) );

	// Allocate space for ALL scattering matrix ptrs
	double ** s_matrix_ptrs = (double **) malloc( n_xs_regions * I.n_egroups * sizeof(double *));

	// Allocate space for ALL scattering data
	double * s_matrix_data = (double *) malloc( n_xs_regions * I.n_egroups * I.n_egroups * 
			sizeof(double));

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
	 * [ ..., Total_XS, nu*SigmaF, Chi, ...]
	 */


	// Allocate space for XS ptrs (by region)
	double *** XS = (double ***) malloc( n_xs_regions * sizeof(double **) );

	// Allocate space for each XS type of interest (total, nu*SigmaF, and chi)
	double ** XS_arrays = (double **) malloc (n_xs_regions * I.n_egroups * sizeof(double *) );

	// Allocate space for total XS data
	double * XS_data = (double *) malloc( n_xs_regions * I.n_egroups * 3 * sizeof(double) );

	// stitch allocation ptrs together for XS data
	for( long i = 0; i < n_xs_regions; i++)
		XS[i] = &XS_arrays[i * n_xs_regions * I.n_egroups];

	for( long i = 0; i < n_xs_regions; i++)
		for(long j = 0; j < I.n_egroups; j++)
			XS[i][j] = &XS_data[i * I.n_egroups * 3 + j * 3];

	// Initialize XS data
	for( long i = 0; i < n_xs_regions; i++)
		for( int j = 0; j < I.n_egroups; j++)
			for( int k = 0; k < 3; k++)
				XS[i][j][k] = urand();

	// Allocate & Initialize Fluxes
	double * Flux = (double *) malloc( I.n_source_regions_per_node * I.n_egroups 
			* sizeof(double));
	for( long i = 0; i < I.n_source_regions_per_node * I.n_egroups; i++ )
		Flux[i] = 1.0;

	// Allocate & Inititalize "source" values
	double * Source = (double *) malloc( I.n_source_regions_per_node * I.n_egroups
			* sizeof(double));
	for( long i = 0; i < I.n_source_regions_per_node; i++ )
		Source[i] = urand();
	
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
		sources[i].flux = &Flux[i * I.n_egroups];
		sources[i].source = &Source[i * I.n_egroups]; 
	}

	// initialize FSR volume
	sources[i].vol = urand();

	// free memory of temporary variables
	free( s_matrices );
	free( XS );

	return sources;
}

void free_sources( Input I, Source * sources )
{
	// Free XS's
	free( sources[0].XS );
	// Free Flux's
	free( sources[0].flux );
	// Free scattering matrices
	free( sources[0].scattering_matrix[0] );
	free( sources[0].scattering_matrix );
	// Free source values
	free( sources[0].source );
}
