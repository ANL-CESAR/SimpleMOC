#include"SimpleMOC_header.h"

Source * initialize_sources( Input I )
{
	// Allocate space
	Source * sources = (Source *) malloc( I.num_source_regions_per_assembly * sizeof(Source) );

	// determine number of cross section regions
	long n_xs_regions = I.num_source_regions_per_assembly / I.n_azimuthal;

	// Allocate scattering matrix matrix ptrs
	double *** matrices = (double ***) malloc( n_xs_regions * sizeof(double**) );

	// Allocate space for ALL scattering matrix ptrs
	double ** matrix_ptrs = (double **) malloc( n_xs_regions * I.n_egroups * sizeof(double *));

	// Allocate space for ALL scattering data
	double * data = (double *) malloc( n_xs_regions * I.n_egroups * I.n_egroups * sizeof(double));

	// Stitch allocation ptrs together
	for( long i = 0; i < n_xs_regions; i++ )
		matrices[i] = &matrix_ptrs[I.n_egroups];

	for( long i = 0; i < n_xs_regions; i++ )
		for( long j = 0; j < I.n_egroups; j++ )
			matrices[i][j] = &data[i * I.n_egroups * I.n_egroups + j * n_egroups];

	// Assign to source regions
	
	for( long i = 0; i < I.num_source_regions_per_assembly; i++ )
	{
		sources[i].scattering_matrix = 
	}

}
