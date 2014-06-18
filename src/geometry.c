#include"SimpleMOC_header.h"

long get_domain_id( double x, double y, double z, Input input, Reactor reactor )
{
	// find x-assembly
	
	return 0;
}

// Calculates constant volume radii
double * determine_radii( Reactor reactor )
{
	// allocate space
	double * radii = (double *) malloc( reactor.n_radial_regions * sizeof(double) ); 

	// Calculate constant volume radii for all radial regions
	for( int i = 0; i < reactor.n_radial_regions; i++ )
		radii[i] = sqrt( reactor.pin_radius * reactor.pin_radius * i / reactor.n_radial_regions ); 

	return radii;
}

// Sets the parameters for the reactor and set the radii
void reactor_init( Reactor * reactor )
{
	reactor->pin_cell_width = 1.26;
	reactor->pin_radius = 0.46;
	reactor->assembly_width = 1.26 * 17;
	reactor->n_radial_regions = 10;
	reactor->n_azimuthal_regions = 8;
	reactor->radii = determine_radii( *reactor );
}
