#include"SimpleMOC_header.h"

// Finds the "localized" region ID of a 3D cartesian coordinate
RegionID get_region_id( double x, double y, double z, Input input, Reactor reactor )
{
	RegionID id;

	// find assembly
	id.assembly = 0;
	
	// find pin cell
	id.pin = 0;
	
	// find the zone inside pin cell
	id.zone = 0;

	return id;
}

// Finds the serialized "global" index of a region 
long get_region_index( RegionID id, Input input, Reactor reactor )
{
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
