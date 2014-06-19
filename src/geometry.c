#include"SimpleMOC_header.h"

// Finds the "localized" region ID of a 3D cartesian coordinate
RegionID get_region_id( double x, double y, double z, Input input, Reactor reactor )
{
	RegionID id;

	// find assembly
	int x_assembly_id = (int) x / reactor.assembly_width;
	int y_assembly_id = (int) y / reactor.assembly_width;
	id.assembly = x_assembly_id + input.x_assemblies * y_assembly_id;
	
	// find pin cell
	id.pin = 0;

	// calculate the (x,y) coordinates in the assembly
	double assembly_x = (double) x - x_assembly_id * reactor.assembly_width;
	double assembly_y = (double) y - y_assembly_id * reactor.assembly_width;

	// find pin cell in the assembly (assuming 17 x 17 assemblies)
    int x_pin_id = (int) assembly_x / reactor.pin_cell_width;
	int y_pin_id = (int) assembly_y / reactor.pin_cell_width;
	id.pin = x_pin_id + 17 * y_pin_id;

	// TODO: Add Z component (3D)

	// calculate the (x,y) coordinates in the pin
	double pin_x = (double) assembly_x - x_pin_id * reactor.pin_cell_width;
	pin_x -= reactor.pin_cell_width/2; // take x coordinate from center of pin
	double pin_y = (double) assembly_y - y_pin_id * reactor.pin_cell_width;
	pin_y -= reactor.pin_cell_width/2; // take y coordinate from center of pin
	double radius = sqrt( pin_x * pin_x + pin_y * pin_y );

	// find the ring inside pin cell
	int ring_id = 0;
    for(int i=0; i < reactor.n_radial_regions; i++)
	{
		if(radius > reactor.radii[i])
			ring_id++;
		else
			break;
	}

	// find the azimuthal region within the ring
	double theta = atan2(pin_y, pin_x);
	double pi = 3.14159265358979323846264338327950;
	double azimuthal_interval = 2 * pi / reactor.n_azimuthal_regions;
	int azimuthal_id = (int) theta / azimuthal_interval;
	
	// compute zone id
	id.zone = azimuthal_id + reactor.n_azimuthal_regions * ring_id;

	return id;
}

// Finds the serialized "global" index of a region 
long get_region_index( RegionID id, Input input, Reactor reactor )
{
	// calculate number of zones in a pin cell (add 1 for outside region)
	int zones_per_pin = (reactor.n_radial_regions + 1) * reactor.n_azimuthal_regions;

	// define number of pins per assembly (assuming 17 x 17 assemblies)
	int pins_per_assembly = 289;

	// TODO: Add Z component (3D)

	// calculate the index number and return
	return id.zone + zones_per_pin * ( id.pin + pins_per_assembly * id.assembly);
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
