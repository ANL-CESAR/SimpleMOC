#include"SimpleMOC_header.h"

// Finds the "localized" region ID of a 3D cartesian coordinate
RegionID get_region_id( double x, double y, double z, Input input, Reactor reactor )
{
	RegionID id;

	// find assembly by x and y coordinates
	int x_assembly_id = (int) x / reactor.assembly_width;
	int y_assembly_id = (int) y / reactor.assembly_width;
	id.assembly = x_assembly_id + input.x_assemblies * y_assembly_id;
	
	// calculate the (x,y) coordinates in the assembly
	double assembly_x = (double) x - x_assembly_id * reactor.assembly_width;
	double assembly_y = (double) y - y_assembly_id * reactor.assembly_width;

	// TODO: Remove assembly geometry details

	// find pin cell in the assembly (assuming 17 x 17 assemblies)
    int x_pin_id = assembly_x / reactor.pin_cell_width;
	int y_pin_id = assembly_y / reactor.pin_cell_width;
	id.pin = x_pin_id + 17 * y_pin_id;

	// calculate the (x,y) coordinates in the pin
	double pin_x = (double) assembly_x - x_pin_id * reactor.pin_cell_width;
	pin_x -= reactor.pin_cell_width / 2; // take x coordinate from center of pin
	double pin_y = (double) assembly_y - y_pin_id * reactor.pin_cell_width;
	pin_y -= reactor.pin_cell_width / 2; // take y coordinate from center of pin
	double radius = sqrt( pin_x * pin_x + pin_y * pin_y );

	// find the ring inside pin cell
	int ring_id = 0;
    for(int i=1; i < reactor.n_radial_regions; i++)
	{
		if(radius >= reactor.radii[i])
			ring_id++;
		else
			break;
	}

	// find the azimuthal region within the ring
	double theta = atan2(pin_y, pin_x) + M_PI;
	double azimuthal_interval = 2 * M_PI / reactor.n_azimuthal_regions;
	int azimuthal_id = theta / azimuthal_interval;
	
	// compute 2D zone id
	int zone_2D = azimuthal_id + reactor.n_azimuthal_regions * ring_id;

	// compute axial layer
	int n_coarse_axial_layers = 400 / input.cai;
	int n_fine_axial_layers = n_coarse_axial_layers * input.fai;
 	int fine_axial_length = input.cai / input.fai;
	int fine_axial_layer = z / fine_axial_length;

	// compute 3D zone id
	id.zone = fine_axial_layer + zone_2D * n_fine_axial_layers;

	return id;
}

// Finds the serialized "global" index of a region 
long get_region_index( RegionID id, Input input, Reactor reactor )
{
	// calculate number of zones in a pin cell
	int n_coarse_axial_layers = 400 / input.cai;
	int zones_per_pin = reactor.n_radial_regions * reactor.n_azimuthal_regions *
		input.fai * n_coarse_axial_layers;

	// define number of pins per assembly (assuming 17 x 17 assemblies)
	int pins_per_assembly = 289;

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
