#include"SimpleMOC_header.h"

// Gets I from user and sets defaults
Input get_input( void )
{
	Input I;

	I.x_assemblies = 15;        // Number of assemblies in the x-axis of the reactor
	I.y_assemblies = 15;        // Number of assemblies in the y-axis of the reactor
	I.cai = 27;                 // This is the number of course axial intervals
	I.fai = 5;                  // This is the number of fine axial intervals per course axial interval
	I.axial_exp = 2;            // Axial source expansion order
	I.radial_ray_sep = 0.1;     // Radial ray separation
	I.axial_z_sep = 0.2;        // Axial stacked z-ray separation
	I.n_azimuthal = 32;         // Number of azimuthal angles
	I.n_polar_angles = 10;      // Number of polar angles
	I.n_egroups = 100;          // Number of energy groups
	I.decompose = 0;            // Turn decomposition on/off (1 on, 0 off)
	I.decomp_assemblies_ax = 1; // Number of assemblies per sub-domain (axially)
	I.assembly_width = 1.26*17; // Width of an assembly - 1.26 x 17 cm
	I.height = 400.0;           // Height of the reactor - 400 cm

	// TODO: Add file/CLI user input

	return I;
}
