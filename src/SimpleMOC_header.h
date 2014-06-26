#ifndef __SimpleMOC_header
#define __SimpleMOC_header

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<rand.h>

// User inputs
typedef struct{
	int x_assemblies;          // Number of assemblies in the x-axis of the reactor
	int y_assemblies;          // Number of assemblies in the y-axis of the reactor
	int cai;                   // This is the number of course axial intervals
	int fai;                   // This is the number of fine axial intervals per course axial interval
	int axial_exp;             // Axial source expansion order
	double radial_ray_sep;     // Radial ray separation
	double axial_z_sep;        // Axial stacked z-ray separation
	int n_azimuthal;           // Number of azimuthal angles
	int n_polar_angles;        // Number of polar angles
	int n_egroups;             // Number of energy groups
	int decompose;             // Turn decomposition on/off (1 on, 0 off)
	int decomp_assemblies_ax;  // Number of assemblies per sub-domain (axially)
	long segments_per_track;   // Average number of segments per track
} Input;

//  Reactor definition
typedef struct{
	double assembly_width;     // Width of an assembly - 1.26 x 17 cm
	double height;             // Height of the reactor - 400 cm
	int n_radial_regions;      // Number of radial regions - default 10
	int n_azimuthal_regions;   // Number of azimuthal regions - default 8
	double * radii;            // Stores the radii of the radial regions
} Reactor;

// Localized geometrical region ID
typedef struct{
	long assembly;             // Assembly ID
	long pin;                  // Pin Cell ID
	long zone;                 // Zone (inside pin cell) ID
} RegionID;

// Cartesian Coordinate Struct
typedef struct{
	double x;
	double y;
	double z;
} Coord;

// Segment Structure
typedef struct{
	double dist;
	long source_id;
} Segment;

// Track2D Structure
typedef struct{
	double az_weight;          // Azimuthal Quadrature Weight (rand)
	double p_weight;           // Polar Quadrature Weight     (rand)
	long n_segments;           // Number of Segments (gaussian)
	Segment * segments;        // Array of Segments
} Track2D;

// Track Structure
typedef struct{
	long track2D_id;           // Link into 2D geometry Track ID
	double p_angle;            // Polar Angle
	double z_height;           // Z-height
	double start_flux;         // Starting (input) flux received from inputting neighbor
	long rank_in;              // MPI rank to receive from
	double end_flux;           // Attenuated (output) flux to send to output neighbor
	long rank_out;             // MPI rank to send to
} Track;

// 11 doubles

// goemetry.c
RegionID get_region_id( double x, double y, double z, Input input, Reactor reactor );
long get_region_index( RegionID id, Input input, Reactor reactor );
double * determine_radii( Reactor reactor );

// init.c
Input get_input( void );
Reactor reactor_init( void );

// io.c
void logo(int version);
void center_print(const char *s, int width);
void border_print(void);
void fancy_int( int a );
void print_input_summary(Input input);

// tester.c
void generate_2D_zone_points(Input input, Reactor reactor, int n_pts);
double urand(void);

#endif
