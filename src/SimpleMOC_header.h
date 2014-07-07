#ifndef __SimpleMOC_header
#define __SimpleMOC_header

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<stdbool.h>
		
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
	double assembly_width;     // Width of an assembly - 1.26 x 17 cm
	double height;             // Height of the reactor - 400 cm
	long n_2D_source_regions_per_assembly; // 3M source regions per assembly (estimate)
	long n_source_regions_per_node; // Number of source regions in a given node
} Input;

// Localized geometrical region ID
typedef struct{
	long assembly;             // Assembly ID
	long pin;                  // Pin Cell ID
	long zone;                 // Zone (inside pin cell) ID
} RegionID;

// Segment Structure
typedef struct{
	double length;
	long source_id;
} Segment;

// Track2D Structure
typedef struct{
	double az_weight;          // Azimuthal Quadrature Weight (rand)
	long n_segments;           // Number of Segments (gaussian)
	Segment * segments;        // Array of Segments
} Track2D;

// Track Structure
typedef struct{
	long track2D_id;           // Link into 2D geometry Track ID
	double p_angle;            // Polar Angle
	double p_weight;           // Polar Quadrature Weight     (rand)
	double z_height;           // Z-height
	double * start_flux;       // Starting (input) flux array received from inputting neighbor
	long rank_in;              // MPI rank to receive from
	double * end_flux;         // Attenuated (output) flux array to send to output neighbor
	long rank_out;             // MPI rank to send to
} Track;

// Flat Source Region Structure
typedef struct{
	double ** XS;
	double ** scattering_matrix;
	double * flux;
	double * source;
	double vol;
} Source;

// Params Structure for easier data pointer passing
typedef struct{
	Track2D * tracks_2D;
	Track * tracks;
	Source * sources; 
	// TODO: need material, XS data etc. separate from source data
} Params;

// init.c
Input get_input( void );
Params build_tracks( Input I );

// io.c
void logo(int version);
void center_print(const char *s, int width);
void border_print(void);
void fancy_int( int a );
void print_input_summary(Input input);

// tracks.c
Track2D * generate_2D_tracks( Input input );
void generate_2D_segments( Input input, Track2D * tracks, long ntracks );
void free_2D_tracks( Track2D * tracks );
Track * generate_tracks(Input input, Track2D * tracks_2D);
void free_tracks( Track * tracks );
long segments_per_2D_track_distribution( Input I );

// utils.c
double urand(void);
double nrand(double mean, double sigma);

// source.c
Source * initialize_sources( Input I );
void free_sources( Input I, Source * sources );

// solver.c
double transport_sweep( Params params, Input I );


// test.c

#endif
