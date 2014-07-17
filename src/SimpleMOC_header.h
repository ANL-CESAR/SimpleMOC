#ifndef __SimpleMOC_header
#define __SimpleMOC_header

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<stdbool.h>
#include<limits.h>
#include<assert.h>

#ifdef MPI
#include<mpi.h>
#endif

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
	bool decompose;            // Turn decomposition on/off
	int decomp_assemblies_ax;  // Number of assemblies per sub-domain (axially)
	long segments_per_track;   // Average number of segments per track
	double assembly_width;     // Width of an assembly - 1.26 x 17 cm
	double height;             // Height of the reactor - 400 cm
	long n_2D_source_regions_per_assembly; // 3M source regions per assembly (estimate)
	long n_source_regions_per_node; // Number of source regions in a given node
	long mype;                 // MPI Rank
	long ntracks_2D;           // Number of 2D tracks (derived)
	int z_stacked;             // Number of z rays (derived)
	long ntracks;              // Total number of 3D tracks per assembly (derived)
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
	double p_weight;           // Polar Quadrature Weight     (rand)
	double z_height;           // Z-height
	double * start_flux;       // Starting (input) flux array received from inputting neighbor
	long rank_in;              // MPI rank to receive from
	double * end_flux;         // Attenuated (output) flux array to send to output neighbor
	long rank_out;             // MPI rank to send to
	double * psi;			   // current angular flux along track
} Track;

// Source Region Structure
typedef struct{
	double ** fine_flux;
	double ** source_params;
	double vol;
	double ** XS;
	double ** scattering_matrix;
} Source;

// Params Structure for easier data pointer passing
typedef struct{
	Track2D * tracks_2D;
	Track *** tracks;
	Source * sources;
   	double * polar_angles;	
} Params;

// MPI 3D Grid info
typedef struct{
	MPI_Comm cart_comm_3d;
	MPI_Datatype Flux_Array;
	int x_pos_src;
	int x_pos_dest;
	int x_neg_src;
	int x_neg_dest;
	int y_pos_src;
	int y_pos_dest;
	int y_neg_src;
	int y_neg_dest;
	int z_pos_src;
	int z_pos_dest;
	int z_neg_src;
	int z_neg_dest;
} CommGrid;

// init.c
Input get_input( void );
Params build_tracks( Input I );
CommGrid init_mpi_grid( Input I );

// io.c
void logo(int version);
void center_print(const char *s, int width);
void border_print(void);
void fancy_int( int a );
void print_input_summary(Input input);

// tracks.c
Track2D * generate_2D_tracks( Input input, size_t * nbytes );
void generate_2D_segments( Input input, Track2D * tracks, long ntracks, size_t * nbytes );
void free_2D_tracks( Track2D * tracks );
Track *** generate_tracks(Input input, Track2D * tracks_2D, size_t * nbytes);
void free_tracks( Track *** tracks );
long segments_per_2D_track_distribution( Input I );
double * generate_polar_angles( Input I );

// utils.c
double urand(void);
double nrand(double mean, double sigma);
double pairwise_sum(double * vector, long size);
double * quadratic_fit(double * data, double xlen, int nx);

// source.c
Source * initialize_sources( Input I, size_t * nbytes );
void free_sources( Input I, Source * sources );

// solver.c
double transport_sweep( Params params, Input I );
void attenuate_fluxes( Track * track, Source * QSR, int fine_id, double ds, int groups, double mu ); 
void renormalize_flux( Params params, Input I );
double update_sources( Params params, Input I, double keff );

// test.c
void gen_norm_pts(double mean, double sigma, int n_pts);

// comms.c
void transfer_boundary_fluxes( Params params, Input I, CommGrid grid);

#endif
