#include"SimpleMOC_header.h"

// Allocates and initializes an array of 2D tracks
Track2D * generate_2D_tracks( Input I, size_t * nbytes )
{
	// Determine number of 2D tracks, a conservative estimate
	long ntracks = I.n_azimuthal * (I.assembly_width * sqrt(2) / I.radial_ray_sep);
	
	// Allocate space for 2D tracks
	Track2D * tracks = (Track2D *) malloc( ntracks * sizeof(Track2D));
	*nbytes += ntracks * sizeof(Track2D);

	// Fill weights with randomized data
	for( int i = 0; i < ntracks; i++ )
		tracks[i].az_weight = urand();

	// Allocate and randomize segments
	generate_2D_segments( I, tracks, ntracks, nbytes );

	return tracks;
}

// Allocates and initializes all segments
void generate_2D_segments( Input I, Track2D * tracks, long ntracks, size_t * nbytes )
{
	// Randomize Number of segments per track, and accumulate total 2D segments in assembly
	long total_segments = 0;
	for( long i = 0; i < ntracks; i++ )
	{
		tracks[i].n_segments = segments_per_2D_track_distribution( I );
		total_segments += tracks[i].n_segments;
	}
	
	// Allocate contiguous space for segments
	Segment * contiguous_segments = (Segment *) malloc( total_segments * sizeof(Segment));
	*nbytes += total_segments * sizeof(Segment);

	// Set segments arrays to correct locations within contiguous allocation
	long idx = 0;
	for( long i = 0; i < ntracks; i++ )
	{
		tracks[i].segments = &contiguous_segments[idx];
		idx += tracks[i].n_segments;
	}

	// Initialize segments to randomized values
	for( long i = 0; i < ntracks; i++ )
	{
		for( long j = 0; j < tracks[i].n_segments; j++ )
		{
			tracks[i].segments[j].length  = urand() * I.assembly_width
				/ tracks[i].n_segments;
			// source ID to be calculated on the fly
		}
	}
}

long segments_per_2D_track_distribution( Input I )
{
	return nrand(I.segments_per_track, sqrt(I.segments_per_track));
}

void free_2D_tracks( Track2D * tracks )
{
	free(tracks[0].segments);
	free(tracks);
}

Track *** generate_tracks(Input I, Track2D * tracks_2D, size_t * nbytes)
{
	// Determine total number of tracks
	long ntracks_2D = I.n_azimuthal * (I.assembly_width * sqrt(2) / I.radial_ray_sep);
	int z_stacked = (int) ( I.height / (I.axial_z_sep * I.decomp_assemblies_ax) );
	long ntracks = ntracks_2D * I.n_polar_angles * z_stacked;  

	// Allocate space for tracks (3D)
	Track *** tracks = (Track ***) malloc( ntracks_2D * sizeof(Track **));
	*nbytes += ntracks_2D * sizeof(Track **);

	// Allocate pointers for tracks associated with a 2D track
	Track ** tracks_in_track2D = (Track **) malloc( ntracks_2D *
		   	I.n_polar_angles * sizeof(Track *));
	*nbytes += ntracks_2D * I.n_polar_angles * sizeof(Track *);

	// Allocate complete array of track data
	Track * track_data = (Track *) malloc( ntracks * sizeof(Track) );
	*nbytes += ntracks * sizeof(Track);
	if(I.mype==0) printf("3D track data requires %zu MB of data...\n", ntracks * sizeof(Track) / 1024 / 1024 );

	// stitch pointers together
	for( long i = 0; i < ntracks_2D; i++ )
		tracks[i] = &tracks_in_track2D[ i * I.n_polar_angles ];

	for( long i = 0; i < ntracks_2D; i++ )
	{
		for( int j = 0; j < I.n_polar_angles; j++ )
		{
			tracks[i][j] = &track_data[ (i * I.n_polar_angles + j) * z_stacked ];
		}
	}

	// Allocate space for Flux Arrays
	size_t flux_bytes_needed = ntracks_2D * I.n_polar_angles * z_stacked * I.n_egroups * 3 * sizeof(double);
	if(I.mype==0) printf("Flux Arrays Require %zu MB of data...\n", flux_bytes_needed / 1024 / 1024);
	double * flux_space = (double *) malloc( flux_bytes_needed );
	*nbytes += flux_bytes_needed;
	size_t flux_idx = 0;

	long offset = ntracks_2D * I.n_polar_angles * z_stacked * I.n_egroups;

	for( long i = 0; i < ntracks_2D; i++ )
	{
		for( int j = 0; j < I.n_polar_angles; j++ )
		{
			for( int k = 0; k < z_stacked; k++ )
			{
				// bottom z heights should only have upward directed polar angles
				// similarly top should only have downward directed polar angles
				if( j < I.n_polar_angles/2 )
					tracks[i][j][k].z_height = I.axial_z_sep * k;
				else
					tracks[i][j][k].z_height = I.axial_z_sep * ( k + 1 );
				
				// set polar weight, NOTE: this is the same for same polar angle
				tracks[i][j][k].p_weight = urand();

				// Allocate start and end flux arrays
				// (Moved allocations out of loop so they are contiguous)
				tracks[i][j][k].start_flux = &flux_space[flux_idx];
				tracks[i][j][k].end_flux = &flux_space[offset + flux_idx];
				tracks[i][j][k].psi = &flux_space[offset*2 + flux_idx];
				flux_idx += I.n_egroups;

				// set incoming flux to 0, perhaps needs to be changed for inner assemblies
				for( int g = 0; g < I.n_egroups; g++)
					tracks[i][j][k].start_flux[g] = 0;
			}
		}
	}

	return tracks;
}

void free_tracks( Track *** tracks )
{
	free(tracks);
}

// assign polar angles
double * generate_polar_angles( Input I )
{
	double * polar_angles = (double *) malloc( I.n_polar_angles * sizeof(double) );
	for( int i = 0; i < I.n_polar_angles; i++)
		polar_angles[i] = M_PI * (i + 0.5) / I.n_polar_angles;
	return polar_angles;
}
