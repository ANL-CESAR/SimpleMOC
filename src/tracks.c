#include"SimpleMOC_header.h"

// Allocates and initializes an array of 2D tracks
Track2D * generate_2D_tracks( Input I )
{
	// Determine number of 2D tracks, a conservative estimate
	long ntracks = I.n_azimuthal * (I.assembly_width * sqrt(2) / I.radial_ray_sep);
	
	// Allocate space for 2D tracks
	Track2D * tracks = (Track2D *) malloc( ntracks * sizeof(Track2D));

	// Fill weights with randomized data
	for( int i = 0; i < ntracks; i++ )
		tracks[i].az_weight = urand();

	// Allocate and randomize segments
	generate_2D_segments( I, tracks, ntracks );

	return tracks;
}

// Allocates and initializes all segments
void generate_2D_segments( Input I, Track2D * tracks, long ntracks )
{
	// Randomize Number of segments per track, and accumulate total 2D segments in assembly
	long total_segments = 0;
	for( long i = 0; i < ntracks; i++ )
	{
		// TODO: Change from even to normal distribution
		tracks[i].n_segments = segments_per_2D_track_distribution( I );
		total_segments += tracks[i].n_segments;
	}
	
	// Allocate contiguous space for segments
	Segment * contiguous_segments = (Segment *) malloc( total_segments * sizeof(Segment));

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
	return rand() % I.segments_per_track;
}

long segments_per_3D_track_distribution( Input I )
{
	return (rand() % I.segments_per_track )* 0.5;
}

void free_2D_tracks( Track2D * tracks )
{
	free(tracks[0].segments);
	free(tracks);
}

Track * generate_tracks(Input I, Track2D * tracks_2D)
{
	// Determine total number of tracks
	long ntracks_2D = I.n_azimuthal * (I.assembly_width * sqrt(2) / I.radial_ray_sep);
	long ntracks = ntracks_2D * (I.n_polar_angles * (int) ( I.height / I.axial_z_sep));  

	// Allocate space for tracks (3D)
	Track * tracks = (Track *) malloc( ntracks * sizeof(Track));

	// Initialize tracks randomly
	// TODO: - a few things left to init regarding domains / MPI
	for( long i = 0; i < ntracks; i++ )
	{
		// TODO: perhaps we might not want to make track2D_id associated randomly
		// since we might parallelize over polar angles (i.e. we can pre-fetch
		// 2D track data)
		tracks[i].track2D_id = rand() % ntracks_2D;
		tracks[i].p_angle = urand() * M_PI;
		tracks[i].p_weight = urand();
		
		// TODO: change these for domain decomposed
		tracks[i].z_height = urand() * I.height;
		tracks[i].start_flux = (double *) malloc( I.n_egroups * sizeof(double) );
		tracks[i].end_flux = (double *) malloc( I.n_egroups * sizeof(double) );

		// set incoming flux to 0
		for( int j = 0; j < I.n_egroups; j++)
			tracks[i].start_flux[j] = 0;
	}

	return tracks;
}

void free_tracks( Track * tracks )
{
	free(tracks);
}
