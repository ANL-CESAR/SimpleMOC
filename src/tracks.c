#include"SimpleMOC_header.h"

// Allocates and initializes an array of 2D tracks
Track2D * generate_2D_tracks( Input I )
{
	// Determine number of 2D tracks
	long ntracks = I.n_azimuthal * (I.assembly_width / I.radial_ray_sep);
	
	// Allocate space for 2D tracks
	Track2D * tracks = (Track2D *) malloc( ntracks * sizeof(Track2D));

	// Fill weights with randomized data
	for( int i = 0; i < ntracks; i++ )
	{
		tracks[i].az_weight = urand();
		tracks[i].p_weight = urand();
	}

	// Allocate and randomize segments
	generate_2D_segments( I, tracks, ntracks );

	return tracks;
}

// Allocates and initializes all segments
void generate_2D_segments( Input I, Track2D * tracks, long ntracks )
{
	// Randomize Number of segments per track, and accumulate total 2D tracks in assembly
	long total_tracks = 0;
	for( long i = 0; i < ntracks; i++ )
	{
		// TODO: Change from even to normal distribution
		tracks[i].n_segments = segments_per_2D_track_distribution( I );
		total_tracks += tracks[i].n_segments;
	}
	
	// Allocate contiguous space for segments
	Segment * contiguous_segments = (Segment *) malloc( total_tracks * sizeof(Segment));

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
			tracks[i].segments[j].length  = urand();
			// TODO: Perhaps this should be assigned in a less random fashion?
			tracks[i].segments[j].source_id = rand() % I.num_source_regions_per_assembly;
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
	long ntracks_2D = I.n_azimuthal * (I.assembly_width / I.radial_ray_sep);
	long ntracks = ntracks_2D * (I.n_polar_angles * ( I.height / I.axial_z_sep));  

	// Allocate space
	Track * tracks = (Track *) malloc( ntracks * sizeof(Track));

	// Initialize tracks randomly
	// TODO: - a few things left to init regarding domains / MPI
	for( long i = 0; i < ntracks; i++ )
	{
		tracks[i].track2D_id = rand() % ntracks_2D;
		tracks[i].p_angle = urand() * 2 * M_PI;
		tracks[i].z_height = I.height;
		tracks[i].start_flux = 0;
	}

	return tracks;
}

void free_tracks( Track * tracks )
{
	free(tracks);
}
