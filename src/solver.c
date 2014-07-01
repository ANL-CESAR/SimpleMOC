#include"SimpleMOC_header.h"

// run one full transport sweep, return k
double transport_sweep( Params params, Input I )
{
	// Determine total number of tracks
	long ntracks_2D = I.n_azimuthal * (I.assembly_width * sqrt(2) / I.radial_ray_sep);
	long ntracks = ntracks_2D * (I.n_polar_angles * (int) ( I.height / I.axial_z_sep));  

	// calculate the height of a node's domain and of each FSR
	double node_delta_z = I.height / I.decomp_assemblies_ax;
	double source_delta_z = I.cai / I.fai;

	// loop over tracks (implicitly azimuthal angles, tracks in azimuthal angles,
	// polar angles, and z stacked rays)
	for (long i = 0; i < ntracks; i++)
	{
		long id = params.tracks[i].track2D_id;

		// initialize z and s 
		// (s = current track length progressed along a 2D segment)
		double z = params.tracks[i].z_height;
		double s;

		// get 2D track segments
		long num_2D_segs = params.tracks_2D[id].n_segments;

		// booleans to determine whether in bounds
		bool in_bounds = true;

		// cycle through all segments in the 2D track laydown
		for (long j =0; j < num_2D_segs; j++)
		{
			bool seg_complete = false;
			s = 0;
			while( !seg_complete )
			{
				// calculate new z coordinate if segment is completed
				double new_z = z + ( params.tracks_2D[id].segments[j].length - s )
				   	/ tan( params.tracks[i].p_angle );

				// check if segment is completed
				if( (int) ( new_z / source_delta_z ) == 
						(int) ( z / source_delta_z) )
					seg_complete = true;

				// otherwise calculate new z coordinate 
				// (hitting the edge of an FSR)
				else
				{
					new_z = (double) ( (int) (z / source_delta_z) );
					if( params.tracks[i].p_angle < M_PI / 2.0)
						new_z += source_delta_z;
					s += ( new_z - z) * tan(params.tracks[i].p_angle);
				}

				// calculate distance traveled in the FSR
				double dist = (new_z - z) / cos(params.tracks[i].p_angle);
					
				// determine if ray is out of bounds
				if( new_z <= 0 || new_z >= node_delta_z)
				{
					in_bounds = false;
					break;
				}

				// TODO: calculate which FSR to tally to
				// TODO: tally to FSR using dist and XS data
				// TODO: update source in FSR for next transport sweep

				// set the new z coordinate
				z = new_z;

			}
			if( !in_bounds )
				break;
		}
	}

	return 0;
}
