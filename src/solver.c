#include"SimpleMOC_header.h"

// run one full transport sweep, return k

// TODO: crank more efficiency and accuracy out of code 
// (i.e. less divisions, pairwise additions, precompute
// values used in many iterations, etc) see OpenMOC

double transport_sweep( Params params, Input I )
{
	printf("Starting transport sweep ...\n");

	// Determine total number of tracks
	long ntracks_2D = I.n_azimuthal * (I.assembly_width * sqrt(2) / I.radial_ray_sep);
	long ntracks = ntracks_2D * (I.n_polar_angles * (int) ( I.height / I.axial_z_sep));  

	// calculate the height of a node's domain and of each FSR
	double node_delta_z = I.height / I.decomp_assemblies_ax;
	double source_delta_z = I.height / (I.cai * I.fai);

	// Start transport sweep

	// loop over tracks (implicitly azimuthal angles, tracks in azimuthal angles,
	// polar angles, and z stacked rays)
	for (long i = 0; i < ntracks; i++)
	{
		// get 2D track ID
		long id = params.tracks[i].track2D_id;

		// initialize z and s 
		// (s = current track length progressed along a 2D segment)
		double z = params.tracks[i].z_height;
		double s;

		// allocate an array to store the temporary flux (psi)
		double * psi = (double *) malloc( I.n_egroups * sizeof(double) );
		for( int k = 0; k < I.n_egroups; k++)
			psi[k] = params.tracks[i].start_flux[k];

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

				// pick a random FSR (cache miss expected)
				long source_id = rand() % I.n_source_regions_per_node;

				// compute weight (azimuthal * polar)
				//TODO: add track weight (area)
				double weight = params.tracks_2D[id].az_weight * params.tracks[i].p_weight;

				// cycle over energy groups
				for( int k = 0; k < I.n_egroups; k++)
				{
					// load XS data
					double sigT = params.sources[source_id].XS[k][0];
					double nuSigF = params.sources[source_id].XS[k][1];
					double chi = params.sources[source_id].XS[k][2];

					// calculate exponential
					// TODO: Maybe compute (1 - exp) ?? (OpenMOC), also use table lookup
					double exponential = exp( - sigT * dist );

					// calculate change in angular flux
					double delta_psi = (psi[k] - params.sources[source_id].source[k]/sigT) *
						(1.0 - exponential);

					// add contribution to new source flux
					params.sources[source_id].flux[k] += delta_psi * weight;
					
					// update angular flux
					psi[k] -= delta_psi;

				}	

				// TODO: initialize source scalar flux with source contribution

				// set the new z coordinate
				z = new_z;

			}
			if( !in_bounds )
				break;
		}
		free(psi);
	}

	// add source contribution to scalar flux in each FSR
	for( int i = 0; i < I.n_source_regions_per_node; i++)
	{
		for( int k = 0; k < I.n_egroups; k++)
		{
			double sigT = params.sources[i].XS[k][0];

			// TODO: determine why this line is here
			params.sources[i].flux[k] *= 0.5;

			// TODO: Use reduced source for computational efficiency
			// ALSO, maybe store 1/volume instead of volume
			params.sources[i].flux[k] = 4 * M_PI * params.sources[i].source[k]
				/ sigT + params.sources[i].flux[k] / (sigT * params.sources[i].vol);
		}
	}

	// tally total fission rate
	// TODO: change to pair-wise summation
	double total_fission_rate = 0;
	for( int i = 0; i < I.n_source_regions_per_node; i++)
	{
		for( int k = 0; k < I.n_egroups; k++)
		{
			total_fission_rate += params.sources[i].flux[k] * params.sources[i].vol
				* params.sources[i].XS[k][1];
		}
	}

	// normalize fluxes by fission reaction rate (TODO: Why by fission rate??)
	double norm_factor = 1.0 / total_fission_rate;
	for( int i = 0; i < I.n_source_regions_per_node; i++)
		for( int k = 0; k < I.n_egroups; k++)
			params.sources[i].flux[k] *= norm_factor;

	// TODO: Normalize boundary fluxes by same factor as well

	// TODO: change k to a real value
	double inverse_k = 2 * urand();

	// calculate new source
	for( int i = 0; i < I.n_source_regions_per_node; i++)
	{
		// allocate new source
		double * new_source = (double * ) malloc(I.n_egroups * sizeof(double));

		// calculate total fission source and scattering source
		double fission_source = 0;
		double scatter_source = 0;

		// TODO: change to pair-wise summantion
		for( int k = 0; k < I.n_egroups; k++ )
		{
			scatter_source = 0;
			for( int k2 = 0; k2 < I.n_egroups; k2++ )
			{
				// compute fission source if not computed yet
				if( k == 0)
					fission_source += params.sources[i].flux[k2] *
						params.sources[i].XS[k2][1];

				// compute scatter source
				// NOTE: this means scatter from k2 -> k
				scatter_source += params.sources[i].scattering_matrix[k][k2] * 
					params.sources[i].flux[k2];
			}

			// normalize fission source by multiplication factor if needed
			if ( k == 0 )
				fission_source *= inverse_k;

			// compuate new total source
			double chi = params.sources[i].XS[k][2];
			new_source[k] = (fission_source * chi + scatter_source) / (4.0 * M_PI);
		}

		// assign new source to the actual source (changing pointers)
		for( int k = 0; k < I.n_egroups; k++ )
			params.sources[i].flux[k] = new_source[k];

		// TODO: free old memory if needed

	}

	// NOTE: See code around line 600 of CPUSolver.cpp in ClosedMOC/ OpenMOC

	// TODO: calculate a real keff
	return 0;
}
