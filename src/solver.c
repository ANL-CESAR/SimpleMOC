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
	int z_stacked = (int) ( I.height / (I.axial_z_sep * I.decomp_assemblies_ax) );
	long ntracks = ntracks_2D * I.n_polar_angles * z_stacked;  

	// calculate the height of a node's domain and of each FSR
	double node_delta_z = I.height / I.decomp_assemblies_ax;
	double fine_delta_z = I.height / (I.cai * I.fai);

	// initialize fluxes in transport sweep
	// TODO: optimize this: it is probably unnecessary depending on how we strcuture
	// communication between nodes; it may be necessary to have this step followed
	// by a barrier
	for( int i = 0; i < ntracks_2D; i++)
		for( int j = 0; j < I.n_polar_angles; j++)
			for( int k = 0; k < z_stacked; k++)
				Track * track = params.tracks[i][j][k];
				for( int g = 0; g < I.n_egroups; g++)
					track->psi[g] = track->start_flux[g];
	
	// Start transport sweep

	// loop over tracks (implicitly azimuthal angles, tracks in azimuthal angles,
	// polar angles, and z stacked rays)
	for (long i = 0; i < ntracks_2D; i++)
	{
		// print progress
		if( i % 50 == 0)
			printf("%s%ld%s%ld\n","2D Tracks Completed = ", i," / ", ntracks_2D);

		// treat positive-z traveling rays first
		for( int j = 0; j < I.n_polar_angles / 2; j++)
		{
			double p_angle = params.polar_angles[j];

			// start with all z stacked rays
			int end_stacked = z_stacked;
			for( int n = 0; n < params.tracks_2D[i].n_segments; n++)
			{
				// calculate distance traveled in cell if segment completed
				double s = params.tracks_2D[i].segments[n].length / sin(p_angle);

				// allocate varaible for distance traveled in an FSR
				double ds;

				// loop over remaining z-stacked rays
				for( int k = 0; k < end_stacked; k++)
				{
					// select current track
					Track * track = &params.tracks[i][j][k];

					// set flag for completeion of segment
					bool seg_complete = false;

					while( !seg_complete )
					{
						// calculate new height based on s (distance traveled in FSR)
						double z = track->z_height + s * cos(p_angle);

						// check if still in same FSR (fine axial interval)
						if( (int) ( track->z_height / fine_delta_z ) == 
								(int) ( z / fine_delta_z ) )
						{
							seg_complete = true;
							ds = s;
						}

						// otherwise, we need to recalculate distances
						else
						{
							// correct z
							int interval = (int) (track->z_height / fine_delta_z);
							z = fine_delta_z * (double) (interval + 1);

							// calculate distance travelled in FSR (ds)
							ds = (z - track->z_height) / cos(p_angle);

							// update track length remaining
							s -= ds;

							// check if out of bounds or track complete
							if( s <= 0 || z >= node_delta_z )
							{
								// mark segment as completed
								seg_complete = true;

								// remember to no longer treat this track
								end_stacked--;

								// reset z height (calculate from k)
								track->z_height = I.axial_z_sep * k;
							}
						}

						// update with new z height or reset if finished
						if( n == params.tracks_2D[i].n_segments - 1 )
							track->z_height = I.axial_z_sep * k;
						else
							track->z_height = z;

						// pick a random FSR (cache miss expected)
						long FSR_id = rand() % I.n_source_regions_per_node;

						// update sources and fluxes from attenuation over FSR
						attenuate_fluxes( track, &params.sources[FSR_id], ds, I.n_egroups );
					}
				}
			}
		}

		// treat negative-z traveling rays next
		for( int j = I.n_polar_angles / 2; j < I.n_polar_angles; j++)
		{
			double p_angle = params.polar_angles[j];
			int begin_stacked = 0;

			for( int n = 0; n < params.tracks_2D[i].n_segments; n++)
			{
				// calculate distance traveled in cell if segment completed
				double s = params.tracks_2D[i].segments[n].length / sin(p_angle);

				// allocate varaible for distance traveled in an FSR
				double ds;

				// loop over all z stacked rays to begin
				for( int k = begin_stacked; k < z_stacked; k++)
				{
					// select current track
					Track * track = &params.tracks[i][j][k];

					// set flag for completeion of segment
					bool seg_complete = false;
					
					while( !seg_complete )
					{
						// calculate new height based on s (distance traveled in FSR)
						double z = track->z_height + s * cos(p_angle);
						
						// check if still in same FSR (fine axial interval)
						// NOTE: a bit of trickery this time using the fact that 
						// 2147483647 is the largest integer value
						int val1 = 2147483647 - (int) (2147483647 - track->z_height
							   	/ fine_delta_z);
						int val2 = 2147483647 - (int) (2147483647 - z / fine_delta_z);
						if( val1 == val2  )
						{
							seg_complete = true;
							ds = s;
						}

						// otherwise, we need to recalculate distances
						else
						{
							// correct z
							int interval = val1 - 1;
							z = fine_delta_z * (double) interval;

							// calculate distance travelled in FSR (ds)
							ds = ( z - track->z_height )	/ cos(p_angle);

							// update track length remaining
							s -= ds;

							// check if out of bounds or track complete
							if( z <= 0 )
							{
								// mark segment as completed
								seg_complete = true;

								// remember to no longer treat this track
								begin_stacked++;

								// reset z height (calculate from k)
								track->z_height = I.axial_z_sep * (k+1);
							}
						}

						// update with new z height or reset if finished
						if( n == params.tracks_2D[i].n_segments - 1 )
							track->z_height = I.axial_z_sep * (k+1); 
						else
							track->z_height = z;

						// pick a random FSR (cache miss expected)
						long FSR_id = rand() % I.n_source_regions_per_node;

						// update sources and fluxes from attenuation over FSR
						attenuate_fluxes( track , &params.sources[FSR_id], ds, I.n_egroups );
					}
				}
			}
		}
	}

	transfer_boundary_fluxes(params);

	// TODO: calculate a real keff, but maybe this can be disregarded?
	return 0;
}

void attenuate_fluxes( Track * track, Source * FSR, double ds, int groups ) 
{
	// compute weight (azimuthal * polar)
	// TODO: add track weight (area), also add azimuthal (tracks_2D[i].az_weight)
	double weight = track->p_weight;

	// cycle over energy groups
	for( int g = 0; g < groups; g++)
	{
		// load XS data
		double sigT = FSR->XS[g][0];
		double nuSigF = FSR->XS[g][1];
		double chi = FSR->XS[g][2];

		// calculate exponential
		// TODO: Maybe compute (1 - exp) ?? (OpenMOC), also use table lookup
		double exponential = exp( - sigT * ds );

		// calculate change in angular flux
		double delta_psi = (track->psi[g] - FSR->source[g]/sigT) *
					(1.0 - exponential);

		// add contribution to new source flux
		FSR->flux[g] += delta_psi * weight;
					
		// update angular flux
		track->psi[g] -= delta_psi;
	}

}	

void transfer_boundary_fluxes( Params params)
{
	return;
}

void add_source_to_flux( Params params, Input I )
{
	// add source contribution to scalar flux in each FSR
	for( int i = 0; i < I.n_source_regions_per_node; i++)
	{
		Source src = params.sources[i];
		for( int k = 0; k < I.n_egroups; k++)
		{
			double sigT = src.XS[k][0];

			// TODO: determine why this line is here
			src.flux[k] *= 0.5;

			// TODO: Use reduced source for computational efficiency
			// ALSO, maybe store 1/volume instead of volume
			src.flux[k] = (4 * M_PI * src.source[k]/ sigT + src.flux[k] / src.vol )
			   	/ sigT;
		}
	}

	// tally total fission rate (pair-wise sum)
	double * fission_rates = malloc( I.n_source_regions_per_node * sizeof(double) );
	for( int i = 0; i < I.n_source_regions_per_node; i++)
	{
		Source src = params.sources[i];
		double * g_fission_rates = malloc( I.n_egroups * sizeof(double) );
		for( int g = 0; g < I.n_egroups; g++)
			g_fission_rates[g] = src.flux[g] * src.vol * src.XS[g][1];
		fission_rates[i] = pairwise_sum( g_fission_rates, I.n_egroups );
	}
	double total_fission_rate = pairwise_sum(fission_rates, 
			I.n_source_regions_per_node);

	// normalize fluxes by fission reaction rate (TODO: Why by fission rate??)
	double norm_factor = 1.0 / total_fission_rate;
	for( int i = 0; i < I.n_source_regions_per_node; i++)
		for( int k = 0; k < I.n_egroups; k++)
			params.sources[i].flux[k] *= norm_factor;

	// TODO: Normalize boundary fluxes by same factor as well for
	// non-vacuum boundary conditions
	return;
}


double update_sources( Params params, Input I, double keff )
{
	// source residual
	double residual;

	// calculate inverse multiplication facotr for efficiency
	double inverse_k = 1.0 / keff;

	// calculate new source
	for( int i = 0; i < I.n_source_regions_per_node; i++)
	{
		Source src = params.sources[i];
		
		// allocate new source
		double * new_source = (double * ) malloc(I.n_egroups * sizeof(double));

		// calculate total fission source and scattering source
		double fission_source;
		double scatter_source;

		// allocate arrays for summation
		double * fission_rates = malloc(I.n_egroups * sizeof(double));
		double * scatter_rates = malloc(I.n_egroups * sizeof(double));
	
		// compute total fission source
		for( int g = 0; g < I.n_egroups; g++ )
			fission_rates[g] = src.flux[g] * src.XS[g][1];
		fission_source = pairwise_sum( fission_rates, (long) I.n_egroups);
		free(fission_rates);

		// normalize fission source by multiplication factor
		fission_source *= inverse_k;

		// compute scattering and new total source for each group
		for( int g = 0; g < I.n_egroups; g++ )
		{
			double * scatter_vector = src.scattering_matrix[g];
			for( int g2 = 0; g2 < I.n_egroups; g2++ )
			{
				// compute scatter source originating from g2 -> g
				scatter_rates[g2] = src.scattering_matrix[g][g2] * 
					src.flux[g2];
			}
			scatter_source = pairwise_sum(scatter_rates, (long) I.n_egroups);

			// compuate new total source
			double chi = src.XS[g][2];
			new_source[g] = (fission_source * chi + scatter_source) / (4.0 * M_PI);
		}

		free(scatter_rates);

		// assign new source to the actual source (changing pointers)
		for( int g = 0; g < I.n_egroups; g++ )
			src.source[g] = new_source[g];

		free(new_source);

	}

	// NOTE: See code around line 600 of CPUSolver.cpp in ClosedMOC/ OpenMOC

	// TODO: calculate real source residual
	residual = 0;
	return residual;
}
