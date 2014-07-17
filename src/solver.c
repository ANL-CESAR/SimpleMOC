#include"SimpleMOC_header.h"

// run one full transport sweep, return k

// TODO: crank more efficiency and accuracy out of code 
// (i.e. less divisions, pairwise additions, precompute
// values used in many iterations, etc) see OpenMOC

double transport_sweep( Params params, Input I )
{
	if(I.mype==0) printf("Starting transport sweep ...\n");

	// Determine total number of tracks
	long ntracks_2D = I.n_azimuthal * (I.assembly_width * sqrt(2) / I.radial_ray_sep);
	int z_stacked = (int) ( I.height / (I.axial_z_sep * I.decomp_assemblies_ax) );
	long ntracks = ntracks_2D * I.n_polar_angles * z_stacked;  

	// calculate the height of a node's domain and of each FSR
	double node_delta_z = I.height / I.decomp_assemblies_ax;
	double fine_delta_z = I.height / (I.cai * I.fai);

	// initialize fluxes in transport sweep
	// TODO: insert a barrier
	for( int i = 0; i < ntracks_2D; i++)
		for( int j = 0; j < I.n_polar_angles; j++)
			for( int k = 0; k < z_stacked; k++)
			{
				Track * track = &params.tracks[i][j][k];
				for( int g = 0; g < I.n_egroups; g++)
					track->psi[g] = track->start_flux[g];
			}

	// Start transport sweep

	// loop over tracks (implicitly azimuthal angles, tracks in azimuthal angles,
	// polar angles, and z stacked rays)
	for (long i = 0; i < ntracks_2D; i++)
	{
		// print progress
		if( i % 50 == 0)
			if(I.mype==0) printf("%s%ld%s%ld\n","2D Tracks Completed = ", i," / ", ntracks_2D);

		// treat positive-z traveling rays first
		for( int j = 0; j < I.n_polar_angles / 2; j++)
		{
			double p_angle = params.polar_angles[j];
			double mu = cos(p_angle);

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
						long QSR_id = rand() % I.n_source_regions_per_node;
						int fine_id = rand() % I.fai;

						// update sources and fluxes from attenuation over FSR
						attenuate_fluxes( track, &params.sources[QSR_id], fine_id, 
								ds, I.n_egroups, mu );
					}
				}
			}
		}

		// treat negative-z traveling rays next
		for( int j = I.n_polar_angles / 2; j < I.n_polar_angles; j++)
		{
			double p_angle = params.polar_angles[j];
			double mu = cos(p_angle);
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
						int val1 = INT_MAX - (int) (INT_MAX - track->z_height
								/ fine_delta_z);
						int val2 = INT_MAX - (int) (INT_MAX - z / fine_delta_z);
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
							ds = ( z - track->z_height ) / cos(p_angle);

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
						int fine_id = rand() % I.fai;

						// update sources and fluxes from attenuation over FSR
						attenuate_fluxes( track , &params.sources[FSR_id], fine_id,
								ds, I.n_egroups, mu );
					}
				}
			}
		}
	}


	// TODO: calculate a real keff
	return 0;
}

void attenuate_fluxes( Track * track, Source * QSR, int fine_id, double ds, int groups, double mu ) 
{
	// compute weight (azimuthal * polar)
	// TODO: add track weight (area), also add azimuthal 
	// (tracks_2D[i].az_weight)
	double weight = track->p_weight;
	double mu2 = mu * mu;

	// load fine source region flux vector
	double * FSR_flux = QSR -> fine_flux[fine_id];

	// cycle over energy groups
	for( int g = 0; g < groups; g++)
	{
		// load total cross section
		double sigT = QSR->XS[g][0];

		// load source components
		double * sourceParams = QSR->source_params[g];
		double q0 = sourceParams[0];
		double q1 = sourceParams[1];
		double q2 = sourceParams[2];

		// calculate exponential
		// TODO: Compute (1 - exp) {OpenMOC} using table lookup
		double expVal = 1.0 - exp( - sigT * ds );

		// add contribution to new source flux
		double tau = sigT * ds;
		double sigT2 = sigT * sigT;
		double flux_integral = (q0 * tau + (sigT * track->psi[g] - q0) * expVal)
			/ sigT2
			+ q1 * mu * (tau * (tau - 2) + 2 * expVal)
			/ (sigT * sigT2)
			+ q2 * mu2 * (tau * (tau * (tau - 3) + 6) - 6 * expVal)
			/ (3 * sigT2 * sigT2);
		FSR_flux[g] += weight * flux_integral;

		// update angular flux
		track->psi[g] = track->psi[g] * (1.0 - expVal) + q0 * expVal / sigT
			+ q1 * mu * (tau - expVal) / sigT2 + q2 * mu2 *
			(tau * (tau - 2) + 2 * expVal) / (sigT2 * sigT);
	}
}	


// renormalize flux for next transport sweep iteration
void renormalize_flux( Params params, Input I )
{
	// tally total fission rate (pair-wise sum)
	double * fission_rates = malloc( I.n_source_regions_per_node * sizeof(double) );
	double * fine_fission_rates = malloc( I.fai * sizeof(double) );
	double * g_fission_rates = malloc( I.n_egroups * sizeof(double) );

	// TODO: Add communication between ranks to accumulate total fission rate
	// (another for loop will be required over ranks)
	for( int i = 0; i < I.n_source_regions_per_node; i++)
	{
		Source src = params.sources[i];
		for( int j = 0; j < I.fai; j++)
		{
			for( int g = 0; g < I.n_egroups; g++)
				g_fission_rates[g] = src.fine_flux[j][g] * src.vol * src.XS[g][1];
			fine_fission_rates[j] = pairwise_sum( g_fission_rates, I.n_egroups );
		}
		fission_rates[i] = pairwise_sum( fine_fission_rates, I.fai );
	}
	double total_fission_rate = pairwise_sum(fission_rates, 
			I.n_source_regions_per_node);

	// free allocated memory
	free(fission_rates);
	free(fine_fission_rates);
	free(g_fission_rates);

	// normalize fluxes by fission reaction rate
	double norm_factor = 1.0 / total_fission_rate;
	for( int i = 0; i < I.n_source_regions_per_node; i++)
	{
		Source * src = &params.sources[i];
		double adjust = norm_factor * 4 * M_PI * I.fai / src->vol;
		for( int k = 0; k < I.fai; k++)
			for( int g = 0; g < I.n_egroups; g++)
				src->fine_flux[k][g] *= adjust;
	}

	// TODO: Normalize boundary fluxes by same factor as well for
	// non-vacuum boundary conditions
	return;
}

// Updates sources for next iteration
double update_sources( Params params, Input I, double keff )
{
	// source residual
	double residual;

	// calculate inverse multiplication facotr for efficiency
	double inverse_k = 1.0 / keff;

	// allocate fine sources
	double ** new_source = (double **) malloc(I.fai * sizeof(double *));
	double * new_source_data = 
		(double *) malloc(I.fai * I.n_egroups * sizeof(double));
	for( int i = 0; i < I.fai; i++)
		new_source[i] = &new_source_data[I.n_egroups * i];
		
	// allocate arrays for summation
	double * fission_rates = malloc(I.n_egroups * sizeof(double));
	double * scatter_rates = malloc(I.n_egroups * sizeof(double));

	// allocate array for data fitting
	double * fit_array = malloc(I.fai * sizeof(double));

	// cycle through all coarse axial intervals to update source
	for( int i = 0; i < I.n_source_regions_per_node; i++)
	{
		Source src = params.sources[i];

		// cycle thorugh all fine axial regions to calculate new source
		for( int j = 0; j < I.fai; j++)
		{
			// calculate total fission source and scattering source
			double fission_source;
			double scatter_source;

			// compute total fission source
			for( int g = 0; g < I.n_egroups; g++ )
				fission_rates[g] = src.fine_flux[j][g] * src.XS[g][1];
			fission_source = pairwise_sum( fission_rates, (long) I.n_egroups);

			// normalize fission source by multiplication factor
			fission_source *= inverse_k;

			// compute scattering and new total source for each group
			for( int g = 0; g < I.n_egroups; g++ )
			{
				for( int g2 = 0; g2 < I.n_egroups; g2++ )
				{
					// compute scatter source originating from g2 -> g
					scatter_rates[g2] = src.scattering_matrix[g][g2] * 
						src.flux[g2];
				}
				scatter_source = pairwise_sum(scatter_rates, (long) I.n_egroups);

				// compuate new total source
				double chi = src.XS[g][2];

				// calculate new source in fine axial interval assuming isotropic
				new_source[j][g] = (fission_source * chi + scatter_source) / (4.0 * M_PI);
			}
		}
		// fit quadratic function to fine sources for each energy group
		for( int g = 0; g < I.n_egroups; g++)
		{
			for( int j = 0; j < I.fai; j++)
				fit_array[j] = new_source[j][g];
			src.source_params[g] = quadratic_fit(fit_array, 
					(int) I.height / (I.fai * I.cai), I.fai);
		}

	}

	// free memory
	free(new_source);
	free(new_source_data);
	free(fission_rates);
	free(scatter_rates);
	free(fit_array);

	// NOTE: See code around line 600 of CPUSolver.cpp in ClosedMOC/ OpenMOC

	// TODO: calculate real source residual
	residual = 0;
	return residual;
}

