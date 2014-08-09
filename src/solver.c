#include"SimpleMOC_header.h"

void attenuate_fluxes( Track * track, Source * QSR, Input * I_in, 
		Params * params_in, float ds, float mu, float az_weight, AttenuateVars * A ) 
{
	Input I = *I_in;
	Params params = *params_in;
	// unload attenuate vars
	float * q0 = A->q0;
	float * q1 = A->q1;
	float * q2 = A->q2;
	float * sigT = A->sigT;
	float * tau = A->tau;
	float * sigT2 = A->sigT2;
	float * expVal = A->expVal;
	float * flux_integral = A->flux_integral;
	float * tally = A->tally;
	float * t1 = A->t1;
	float * t2 = A->t2;
	float * t3 = A->t3;
	float * f1 = A->f1;
	float * f2 = A->f2;
	float * f3 = A->f3;

	// compute fine axial interval spacing
	float dz = I.height / (I.fai * I.decomp_assemblies_ax * I.cai);

	// compute fine axial region ID
	int fine_id = (int) ( I.height / dz ) % I.cai;

	// compute z height in cell
	float zin = track->z_height - dz * ( (int)( track->z_height / dz ) + 0.5 );

	// compute weight (azimuthal * polar)
	// NOTE: real app would also have volume weight component
	float weight = track->p_weight * az_weight;
	float mu2 = mu * mu;

	// load fine source region flux vector
	float * FSR_flux = QSR -> fine_flux[fine_id];

	if( fine_id == 0 )
	{
		// cycle over energy groups
		#pragma simd
		for( int g = 0; g < I.n_egroups; g++)
		{
			// load neighboring sources
			float y2 = QSR->fine_source[fine_id][g];
			float y3 = QSR->fine_source[fine_id+1][g];

			// do linear "fitting"
			float c0 = y2;
			float c1 = (y3 - y2) / dz;

			// calculate q0, q1, q2
			q0[g] = c0 + c1*zin;
			q1[g] = c1;
			q2[g] = 0;
		}
	}
	else if ( fine_id == I.fai - 1 )
	{
		// cycle over energy groups
		#pragma simd
		for( int g = 0; g < I.n_egroups; g++)
		{
			// load neighboring sources
			float y1 = QSR->fine_source[fine_id-1][g];
			float y2 = QSR->fine_source[fine_id][g];

			// do linear "fitting"
			float c0 = y2;
			float c1 = (y2 - y1) / dz;

			// calculate q0, q1, q2
			q0[g] = c0 + c1*zin;
			q1[g] = c1;
			q2[g] = 0;
		}
	}
	else
	{
		// cycle over energy groups
		#pragma simd
		for( int g = 0; g < I.n_egroups; g++)
		{
			// load neighboring sources
			float y1 = QSR->fine_source[fine_id-1][g];
			float y2 = QSR->fine_source[fine_id][g];
			float y3 = QSR->fine_source[fine_id+1][g];

			// do quadratic "fitting"
			float c0 = y2;
			float c1 = (y1 - y3) / (2*dz);
			float c2 = (y1 - 2*y2 + y3) / (2*dz*dz);

			// calculate q0, q1, q2
			q0[g] = c0 + c1*zin + c2*zin*zin;
			q1[g] = c1 + 2*c2*zin;
			q2[g] = c2;
		}
	}


	// cycle over energy groups
	#pragma simd
	for( int g = 0; g < I.n_egroups; g++)
	{
		// load total cross section
		sigT[g] = QSR->sigT[g];

		// calculate common values for efficiency
		tau[g] = sigT[g] * ds;
		sigT2[g] = sigT[g] * sigT[g];
	}

	// cycle over energy groups
	#pragma simd
	for( int g = 0; g < I.n_egroups; g++)
		expVal[g] = interpolateTable( params.expTable, tau[g] );  

	// Flux Integral
	
	// Term 1
	#pragma simd
	for( int g = 0; g < I.n_egroups; g++)
	{
		f1[g] = (q0[g] * tau[g] + (sigT[g] * track->psi[g] - q0[g]) * expVal[g]) / sigT2[g]; 
	}
	// Term 2
	#pragma simd
	for( int g = 0; g < I.n_egroups; g++)
	{
		f2[g] = q1[g] * mu * (tau[g] * (tau[g] - 2.f) + 2.f * expVal[g]) / (sigT[g] * sigT2[g]); 
	}
	// Term 3
	#pragma simd
	for( int g = 0; g < I.n_egroups; g++)
	{
		f3[g] = q2[g] * mu2 * (tau[g] * (tau[g] * (tau[g] - 3.f) + 6.f) - 6.f * expVal[g]) / (3.f * sigT2[g] * sigT2[g]);
	}
	// Total
	#pragma simd
	for( int g = 0; g < I.n_egroups; g++)
	{
		flux_integral[g] = f1[g] + f2[g] + f3[g];
	}

	/*
	   #pragma simd
	   for( int g = 0; g < I.n_egroups; g++)
	   {
	// add contribution to new source flux
	flux_integral[g] = (q0[g] * tau[g] + (sigT[g] * track->psi[g] - q0[g]) * expVal[g]) / sigT2[g] + q1[g] * mu * (tau[g] * (tau[g] - 2.f) + 2.f * expVal[g]) / (sigT[g] * sigT2[g]) + q2[g] * mu2 * (tau[g] * (tau[g] * (tau[g] - 3.f) + 6.f) - 6.f * expVal[g]) / (3.f * sigT2[g] * sigT2[g]);
	}
	*/

	#pragma simd
	for( int g = 0; g < I.n_egroups; g++)
	{
		// Prepare tally
		tally[g] = weight * flux_integral[g];
	}

	#ifdef OPENMP
	omp_set_lock(QSR->locks + fine_id);
	#endif

	#pragma simd
	for( int g = 0; g < I.n_egroups; g++)
	{
		FSR_flux[g] += tally[g];
	}

	#ifdef OPENMP
	omp_unset_lock(QSR->locks + fine_id);
	#endif


	#pragma simd
	for( int g = 0; g < I.n_egroups; g++)
	{
		t1[g] = q0[g] * expVal[g] / sigT[g]; 
	}
	#pragma simd
	for( int g = 0; g < I.n_egroups; g++)
	{
		t2[g] = q1[g] * mu * (tau[g] - expVal[g]) / sigT2[g]; 
	}
	#pragma simd
	for( int g = 0; g < I.n_egroups; g++)
	{
		t3[g] =	q2[g] * mu2 * (tau[g] * (tau[g] - 2.f) + 2.f * expVal[g]) / (sigT2[g] * sigT[g]);
	}
	#pragma simd
	for( int g = 0; g < I.n_egroups; g++)
	{
		track->psi[g] = track->psi[g] * (1.f - expVal[g]) + t1[g] + t2[g] + t3[g];
	}

	/*
	   #pragma simd
	   for( int g = 0; g < I.n_egroups; g++)
	   {
	// update angular flux
	track->psi[g] = track->psi[g] * (1.f - expVal[g]) + q0[g] * expVal[g] / sigT[g] + q1[g] * mu * (tau[g] - expVal[g]) / sigT2[g] + q2[g] * mu2 * (tau[g] * (tau[g] - 2.f) + 2.f * expVal[g]) / (sigT2[g] * sigT[g]);
	}
	*/
}	

// run one full transport sweep, return k
void transport_sweep( Params params, Input I )
{
	if(I.mype==0) printf("Starting transport sweep ...\n");

	// calculate the height of a node's domain and of each FSR
	double node_delta_z = I.height / I.decomp_assemblies_ax;
	double fine_delta_z = node_delta_z / (I.cai * I.fai);

	// initialize fluxes in transport sweep
	for( int i = 0; i < I.ntracks_2D; i++)
		for( int j = 0; j < I.n_polar_angles; j++)
			for( int k = 0; k < I.z_stacked; k++)
			{
				Track * track = &params.tracks[i][j][k];
				for( int g = 0; g < I.n_egroups; g++)
					track->psi[g] = track->start_flux[g];
			}

	/* loop over tracks (implicitly azimuthal angles, tracks in azimuthal 
	 * angles, polar angles, and z stacked rays) */

	//print_Input_struct( I );

	#pragma omp parallel default(none) \
	shared( I, params, node_delta_z, fine_delta_z ) 
	{
		#ifdef OPENMP
		int thread = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		unsigned int seed = time(NULL) * (thread+1);
		#endif
	//print_Input_struct( I );

		#ifdef PAPI
		int eventset = PAPI_NULL;
		int num_papi_events;
		#pragma omp critical
		{
			counter_init(&eventset, &num_papi_events, I);
		}
		#endif

		AttenuateVars A;
		A.q0 = (float *) malloc( I.n_egroups * sizeof(float));
		A.q1 = (float *) malloc( I.n_egroups * sizeof(float));
		A.q2 = (float *) malloc( I.n_egroups * sizeof(float));
		A.sigT = (float *) malloc( I.n_egroups * sizeof(float));
		A.tau = (float *) malloc( I.n_egroups * sizeof(float));
		A.sigT2 = (float *) malloc( I.n_egroups * sizeof(float));
		A.expVal = (float *) malloc( I.n_egroups * sizeof(float));
		A.flux_integral = (float *) malloc( I.n_egroups * sizeof(float));
		A.tally = (float *) malloc( I.n_egroups * sizeof(float));
		A.t1 = (float *) malloc( I.n_egroups * sizeof(float));
		A.t2 = (float *) malloc( I.n_egroups * sizeof(float));
		A.t3 = (float *) malloc( I.n_egroups * sizeof(float));
		A.f1 = (float *) malloc( I.n_egroups * sizeof(float));
		A.f2 = (float *) malloc( I.n_egroups * sizeof(float));
		A.f3 = (float *) malloc( I.n_egroups * sizeof(float));

		#pragma omp for schedule( dynamic ) 
		for (long i = 0; i < I.ntracks_2D; i++)
		{
			// print progress
			#ifdef OPENMP
			if(I.mype==0 && thread == 0)
			{
				printf("\rAttenuating Tracks... (%.0lf%% completed)",
						(i / ( (double)I.ntracks_2D / (double) nthreads ))
						/ (double) nthreads * 100.0);
			}
			#else
			if( i % 50 == 0)
				if(I.mype==0)
					printf("%s%ld%s%ld\n","2D Tracks Completed = ", i," / ", 
							I.ntracks_2D );
			#endif


			// treat positive-z traveling rays first
			bool pos_z_dir = true;
			for( int j = 0; j < I.n_polar_angles; j++)
			{
				if( j == I.n_polar_angles / 2 )
					pos_z_dir = false;
				float p_angle = params.polar_angles[j];
				float mu = cos(p_angle);

				// start with all z stacked rays
				int begin_stacked = 0;
				int end_stacked = I.z_stacked;

				for( int n = 0; n < params.tracks_2D[i].n_segments; n++)
				{
					// calculate distance traveled in cell if segment completed
					float s_full = params.tracks_2D[i].segments[n].length 
						/ sin(p_angle);

					// allocate varaible for distance traveled in an FSR
					float ds = 0;

					// loop over remaining z-stacked rays
					for( int k = begin_stacked; k < end_stacked; k++)
					{
						// initialize s to full length
						float s = s_full;

						// select current track
						Track * track = &params.tracks[i][j][k];

						// set flag for completeion of segment
						bool seg_complete = false;

						// calculate interval
						int curr_interval;
						if( pos_z_dir)
							curr_interval = get_pos_interval(track->z_height, 
									fine_delta_z);
						else
							curr_interval = get_neg_interval(track->z_height, 
									fine_delta_z);

						while( !seg_complete )
						{
							// flag to reset z position
							bool reset = false;

							/* calculate new height based on s 
							 * (distance traveled in FSR) */
							float z = track->z_height + s * cos(p_angle);

							// check if still in same FSR (fine axial interval)
							int new_interval;
							if( pos_z_dir )
								new_interval = get_pos_interval(z, 
										fine_delta_z);
							else
								new_interval = get_neg_interval(z,
										fine_delta_z);

							if( new_interval == curr_interval )
							{
								seg_complete = true;
								ds = s;
							}

							// otherwise, we need to recalculate distances
							else
							{
								// correct z
								if( pos_z_dir )
								{
									curr_interval++;
									z = fine_delta_z * (float) curr_interval;
								}
								else{
									curr_interval--;
									z = fine_delta_z * (float) curr_interval;
								}

								// calculate distance travelled in FSR (ds)
								ds = (z - track->z_height) / cos(p_angle);

								// update track length remaining
								s -= ds;

								/* check remaining track length to protect
								 * against potential roundoff errors */
								if( s <= 0 )
									seg_complete = true;

								// check if out of bounds or track complete
								if( z <= 0 || z >= node_delta_z )
								{
									// mark segment as completed
									seg_complete = true;

									// remember to no longer treat this track
									if ( pos_z_dir )
										end_stacked--;
									else
										begin_stacked++;

									// reset z height
									reset = true;
								}
							}

							// pick a random FSR (cache miss expected)
							#ifdef OPENMP
							long QSR_id = rand_r(&seed) % 
								I.n_source_regions_per_node;
							#else
							long QSR_id = rand() % 
								I.n_source_regions_per_node;
							#endif

							/* update sources and fluxes from attenuation 
							 * over FSR */
							attenuate_fluxes( track, params.sources +QSR_id, &I, &params, ds, mu, params.tracks_2D[i].az_weight, &A );

							// update with new z height or reset if finished
							if( n == params.tracks_2D[i].n_segments - 1  
									|| reset)
							{
								if( pos_z_dir)
									track->z_height = I.axial_z_sep * k;
								else
									track->z_height = I.axial_z_sep * (k+1);
							}
							else
								track->z_height = z;

						}
					}
				}
			}
		}
		#ifdef OPENMP
		if(thread == 0 && I.mype==0) printf("\n");
		#endif

		#ifdef PAPI
        if( thread == 0 )
        {
            printf("\n");
            border_print();
            center_print("PAPI COUNTER RESULTS", 79);
            border_print();
            printf("Count          \tSmybol      \tDescription\n");
        }
        {
        #pragma omp barrier
        }
        counter_stop(&eventset, num_papi_events, &I);
        #endif
	}

	return;
}

/* returns integer number for axial interval for tracks traveling in the
 *  positive direction */
int get_pos_interval( float z, float dz)
{
	int interval = (int) (z/dz);
	return interval;
}

/* returns integer number for axial interval for tracks traveling in the 
 * negative direction */
int get_neg_interval( float z, float dz)
{
	// NOTE: a bit of trickery using floors to obtain ceils 
	int interval = INT_MAX - (int) ( INT_MAX - ( z / dz ) );
	return interval;
}


void alt_attenuate_fluxes( Track * track, Source * QSR, Input I, 
		Params params, float ds, float mu, float az_weight ) 
{
	// compute fine axial interval spacing
	float dz = I.height / (I.fai * I.decomp_assemblies_ax * I.cai);

	// compute fine axial region ID
	int fine_id = (int) ( I.height / dz ) % I.cai;

	// compute z height in cell
	float zin = track->z_height - dz * ( (int)( track->z_height / dz ) + 0.5 );

	// compute weight (azimuthal * polar)
	// NOTE: real app would also have volume weight component
	float weight = track->p_weight * az_weight;
	float mu2 = mu * mu;

	// load fine source region flux vector
	float * FSR_flux = QSR -> fine_flux[fine_id];

	// cycle over energy groups
	for( int g = 0; g < I.n_egroups; g++)
	{
		// load total cross section
		float sigT = QSR->sigT[g];

		// define source parameters
		float q0, q1, q2;

		// calculate source components
		if( fine_id == 0 )
		{
			// load neighboring sources
			float y2 = QSR->fine_source[fine_id][g];
			float y3 = QSR->fine_source[fine_id+1][g];

			// do linear "fitting"
			float c0 = y2;
			float c1 = (y3 - y2) / dz;

			// calculate q0, q1, q2
			q0 = c0 + c1*zin;
			q1 = c1;
			q2 = 0;
		}
		else if( fine_id == I.fai - 1 )
		{
			// load neighboring sources
			float y1 = QSR->fine_source[fine_id-1][g];
			float y2 = QSR->fine_source[fine_id][g];

			// do linear "fitting"
			float c0 = y2;
			float c1 = (y2 - y1) / dz;

			// calculate q0, q1, q2
			q0 = c0 + c1*zin;
			q1 = c1;
			q2 = 0;
		}		
		else
		{
			// load neighboring sources
			float y1 = QSR->fine_source[fine_id-1][g];
			float y2 = QSR->fine_source[fine_id][g];
			float y3 = QSR->fine_source[fine_id+1][g];

			// do quadratic "fitting"
			float c0 = y2;
			float c1 = (y1 - y3) / (2*dz);
			float c2 = (y1 - 2*y2 + y3) / (2*dz*dz);

			// calculate q0, q1, q2
			q0 = c0 + c1*zin + c2*zin*zin;
			q1 = c1 + 2*c2*zin;
			q2 = c2;
		}

		// calculate common values for efficiency
		float tau = sigT * ds;
		float sigT2 = sigT * sigT;

		// compute exponential ( 1 - exp(-x) ) using table lookup
		float expVal = interpolateTable( params.expTable, tau );  

		// add contribution to new source flux
		float flux_integral = (q0 * tau + (sigT * track->psi[g] - q0) * expVal)
			/ sigT2
			+ q1 * mu * (tau * (tau - 2) + 2 * expVal)
			/ (sigT * sigT2)
			+ q2 * mu2 * (tau * (tau * (tau - 3) + 6) - 6 * expVal)
			/ (3 * sigT2 * sigT2);

		#pragma omp atomic
		FSR_flux[g] += weight * flux_integral;

		// update angular flux
		track->psi[g] = track->psi[g] * (1.0 - expVal) + q0 * expVal / sigT
			+ q1 * mu * (tau - expVal) / sigT2 + q2 * mu2 *
			(tau * (tau - 2) + 2 * expVal) / (sigT2 * sigT);
	}
}

// renormalize flux for next transport sweep iteration
void renormalize_flux( Params params, Input I, CommGrid grid )
{
	if( I.mype == 0 ) printf("Renormalizing Flux...\n");
	// tally total fission rate (pair-wise sum)
	float * fission_rates = malloc( I.n_source_regions_per_node 
			* sizeof(float) );

	float * fine_fission_rates = malloc( I.fai * sizeof(float) );
	float * g_fission_rates = malloc( I.n_egroups * sizeof(float) );

	// accumulate total fission rate on node domain
	for( int i = 0; i < I.n_source_regions_per_node; i++)
	{
		Source src = params.sources[i];
		for( int j = 0; j < I.fai; j++)
		{
			for( int g = 0; g < I.n_egroups; g++)
				g_fission_rates[g] = src.fine_flux[j][g] * src.vol 
					* src.XS[g][0];
			fine_fission_rates[j] = pairwise_sum( g_fission_rates, 
					I.n_egroups );
		}
		fission_rates[i] = pairwise_sum( fine_fission_rates, I.fai );
	}
	float node_fission_rate = pairwise_sum(fission_rates, 
			I.n_source_regions_per_node);

	#ifdef MPI	
	// accumulate total fission rate by MPI Allreduce
	float total_fission_rate = 0;
	MPI_Barrier(grid.cart_comm_3d);
	MPI_Allreduce( &node_fission_rate, // Send Buffer
			&total_fission_rate,    // Receive Buffer
			1,                    	// Element Count
			MPI_FLOAT,           	// Element Type
			MPI_SUM,              	// Reduciton Operation Type
			grid.cart_comm_3d );  	// MPI Communicator
	MPI_Barrier(grid.cart_comm_3d);
	#else
	float total_fission_rate = node_fission_rate;
	#endif

	// free allocated memory
	free(fission_rates);
	free(fine_fission_rates);
	free(g_fission_rates);

	// normalize fluxes by fission reaction rate
	float norm_factor = 1.0 / total_fission_rate;
	for( int i = 0; i < I.n_source_regions_per_node; i++)
	{
		Source * src = &params.sources[i];
		float adjust = norm_factor * 4 * M_PI * I.fai / src->vol;
		for( int k = 0; k < I.fai; k++)
			for( int g = 0; g < I.n_egroups; g++)
				src->fine_flux[k][g] *= adjust;
	}

	// normalize boundary fluxes by same factor
	for( long i = 0; i < I.ntracks_2D; i++)
		for( int j = 0; j < I.n_polar_angles; j++)
			for( int k = 0; k < I.z_stacked; k++)
				for( int g = 0; g < I.n_egroups; g++)
					params.tracks[i][j][k].start_flux[g] *= norm_factor;

	if( I.mype == 0 ) printf("Renormalizing Flux Complete.\n");
	return;
}

// Updates sources for next iteration
float update_sources( Params params, Input I, float keff )
{
	// source residual
	float residual;

	// calculate inverse multiplication facotr for efficiency
	float inverse_k = 1.0 / keff;

	// allocate residual arrays
	float * group_res = (float *) malloc(I.n_egroups * sizeof(float));
	float * fine_res = (float *) malloc(I.fai * sizeof(float));
	float * residuals = (float *) malloc(I.n_source_regions_per_node 
			* sizeof(float));

	// allocate arrays for summation
	float * fission_rates = malloc(I.n_egroups * sizeof(float));
	float * scatter_rates = malloc(I.n_egroups * sizeof(float));

	// cycle through all coarse axial intervals to update source
	for( long i = 0; i < I.n_source_regions_per_node; i++)
	{
		Source src = params.sources[i];

		// cycle thorugh all fine axial regions to calculate new source
		for( int j = 0; j < I.fai; j++)
		{
			// calculate total fission source and scattering source
			float fission_source;
			float scatter_source;

			// compute total fission source
			for( int g = 0; g < I.n_egroups; g++ )
				fission_rates[g] = src.fine_flux[j][g] * src.XS[g][0];
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
						src.fine_flux[j][g2];
				}
				scatter_source = pairwise_sum(scatter_rates, 
						(long) I.n_egroups);

				// compuate new total source
				float chi = src.XS[g][2];

				// calculate new fine source
				float newSrc = (fission_source * chi + scatter_source) 
					/ (4.0 * M_PI);

				// calculate residual
				float oldSrc = src.fine_source[j][g];
				group_res[g] = (newSrc - oldSrc) * (newSrc - oldSrc)
					/ (oldSrc * oldSrc);

				/* calculate new source in fine axial interval assuming 
				 * isotropic source components */
				src.fine_source[j][g] = newSrc;
			}
			fine_res[j] = pairwise_sum(group_res, (long) I.n_egroups);
		}
		residuals[i] = pairwise_sum(fine_res, (long) I.fai);
	}

	// calculate source residual
	residual = pairwise_sum(residuals, I.n_source_regions_per_node);

	// free memory
	free(fission_rates);
	free(scatter_rates);
	free(group_res);
	free(fine_res);
	free(residuals);


	// NOTE: See code around line 600 of CPUSolver.cpp in ClosedMOC/ OpenMOC

	return residual;
}

float compute_keff(Params params, Input I, CommGrid grid)
{
	// allocate temporary memory
	float * sigma = malloc( I.n_egroups * sizeof(float) );
	float * group_rates = malloc( I.n_egroups * sizeof(float) );
	float * fine_rates = malloc( I.fai * sizeof(float) );
	float * QSR_rates = malloc( I.n_source_regions_per_node * sizeof(float) );

	///////////////////////////////////////////////////////////////////////////

	// compute total absorption rate, looping over source regions
	for( long i = 0; i < I.n_source_regions_per_node; i++)
	{
		// load absorption XS data
		Source src = params.sources[i];
		for( int g = 0; g < I.n_egroups; g++)	
			sigma[g] = src.XS[g][1];

		for( int j = 0; j < I.fai; j++ )
		{
			// calculate absorption rates
			float * fine_flux = src.fine_flux[j];
			for( int g = 0; g < I.n_egroups; g++)
				group_rates[g] = sigma[g] * fine_flux[g];

			// sum absorption over all energy groups
			fine_rates[j] = pairwise_sum( group_rates, (long) I.n_egroups );
		}
		// sum absorption over all fine axial intervals
		QSR_rates[i] = pairwise_sum( fine_rates, (long) I.fai );
	}
	// sum absorption over all source regions in a node
	float node_abs = pairwise_sum( QSR_rates, I.n_source_regions_per_node);

	///////////////////////////////////////////////////////////////////////////

	// compute total absorption rate, looping over source regions
	for( long i = 0; i < I.n_source_regions_per_node; i++)
	{
		// load nuSigmaF XS data
		Source src = params.sources[i];
		for( int g = 0; g < I.n_egroups; g++)	
			sigma[g] = src.XS[g][0];

		for( int j = 0; j < I.fai; j++ )
		{
			// calculate absorption rates
			float * fine_flux = src.fine_flux[j];
			for( int g = 0; g < I.n_egroups; g++)
				group_rates[g] = sigma[g] * fine_flux[g];

			// sum fission over all energy groups
			fine_rates[j] = pairwise_sum( group_rates, (long) I.n_egroups );
		}
		// sum fission over all fine axial intervals
		QSR_rates[i] = pairwise_sum( fine_rates, (long) I.fai );
	}
	// sum fission over all source regions in a node
	float node_fission = pairwise_sum( QSR_rates, I.n_source_regions_per_node);

	///////////////////////////////////////////////////////////////////////////

	// MPi Reduction
	float tot_abs = 0;
	float tot_fission = 0;
	float leakage = 0;

	#ifdef MPI

	// Total Absorption Reduction
	MPI_Reduce( &node_abs,    		// Send Buffer
			&tot_abs,      			// Receive Buffer
			1,                  	// Element Count
			MPI_FLOAT,          	// Element Type
			MPI_SUM,            	// Reduciton Operation Type
			0,                  	// Master Rank
			grid.cart_comm_3d );	// MPI Communicator

	// Total Fission Reduction
	MPI_Reduce( &node_fission,     	// Send Buffer
			&tot_fission,  			// Receive Buffer
			1,                    	// Element Count
			MPI_FLOAT,           	// Element Type
			MPI_SUM,              	// Reduciton Operation Type
			0,                    	// Master Rank
			grid.cart_comm_3d );  	// MPI Communicator

	// Total Leakage Reduction
	MPI_Reduce( params.leakage,  	// Send Buffer
			&leakage,      			// Receive Buffer
			1,                    	// Element Count
			MPI_FLOAT,           	// Element Type
			MPI_SUM,              	// Reduciton Operation Type
			0,                    	// Master Rank
			grid.cart_comm_3d );  	// MPI Communicator

	MPI_Barrier(grid.cart_comm_3d);

	// calculate keff
	float keff = tot_fission/ (tot_abs + leakage);
	#else
	float keff = node_fission / (node_abs + *params.leakage);
	#endif

	///////////////////////////////////////////////////////////////////////////

	// free memory
	free(sigma);
	free(group_rates);
	free(fine_rates);
	free(QSR_rates);

	return keff;
}

/* Interpolates a formed exponential table to compute ( 1- exp(-x) )
 *  at the desired x value */
float interpolateTable( Table table, float x)
{
	// check to ensure value is in domain
	if( x > table.maxVal )
		return 1.0;
	else
	{
		int interval = (int) ( x / table.dx + 0.5 * table.dx );
		/*
		if( interval >= table.N || interval < 0)
		{
			printf( "Interval = %d\n", interval);
			printf( "N = %d\n", table.N);
			printf( "x = %f\n", x);
			printf( "dx = %f\n", table.dx);
			exit(1);
		}
		*/
		float slope = table.values[ 2 * interval ];
		float intercept = table.values[ 2 * interval + 1 ];
		float val = slope * x + intercept;
		return val;
	}
}
