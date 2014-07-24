#include"SimpleMOC_header.h"

#ifdef MPI
// Faster Transfer information between nodes (angular fluxes)
void fast_transfer_boundary_fluxes( Params params, Input I, CommGrid grid)
{
	MPI_Barrier(grid.cart_comm_3d);
	if(I.mype==0) printf("Beginning Inter-Node Border Flux Transfer...\n");

	float h = I.domain_height;
	float x = I.assembly_width;

	// calculate number of tracks to each surface
	long ntracks_per_axial_direction  = I.ntracks * x / (2*x + 4*h);
	long ntracks_per_radial_direction = I.ntracks * h / (2*x + 4*h);

	// correct so that all tracks are used and are symmetric
	long remaining_tracks = I.ntracks - 2 * ntracks_per_axial_direction
	   - 4 * ntracks_per_radial_direction;

	long add_radial = remaining_tracks * ( 4*h / (2*x + 4*h) );
	add_radial = 4 * (add_radial / 4);
	ntracks_per_raidal_direction += add_radial / 4;
	
	long add_axial = I.ntracks - add_radial;
	ntracks_per_axial_direction += add_axial / 2;

	// Calculate all requests needed
	long max_requests = ntracks_per_radial_direction / 100;
	max_requests *= 4;
	max_requests += 2 * (ntracks_per_axial_direction / 100 );

	// One for send , one for received
	max_requests *= 2;

	long send_idx = 0;
	MPI_Status stat;
	MPI_Request *request = (MPI_Request *) malloc( max_requests * 
			sizeof(MPI_Request));

	long n_requests = 0;

	// Use knowledge of underlying flux structure for efficiency
	float * flux_array = params.tracks[0][0][0].start_flux;

	long msg_id = 0;
	long req_id = 0;

	/////////////////// Launch All Sends ////////////////////////

	// make an array of radial receiving sources
	int send_dest[6] = 
	{
		grid.x_pos_dest,
		grid.x_neg_dest,
		grid.y_pos_dest,
		grid.y_neg_dest,
		grid.z_pos_dest,
		grid.z_neg_dest
	};

	// make an array of number of messages
	long num_messages[6] =
	{
		ntracks_per_radial_direction/100,
		ntracks_per_radial_direction/100,
		ntracks_per_radial_direction/100,
		ntracks_per_radial_direction/100,
		ntracks_per_axial_direction/100,
		ntracks_per_axial_direction/100
	};


	// loop over all rectangular surfaces
	for( int i = 0; i < 6; i++ )
	{
		// check if border assembly
		if( send_dest[i] == -1 )
		{
			* params.leakage += pairwise_sum( &flux_array[send_idx],
					num_messages[i] * I.n_egroups * 100 );

			send_idx += num_messages[i] * I.n_egroups * 100;
		}
		else
		{
			// Launch  Sends
			for( long j = 0; j < num_messages[i]; j++ )
			{
				MPI_Isend(
					&flux_array[send_idx],   // Send Buffer
					100,                     // Number of Elements
					grid.Flux_Array,         /* Type of element 
											    (all energy group array) */
					send_dest[i],	         // Destination MPI rank
					msg_id,                  // Message ID
					grid.cart_comm_3d,       // MPI Communicator
					&request[req_id] );      /* MPI Request (to monitor 
												when call finishes) */
				send_idx += (long) I.n_egroups*100;
				msg_id++;
				req_id++;
				n_requests++;
			}
		}
	}

	/////////////////// Launch All Receives ////////////////////////

	msg_id = 0;
	long recv_idx = send_idx;

	// make an array of radial receiving sources
	int rec_sources[6] = 
	{
		grid.x_pos_src,
		grid.x_neg_src,
		grid.y_pos_src,
		grid.y_neg_src,
		grid.z_pos_src,
		grid.z_neg_src
	};

	// loop over all rectangular surfaces
	for( int i = 0; i < 6; i++ )
	{
		// check if border assembly
		if( rec_sources[i] == -1)
		{
			long dim = num_messages[i] * I.n_egroups * 100;
			for( long j =0; j < dim; j++)
				flux_array[recv_idx + j] = 0;
			recv_idx += dim;
		}
		else
		{
			// Launch Receives
			for( long j = 0; j < num_messages[i]; j++ )
			{
				MPI_Irecv(
					&flux_array[recv_idx],   // Recv Buffer
					100,                     // Number of Elements
					grid.Flux_Array,         /* Type of element 
												(all energy group array) */
					rec_sources[i],          // MPI rank to Receive From
					msg_id,                  // Message ID
					grid.cart_comm_3d,       // MPI Communicator
					&request[req_id] );      /* MPI Request (to monitor 
												when call finishes) */

				recv_idx += (long) I.n_egroups*100;
				msg_id++;
				req_id++;
				n_requests++;
			}
		}
	}

	
	/////////////////////////////////////////////////////////////////////////

	// Wait for all Communication to Complete
	for( long i = 0; i < n_requests; i++ )
	{
		MPI_Wait( &request[i], &stat );
	}

	MPI_Barrier( grid.cart_comm_3d );

	if(I.mype==0) printf("Finished Inter-Node Border Flux Transfer.\n");

	free(request);
}
#endif

/*
#ifdef MPI
// Transfer information between nodes (angular fluxes)
// TODO: Need to handle border cases!
void transfer_boundary_fluxes( Params params, Input I, CommGrid grid)
{
if(I.mype==0) printf("Beginning Inter-Node Border Flux Transfer...\n");
// MPI code is written assuming that the final app has a very orderly
// data structure, with all border fluxes organized in an optimal
// contiguous manner, i.e.,:
//
// Track 1: Input Flux 1
//          Input Flux 2
//          ...
//          Input Flux n
// Track 2: Input Flux 1
//          Input Flux 2
//          ...
//          Input Flux n
// ...
// Track M: input Flux 1
//          Input Flux 2
//          ...
//          Input Flux n
// Track 1: Output Flux 1
//          Output Flux 2
//          ...
//          Outpuf Flux n
// Track 2: Output Flux 1
//          ...
//          etc etc
//
// In this manner, we can get by with only 6 large sendrecvs and no
// manual buffer packing. It doesn't actually matter if our app
// isn't quite laid out like this, since it's not a matter of
// performance just how elaborate our initialization code is

// Note we are introducing a new input restriction for convenience,
// namely, that the total number of tracks is divisible by 6
int elements = I.ntracks / 6;
//TODO: Proper checking of this -- right now we aren't actually send all data
//assert( I.ntracks % 6 == 0 );
long send_idx;
long recv_idx;
MPI_Status stat;

// Use knowledge of underlying flux structure for efficiency
float * flux_array = params.tracks[0][0][0].start_flux;

// X Positive Direction
send_idx = 0;
long dim = (long) elements * (long) I.n_egroups;
recv_idx = send_idx + dim;

// TODO: DO THIS FOR ALL TRANSFERS
// FIXME: determine vacant neighbors using -1 or something else 
if( grid.x_pos_dest == -1)
{
leakage += pairwise_sum( &flux_array[send_idx], dim);
// FIXME: Do simple receive
}
else if( grid.x_pos_src == -1 )
{
// FIXME: Do simple send
for( long i = 0; i < dim; i++)
flux_array[recv_idx + i] = 0;
}
else{
MPI_Sendrecv(
&flux_array[send_idx], // send buffer
elements,              // Number of elements to send
grid.Flux_Array,       // Type of element
	grid.x_pos_dest,       // Destination MPI Rank
	0,                     // Send Tag
	&flux_array[recv_idx], // Receive Buffer
	elements,              // Number of elements to receive
	grid.Flux_Array,       // Type of element
	grid.x_pos_src,        // Source MPI Rank
	MPI_ANY_TAG,           // Receive Tag
	grid.cart_comm_3d,     // Communicator (3d cart comm)
	&stat                  // Status variable (non-blocking)
	);
	}

// X Negative Direction
send_idx = recv_idx + dim;
recv_idx = send_idx + dim;
MPI_Sendrecv(
		&flux_array[send_idx], // send buffer
		elements,              // Number of elements to send
		grid.Flux_Array,       // Type of element
		grid.x_neg_dest,       // Destination MPI Rank
		0,                     // Send Tag
		&flux_array[recv_idx], // Receive Buffer
		elements,              // Number of elements to receive
		grid.Flux_Array,       // Type of element
		grid.x_neg_src,        // Source MPI Rank
		MPI_ANY_TAG,           // Receive Tag
		grid.cart_comm_3d,     // Communicator (3d cart comm)
		&stat                  // Status variable (non-blocking)
		);

// Y Positive Direction
send_idx = recv_idx + dim;
recv_idx = send_idx + dim;
MPI_Sendrecv(
		&flux_array[send_idx], // send buffer
		elements,              // Number of elements to send
		grid.Flux_Array,       // Type of element
		grid.y_pos_dest,       // Destination MPI Rank
		0,                     // Send Tag
		&flux_array[recv_idx], // Receive Buffer
		elements,              // Number of elements to receive
		grid.Flux_Array,       // Type of element
		grid.y_pos_src,        // Source MPI Rank
		MPI_ANY_TAG,           // Receive Tag
		grid.cart_comm_3d,     // Communicator (3d cart comm)
		&stat                  // Status variable (non-blocking)
		);

// Y Negative Direction
send_idx = recv_idx + dim;
recv_idx = send_idx + dim;
MPI_Sendrecv(
		&flux_array[send_idx], // send buffer
		elements,              // Number of elements to send
		grid.Flux_Array,       // Type of element
		grid.y_neg_dest,       // Destination MPI Rank
		0,                     // Send Tag
		&flux_array[recv_idx], // Receive Buffer
		elements,              // Number of elements to receive
		grid.Flux_Array,       // Type of element
		grid.y_neg_src,        // Source MPI Rank
		MPI_ANY_TAG,           // Receive Tag
		grid.cart_comm_3d,     // Communicator (3d cart comm)
		&stat                  // Status variable (non-blocking)
		);

// Z Positive Direction
send_idx = recv_idx + dim;
recv_idx = send_idx + dim;
MPI_Sendrecv(
		&flux_array[send_idx], // send buffer
		elements,              // Number of elements to send
		grid.Flux_Array,       // Type of element
		grid.z_pos_dest,       // Destination MPI Rank
		0,                     // Send Tag
		&flux_array[recv_idx], // Receive Buffer
		elements,              // Number of elements to receive
		grid.Flux_Array,       // Type of element
		grid.z_pos_src,        // Source MPI Rank
		MPI_ANY_TAG,           // Receive Tag
		grid.cart_comm_3d,     // Communicator (3d cart comm)
		&stat                  // Status variable (non-blocking)
		);

// Z Negative Direction
send_idx = recv_idx + dim;
recv_idx = send_idx + dim;
MPI_Sendrecv(
		&flux_array[send_idx], // send buffer
		elements,              // Number of elements to send
		grid.Flux_Array,       // Type of element
		grid.z_neg_dest,       // Destination MPI Rank
		0,                     // Send Tag
		&flux_array[recv_idx], // Receive Buffer
		elements,              // Number of elements to receive
		grid.Flux_Array,       // Type of element
		grid.z_neg_src,        // Source MPI Rank
		MPI_ANY_TAG,           // Receive Tag
		grid.cart_comm_3d,     // Communicator (3d cart comm)
		&stat                  // Status variable (non-blocking)
		);

MPI_Barrier( grid.cart_comm_3d );

if(I.mype==0) printf("Finished Inter-Node Border Flux Transfer.\n");
}
#endif
*/
