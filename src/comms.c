#include"SimpleMOC_header.h"

#ifdef MPI
// Faster Transfer information between nodes (angular fluxes)
void fast_transfer_boundary_fluxes( Params params, Input I, CommGrid grid)
{
	if(I.mype==0) printf("Beginning Inter-Node Border Flux Transfer...\n");

	// Note we are introducing a new input restriction for convenience,
	// namely, that the total number of tracks is divisible by 6
	int elements = I.ntracks / 6;

	double h = I.domain_height;
	double x = I.assembly_width;

	long ntracks_per_axial_direction  = I.ntracks * x / (2*x + 4*h);
	long ntracks_per_radial_direction = I.ntracks * h / (2*x + 4*h);

	long send_idx = 0;
	MPI_Status stat;
	MPI_Request *request = (MPI_Request *) malloc( 2 * I.ntracks * sizeof(MPI_Request));

	// Use knowledge of underlying flux structure for efficiency
	double * flux_array = params.tracks[0][0][0].start_flux;

	long dim = (long) elements * (long) I.n_egroups;
	long msg_id = 0;

	/////////////////// Launch All Sends ////////////////////////

	// check if border assembly
	if( grid.x_pos_dest == -1 )
	{
		params.leakage += pairwise_sum( &flux_array[send_idx], dim );
		send_idx += dim;
	}
	else
	{
		// Launch X positive Sends
		for( long i = 0; i < ntracks_per_radial_direction; i++ )
		{
			MPI_Isend(
					&flux_array[send_idx],   // Send Buffer
					1,                       // Number of Elements
					grid.Flux_Array,         // Type of element (all energy group array)
					grid.x_pos_dest,         // Destination MPI rank
					msg_id,                  // Message ID
					grid.cart_comm_3d,       // MPI Communicator
					&request[msg_id] );      // MPI Request (to monitor when call finishes)
			send_idx += (long) I.n_egroups;
			msg_id++;
		}
	}

	// check if border assembly
	if( grid.x_neg_dest == -1 )
	{
		params.leakage += pairwise_sum( &flux_array[send_idx], dim );
		send_idx += dim;
	}
	else
	{
		// Launch X negative Sends
		for( long i = 0; i < ntracks_per_radial_direction; i++ )
		{
			MPI_Isend(
					&flux_array[send_idx],   // Send Buffer
					1,                       // Number of Elements
					grid.Flux_Array,         // Type of element (all energy group array)
					grid.x_neg_dest,         // Destination MPI rank
					msg_id,                  // Message ID
					grid.cart_comm_3d,       // MPI Communicator
					&request[msg_id] );      // MPI Request (to monitor when call finishes)
			send_idx += (long) I.n_egroups;
			msg_id++;
		}
	}

	// check if border assembly
	if( grid.y_pos_dest == -1 )
	{
		params.leakage += pairwise_sum( &flux_array[send_idx], dim );
		send_idx += dim;
	}
	else
	{	
		// Launch Y positive Sends
		for( long i = 0; i < ntracks_per_radial_direction; i++ )
		{
			MPI_Isend(
					&flux_array[send_idx],   // Send Buffer
					1,                       // Number of Elements
					grid.Flux_Array,         // Type of element (all energy group array)
					grid.y_pos_dest,         // Destination MPI rank
					msg_id,                  // Message ID
					grid.cart_comm_3d,       // MPI Communicator
					&request[msg_id] );      // MPI Request (to monitor when call finishes)
			send_idx += (long) I.n_egroups;
			msg_id++;
		}
	}
	// check if border assembly
	if( grid.y_neg_dest == -1 )
	{
		params.leakage += pairwise_sum( &flux_array[send_idx], dim );
		send_idx += dim;
	}
	else
	{	
		// Launch Y negative Sends
		for( long i = 0; i < ntracks_per_radial_direction; i++ )
		{
			MPI_Isend(
					&flux_array[send_idx],   // Send Buffer
					1,                       // Number of Elements
					grid.Flux_Array,         // Type of element (all energy group array)
					grid.y_neg_dest,         // Destination MPI rank
					msg_id,                  // Message ID
					grid.cart_comm_3d,       // MPI Communicator
					&request[msg_id] );      // MPI Request (to monitor when call finishes)
			send_idx += (long) I.n_egroups;
			msg_id++;
		}
	}

	// check if border assembly
	if( grid.z_pos_dest == -1 )
	{
		params.leakage += pairwise_sum( &flux_array[send_idx], dim );
		send_idx += dim;
	}
	else
	{	
		// Launch Z positive Sends
		for( long i = 0; i < ntracks_per_axial_direction; i++ )
		{
			MPI_Isend(
					&flux_array[send_idx],   // Send Buffer
					1,                       // Number of Elements
					grid.Flux_Array,         // Type of element (all energy group array)
					grid.z_pos_dest,         // Destination MPI rank
					msg_id,                  // Message ID
					grid.cart_comm_3d,       // MPI Communicator
					&request[msg_id] );      // MPI Request (to monitor when call finishes)
			send_idx += (long) I.n_egroups;
			msg_id++;
		}
	}

	// check if border assembly
	if( grid.z_neg_dest == -1 )
	{
		params.leakage += pairwise_sum( &flux_array[send_idx], dim );
		send_idx += dim;
	}
	else
	{	
		// Launch Z negative Sends
		for( long i = 0; i < ntracks_per_axial_direction; i++ )
		{
			MPI_Isend(
					&flux_array[send_idx],   // Send Buffer
					1,                       // Number of Elements
					grid.Flux_Array,         // Type of element (all energy group array)
					grid.z_neg_dest,         // Destination MPI rank
					msg_id,                  // Message ID
					grid.cart_comm_3d,       // MPI Communicator
					&request[msg_id] );      // MPI Request (to monitor when call finishes)
			send_idx += (long) I.n_egroups;
			msg_id++;
		}
	}

	/////////////////// Launch All Receives ////////////////////////

	msg_id = 0;
	long recv_idx = send_idx;

	// check if border assembly
	if( grid.x_pos_src == -1)
	{
		for( long i =0; i < dim; i++)
			flux_array[recv_idx + i] = 0;
		recv_idx += dim;
	}
	else
	{
		// Launch X positive Receives
		for( long i = 0; i < ntracks_per_radial_direction; i++ )
		{
			MPI_Irecv(
					&flux_array[recv_idx],   // Recv Buffer
					1,                       // Number of Elements
					grid.Flux_Array,         // Type of element (all energy group array)
					grid.x_pos_src,          // MPI rank to Receive From
					msg_id,                  // Message ID
					grid.cart_comm_3d,       // MPI Communicator
					&request[msg_id] );      // MPI Request (to monitor when call finishes)
			recv_idx += (long) I.n_egroups;
			msg_id++;
		}
	}

	// check if border assembly
	if( grid.x_pos_src == -1)
	{
		for( long i =0; i < dim; i++)
			flux_array[recv_idx + i] = 0;
		recv_idx += dim;
	}
	else
	{
		// Launch X negative Receives
		for( long i = 0; i < ntracks_per_radial_direction; i++ )
		{
			MPI_Irecv(
					&flux_array[recv_idx],   // Recv Buffer
					1,                       // Number of Elements
					grid.Flux_Array,         // Type of element (all energy group array)
					grid.x_neg_src,          // MPI rank to Receive From
					msg_id,                  // Message ID
					grid.cart_comm_3d,       // MPI Communicator
					&request[msg_id] );      // MPI Request (to monitor when call finishes)
			recv_idx += (long) I.n_egroups;
			msg_id++;
		}
	}

	// check if border assembly
	if( grid.x_pos_src == -1)
	{
		for( long i =0; i < dim; i++)
			flux_array[recv_idx + i] = 0;
		recv_idx += dim;
	}
	else
	{
		// Launch Y positive Receives
		for( long i = 0; i < ntracks_per_radial_direction; i++ )
		{
			MPI_Irecv(
					&flux_array[recv_idx],   // Recv Buffer
					1,                       // Number of Elements
					grid.Flux_Array,         // Type of element (all energy group array)
					grid.y_pos_src,          // MPI rank to Receive From
					msg_id,                  // Message ID
					grid.cart_comm_3d,       // MPI Communicator
					&request[msg_id] );      // MPI Request (to monitor when call finishes)
			recv_idx += (long) I.n_egroups;
			msg_id++;
		}
	}

	// check if border assembly
	if( grid.x_pos_src == -1)
	{
		for( long i =0; i < dim; i++)
			flux_array[recv_idx + i] = 0;
		recv_idx += dim;
	}
	else
	{
		// Launch Y negative Receives
		for( long i = 0; i < ntracks_per_radial_direction; i++ )
		{
			MPI_Irecv(
					&flux_array[recv_idx],   // Recv Buffer
					1,                       // Number of Elements
					grid.Flux_Array,         // Type of element (all energy group array)
					grid.y_neg_src,          // MPI rank to Receive From
					msg_id,                  // Message ID
					grid.cart_comm_3d,       // MPI Communicator
					&request[msg_id] );      // MPI Request (to monitor when call finishes)
			recv_idx += (long) I.n_egroups;
			msg_id++;
		}
	}

	// check if border assembly
	if( grid.x_pos_src == -1)
	{
		for( long i =0; i < dim; i++)
			flux_array[recv_idx + i] = 0;
		recv_idx += dim;
	}
	else
	{
		// Launch Z positive Receives
		for( long i = 0; i < ntracks_per_axial_direction; i++ )
		{
			MPI_Irecv(
					&flux_array[recv_idx],   // Recv Buffer
					1,                       // Number of Elements
					grid.Flux_Array,         // Type of element (all energy group array)
					grid.z_pos_src,          // MPI rank to Receive From
					msg_id,                  // Message ID
					grid.cart_comm_3d,       // MPI Communicator
					&request[msg_id] );      // MPI Request (to monitor when call finishes)
			recv_idx += (long) I.n_egroups;
			msg_id++;
		}
	}

	// check if border assembly
	if( grid.x_pos_src == -1)
	{
		for( long i =0; i < dim; i++)
			flux_array[recv_idx + i] = 0;
		recv_idx += dim;
	}
	else
	{
		// Launch Z negative Receives
		for( long i = 0; i < ntracks_per_axial_direction; i++ )
		{
			MPI_Irecv(
					&flux_array[recv_idx],   // Recv Buffer
					1,                       // Number of Elements
					grid.Flux_Array,         // Type of element (all energy group array)
					grid.z_neg_src,          // MPI rank to Receive From
					msg_id,                  // Message ID
					grid.cart_comm_3d,       // MPI Communicator
					&request[msg_id] );      // MPI Request (to monitor when call finishes)
			recv_idx += (long) I.n_egroups;
			msg_id++;
		}
	}

	/////////////////////////////////////////////////////////////////////////

	// Wait for all Communication to Complete
	for( long i = 0; i < I.ntracks * 2; i++ )
		MPI_Wait( &request[i], &stat );

	MPI_Barrier( grid.cart_comm_3d );

	if(I.mype==0) printf("Finished Inter-Node Border Flux Transfer.\n");
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
double * flux_array = params.tracks[0][0][0].start_flux;

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
