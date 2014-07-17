#include"SimpleMOC_header.h"


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

	double * flux_array = params.tracks[0][0][0].start_flux;

	// X Positive Direction
	send_idx = 0;
	recv_idx = send_idx + (long) elements * (long) I.n_egroups;
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

	// X Negative Direction
	send_idx = recv_idx + (long) elements * (long) I.n_egroups;
	recv_idx = send_idx + (long) elements * (long) I.n_egroups;
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
	send_idx = recv_idx + (long) elements * (long) I.n_egroups;
	recv_idx = send_idx + (long) elements * (long) I.n_egroups;
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
	send_idx = recv_idx + (long) elements * (long) I.n_egroups;
	recv_idx = send_idx + (long) elements * (long) I.n_egroups;
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
	send_idx = recv_idx + (long) elements * (long) I.n_egroups;
	recv_idx = send_idx + (long) elements * (long) I.n_egroups;
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
	send_idx = recv_idx + (long) elements * (long) I.n_egroups;
	recv_idx = send_idx + (long) elements * (long) I.n_egroups;
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
