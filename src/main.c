#include"SimpleMOC_header.h"

int main( int argc, char * argv[] )
{
	int version = 0;
	int mype = 0;
	int nranks;
	//test

	#ifdef MPI
	MPI_Status stat;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif

	if( mype == 0 )
		logo(version);

	srand(time(NULL) * (mype+1));

	Input input = get_input();
	
	CommGrid grid = init_mpi_grid( input );

	if( mype == 0 )
		print_input_summary(input);

	Params params = build_tracks( input );

	double res;
	double keff = 1.0;
	int num_iters = 1;

	for( int i = 0; i < num_iters; i++)
	{
		transport_sweep(params, input);                // Local
		#ifdef MPI
		fast_transfer_boundary_fluxes(params, input, grid); // MPI Caretesian Shift Comms
		#endif
		renormalize_flux(params,input, grid);                // MPI Global Accumulate
		res = update_sources(params, input, keff);     // Local
		keff = compute_keff(params, input, grid);            // MPI Global Accumulate
		if( mype == 0 )
			printf("keff = %lf\n", keff);
	}

	free_2D_tracks( params.tracks_2D );
	free_tracks( params.tracks );
	
	#ifdef MPI
	MPI_Finalize();
	#endif

	return 0;
}
