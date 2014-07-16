#include"SimpleMOC_header.h"

#ifdef MPI
#include<mpi.h>
#endif

int main( int argc, char * argv[] )
{
	int version = 0;
	int mype = 0;
	int nranks;

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

	if( mype == 0 )
		print_input_summary(input);

	Params params = build_tracks( input );

	double res;
	double keff = 1.0;
	int num_iters = 100;

	for( int i = 0; i < num_iters; i++)
	{
		keff = transport_sweep(params, input);
		renormalize_flux(params,input);
		//res = update_sources(params, input, keff);
	}

	free_2D_tracks( params.tracks_2D );
	free_tracks( params.tracks );
	
	#ifdef MPI
	MPI_Finalize();
	#endif

	return 0;
}
