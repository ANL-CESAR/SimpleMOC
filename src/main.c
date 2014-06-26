#include"SimpleMOC_header.h"

int main( int argc, char * argv[] )
{
	int version = 0;
	logo(version);
	srand(time(NULL));

	Reactor reactor = reactor_init();
	Input input = get_input();

	// RegionID test
	RegionID id = get_region_id( 1.23, 54.32, 4.9, input, reactor);	
	printf("Region ID:\n\tAssembly: %ld\n\tPin Cell: %ld\n\tZone: %ld\n",
			id.assembly, id.pin, id.zone);

	// Serialized ID test
	long index = get_region_index( id, input, reactor );
	printf("Serialized Index: %ld\n", index);


	// Run 2D diagnostic test
	printf("Running 2D diagnostic test\n");
	generate_2D_zone_points(input, reactor, 1000000);

	return 0;
}
