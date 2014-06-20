#include"SimpleMOC_header.h"

int main( int argc, char * argv[] )
{
	int version = 0;
	logo(version);

	Reactor reactor = reactor_init();
	Input input = get_input();

	// RegionID test
	RegionID id = get_region_id( 1.23, 54.32, 4.9, input, reactor);	
	printf("Region ID:\n\tAssembly: %d\n\tPin Cell: %d\n\tZone: %d\n",
			id.assembly, id.pin, id.zone);

	// Serialized ID test
	long index = get_region_index( id, input, reactor );
	printf("Serialized Index: %ld\n", index);


	return 0;
}
