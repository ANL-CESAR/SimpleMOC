#include"SimpleMOC_header.h"

int main( int argc, char * argv[] )
{
	int version = 0;
	logo(version);
	srand(time(NULL));

	Input input = get_input();

	Params params = build_tracks( input );

	return 0;
}
