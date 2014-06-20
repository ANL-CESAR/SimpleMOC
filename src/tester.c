#include"SimpleMOC_header.h"

// function to test 2D zone finding
// creates points that "draw" the 2D projection of the geometry
// points are output to a text file to be parsed and plotted by
// an external plotter such as Python's matplotlib
void generate_2D_zone_points(Input input, Reactor reactor, int n_pts)
{
	// generate output file
	FILE * out;
	out = fopen("gen_points_2D.txt","w");
	fprintf(out, "X\tY\tZ\tAssembly\tPin\tZone\tIndex");

	// z has range (0,400) - hard coded
	double zmax = 400;

	// calculate x range
	double xmax = input.x_assemblies * reactor.assembly_width;

	// calculate y range
	double ymax = input.y_assemblies * reactor.assembly_width;

	// generate z value which will be common for all points
	double zpt = zmax * urand();

	// define xpt, ypt
	double xpt;
	double ypt;
	RegionID id;
	long index;


	// generate points
	for(int i = 0; i < n_pts; i++)
	{
		xpt = xmax * urand();
		ypt = ymax * urand();
		id = get_region_id(xpt, ypt, zpt, input, reactor);	
		index = get_region_index(id, input, reactor);
		fprintf(out, "%f\t%f\t%f\t", xpt, ypt, zpt);
		fprintf(out, "%ld\t%ld\t%ld\t%ld\n", id.assembly, id.pin, id.zone, index);
	}	

	// close stream	
	fclose(out);
	
	return;
}

// generates a random number between 0 and 1
double urand(void)
{
	return (double)rand() / (double)RAND_MAX;
}
