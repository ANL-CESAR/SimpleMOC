#include"Simple_MOC.h"

track* generate_tracks(Input input, Reactor reactor)
{
	int n_azim = input.n_azimuthal / 2;

	// allocate space for variables
	double * phi_eff = (double *) malloc(input.n_azimuthal * sizeof(double));
	double * dx_eff = (double *) malloc(input.n_azimuthal * sizeof(double));
	double * dy_eff = (double *) malloc(input.n_azimuthal * sizeof(double));
	double * ds_eff = (double *) malloc(input.n_azimuthal * sizeof(double));
	double * azim_weight = (double *) malloc(input.n_azimuthal * sizeof(double));
	int * n_x = (int *) malloc(input.n_azimuthal * sizeof(int));
	int * n_y = (int *) malloc(input.n_azimuthal * sizeof(int));
	int * n_tracks = (int *) malloc(input.n_azimuthal * sizeof(int));

	for(int i = 0; i < n_azim; i++)
	{
		// calculate desired azimuthal angle
		double phi = M_PI/n_azim * (i + 0.5);

		// calculate integer number of intersections with x and y axes
		n_x[i] = (int) ( reactor.assembly_width / input.radial_ray_sep 
				* sin(phi) ) + 1;
		n_y[i] = (int) ( reactor.assembly_width / input.radial_ray_sep
				* cos(phi) ) + 1;
		n_tracks[i] = n_x[i] + n_y[i]

		// compute new azimuthal angle such that there are an integer number
		// of radial ray seperations accross the assembly. The new angle will
		// be close to the desired angle

		phi_eff[i] = atan( (double) n_x[i] / (double) n_y[i]);

		// correct angle to correct zone if necessary (i.e. if angle was
		// originally in quadrant II)

		if(phi > M_PI/2)
			phi_eff[i] = M_PI - phi_eff[i];
		


	    // calculate new effective track spacing
		dx_eff[i] = reactor.assembly_width / n_x[i];
		dy_eff[i] = reactor.assembly_width / n_y[i];
		ds_eff[i] = dx_eff[i] * sin(phi_eff[i]);
	}

	// compute azimuthal angle quadrature weights
	double dphi_1, dphi_2;
	for (int i = 0; i < n_azim; i++)
	{
		// compute azimuthal space "owned" by angle in the clockwise direction
		if(i == 1)
			dphi_1 = phi_eff[i];
		else
			dphi_1 = 0.5 * (phi_eff[i] - phi_eff[i-1]);

		// compute azimuthal space "owned" by angle in the 
		// counter clockwise direction
		if(i == n_azim - 1)
			dphi_2 = M_PI - phi_eff[i];
		else
			dphi_2 = 0.5 * (phi_eff[i+1] - phi_eff[i]);

		azim_weight[i] = ds_eff[i] * (dphi_1 + dphi_2) / M_PI;
	}

	// compute track starting and ending points and create tracks
	double x1, x2, y1, y2, S;
	for (int i = 0; i < n_azim; i++)
	{
		// TODO: Determine Z-coordinate start

		// allocate memory for tracks in each azimuthal angle
		track * azim_tracks = (track *) malloc(n_tracks[i] * sizeof(track));
		for (int j = 0; j < n_tracks[i]; j++)
		{
			track new_track;
			azim_tracks[j] = new_track;
		}

		// compute starting points for tracks originating on the x axis
		for (int j = 0; j < n_x[i]; j++)
		{
			x1 = dx_eff[i] * (0.5 + j);
			y1 = 0;
			azim_tracks[j].set_start(x1,y1);
		}

		// compute starting points for tracks on y axis
		for (int j = 0; j < n_y[i]; j++)
		{
			// tracks starting on left edge, pointing diagonally right
			if( phi_eff[i] < M_PI/2)
				x1 = 0;
			// tracks starting on right edge pointing diagonally left
			else
				x1 = reactor.assembly_width;

			y1 = dy_eff[i] * (0.5 + j);
			azim_tracks[j].set_start(x1,y1);
		}

		// for all tracks, compute the ending points
		for (int j = 0; j < n_tracks[i]; j++)
		{
			// try hitting a y axis
			x2 = something;
			y2 = something;
			if(y2 > reactor.assembly_width)
			{
				// try hitting an x axis
				y2 = 

			S = reactor.assembly_width/sin(phi


		// TODO: Generate a synthetic number of segments based off of azimuthal
		// angle and starting location (track length)

		// TODO: Determine a better way of organizing tracks
