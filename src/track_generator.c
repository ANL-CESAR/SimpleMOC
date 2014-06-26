#include"SimpleMOC_header.h"

long calculate_num_2D_tracks( Input input, Reactor reactor )
{
	long tracks = input.n_azimuthal * (reactor.assembly_width / input.radial_ray_sep);
	return tracks;
}

// Generates an array of track structures.
// Each track structure has:
//   - Starting cartesian coord
//   - End cartesian coord
//   - Track length
//   - Azimuthal Angle
//   - Polar Angle
//   - Quadrature Weights
Track * generate_tracks(Input input, Reactor reactor)
{
	int n_azim = input.n_azimuthal / 2;

	// Allocate space for tracks
	//long num_tracks = input.n_azimuthal * input.n_polar_angles;
	//Track * tracks = (Track *) malloc( 2 * num_tracks * sizeof(Track)); 
	//
	
	// We want to determine total number of tracks in the entire system
	

	// allocate space for variables
	double * phi_eff = (double *) malloc(input.n_azimuthal * sizeof(double));
	double * dx_eff = (double *) malloc(input.n_azimuthal * sizeof(double));
	double * dy_eff = (double *) malloc(input.n_azimuthal * sizeof(double));
	double * ds_eff = (double *) malloc(input.n_azimuthal * sizeof(double));
	double * azim_weight = (double *) malloc(input.n_azimuthal * sizeof(double));
	int * n_x = (int *) malloc(input.n_azimuthal * input.n_polar_angles * sizeof(int));
	int * n_y = (int *) malloc(input.n_azimuthal * input.n_polar_angles * sizeof(int));
	int * n_tracks = (int *) malloc(input.n_azimuthal * input.n_polar_angles * sizeof(int));

	// number of axial tracks per polar angle
	int n_z = (int) ( reactor.height / input.axial_z_sep );

	for(int i = 0; i < n_azim; i++)
	{
		// calculate desired azimuthal angle
		double phi = M_PI/n_azim * (i + 0.5);

		// calculate integer number of intersections with x and y axes
		n_x[i] = (int) ( reactor.assembly_width / input.radial_ray_sep 
				* sin(phi) ) + 1;
		n_y[i] = (int) ( reactor.assembly_width / input.radial_ray_sep
				* cos(phi) ) + 1;
		n_tracks[i] = (n_x[i] + n_y[i]) * n_z * input.n_polar_angles;

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
	double x1, x2, y1, y2, z1, z2, S;
	for (int i = 0; i < n_azim; i++)
	{
		// TODO: Determine Z-coordinate start
		z1 = 0;

		// allocate memory for tracks in each azimuthal angle
		track * azim_tracks = (track *) malloc(n_tracks[i] * sizeof(track));

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
			// load in assembly width for reference
			double w = reactor.assembly_width;

			// define all possible intersection points (each surface)
			double * points_x = (double *) malloc(6 * sizeof(double));
			double * points_y = (double *) malloc(6 * sizeof(double));
			double * points_z = (double *) malloc(6 * sizeof(double));

			// left surface (x = 0)
			points_x[0] = 0;
			// delta-y = delta-x * tan(phi)
			points_y[0] = y1 - x1 * tan(phi_eff[i]);
			// delta-z = delta-x / (cos(phi) * tan(theta))
			points_z[0] = z1 - x1 / (cos(phi_eff[i]) * tan(polar[i]));

			// right surface (x = assembly width)
			points_x[1] = w;
			// delta-y = delta-x * tan(phi)
			points_y[1] = y1 + (w - x1) * tan(phi_eff[i]);
			// delta-z = delta-x / (cos(phi) * tan(theta))
			points_z[1] = z1 + (w - x1) / (cos(phi_eff[i]) * tan(polar[i]));
			
			// bottom surface (y = 0)
			points_y[2] = 0;
			// delta-x = delta-y / tan(phi)
			points_x[2] = x1 - y1 / tan(phi_eff[i]);
			// delta-z = delta-y / (sin(phi) * tan(theta))
			points_z[2] = z1 - y1 / (sin(phi_eff[i]) * tan(polar[i]));

			// top surface (y = assembly width)
			points_y[3] = w;
			// delta-x = delta-y / tan(phi)
			points_x[3] = x1 + (w - y1) / tan(phi_eff[i]);
			// delta-z = delta-y / (sin(phi) * tan(theta))
			points_z[3] = z1 + (w - y1) / (sin(phi_eff[i]) * tan(polar[i]));

			// down surface (z = 0)
			points_z[4] = 0;
			// delta-x = delta-z * cos(phi) * tan(theta)
			points_x[4] = x1 - z1 * cos(phi_eff[i]) * tan(polar[i]);
			// delta-y = delta-z * sin(phi) * tan(theta)
			points_y[4] = y1 - z1 * sin(phi_eff[i]) * tan(polar[i]);

			// up surface (z = domain height)
			// TODO: Calculate domain height
			double h = 10;
			points_z[4] = h;
			// delta-x = delta-z * cos(phi) * tan(theta)
			points_x[4] = x1 + (h - z1) * cos(phi_eff[i]) * tan(polar[i]);
			// delta-y = delta-z * sin(phi) * tan(theta)
			points_y[4] = y1 + (h - z1) * sin(phi_eff[i]) * tan(polar[i]);

			for (int k = 0; k < 6; k++)
			{
				// try each plane
				x2 = points_x[k];
				y2 = points_y[k];
				z2 = points_z[k];

				// test to make sure plane is valid
				if( x2 >= 0 && x2 <= w &&
						y2 >= 0 && y2 <= w &&
						z2 >= 0 && z2 <= h)
					break;
				else if (k == 5)
					// could not find plane
					//TODO: report error message
			}


		// TODO: Generate a synthetic number of segments based off of azimuthal
		// angle and starting location (track length)

		// TODO: Determine a better way of organizing tracks
		/* TODO: Find connecting track in other assembly/region: this will
		 * be done by noting which plane we intersect and thus determining 
		 * which domain needs to be passed information. We then correct 
		 * the coordinate so it is the correct coordinates within that
		 * new assembly/region. Ex: if we connect with the assebmly/region
		 * above we would subtract the domain height fromt he z coordinate */
