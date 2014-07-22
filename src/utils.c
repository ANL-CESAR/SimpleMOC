#include"SimpleMOC_header.h"

// generates a random number between 0 and 1
double urand(void)
{
	return (double)rand() / (double)RAND_MAX;
}

// generates a random number from a normal distribution
// of given mean and standard deviation (Box-Muller)
double nrand(double mean, double sigma)
{
	// generate two random numbers
	double rand1 = urand();
	double rand2 = urand();

	// use first for "exponential part" and second for
	// "azimuthal part" to create a standard normal random number
	double x = sqrt( -2 * log(rand1) ) * cos( 2 * M_PI * rand2);

	// NOTE: We can generate a second random number if
	// neeeded
	// y = sqrt(-2*log(rand1)) * sin( 2 * M_PI * rand2);
	
	// shift random number to appropriate normal distribution and return
	return x * sigma + mean;
}

// pairwise summation for long arrays
// Note: should be templated, but this is C
double pairwise_sum( double * vector, long size ){
	double sum = 0;

	// Base case: perform summation if size <= 16
	if( size <= 16)
		for( int i = 0; i < size; i++ )
			sum += vector[i];

	else
	{
		// otherwise, split
		sum = pairwise_sum( &vector[0], size/2 ) +
			pairwise_sum( &vector[size/2], size - size/2 );
	}
	
	return sum;
}

// Builds a table of exponential values for linear interpolation
Table buildExponentialTable( double precision, double maxVal )
{
	// define table
	Table table;

	// compute number of arry values
	int N = (int) ( maxVal * sqrt(1.0 / ( 8.0 * precision * 0.01 ) ) ); 

	// compute spacing
	double dx = maxVal / (double) N;

	// allocate an array to store information
	double * tableVals = malloc( 2 * N * sizeof(double) );

	// store linear segment information (slope and y-intercept)
	for( int n = 0; n < N; n++ )
	{
		// compute slope and y-intercept for ( 1 - exp(-x) )
		double exponential = exp( - n * dx );
		tableVals[ 2*n ] = - exponential;
		tableVals[ 2*n + 1 ] = 1 + ( n * dx - 1 ) * exponential;
	}

	// assign data to table
	table.dx = dx;
	table.values = tableVals;
	table.maxVal = maxVal;

	return table;
}

// Interpolates a formed exponential table to compute ( 1- exp(-x) )
// at the desired x value
double interpolateTable( Table table, double x)
{
	// check to ensure value is in domain
	if( x > table.maxVal )
		return 1.0;
	else
	{
		int interval = (int) ( x / table.dx + 0.5 * table.dx );
		double slope = table.values[ 2 * interval ];
		double intercept = table.values[ 2 * interval + 1 ];
		double val = slope * x + intercept;
		return val;
	}
}

