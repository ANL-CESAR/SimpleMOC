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
