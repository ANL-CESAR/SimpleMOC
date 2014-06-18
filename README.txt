===============================================================================
    
              _____ _                 _      __  __  ____   _____ 
             / ____(_)               | |    |  \/  |/ __ \ / ____|
            | (___  _ _ __ ___  _ __ | | ___| \  / | |  | | |     
             \___ \| | '_ ` _ \| '_ \| |/ _ \ |\/| | |  | | |     
             ____) | | | | | | | |_) | |  __/ |  | | |__| | |____ 
            |_____/|_|_| |_| |_| .__/|_|\___|_|  |_|\____/ \_____|
                               | |                                
                               |_|                                
    
===============================================================================
    
                                  Developed at
                   The Massachusetts Institute of Technology
                                      and
                          Argonne National Laboratory
    
===============================================================================
Purpose
===============================================================================

The purpose of this mini-app is to demonstrate the performance
characterterics and viability of the Method of Characteristics (MOC)
for 3D neutron transport calculations in the context of full scale
light water reactor simulation.

===============================================================================
Strawman Reactor Defintion
===============================================================================

For the purposes of simplicity this mini-app uses a conservative "strawman"
reactor model to represent a good target problem for full core reactor
simualations to be run on exascale class supercomputers. Arbitrary
user-defined geometries are not supported.

The reactor model used by this application is detailed as follows:

Fuel Pin:
	- Composed of 10 radial rings, representing materials similar to physical
	  geometry (i.e., fuel, He gap, cladding, water). Slices are taken
	  as volumetrically equal for the purposes of simplicity.
	- This represents an extremely conservative model with hundreds of fuel
	  nuclides in all pins (i.e., not fresh fuel).

Assembly:
	- 17 x 17 lattice cells, each containing a fuel pin surrounded by
	  borated water.

Reactor:
	- The assemblies are arranged in the manner of the MIT BEAVRS reactor
	  benchmark. 
