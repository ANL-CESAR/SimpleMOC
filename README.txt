===============================================================================
    
              _____ _                 _      __  __  ____   _____ 
             / ____(_)               | |    |  \/  |/ __ \ / ____|
            | (___  _ _ __ ___  _ __ | | ___| \  / | |  | | |     
             \___ \| | '_ ` _ \| '_ \| |/ _ \ |\/| | |  | | |     
             ____) | | | | | | | |_) | |  __/ |  | | |__| | |____ 
            |_____/|_|_| |_| |_| .__/|_|\___|_|  |_|\____/ \_____|
                               | |                                
                               |_|                                
    
==============================================================================
Contact Information
==============================================================================

Organizations:     Computational Reactor Physics Group
                   Massachusetts Institute of Technology

                   Center for Exascale Simulation of Advanced Reactors (CESAR)
                   Argonne National Laboratory

Development Leads: Geoffrey Gunow <geogunow@mit.edu>
                   John Tramm     <jtramm@mcs.anl.gov>
    
===============================================================================
What is SimpleMOC?
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
	- The assemblies are arranged in a rectangular shape, as defined by
	  the user. 

===============================================================================
Input Variables
===============================================================================

x_assemblies
	Number of assemblies in the x-axis of the reactor

y_assemblies
	Number of assemblies in the y-axis of the reactor

cai
	This is the number of course axial intervals

fai
	This is the number of fine axial intervals per course axial interval

axial_exp
	Axial source expansion order

radial_ray_sep
	Radial ray separation

axial_z_sep
	Axial stacked z-ray separation

n_azimuthal
	Number of azimuthal angles

n_polar_angles
	Number of polar angles

n_egroups
	Number of energy groups

decompose
	Turn decomposition on/off

decomp_assemblies_ax
	Number of assemblies per sub-domain (axially)
