===============================================================================
    
              _____ _                 _      __  __  ____   _____ 
             / ____(_)               | |    |  \/  |/ __ \ / ____|
            | (___  _ _ __ ___  _ __ | | ___| \  / | |  | | |     
             \___ \| | '_ ` _ \| '_ \| |/ _ \ |\/| | |  | | |     
             ____) | | | | | | | |_) | |  __/ |  | | |__| | |____ 
            |_____/|_|_| |_| |_| .__/|_|\___|_|  |_|\____/ \_____|
                               | |                                
                               |_|                                

                                   Version 4
    
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

More information on SimpleMOC can be found in the following publication:

Geoffrey Gunow, John Tramm, Benoit Forget, Kord Smith, and Tim He. SimpleMOC
– A performance abstraction for 3D MOC. In ANS & M&C 2015 - Joint
International Conference on Mathematics and Computation (M&C), Supercomputing
in Nuclear Applications (SNA) and the Monte Carlo (MC) Method, 2015.
http://www.mcs.anl.gov/publication/simplemoc-performance-abstraction-3d-moc

==============================================================================
Quick Start Guide
==============================================================================

Download----------------------------------------------------------------------

	For the most up-to-date version of SimpleMOC, we recommend that you
	download from our git repository. This can be accomplished via
	cloning the repository from the command line, or by downloading a zip
	from our github page.

	Git Repository Clone:
		
		Use the following command to clone SimpleMOC to your machine:

		>$ git clone https://github.com/ANL-CESAR/SimpleMOC.git

		Once cloned, you can update the code to the newest version
		using the following command (when in the SimpleMOC directory):

		>$ git pull

Compilation-------------------------------------------------------------------

	To compile XSBench with default (serial mode) settings, use the following
	command:

	>$ make

	To enable shared memory (OpenMP) or distributed memory (MPI) paralleism,
	set the OpenMP and/or MPI flags to "yes" in the makefile before building.
	See below for more details regarding advanced compilation options.

Running SimpleMOC-------------------------------------------------------------

	To run SimpleMOC with default settings, use the following command:

	>$ ./SimpleMOC

	For non-default settings, SimpleMOC supports the following command line
	options:	

	Usage: ./SimpleMOC <options>
	Options include:
	  -t <threads>     Number of OpenMP threads to run
	  -i <filename>    Input file name with customized parameters.
      -p <PAPI event>  PAPI event name to count (1 only)
      -s               Small problem flag
	Default (no arguments given) runs a small model using baked-in parameters.

==============================================================================
Advanced Compilation, Debugging, Optimization, and Profiling
==============================================================================

There are a number of switches that can be set at the top of the makefile
to enable MPI and OpenMP parallelism, along with more advanced compilation
features.

Here is a sample of the control panel at the top of the makefile:

COMPILER    = gnu
MPI         = no
OPENMP      = no
OPTIMIZE    = yes
DEBUG       = no
PROFILE     = no
PAPI        = no

Explanation of Flags:

COMPILER <gnu, intel, ibm> - This selects your compiler.

MPI - Enables MPI support in the code. Distribution of the reactor domain
      accross ranks should be set in the input file.

OpenMP - Enables OpenMP support in the code. By default, the code will
         run using the maximum number of threads on the system, unless
         otherwise specified with the "-t" command line argument.

OPTIMIZE - Adds compiler optimization flag "-O3".

DEBUG - Adds the compiler flag "-g".

PROFILE - Adds the compiler flag "-pg".

PAPI - Enables PAPI support in the code. You may need to alter the makefile
       or your environment to ensure proper linking with the PAPI library.
       See PAPI section below for more details.

===============================================================================
SimpleMOC Strawman Reactor Defintion
===============================================================================

For the purposes of simplicity this mini-app uses a conservative "strawman"
reactor model to represent a good target problem for full core reactor
simualations to be run on exascale class supercomputers. Arbitrary
user-defined geometries are not supported.

===============================================================================
Input Variables
===============================================================================

By default, the program will run with a default set of inputs
so that around 13 GB of memory is used on a node.

A smaller problem, using only 1 GB of memory, can be run using the "-s"
command line option. This problem is quite small, but may be useful for
error checking on smaller systems with limited memory capacities.

Custom inputs can be defined via input file using the "-i" command line
argument.

A sample input file is given in "default.in". The variables the user can
define are given below:

x_assemblies   : Number of assemblies in the x-axis of the reactor
y_assemblies   : Number of assemblies in the y-axis of the reactor
cai            : Number of coarse axial intervals
fai            : Number of fine axial intervals per coarse  axial interval
axial_exp      : Axial source expansion order
radial_ray_sep :  Radial ray separation
axial_z_sep    :  Axial stacked z-ray separation
n_azimuthal    :  Number of azimuthal angles (should be 32)
n_polar_angles :  Number of polar angles
n_egroups      :  Number of energy groups
decompose      : Turn decomposition on/off (true = 1, flase = 0) 
decomp_assemblies_ax : Number of assemblies per sub-domain (axially) 
segments_per_track : Average number of segments per track
assembly_width : Width of an assembly - 1.26 x 17 cm
height         : Height of the reactor
precision      : precision for source convergence
n_2D_source_regions_per_assembly : 2D src regions per assembly 
papi_event_set : PAPI Event Set Choice 

===============================================================================
Citing SimpleMOC
===============================================================================

Papers that cite SimpleMOC should in general refer to the following
publication:

Geoffrey Gunow, John Tramm, Benoit Forget, Kord Smith, and Tim He. SimpleMOC
– A performance abstraction for 3D MOC. In ANS & M&C 2015 - Joint
International Conference on Mathematics and Computation (M&C), Supercomputing
in Nuclear Applications (SNA) and the Monte Carlo (MC) Method, 2015.
http://www.mcs.anl.gov/publication/simplemoc-performance-abstraction-3d-moc

The bibtext entry for this paper is included below:

@inproceedings{Gunow2015,
author = {Gunow, Geoffrey and Tramm, John and Forget, Benoit and Smith, Kord and He, Tim},
keywords = {high performance computing,method of characteristics},
booktitle = {ANS \& M\&C 2015 - Joint International Conference on Mathematics and Computation (M\&C), Supercomputing in Nuclear Applications (SNA) and the Monte Carlo (MC) Method},
title = {{SimpleMOC} -- A PERFORMANCE ABSTRACTION FOR {3D MOC}},
year = {2015}
}
===============================================================================
