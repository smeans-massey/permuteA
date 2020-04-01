# permuteA
builds adjacency matrices for networks

Suite of files for demonstrating usage of Permute Assembly method matlab function. You will need matlab to fully utilise these files. Octave may substitute, but no guarantees as these are not tested within that platform.

2020 - Shawn Means
School of Natural and Computational Sciences
Massey University, New Zealand
s.means@massey.ac.nz

These are presented as supplementary material for the article "A Permutation Method for Network Assembly" (in submission, 2020). If you have any problems with these items feel free to contact at email above.

Files include the following items:

File / Description
-----------------------------------------------------------------------
permuteA_beta3.m
			
	Matlab function for assembly of adjacency matrices utilising the 'Permutation Method'

permute_assby_demo_control1p1.m	
	
	Control script for demonstrating usage of above function file; intended to run within a matlab gui environment

permute_assby_demo_command_line_example1.m

	Script demonstrating command line or 'batch' mode execution of above function -- suited for larger network assemblies

permute_assby_sample_weight_figure_generation1.m

	Script for generation of sample weights and estimated metric of sample space for test network realisations. This re-generates Figures 5 and 7 from the article (progression of spectral radius estimate and histogram distribution of final matrix realisations).

Below are sample sequence datasets, utilised in above scripts:

N20_kin_5-19_kout_5-19_rho0.54.mat
N100_kin_20-50_kout_20-50_rho0.5.mat
N1000_kin_251-999_kout_251-999_rho0.46.mat
N2000_kin_1000-1999_kout_1000-1999_rho0.017.mat
N5000_kin_750-2000_kout_750-2000_rho0.mat
N10000_kin_5000-9900_kout_5000-9900_rho0.5.mat
N10000_kin_7000-9900_kout_7000-9900_rho0.5.mat

Below are sample generations of adjacency matrices, produced by the 'command_line_example' above:

permute_batch_example_data_N10000_5000-9900.log
permute_batch_example_data_N10000_5000-9900.mat
permute_batch_example_data_N10000_7000-9900.log
permute_batch_example_data_N10000_7000-9900.mat

Below is set of initial permuted matrix types ('A1') utilised in the 'sample_weight' generation example above:

permute_assby_3x3_A1_types_wselfloop.mat

See each respective .m script file above for more on these datasets.
