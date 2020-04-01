%
%Control script demonstrating usage of permuteA function
%
%2020 - Shawn Means
%School of Natural and Computational Sciences
%Massey University, New Zealand
%s.means@massey.ac.nz
%
%this .m file is part of a .zip archive providing access to functional Matlab implementation of the
%'Permutation Method' as described in "A Permutation Method for Network Assembly" (in submission,
%2020). If you have any problems with these items feel free to contact at email above.

%
%.zip file should include following demonstration in and out degree sequence datasets;
%Each section should be independently executable in the matlab gui via 'Run Section' command
%
%Note, data filename indicates the size of system (N=...) and range of in- and out- 
%degree sequences (kin_... kout...) and any correlation between them (rho...)
%of course, you are quite welcome to try any in- and out- sequences you generate or import
%but the vectors kin and kout passed to the assembly function need to be column vectors of length N:

%load ./N100_kin_20-50_kout_20-50_rho0.5.mat;
%load ./N1000_kin_251-999_kout_251-999_rho0.46.mat;
%load ./N2000_kin_1000-1999_kout_1000-1999_rho0.017.mat;
%load ./N5000_kin_750-2000_kout_750-2000_rho0.mat;

%below sequences take longer...
%load ./N10000_kin_7000-9900_kout_7000-9900_rho0.5.mat;
%load ./N10000_kin_5000-9900_kout_5000-9900_rho0.5.mat;

%Each section is 'clickable' in the matlab gui; select the section you want to try and click 'Run
%Section' in the Matlab toobar.

%%
%Simple Network N=100 Example  (Click here and 'Run Section' to try it out)
%
clear;close all;
load ./N100_kin_20-50_kout_20-50_rho0.5.mat;

%call the actual function. By default it builds a 'simple' network: no multiple edges or self-loops.
A = permuteA_beta3(kin,kout);

%this matrix should still satisfy both in- and out- degree constraints -- if the assembly was
%successful
fprintf('\n\tTesting "A" for conformation to degree sequences...\n')
nnz_kin = nnz(sum(A,2) - kin)
nnz_kout = nnz(sum(A,1) - kout')

%%


%%
%Multi-Edge Example - N=100 (Click here and 'Run Section' to try it out)
%
%if multi-edges are desired at a specific fraction, pass this option through either "long" option
%string 'multi_edge_target_proportion' or abbreviated 'targ_prop', here for 10% multi-edges:
clear;close all;
load ./N100_kin_20-50_kout_20-50_rho0.5.mat;

A = permuteA_beta3(kin,kout,'multi_edge_target_proportion',0.1);

%this matrix should roughly satisfy the 10% target proportion of multi-edges
fprintf('\n\tTesting "A" for target proportion multi-edges...\n')
num_multi_edges = sum(A(A>1)-1);
tot_edges = sum(sum(A));
pct_multi_edges = num_multi_edges/tot_edges

%%

%%
%Multi-Edge Example - very high target percentage
% 
%very high proportion requested may not be obtainable
clear;close all;
load ./N100_kin_20-50_kout_20-50_rho0.5.mat;
A = permuteA_beta3(kin,kout,'targ_prop',0.999);

%the above matrix may not be a successful assembly:
fprintf('\n\tTesting very high multi-edge target proportion "A"...\n')
nnz_kin = nnz(sum(A,2) - kin)
nnz_kout = nnz(sum(A,1) - kout')

%and may not hit the requested target percentage:
num_multi_edges = sum(A(A>1)-1);
tot_edges = sum(sum(A));
pct_multi_edges = num_multi_edges/tot_edges

%%
%%Multi-Edge Example - capturing return success flag
%
%Instead of testing the resulting A for successfull assembly, capture status with additional output variables:
clear;close all;
load ./N100_kin_20-50_kout_20-50_rho0.5.mat;
[A, success] = permuteA_beta3(kin,kout,'targ_prop',0.999);

if success
    fprintf('\n\tAssembly success!\n')
else
    fprintf('\n\tAssembly failure!\n')
end

fprintf('\n\tTesting failed/successful "A"...\n')
nnz_kin = nnz(sum(A,2) - kin)
nnz_kout = nnz(sum(A,1) - kout')

%Note, although we may obtain a successful assembly, where A satisfies both kin and kout, it may not
%present the targeted percentage of multi-edges:
num_multi_edges = sum(A(A>1)-1);
tot_edges = sum(sum(A));
pct_multi_edges = num_multi_edges/tot_edges

%%

%%
%Multi-Edge Example - larger N gives better chance hitting target percentage
%(or just higher density network...). Note, such an extreme target here of 99.9% is likely to
%fail regardless. This may take several attempts...
clear;close all;
load ./N1000_kin_251-999_kout_251-999_rho0.46.mat;
[A, success] = permuteA_beta3(kin,kout,'targ_prop',0.999);
if success
    fprintf('\n\tAssembly success!\n')
else
    fprintf('\n\tAssembly failure!\n')
end

fprintf('\n\tTesting failed/successful "A"...\n')
nnz_kin = nnz(sum(A,2) - kin)
nnz_kout = nnz(sum(A,1) - kout')

num_multi_edges = sum(A(A>1)-1);
tot_edges = sum(sum(A));
pct_multi_edges = num_multi_edges/tot_edges
%%

%%
%Sample Weights Example - N=100
%
%finally, if sampling weights for returned A are needed, pass the following option
%'calculate_sampling_weights' or the abbreviated 'weights' as a '1' to enable, and don't forget to
%include additional output variable to capture the data
%Note, this will slow down the procedure! Large networks will fast become unwieldy.
%i.e., the N=100 sample sequence:
clear;close all;
load ./N100_kin_20-50_kout_20-50_rho0.5.mat;
[A, success, weight] = permuteA_beta3(kin,kout,'calculate_sampling_weights',1);

if success
    fprintf('\n\tSampling weight: %d...\n',weight)
end

%for even a relatively small network of N=100, the sampling space can be enormous and the weights grow
%quite fast. Next example shows how to capture internal data structures from the function to observe
%this.


%%
%Sample Weights Example - Internal Data Structures 
%we can capture all the internal data structs of the permuteA function with one more output
%variable:
clear;close all;
load ./N100_kin_20-50_kout_20-50_rho0.5.mat;
[A, success, weight, internal_data] = permuteA_beta3(kin,kout,'weights',1);

%internal_data contains numerous structures and sub-structures (mostly for debugging and such). 
%Of interest for the sample weights are the fields:
%internal_data.open_edge_data.weight_num_open_rows and 
%internal_data.open_edge_data.weight_num_open_cols
%these indicate the possible row and column selections available at every edge swap
%The final weights are computed thus:
%weight = Product[weight_num_open_rows.*weight_num_open_cols]
%which may be too large a number for the range of available integers in matlab; observe
prod_cols_rows = internal_data.open_edge_data.weight_num_open_cols.*...
    internal_data.open_edge_data.weight_num_open_rows;
for cur_prod = 1:length(prod_cols_rows)
    prod_trunc(cur_prod) = prod(prod_cols_rows(1:cur_prod));
end

%and for the N=100 network the size of these weights gets enormous - order 10^300 
semilogy(prod_trunc,'k+')

%%

%%
%Sample Weights Example - Smaller sequence
%hence, try the following sequence instead
clear;close all;
load ./N20_kin_5-19_kout_5-19_rho0.54.mat;
[A, success, weight, internal_data] = permuteA_beta3(kin,kout,'weights',1);

%the above small network exhibits a high failure rate -- only so many edges to work with.
%%

%%
%Sample Weights Example - Changing internal parameters

%The default parameters in the function may not be suitable for this network. If you haven't
%already, open the function 'permuteA_beta3.m' and search for the string
%'DEFAULT_MAX_DONOR_LOOPS'. It's default value is simply 2*N, or here 40 loops.
%Try increasing this to 4*N, say, and attempt another assembly:

%DEFAULT_MAX_DONOR_LOOPS = 4*N should prove more successful with the N=20 sequence
clear;close all;
load ./N20_kin_5-19_kout_5-19_rho0.54.mat;
[A, success, weight, internal_data] = permuteA_beta3(kin,kout,'calculate_sampling_weights',1);

prod_cols_rows = internal_data.open_edge_data.weight_num_open_cols.*...
    internal_data.open_edge_data.weight_num_open_rows;
for cur_prod = 1:length(prod_cols_rows)
    prod_trunc(cur_prod) = prod(prod_cols_rows(1:cur_prod));
end

%for the N=20 network, the weights still grow quickly but are somewhat more manageable than the
%N=100 case
semilogy(prod_trunc,'k+')

%this is partly why the companion paper for this method presents results for a N=3 system -- much
%more manageable!  The script permute_assby_sample_weight_figure_generation1.m demonstrates 
%how the paper produced the figures showing weighted average estimates of spectral radius for the given
%system. 

%%

%%
%More Internal Parameters and Run-time info
%

%Two other parameters in the function can dramatically affect performance of assemblies:
%DEFAULT_INERT_SWAP_ACTIVATE_RATIO and DEFAULT_OPEN_EDGE_SEARCH_SWITCH_THRESH. Changing these values
%affects how often inert swaps are performed (when the code activates inert swapping) and when to
%change methods for finding suitable edges to swap.  Search for these if you are interested and try
%manipulating these values -- your mileage may vary.
%
%Additionally, you can see much more output from the function as it runs by setting the variable
%'verbosity_level' (at the top of the function) to '1' to enable progression displays. 
%Note, this will display significant amounts of information and may be more suitable for command-line usage. 
%As an example observe:

%with verbosity_level = 1
%
clear;close all;
load ./N2000_kin_1000-1999_kout_1000-1999_rho0.017.mat;
[A, success] = permuteA_beta3(kin,kout);

%running this in the Matlab gui will dump a ton of printouts into the command window, such as
%progression of the method (Distance from target kout), number of swaps performed (out of total
%possible), whether inert shuffling occurred, time elapsed (this round, and overall), etc.

%For large networks, running in a batch mode from a terminal command line (not within matlab GUI) is
%advised. See examples of batch commands from a terminal command line in the file 
%permute_assby_demo_command_line_example1.m



