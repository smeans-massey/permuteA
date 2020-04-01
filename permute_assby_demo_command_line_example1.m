%
%Control script demonstrating command-line usage of permuteA function
%
%2020 - Shawn Means
%School of Natural and Computational Sciences
%Massey University, New Zealand
%s.means@massey.ac.nz
%
%this .m file is part of a .zip archive providing access to functional Matlab implementation of the
%'Permutation Method' as described in "A Permutation Method for Network Assembly" (in submission,
%2020). If you have any problems with these items feel free to contact at email above.

%If you aren't familiar with terminal commands (such as ls, cd, grep) this may not be for
%you! Further, this assumes you're on a macos or linux system. Windows users you're on your own -- sorry!

%This particular script is a batch-control file for loading a sequence of in and out degrees and
%calling the permute assembly function without all the matlab GUI overhead. Further, by redirecting
%output of the command, the method's progression can be more easily observed and analysed for
%possible parameter variations or suitability of user-requested items, e.g., multi-edge proportions
%are just too high for the network, etc.

%note: you must have matlab in your PATH variable for this to work, e.g., the command
% %which matlab
%should give a response. Matlab is aliased in this configuration thus:
% matlab:          aliased to open -a MATLAB_R2017a

%or you can use a full pathway thus, with appended switches for no GUI stuff:
% %/Applications/MATLAB_R2017a.app/bin/matlab -nodisplay -nojvm
%will launch a terminal session of matlab.

%For this example script 'cd' into the directory with these files (unzipped to
%'permute_assby_tutorial'):
% %cd permute_assby_tutorial

%Then the below 'batch' terminal command should work: 
%it invokes matlab and redirects this script to matlab and
%directs matlab output to the file 'outputlog', and runs the process in the background with the '&'
%at the end:
% %/Applications/MATLAB_R2017a.app/bin/matlab -nodisplay -nojvm < ./permute_assby_demo_command_line_example1.m >& outputlog &
%note: you can observe progress of the job at the command line with 
% %tail -f outputlog
%(if you launched the job in the background with the ampersand at the end of above command)

fprintf('\n\tExecuting example script for permute assembly method...\n')

%%matlab section
clear;close all;
load ./N10000_kin_7000-9900_kout_7000-9900_rho0.5.mat;
[A, success] = permuteA_beta3(kin,kout);

%save this result
save('permute_batch_example_data1.mat')
%%

%if all goes well, you should see the above datafile in current directory loadable for your matlab
%pleasure. 

%An example of the above assembly is given in the following .mat & .log file
% that took 1536.432 secs (~26 mins) on a macbook pro.
%datafile: permute_batch_example_data_N10000_7000-9900.mat
%logfile:  permute_batch_example_data_N10000_7000-9900.log

%next example is a bit harder to assemble - there seems to be greater difficulty in assembly of
%sequences with greater ranges over minimum and maximum degree; comment out the above matlab commands and 
%uncomment the below to try this sequence with exact same terminal command as above:
% %/Applications/MATLAB_R2017a.app/bin/matlab -nodisplay -nojvm < ./permute_assby_demo_command_line_example1.m >& outputlog &

%% matlab section
% clear;close all;
% load ./N10000_kin_5000-9900_kout_5000-9900_rho0.5.mat;
% [A2, success2] = permuteA_beta3(kin,kout);
% 
% %save this result
% save('permute_batch_example_data2.mat')
%%

%An example of the above assembly is given in the following .mat & .log file
% that took 2812.661 secs (~47 mins) on a macbook pro.
%datafile: permute_batch_example_data_N10000_5000-9900.mat
%logfile:  permute_batch_example_data_N10000_5000-9900.log

