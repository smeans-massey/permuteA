function [A_out, varargout] = permuteA_beta3(kin,kout,varargin)
%[A, <opt success?>, <opt sampling weight> <opt debug data>] = permuteA_beta3(kin,kout,varargin)
%
%Function builds adjacency network matrix 'A' from input of degree sequences 'kin' and 'kout' - 
%both which must be sized (1,N).
%
%By default, A is built without multi-edges or self-loops - a 'simple' graph' - 
%if you desire otherwise see below for additional options.
%
%Returned value is adjacency matrix 'A' whether successful or a failed effort. Additional returned
%variables provided in this order if requested at function call:
%
%   success_flag (0 fail / 1 success), sampling_weight, assembly_data
%
%assembly_data includes a lot of structures internal to this function for debugging, etc.
%
%This function performs tests on passed kin, kout sequences ensuring graphicality, e.g., 
%
%   (1) sum(kin) = sum(kout) 
%
%   (2) aka the Gale-Ryser Inequalities
%
%       Sum_i=1...n [min(kout(i),j]  >=   Sum_i=1...j [kin(i)]
%
%       for all j in [1,...,n-1]
%
%additional options may be passed including the following:
%   Text String (abbrev)            Value(s)    Description
%   multi_edge_target_proportion    0 - 1.0     Proportion of multi-edges in matrix (Default: 0 or none)
%       (targ_prop)
%   auto_connects_permitted         0 / 1       Disallow / Permit self-loops (Default: Disallow)
%       (self_loops)
%   calculate_sampling_weights      0 / 1       Disable / Enable calculation sample weight for
%       (weights)                               returned A (Default: Disabled)
%   user_a1                         Matrix      Initial permuted 'A1' provided by user
%       (a1)                                    matrix must be sized NxN and satisfy sum(A1,2) == kin
%                                               (Default: Assemble Internally)
%
%Note - additional options internal to this function in subroutine permA_init_structs; alter them at
%your own risk.
%
%2020 - Shawn Means 
%School of Natural and Computational Sciences
%Massey University, New Zealand
%s.means@massey.ac.nz
%
 

%%%%%%%%%%%%%%%%%%
%global vars
%%%%%%%%%%%%%%%%%%

global N debug_level verbosity_level

debug_level = 0;
%set below to '0' for silent operation
verbosity_level = 0;

if verbosity_level
    fprintf('\n\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n\n')
    fprintf('permute A Adjacency Matrix Assembly called')
    fprintf('\n\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n\n')
    
end

%%%%%%%%%%%%%%%%%%
%parse arguments
%%%%%%%%%%%%%%%%%%

[parse_flag,parse_msg,passed_parms_data] = permA_parse_args(kin,kout,nargin,...
    varargin,nargout);

%%%%%%%%%%%%%%%%%%
%initialise structures
%%%%%%%%%%%%%%%%%%
performance.overall_start = tic;

%call initialisation function - passing the parsed options from user
[init_flag, init_msg, matrix_data, multi_edge_data, auto_edge_data, assby_data, open_edge_data, inert_swap_data, performance] = ...
    permA_init_structs(kin,kout,passed_parms_data,performance);

%assby_data

A_out = [];
data_out = [];

%%%%%%%%%%%%%%%%%%
% error tests
%%%%%%%%%%%%%%%%%%

%ensure kin and kout col vectors
if ~iscolumn(kin) && ~iscolumn(kout)
    error('Input degree sequences must be column vectors\n')
end

%get and check size of kin kout
if matrix_data.Nkin ~= matrix_data.Nkout
    error('Input degree sequences must be same length\n');
else
    N = matrix_data.Nkin;
end

if verbosity_level
    fprintf('\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n')
    fprintf('\n\tSize of network N=%d\n',N);
end

%call tests for graphicality of kin & kout

if verbosity_level
    fprintf('\tGraphicality tests...')
end

[graphic_test_flag, graphic_test_msg] = permA_test_kin_kout(kin,kout,matrix_data);

if ~graphic_test_flag
    error(graphic_test_msg)
else
    
    if verbosity_level

        fprintf('passed\n\n')
    end
end


%%%%%%%%%%%%%%%%%%
%options header
%%%%%%%%%%%%%%%%%%
if verbosity_level

    fprintf('\tAssembling network with following options:\n\n');


    %if multi-edges flagged, display, else print default setting
    if multi_edge_data.flag
        fprintf('\t\t** Multi-edges requested\n\t\t\ttarget proportion: %3.2f\n',...
            multi_edge_data.targ_prop);
        fprintf('\t\t\tacceptable target deviation: %3.2f\n',...
            multi_edge_data.epsilon);

        %if other options passed display passed value or default
        if passed_parms_data.opts_flag(2)
           fprintf('\t\t\tMax attempts target pct: %d\n',...
               passed_parms_data.val{2});
        %else

        end


    else
        fprintf('\t\t** Multi-edges disallowed\n');

    end

    %if auto-edges flagged, display, else print default
    if auto_edge_data.flag
        fprintf('\t\t** Auto-edges/self-loops permitted\n');

    else
        fprintf('\t\t** Auto-edges/self-loops disallowed\n');

    end
    
    %if user requested weight calculations and return 
    if open_edge_data.calc_weights
       fprintf('\t\t** Sampling Weight Calculation requested\n') 
        
    end

end

%if user-supplied A1 matrix provided, notify
if matrix_data.user_A1_flag
    
    if verbosity_level
    
        fprintf('\n\t\t** User-supplied initial permutation matrix, A1, provided...')
    end
    
    %test for suitable size of A1
    %probably want to test for NAN's, real-valued, integers...
    if size(matrix_data.A1,1) == N && size(matrix_data.A1,2) == N
        if verbosity_level

            fprintf('proceeding\n')
        end
        
    else
        error('A1 must be sized (%d,%d)',N,N)
        
    end
    
%else
    
    
end

if verbosity_level

    fprintf('\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
end





%%%%%%%%%%%%%%%%%%
% preallocation
%%%%%%%%%%%%%%%%%%

matrix_data.A0 = zeros(N,N);    %initial pre-permuted A, violating kout

%if user-supplied A1, don't zero it out...
if ~matrix_data.user_A1_flag
    matrix_data.A1 = zeros(N,N);    %secondary, permuted A, violating kout
end
matrix_data.Ai = zeros(N,N);    %interim A, gradually shifted to kout
%matrix_data.A02 = zeros(N,N);
matrix_data.node_idx = [1:N];
% matrix_data.node_edges = cell(1,N); 


%%%%%%%%%%%%%%%%%%
% random seeding
%%%%%%%%%%%%%%%%%%

%default behaviour: seed from time for unique sequences
rng('shuffle') 

%%%%%%%%%%%%%%%%%%
% A0 assembly
%%%%%%%%%%%%%%%%%%

%fprintf('\n\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n\n')

%if A1 provided, skip A0 assembly
if ~matrix_data.user_A1_flag
            
    if verbosity_level

        fprintf('\n\n\tAssembling "A0" initial matrix...\n');
    end

    %if multi-edges flagged, we make max_attempts at hitting target percentage, adjusting the
    %distribution of multi-edges for selection until we reach within epsilon

    if multi_edge_data.flag

        if debug_level >= 1

            fprintf('\t\tMulti-Edges Requested------\n\t\t\tTarget Percentage: %2.3f\n',multi_edge_data.targ_prop);
            fprintf('\t\t\tmaximum attempts per node: %d\n\t\t\tacceptable deviation: %2.3f\n',...
                multi_edge_data.max_pick_tries,multi_edge_data.epsilon);
            fprintf('\t\t\tmulti-edge distn parm, lambda: %2.3f\n',multi_edge_data.lambda);

        end

    end

    performance.A0_start = tic;

    %assemble options struct here...

    num_A0_attempts = 1;
    %targ_pct_hit = 0;
    exit_A0_build = 0;

    %initialise A0 history struct (may move this to init function)
    A0_hist.A0 = zeros(N,N); %for saving previous A0 in case of error / rejection
    A0_hist.kin = zeros(N,1);
    A0_hist.A0_kin_diff = zeros(N,1);


    %loop until we hit target percentage or exceed max attempts
    while ~exit_A0_build && (num_A0_attempts <= multi_edge_data.max_pct_attempts)

        if verbosity_level

            fprintf('\n\tA0 assembly attempt #%d/%d...\n',num_A0_attempts,multi_edge_data.max_pct_attempts)
        end

        performance.A0_attempt_start = tic;

        if debug_level >= 1
            fprintf('\n\tnnz(A0): %d\n',nnz(matrix_data.A0));
        end

        %call internal A0 assembler
        [A0_flag, A0_msg, matrix_data, multi_edge_data] = ...
            permA_assemble_A0(kin,kout,matrix_data,multi_edge_data,auto_edge_data);


        %test for target percentage hit -- if multi-edges requested
        if multi_edge_data.flag

            %compare resulting multi-edge proportion with target, if within epsilon flag success
            %cur_multi_pct = matrix_data.A0_num_multi/matrix_data.sum_kin;

            if verbosity_level

                fprintf('\n\t\tCurrent returned multi-edge pct: %2.3f...\n',multi_edge_data.A0_pct_multi)
            end    

            multi_edge_data.multi_err_hist(num_A0_attempts) = abs(multi_edge_data.A0_pct_multi - multi_edge_data.targ_prop);
    %         cur_err = abs(multi_edge_data.A0_pct_multi - multi_edge_data.targ_prop);
    %         multi_edge_err_hist(num_A0_attempts) = cur_err;

            if debug_level >= 1
               fprintf('\t\t\t(deviation from target: %e)\n',multi_edge_data.multi_err_hist(num_A0_attempts))

            end

            %flags for error / error types / exit condition testing
            error_exit_type = 0;


            if multi_edge_data.multi_err_hist(num_A0_attempts) <= multi_edge_data.epsilon
    %             fprintf('\t\tMulti-edge proportion within tolerance!\n\t\t\taccepting "A0"...\n',multi_edge_data.A0_pct_multi)
                error_exit_type = 1; 
                %targ_pct_hit = 1;

            elseif (num_A0_attempts>1) &&  ...
                    (multi_edge_data.multi_err_hist(num_A0_attempts) > multi_edge_data.multi_err_hist(num_A0_attempts-1)) 
                error_exit_type = 2;
    %             fprintf('\t\tMulti-edge error deterioriating; previous: %2.3f, current: %2.3f...\n',...
    %                 multi_edge_data.multi_err_hist(num_A0_attempts),multi_edge_data.multi_err_hist(num_A0_attempts-1))

            elseif (num_A0_attempts>1) && (matrix_data.A0_kin_errnum > 0)
                error_exit_type = 3;
    %             fprintf('\t\tA0 build violates prescribed kin (%d nodes incorrect degree); rejecting...\n',matrix_data.A0_kin_errnum)

            elseif (num_A0_attempts + 1) > multi_edge_data.max_pct_attempts
                error_exit_type = 4;
    %             fprintf('\t\tMaximum attempts (%d) reached!\n\t\t\texiting with current "A0"...\n',multi_edge_data.max_pct_attempts)
    %             fprintf('\t\t\tIncrease "max_pct_attempts" to improve result...\n')
            %else
                  %This state handler shifted to switch-case = 0 below
    %             %we failed to hit target percentage within epsilon; only option is adjust the
    %             %exponential distribution for multi-edges via lambda:
    %             multi_edge_data.lambda = multi_edge_data.lambda*2;
    %             fprintf('\t\tIncreasing exponential distribution parameter to %2.3f...\n',multi_edge_data.lambda)
    %             
    %             %clear appended / modified fields to matrix_data and multi_edge_data from previous
    %             %effort
    % %             rmfield(multi_edge_data,{'kin_mult','skip_multi','A0_pct_multi'});
    % %             rmfield(matrix_data,{'A0','A0_kin','A0_kin_diff','A0_kin_errnum','A0_num_multi'});
    %             
    %             multi_edge_data.kin_mult = [];
    %             multi_edge_data.skip_multi = [];
    %             multi_edge_data.A0_pct_multi = [];
    %             
    %             matrix_data.A0 = zeros(N,N);
    %             %matrix_data.A02 = zeros(N,N);
    %             matrix_data.A0_kin = [];
    %             matrix_data.A0_kin_diff = [];
    %             matrix_data.A0_kin_errnum = [];
    %             matrix_data.A0_num_multi = [];


            end

            %depending on error_exit type, determine action
            switch error_exit_type

                %continue A0 attempt
                case 0
                    %we failed to hit target percentage within epsilon; only option is adjust the
                    %exponential distribution for multi-edges via lambda:
                    multi_edge_data.lambda = multi_edge_data.lambda*2;

                    if verbosity_level

                        fprintf('\t\tIncreasing multi-edge distribution parameter to %2.3f...\n',multi_edge_data.lambda)
                    end

                    %save the previous A0 build 
                    A0_hist.A0 = matrix_data.A0;
                    A0_hist.A0_kin = matrix_data.A0_kin;
                    A0_hist.A0_kin_diff = matrix_data.A0_kin_diff;
                    A0_hist.A0_kin_errnum = matrix_data.A0_kin_errnum;
                    A0_hist.A0_num_multi = matrix_data.A0_num_multi;

                    A0_hist.A0_pct_multi = multi_edge_data.A0_pct_multi;


                    %clear appended / modified fields to matrix_data and multi_edge_data from previous
                    %effort
                    multi_edge_data.kin_mult = [];
                    multi_edge_data.skip_multi = [];
                    multi_edge_data.A0_pct_multi = [];

                    matrix_data.A0 = zeros(N,N);
                    %matrix_data.A02 = zeros(N,N);
                    matrix_data.A0_kin = [];
                    matrix_data.A0_kin_diff = [];
                    matrix_data.A0_kin_errnum = [];
                    matrix_data.A0_num_multi = [];


                %success!     
                case 1
                    if verbosity_level

                        fprintf('\t\tMulti-edge proportion (%2.3f) within tolerance!\n\t\t\taccepting "A0"...\n',multi_edge_data.A0_pct_multi)
                    end
                    exit_A0_build = 1;

                %error type: edge proportion deviation deterioriating
                case 2
                    if verbosity_level

                        fprintf('\t\tMulti-edge error deterioriating; previous: %2.3f, current: %2.3f...\n',...
                            multi_edge_data.multi_err_hist(num_A0_attempts),multi_edge_data.multi_err_hist(num_A0_attempts-1))
                    end        

                    exit_A0_build = 1;

                %error type: catastrophic; A0 violates kin    
                case 3
                    if verbosity_level

                        fprintf('\t\tA0 build violates prescribed kin (%d nodes incorrect degree); rejecting...\n',matrix_data.A0_kin_errnum)
                    end

                    exit_A0_build = 1;


                %maximum attempts hit; exit
                case 4
                    if verbosity_level

                        fprintf('\t\tMaximum attempts (%d) reached!\n\t\t\texiting with current "A0"...\n',multi_edge_data.max_pct_attempts)
                    end       
                    exit_A0_build = 1;

            end

        else
            %simply ensure we have _no_ multi-edges since not requested
            if matrix_data.A0_num_multi ~=0 
                error('\n\tA0 constructor returned with multi-edges when not requested...exiting...\n')
            else
                %flags for error / error types / exit condition testing - this is NULL for current
                %state, but set to '0'
                error_exit_type = 0;

                exit_A0_build = 1;

            end

        end

        performance.A0_attempt_finish(num_A0_attempts) = toc(performance.A0_attempt_start);
        if verbosity_level

            fprintf('\n\t\t%2.3fs elapsed this attempt\n',performance.A0_attempt_finish(num_A0_attempts))
        end

        num_A0_attempts = num_A0_attempts + 1;


    end


    %if last A0 rejected, restore previous A0 and display statistics 
    if error_exit_type == 2 || error_exit_type == 3
        %restore A0 from previous attempt

        if verbosity_level

            fprintf('\n\tRestoring A0 build from previous attempt:\n')
        end

        %matrix data struct
        matrix_data.A0 = A0_hist.A0;
        matrix_data.A0_kin = A0_hist.A0_kin;
        matrix_data.A0_kin_diff = A0_hist.A0_kin_diff;
        matrix_data.A0_kin_errnum = A0_hist.A0_kin_errnum;
        matrix_data.A0_num_multi = A0_hist.A0_num_multi;

        %multi_edge data struct
        multi_edge_data.A0_pct_multi = A0_hist.A0_pct_multi;
        multi_edge_data.multi_err_hist = multi_edge_data.multi_err_hist(1:end-1);


        if verbosity_level

            fprintf('\n\t\tFinal multi-edge pct: %2.3f...\n',multi_edge_data.A0_pct_multi)
            fprintf('\t\t\t(deviation from target: %e)\n',multi_edge_data.multi_err_hist(end))

        end    
    end

    performance.A0_finish = toc(performance.A0_start);

    if verbosity_level

        fprintf('\n\tA0 Assembly complete in %2.3fs',performance.A0_finish)
        fprintf('\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n')
    end

else
    if verbosity_level
    
        fprintf('\n\t\t** User-supplied initial permutation matrix, A1, provided; skipping A0 generation...')
    end
    
    
end

%%%%%%%%%%%%%%%%%%
% A1 permutation
%%%%%%%%%%%%%%%%%%

%if not user-supplied A1, get permuted A1, otherwise proceed
if ~matrix_data.user_A1_flag


    %fprintf('\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n')
    %fprintf('\n\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n')
    
    if verbosity_level
    
        fprintf('\n\tPermuting "A0" into "A1"\n')
        if auto_edge_data.flag
           fprintf('\t\tAuto-edges permitted...\n')
        else
            fprintf('\t\tAuto-edges disallowed...\n')

        end
    end
    performance.A1_start = tic;

    %call internal A0->A1 permuter;
    [permuteA1_flag, parmuteA1_msg, matrix_data] = permA_permute_A1(kin,matrix_data,auto_edge_data);


    performance.A1_finish = toc(performance.A1_start);
    if verbosity_level

        fprintf('\tcomplete in %2.3fs\n',performance.A1_finish)
    end
else
    if verbosity_level

        fprintf('\n\tUser-supplied "A1"; skipping A0 permutation...\n')
    end
    
end

if verbosity_level

    fprintf('\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n')
end
%fprintf('\n\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n')

%%%%%%%%%%%%%%%%%%
% Ai assembly -> Af
%%%%%%%%%%%%%%%%%%

if verbosity_level

    fprintf('\n\t\t///// Ai assembly sequence underway ///// \n');
    fprintf('\n\tShifting edges from "donor" to "recipient" nodes...\n')
end

performance.Ai_start = tic;

[assemble_Af_flag, assemble_Af_msg, matrix_data, assby_data, open_edge_data, inert_swap_data, performance] = ...
    permA_assemble_Af(kin,kout,matrix_data,auto_edge_data,assby_data,open_edge_data,inert_swap_data,performance);


performance.Ai_finish = toc(performance.Ai_start);

%
%print exit header
%

if verbosity_level

    fprintf('\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n')
    fprintf('\tPermutation Method Complete:\n')

    if assby_data.success
       fprintf('\n\tAssembly Successful in %2.3f secs',performance.Ai_finish) 

    else
       fprintf('\n\tAssembly FAILED after %2.3f secs',performance.Ai_finish)


    end

    fprintf('\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n')
end

performance.overall_finish = toc(performance.overall_start);


if verbosity_level

    fprintf('\n\t\tSwaps performed: %d; Inert Edge Shuffles: %d\n',assby_data.tot_swaps,inert_swap_data.tot_num_swaps)
    fprintf('\n\t\tOperation completed in %2.3f secs (%2.3f ave each loop)\n',...
        performance.overall_finish,performance.Ai_finish/assby_data.num_donor_loops);

    fprintf('\n\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n\n')
    fprintf('>> permute A Adjacency Matrix Assembly finished <<')
    fprintf('\n\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n\n')
end

%return final matrix A if successful:
if assby_data.success
    A_out = matrix_data.Ai; 
    
    
else
    %if assembly failed return A...?
    A_out = matrix_data.Ai;
    
    
end


%if user-requested additional data assemble here

%if user requested weight calculations, compute and return those here
if open_edge_data.calc_weights

   if debug_level >= 2
      fprintf('\n\tComputing weights for final matrix A...\n')
      
   end

  %inert swapping may introduce zeros into these available open rows / cols
  %squeeze out zero entries then compute
  open_edge_data.weight_num_open_rows = open_edge_data.weight_num_open_rows(open_edge_data.weight_num_open_rows~=0);
  open_edge_data.weight_num_open_cols = open_edge_data.weight_num_open_cols(open_edge_data.weight_num_open_cols~=0);
  open_edge_data.weights = prod(open_edge_data.weight_num_open_rows.*...
      open_edge_data.weight_num_open_cols);

end

%fprintf('\n\tnargout: %d\n',nargout);

%parse out from number of arguments out (return vars) and internal switches for what data to
%return...
switch nargout
    
    %if only 1 additional output var, simply return success flag
    case 2
        %fprintf('\n\treturning success flag %d\n',assby_data.success);
        %default: 2nd output variable is success flag
        varargout{1} = assby_data.success;
        

    %if 2 additional, return success flag, and either sample weights (if requested) or full
    %debugging structures (if no sample weights)
    case 3
        %fprintf('\n\treturning success flag %d\n',assby_data.success);
        %default: 2nd output variable is success flag
        varargout{1} = assby_data.success;
    
        %test if sample weights requested if not then return debug data only
        if open_edge_data.calc_weights
            
           %fprintf('\n\treturning sampling weight %d\n',open_edge_data.weights);
           %assign 3rd output var with weight calculations -- if requested
           varargout{2} = open_edge_data.weights;
            
        else
            
            %fprintf('\n\treturning debug structures\n')
            %assemble data_out -- debugging structures
            varargout{2}.matrix_data = matrix_data;
            varargout{2}.multi_edge_data = multi_edge_data;
            varargout{2}.auto_edge_data = auto_edge_data;
            varargout{2}.passed_parms_data = passed_parms_data;
            %if no A1 supplied, return A0 data
            if ~matrix_data.user_A1_flag
                varargout{2}.A0_hist = A0_hist;
            end
            varargout{2}.assby_data = assby_data;
            varargout{2}.open_edge_data = open_edge_data;  
            varargout{2}.inert_swap_data = inert_swap_data;
            varargout{2}.performance = performance;
            
        end
    
    %if 3 additiona output vars, return everything (so far??)
    case 4
        %fprintf('\n\treturning success flag %d\n',assby_data.success);
        %default: 2nd output variable is success flag
        varargout{1} = assby_data.success;

        %fprintf('\n\treturning sampling weight %d\n',open_edge_data.weights);
        %assign 3rd output var with weight calculations -- if requested
        varargout{2} = open_edge_data.weights;
        
        %fprintf('\n\treturning debug structures\n')
        %assemble data_out -- debugging structures
        varargout{3}.matrix_data = matrix_data;
        varargout{3}.multi_edge_data = multi_edge_data;
        varargout{3}.auto_edge_data = auto_edge_data;
        varargout{3}.passed_parms_data = passed_parms_data;
        %if no A1 supplied, return A0 data
        if ~matrix_data.user_A1_flag
            varargout{3}.A0_hist = A0_hist;
        end
        varargout{3}.assby_data = assby_data;
        varargout{3}.open_edge_data = open_edge_data;  
        varargout{3}.inert_swap_data = inert_swap_data;
        varargout{3}.performance = performance;
        
        
    
end

% if nargout == 2
%     fprintf('\n\treturning success flag %d\n',assby_data.success);
%     %default: 2nd output variable is success flag
%     varargout{1} = assby_data.success;
%     
%     
% elseif nargout == 3 && open_edge_data.calc_weights
%     fprintf('\n\treturning success flag %d\n',assby_data.success);
%     %default: 2nd output variable is success flag
%     varargout{1} = assby_data.success;
%     
%    fprintf('\n\treturning sampling weight %d\n',open_edge_data.weights);
%    %assign 3rd output var with weight calculations -- if requested
%    varargout{2} = open_edge_data.weights;
%     
% elseif (nargout == 3 && ~open_edge_data.calc_weights) || nargout > 3
%     fprintf('\n\treturning debug structures\n')
%     %assemble data_out -- debugging structures
%     varargout{2}.matrix_data = matrix_data;
%     varargout{2}.multi_edge_data = multi_edge_data;
%     varargout{2}.auto_edge_data = auto_edge_data;
%     varargout{2}.passed_parms_data = passed_parms_data;
%     %if no A1 supplied, return A0 data
%     if ~matrix_data.user_A1_flag
%         varargout{2}.A0_hist = A0_hist;
%     end
%     varargout{2}.assby_data = assby_data;
%     varargout{2}.open_edge_data = open_edge_data;  
%     varargout{2}.inert_swap_data = inert_swap_data;
%     varargout{2}.performance = performance;
% 
% %     %if user requested weight calculations, compute and return those here
% %     if open_edge_data.calc_weights
% % 
% %        if debug_level >= 2
% %           fprintf('\n\tComputing weights for final matrix A...\n')
% % 
% %        end
% % 
% %       open_edge_data.weights = prod(open_edge_data.weight_num_open_rows.*...
% %           open_edge_data.weight_num_open_cols);
% %     end
% 
%     %if all data structs requested, ensure above weight calc stuff included...
%     %data_out.open_edge_data = open_edge_data;
%     
% end

%%%%%%%%%%%%%%%%%%
%end permuteA
%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%
% internal functions
%%%%%%%%%%%%%%%%%

function [out_flag, out_msg, matrix_data, assby_data, open_edge_data, inert_swap_data, performance] = ...
    permA_assemble_Af(kin,kout,matrix_data,auto_edge_data,assby_data,open_edge_data,inert_swap_data,performance)
%internal function taking A1 as input and swapping edges between 'donor' and 'recipient' nodes until
%sum(Ai,1) == kout

    global N debug_level verbosity_level
    
    out_flag = 0;
    out_msg = 'permA_assemble_Af';
    out_data = [];

    if debug_level >= 1
        fprintf('permA_assemble_Af called...\n');
    end
    
    %initially, copy over A1 to Ai
    matrix_data.Ai = matrix_data.A1;
    
    %get current 'state' of column sum for Ai, or deviation of kout^i from target kout
    state = 1; %initialise
    [out_flag, out_msg, assby_data, inert_swap_data, matrix_data] = permA_update_state_Ai(kout,matrix_data,assby_data,inert_swap_data,state);

    if debug_level >= 2 && matrix_data.debug_tracking
        fprintf('\n\tPost-initialisation matrix_data.tracking{%d}...\n\n',assby_data.num_updates)
        matrix_data.tracking{assby_data.num_updates}
    end
    
    
    %ensure no fatal violations occurred so far (test kin = sum(A,2), sum(col_sum_A) = 0, etc)
    [out_flag, out_msg] = permA_errcheck_Ai(kin,matrix_data,assby_data);
    
    if out_flag < 0
       error(out_msg) 
    end
   
    if verbosity_level

        %fprintf('\n\tInitial distance from target kout: %2.3f...\n',assby_data.norm_targ_dist(assby_data.num_updates))
        fprintf('\tInitial Distance from Target kout\n\t\t\t\t\t>>>> %2.3f\n',assby_data.norm_targ_dist(assby_data.num_updates))

        %fprintf('\t\t\tNumber donors / recipients: %d / %d\n',assby_data.num_donors,assby_data.num_recips)
        fprintf('\t\t\tNumber donors / recipients / inerts: %d / %d / %d\n',...
            assby_data.num_donors,assby_data.num_recips,assby_data.num_inerts)
    end    
    
    %test if A1 permutation actually hit target Af
    if assby_data.norm_targ_dist(assby_data.num_updates) == 0
        
        if verbosity_level
            
           fprintf('\n\t!!! Target kout hit with A1: Its Your Lucky Day !!!\n\n')
            
            
        end
        
        %if user requested weight calculations return '1' since this is the only option!
        if open_edge_data.calc_weights
            
            %bump counter
            open_edge_data.weight_count = open_edge_data.weight_count + 1;

            %simply return '1' for both row and col
            open_edge_data.weight_num_open_rows(open_edge_data.weight_count) = 1;
            open_edge_data.weight_num_open_cols(open_edge_data.weight_count) = 1;
            
        end
        
        %flag this success type
        %success types: 0 <-> no success, 1 <-> full assembly (A0 -> A1 -> Ai -> Af)
        %2 <-> A0 -> A1 == Af
        assby_data.success = 1;
        assby_data.success_type = 2;
        assby_data.num_donor_loops = 0;        
        
    else
        assby_data.success = 0;
        assby_data.success_type = 0;
        assby_data.num_donor_loops = 1;
        
    end
    
    
    if debug_level >=3
       fprintf('assby_data.success: %d, assby_data.num_donor_loops: %d, assby_data.max_donor_loops: %d\n',...
           assby_data.success,assby_data.num_donor_loops,assby_data.max_donor_loops)
        
    end
    
    
    %main assembly loop
    while ~assby_data.success && (assby_data.num_donor_loops <= assby_data.max_donor_loops)
        
        performance.last_donor_loop_start = tic;
        
        if (assby_data.num_donor_loops > 1)
            
            if verbosity_level
            
            
                fprintf('\tCurrent Distance from Target kout\n\t\t\t\t\t>>>> %2.3f\n',assby_data.norm_targ_dist(assby_data.num_updates))
                fprintf('\tnumber swaps####\n\t\t\t\tlast loop:\t%d/%d possible\n\t\t\t\ttotal overall:\t%d\n',...
                    assby_data.num_swaps_this_loop,assby_data.prev_num_donors,assby_data.tot_swaps);
                fprintf('\tTime elapsed----\n\t\t\tlast donor loop: %e secs;\n',performance.time_elapsed_last_donor_loop)
                fprintf('\t\t\t\toverall: %e\n',toc(performance.Ai_start))
    %             if inert_swap_prev_loop
    %                 fprintf('\n\tinert swaps previous loop...\n')
    %             end        

%                 fprintf('\n  Donor Loop %d/%d; (#donors/recips: %d/%d)\n',...
%                     assby_data.num_donor_loops,assby_data.max_donor_loops,assby_data.num_donors,assby_data.num_recips)

                %below just print calculated #inerts instead of the assby_data.num_inerts field;
                %it's not intended for progression statistics...
%                 fprintf('\n  Donor Loop %d/%d; (#donors/recips/inerts: %d/%d/%d)\n',...
%                     assby_data.num_donor_loops,assby_data.max_donor_loops,...
%                     assby_data.num_donors,assby_data.num_recips,assby_data.num_inerts)
                
                fprintf('\n  Donor Loop %d/%d; (#donors/recips/inerts: %d/%d/%d)\n',...
                    assby_data.num_donor_loops,assby_data.max_donor_loops,...
                    assby_data.num_donors,assby_data.num_recips,(N - (assby_data.num_donors + assby_data.num_recips)))
                
                if open_edge_data.search_mode==1
%                     fprintf('\t\t("Exhaustive" Open Edge Search (mean ratio: %2.3f))\n',...
%                         open_edge_data.ave_recent_search_ratio)
                    fprintf('\t\t("Exhaustive" Open Edge Search (mean ratio: %2.3f))\n',...
                        mean(open_edge_data.search_history))
                    
                else
                    fprintf('\t\t("Quick" Open Edge Search)\n')
                end
                
            end

            %test for multi-edge count change...
            if debug_level >= 2
                %calculate number of multiedges before processing
                num_multiedge_cur = sum(matrix_data.Ai(matrix_data.Ai>1)) - nnz(matrix_data.Ai>1);
                %test for deviation
                if abs(num_multiedge_cur - matrix_data.num_multi) > 0
                      fprintf('\n\tWARNING! change in multi-edge number from pre (%d) to post (%d)...\n',...
                        matrix_data.num_multi,num_multiedge_cur);
                end
            end

%         elseif (assby_data.num_donor_loops == 2)
% 
%             fprintf('\tCurrent Distance from Target kout\n\t\t\t\t\t>>>> %2.3f\n',assby_data.norm_targ_dist(assby_data.num_updates))
%             fprintf('\tnumber swaps####\n\t\t\t\tlast loop:\t%d/%d possible\n\t\t\t\ttotal overall:\t%d\n',...
%                 assby_data.num_swaps_this_loop,assby_data.prev_num_donors,assby_data.tot_swaps);
%             fprintf('\tTime elapsed----\n\t\t\tlast donor loop: %e secs;\n',performance.time_elapsed_last_donor_loop)

        elseif (assby_data.num_donor_loops == 1)
            
            if verbosity_level
                fprintf('\nFirst Loop over "Donor" Nodes; Maximum #Loops %d\n',...
                    assby_data.max_donor_loops)
                %fprintf('\t\tDistance from Target kout\n\t\t\t\t\t>>>> %e\n',assby_data.norm_targ_dist(assby_data.num_updates))
            end
        end

        if debug_level >=1
            fprintf('\tdonor list:\n')
            assby_data.donors(1:assby_data.num_donors)
        end

        
        assby_data.swaps_this_loop = 0;
        assby_data.num_swaps_this_loop = 0;
        
        %loop over all donors
        for cur_donor = 1:assby_data.num_donors
        %for cur_donor = 1:9
            
            assby_data.donor_col = assby_data.donors(cur_donor);
            if debug_level >=1
                if (assby_data.num_donor_loops ==1)
                    fprintf('\n\tprocessing donor %d/%d in col# %d...\n',cur_donor,assby_data.num_donors,assby_data.donor_col);
                else
                    fprintf('\n\t(loop %d/%d) processing donor %d/%d in col# %d...\n',...
                        assby_data.num_donor_loops,assby_data.max_donor_loops,cur_donor,assby_data.num_donors,assby_data.donor_col);
                end
            end
            
            perform_swap = 0;
            
           %find rows with edges in this donor col
           %open_edge_data.donor_edge_rows = matrix_data.node_idx(logical(matrix_data.Ai(:,assby_data.donor_col)));
           %num_donor_edge_rows = length(donor_edge_rows);

%            if debug_level >= 2
%                fprintf('\n\trows with edges: \n');
%                donor_edge_rows
%            end
           
           %examine history of open-edge search attemtps; if too many failures over set memory
           %length, switch to 'quick' mode
           %open_edge_data.search_mode = 1;
%            %open_edge_data.search_history
%            if nnz(open_edge_data.search_history(1:open_edge_data.search_history_test_length) <= ...
%                    open_edge_data.search_switch_thresh)
%                open_edge_data.search_mode = 2;
%                
%                %if debug_level >=1
%                    fprintf('open edge failure history (%d) exceeds threshold (%d); switching to "quick" search\n',...
%                        open_edge_data.search_history_test_length,open_edge_data.search_switch_thresh)
%                %end
%                
%            end
           
           %determine if among these rows, we have suitable 'open' edge
           
%            %if in default search mode (row-by-row), test history of attempted row-searches:
%            if open_edge_data.search_mode == 1
%                
%                 %examine history last 10 searches; if ratios exceed switch threshold, change from
%                 %exhaustive row-by-row search to alternate, sub-matrix search
%                 ave_recent_search_ratio = mean(open_edge_data.search_history);
%                 
%                 if ave_recent_search_ratio > open_edge_data.search_switch_thresh
%                     %if debug_level >= 2
%                        fprintf('average recent search ratio (%2.3f) exceeds switch threshold (%2.3f)...\n',...
%                            ave_recent_search_ratio,open_edge_data.search_switch_thresh)
% 
%                    %end
%                    
%                    open_edge_data.search_mode = 2;
%                    
%                 end              
%                
%            end
%            
           %for now set 'mode' to exhaustive or '1'  <- this is default now
           %for tracking probabilities - set to '2' or quick since provides matrices of open edges
           %open_edge_data.search_mode = 2;  
           
           
           %if in default search mode (row-by-row), test history of attempted row-searches:
           if open_edge_data.search_mode == 1
               
                %examine history last 10 searches; if ratios exceed switch threshold, change from
                %exhaustive row-by-row search to alternate, sub-matrix search
                ave_recent_search_ratio = mean(open_edge_data.search_history);
                
                if ave_recent_search_ratio > open_edge_data.search_switch_thresh
                    if debug_level >= 2
                       fprintf('average recent search ratio (%2.3f) exceeds switch threshold (%2.3f)...\n',...
                           ave_recent_search_ratio,open_edge_data.search_switch_thresh)

                   end
                   
                   open_edge_data.search_mode = 2;
                   
                end              
               
           end
           
           %for now set 'mode' to exhaustive or '1'  <- this is default now
           %open_edge_data.search_mode = 1;  
           
           [out_flag, out_msg, open_edge_data, assby_data, performance] = ...
               permA_if_open_edge(matrix_data,assby_data,auto_edge_data,open_edge_data,performance);
           
%            %test comparison search method 2 
%            open_edge_data.search_mode = 2;  
%            [out_flag, out_msg, open_edge_data, assby_data, performance] = ...
%                permA_if_open_edge(matrix_data,assby_data,auto_edge_data,open_edge_data,performance);
           
           
                     
           %if open edges...find list of potential open cols among recip list
           if open_edge_data.flag
               %narrow list potentials recip cols to those with suitable difference
               %with donor edge(s): ideally want to move _one_ edge so reduce list to those 1 less than
               %donor edge value; we found these above in 'current_row_select'
               
               %below mapping of open edges moved to if_open_edge
               %recips_trim_open_edges = assby_data.recips_trim(open_edge_data.current_row_select);
               %open_edge_data.recips_trim_open_edges = assby_data.recips_trim(open_edge_data.current_row_select);
                
               if debug_level>=2
                   open_edge_data.recips_trim_open_edges
               end
               assby_data.target_open_edge_col = -1;
               
               %this list should be length >= 1 since open_edge found
               if nnz(open_edge_data.recips_trim_open_edges) == 1
                    if debug_level>=1
                        fprintf('\n\tOne open edge available - selecting %d for swap...\n',open_edge_data.recips_trim_open_edges);
                    end
                   assby_data.target_open_edge_col = open_edge_data.recips_trim_open_edges;
                   perform_swap = 1;

               else
                   %assby_data.target_open_edge_col = open_edge_data.recips_trim_open_edges(randperm(length(open_edge_data.recips_trim_open_edges),1));
                   assby_data.target_open_edge_col = open_edge_data.recips_trim_open_edges(randi(length(open_edge_data.recips_trim_open_edges)));
                   perform_swap = 1;
                    if debug_level>=1
                        fprintf('\n\tMultiple open edges available - randperm: %d selected...\n',assby_data.target_open_edge_col)
                    end
               end

               
                
           %otherwise skip to next donor col...    
           else
                if debug_level>=1
                     fprintf('\n\tRecipients with no suitable open edges available...\n')

                     %>>> should this be '<' num_donors? if have just one donor
                     if cur_donor < assby_data.num_donors
                         fprintf('\n\t\tAlternate Donors remain; skipping...\n');
                     else
                         fprintf('\n\t\tNo Alternate Donors remain; no swap this round...\n');
                     end

                end
               
               
           
           end
           

 
           if perform_swap 
                if debug_level>=1
                   fprintf('\n\tEdge Swap: Row %d, Donor Col %d, Recip Col %d...\n',...
                       open_edge_data.row_select_idx,assby_data.donor_col,assby_data.target_open_edge_col);
                end

                assby_data.swaps_this_loop = 1;
                assby_data.num_swaps_this_loop = assby_data.num_swaps_this_loop + 1;

                if debug_level>=2
                    fprintf('\n\tShifting 1 edge from (%d,%d) to (%d,%d)...\n',...
                        open_edge_data.row_select_idx,assby_data.donor_col,open_edge_data.row_select_idx,assby_data.target_open_edge_col)
                end
                
                %matrix_data.Ai
                %sum(matrix_data.Ai)


                %perform the actual edge swap
                matrix_data.Ai(open_edge_data.row_select_idx,assby_data.donor_col) = ...
                    matrix_data.Ai(open_edge_data.row_select_idx,assby_data.donor_col) - 1;
                matrix_data.Ai(open_edge_data.row_select_idx,assby_data.target_open_edge_col) = ...
                    matrix_data.Ai(open_edge_data.row_select_idx,assby_data.target_open_edge_col) + 1;
                
                %matrix_data.Ai
                %sum(matrix_data.Ai)
                
                %assby_data.A_col_sum


                if debug_level >=2
                   fprintf('\n\tUpdating col_sum_diff: \n');
    %                fprintf('\tdonor_col by %d; target_col by %d...\n',...
    %                    -(source_val - target_val),(source_val - target_val));
                    fprintf('\tdonor_col by -1; target_col by +1...\n')

                end

                %call state update function
                state = 2; %2 <-> 'update col_sum_diff, recip list
                [out_flag, out_msg, assby_data, inert_swap_data, matrix_data] = permA_update_state_Ai(kout,matrix_data,assby_data,inert_swap_data,state);
                
                

           else

               if debug_level>=1
                    fprintf('\n\tNOT performing swap; skipping...\n');
               end

               %skip tracking - may activate only with debuging...
                assby_data.num_donor_skips = assby_data.num_donor_skips + 1;
                assby_data.skip_donor_tracking(assby_data.num_donor_skips) = assby_data.donor_col;

           end
            
           
           %test current distance from target kout; if hit, flag for success!
           if assby_data.norm_targ_dist(assby_data.num_updates) == 0
               if debug_level>=1
                   fprintf('\n\n\t\tTarget kout distance: %2.3f; Success!\n',assby_data.norm_targ_dist(assby_data.num_updates))
               end
               assby_data.success = 1;
               assby_data.success_type = 1;
           end
            
            
        end
        
        assby_data.tot_swaps = assby_data.tot_swaps + assby_data.num_swaps_this_loop;
       
        %calculate ratio of swaps to potential max
        assby_data.swap_ratio = assby_data.num_swaps_this_loop/assby_data.num_donors;
                
        %
        %inert swapping here...
        %
        
        
%         %if no swaps this loop, call inert edge swapper; simply exchanges edges between nodes already at
%         %full kout with those that aren't, hopefully to permit subsequent swaps between donors & recips
%         if ~assby_data.swaps_this_loop

        %test for number of swaps falling below activation threshold - if so, call inert shuffling
        %if we actually have inerts...
        
        %update state for num_inerts needs fixing; for now utilise difference from num_donors and
        %num_recips
%        if (assby_data.swap_ratio <= inert_swap_data.activate_ratio) && (assby_data.num_inerts > 0)
         temp_num_inerts = N - (assby_data.num_donors + assby_data.num_recips);
         if (assby_data.swap_ratio <= inert_swap_data.activate_ratio) && (temp_num_inerts > 0)
           
            %if debug_level >=2
%                fprintf('\n\tdonor->recip swap ratio (%2.3f) below inert shuffle activation threshold (%2.3f)\n',...
%                    assby_data.swap_ratio,inert_swap_data.activate_ratio)
                
            %end

            if verbosity_level

                fprintf('\n\tInert shuffling activated -->\n\t\t(current swap ratio/threshold: %2.3f/%2.3f)\n',...
                   assby_data.swap_ratio,inert_swap_data.activate_ratio)
            end            
            
            %fprintf('\n\tNo swaps performed previous donor loop --\n\t\tshuffling edges with inert nodes...');
            
            %call state update function - only update for inert swaps
            state = 4; %4<-> update inert_swap lists
            [out_flag, out_msg, assby_data, inert_swap_data, matrix_data] = permA_update_state_Ai(kout,matrix_data,assby_data,inert_swap_data,state);
            
            %call inert swap routine
            performance.inert_swap_start = tic;
            [out_flag, out_msg, matrix_data, assby_data, inert_swap_data] = permA_inert_swap(matrix_data,assby_data,inert_swap_data,auto_edge_data);
            performance.inert_swap_finish = toc(performance.inert_swap_start);
            
            if verbosity_level
                fprintf('\t\tcompleted %d shuffles in %2.3f secs\n\n',inert_swap_data.num_swaps,performance.inert_swap_finish)
            end
            
            %test for num_recip error:
            num_recip_col_sum = nnz(assby_data.col_sum_diff < 0);
            if (num_recip_col_sum ~= assby_data.num_recips)
                error('recip number deviation')
            end
            
            
        else
            if debug_level >=1
                fprintf('\n\tNo Inert Shuffling--\n');
                fprintf('\tassby_data.swap_ratio %2.3f, assby_data.num_inerts %d\n',...
                    assby_data.swap_ratio,assby_data.num_inerts);
            end
        end
        
        %if %swaps performed; update the donor list
        if assby_data.num_swaps_this_loop > 0
            %call state update function
            state = 3; %3 <-> 'update donor_list
            [out_flag, out_msg, assby_data, inert_swap_data, matrix_data] = permA_update_state_Ai(kout,matrix_data,assby_data,inert_swap_data,state);
        end
                
        %
        %donor loop data tracking...
        %
        
        assby_data.donor_loop_num_swaps(assby_data.num_donor_loops) = assby_data.num_swaps_this_loop;
        assby_data.donor_loop_norm_tracking(assby_data.num_donor_loops) = assby_data.norm_targ_dist(assby_data.num_updates);
        
        assby_data.num_donor_loops = assby_data.num_donor_loops + 1;

        %update time tracker this last loop
        performance.time_elapsed_last_donor_loop = toc(performance.last_donor_loop_start);
        
        
    end
    
    
    
end

function [out_flag, out_msg, matrix_data, assby_data, inert_swap_data] = permA_inert_swap(matrix_data,assby_data,inert_swap_data,auto_edge_data)
%internal function for swapping out entries between 'inert' and 'donor' or 'recipient' columns of A

    global N debug_level
    
    out_flag = 0;
    out_msg = 'permA_inert_swap';
    
    %flag returns '0' for no swaps performed (or failure), '1' for performed
    assby_data.swapped_flag = 0;
    
    if debug_level >= 1
        fprintf('permA_inert_swap called...\n');
        
    end
    
    if debug_level >= 1
       fprintf('\n\tpre swap: topick: num donors: %d recips: %d inerts: %d...\n',...
           inert_swap_data.num_donors_topick,inert_swap_data.num_recips_topick,inert_swap_data.num_inert_topick)
    end
    
    
    %determine which pool of recips or donors is larger and pick as target for swapping
    if inert_swap_data.num_recips_topick > inert_swap_data.num_donors_topick
        pick_type = 1; %recips
    elseif inert_swap_data.num_recips_topick < inert_swap_data.num_donors_topick
        pick_type = 2; %donors
    %or if equivalent, flip a coin between them
    else
        pick_type = randi(2);
    end
    
    if debug_level >= 1
        fprintf('\n\tpick type %d...\n',pick_type)
    end
    
    %ensure we have at least 1 inert to attempt swapping, and set max_attemtps depending on pick
    %type
    if inert_swap_data.num_inert_topick == 0
        max_attempts = 0;
    else
        if pick_type == 1
            max_attempts = inert_swap_data.num_recips_topick*inert_swap_data.num_inert_topick;
        else
            max_attempts = inert_swap_data.num_donors_topick*inert_swap_data.num_inert_topick;
        end
    end
    
    try_swap = 1;
    num_attempt_swap_rounds = 0;
    inert_swap_data.num_swaps = 0;
    
    %keep attemtping swaps as long as we have >0 inerts to pick from...
    while try_swap && (num_attempt_swap_rounds <= max_attempts)
        
        %select random inert from available set
        cur_inert_idx = inert_swap_data.inert_idx_topick(randi(inert_swap_data.num_inert_topick));
        
        if debug_level>=1
            fprintf('\n\tSelected inert col #%d...\n',cur_inert_idx)
        end
        
       %select target node (either recip / donor as pick_type)
       if pick_type == 1
           swap_pick_idx = inert_swap_data.recip_idx_topick(randi(inert_swap_data.num_recips_topick));
           
       else
           swap_pick_idx = inert_swap_data.donor_idx_topick(randi(inert_swap_data.num_donors_topick));
           
       end
       
%        if debug_level>=1
%           pick_str = [pick_str num2str(swap_pick_idx) ' picked for swapping...'];
%           fprintf('%s\n',pick_str);
%        end
       

       %get indexes rows with open edges or '0'
       %can exploit 'logical' here since all entries should be 0 or 1
       %cur_inert_open_edge_idx = ~logical(matrix_data.Ai(:,cur_inert_idx));
       
%        if debug_level >=2
%            %save copy of Ai for testing auto-edge err
%            fprintf('\n\tsaving copy of Ai to inert swap data for auto-edge debug...\n')
%            inert_swap_data.Ai_copy = matrix_data.Ai;
%        end
       
       inert_open_edge_idx = matrix_data.node_idx(~logical(matrix_data.Ai(:,cur_inert_idx)));
       %omit self-indexed open edges avoiding introducting auto-connects
       inert_open_edge_idx = inert_open_edge_idx(inert_open_edge_idx~=cur_inert_idx);
       %and find the filled edges for this corresponding target node (recip / donor)
       swap_pick_pair_filled_edge_idx = inert_open_edge_idx(logical(matrix_data.Ai(inert_open_edge_idx,swap_pick_idx)));
        
       %now opposite: inert filled edges and swap pick open...
       %inert_filled_edge_idx = node_idx(logical(A(:,cur_inert_idx)));
       
       %below statement with opposite logical filtering from original
       %code...may introduce auto-edges
       %inert_filled_edge_idx = matrix_data.node_idx(~logical(matrix_data.Ai(:,cur_inert_idx)));
       
       %below statement with same logical filtering from original code...
       inert_filled_edge_idx = matrix_data.node_idx(logical(matrix_data.Ai(:,cur_inert_idx)));
       swap_pick_pair_open_edge_idx = inert_filled_edge_idx(~logical(matrix_data.Ai(inert_filled_edge_idx,swap_pick_idx)));
       %omit self-indexed open edges avoiding introducting auto-connects
       swap_pick_pair_open_edge_idx = swap_pick_pair_open_edge_idx(swap_pick_pair_open_edge_idx~=swap_pick_idx);
       
       
       %if we have filled and open pairs of edges to swap, swap!
       if length(swap_pick_pair_filled_edge_idx) >= 1 && length(swap_pick_pair_open_edge_idx) >= 1

        if debug_level>=1
              fprintf('\tsufficient edge pairs (0,X)->(X,0) for swapping...\n');
        end
        
          %of the selected filled and open pairs, randomly pick one of each
          swap_pick_pair_filled_select_idx = ...
              swap_pick_pair_filled_edge_idx(randi(length(swap_pick_pair_filled_edge_idx)));
          swap_pick_pair_open_select_idx = ...
              swap_pick_pair_open_edge_idx(randi(length(swap_pick_pair_open_edge_idx)));
        
        if debug_level>=1
              fprintf('\tswapping "filled" (%d) with "open" (%d)...\n',...
                  swap_pick_pair_filled_select_idx,swap_pick_pair_open_select_idx)
        end       
        
        
          inert_swap_data.num_swaps = inert_swap_data.num_swaps + 1;
          if debug_level>=2
              inert_swap_data.A_inert_copy{inert_swap_data.num_swaps} = matrix_data.Ai;
          end
          
            %save for target tracking
            inert_swap_data.swap_target_tracking(inert_swap_data.num_swaps) = swap_pick_idx;

          %swap 'em - below method may be slower but works with multi-edge entries
          %that are > 1...
          matrix_data.Ai(swap_pick_pair_filled_select_idx,cur_inert_idx) = ...
              matrix_data.Ai(swap_pick_pair_filled_select_idx,cur_inert_idx) + 1;
          matrix_data.Ai(swap_pick_pair_open_select_idx,cur_inert_idx) = ...
              matrix_data.Ai(swap_pick_pair_open_select_idx,cur_inert_idx) - 1;
          
          matrix_data.Ai(swap_pick_pair_filled_select_idx,swap_pick_idx) = ...
              matrix_data.Ai(swap_pick_pair_filled_select_idx,swap_pick_idx) - 1;
          matrix_data.Ai(swap_pick_pair_open_select_idx,swap_pick_idx) = ...
              matrix_data.Ai(swap_pick_pair_open_select_idx,swap_pick_idx) + 1;
          
          %test for introduction auto-edges - if not flagged
          if ~auto_edge_data.flag
          %if debug_level>=2
             num_auto_edges_post_inert_swap =  nnz(diag(matrix_data.Ai));
             if num_auto_edges_post_inert_swap > 0
                 inert_swap_data.auto_edge_err_num = inert_swap_data.auto_edge_err_num + 1;
                 
%                  %save data 
%                  inert_swap_data.cur_inert_idx_copy{inert_swap_data.auto_edge_err_num} = cur_inert_idx;
%                  inert_swap_data.inert_open_edge_idx_copy{inert_swap_data.auto_edge_err_num} = inert_open_edge_idx;
%                  inert_swap_data.swap_pick_idx_copy{inert_swap_data.auto_edge_err_num} = swap_pick_idx;
%                  inert_swap_data.swap_pick_pair_filled_edge_idx_copy{inert_swap_data.auto_edge_err_num} = swap_pick_pair_filled_edge_idx;
%                  inert_swap_data.
                 error('auto edges introduced post-inert swap')
             end
          end
          %end
          
          
          %winnow down the pick lists due to successful swap
            %remove this inert from list
%            inert_swap_data.inert_idx_topick = inert_swap_data.inert_idx_topick(inert_swap_data.inert_idx_topick~=cur_inert_idx);
%             inert_swap_data.num_inert_topick = length(inert_swap_data.inert_idx_topick);
            inert_swap_data.inert_idx_topick(inert_swap_data.inert_idx_topick == cur_inert_idx) = [];
            inert_swap_data.num_inert_topick = inert_swap_data.num_inert_topick - 1;
            
           %winnow out this pick from available pool depending on donor / recip pick

            if (pick_type == 1)
                
%                inert_swap_data.recip_idx_topick = inert_swap_data.recip_idx_topick(inert_swap_data.recip_idx_topick~=swap_pick_idx);
%                inert_swap_data.num_recips_topick = length(inert_swap_data.recip_idx_topick)
                inert_swap_data.recip_idx_topick(inert_swap_data.recip_idx_topick ==  swap_pick_idx) = [];
                inert_swap_data.num_recips_topick = inert_swap_data.num_recips_topick - 1;
                
            else
%                inert_swap_data.donor_idx_topick = inert_swap_data.donor_idx_topick(inert_swap_data.donor_idx_topick~=swap_pick_idx);
%                inert_swap_data.num_donors_topick = length(inert_swap_data.donor_idx_topick)
                inert_swap_data.donor_idx_topick(inert_swap_data.donor_idx_topick == swap_pick_idx) = [];
                inert_swap_data.num_donors_topick = inert_swap_data.num_donors_topick - 1;
                
            end
            
            if debug_level >= 1
               fprintf('\n\tpost swap: topick num donors: %d recips: %d inerts: %d...\n',...
                   inert_swap_data.num_donors_topick,inert_swap_data.num_recips_topick,inert_swap_data.num_inert_topick)
            
               fprintf('\t (logic col_sum_diff) num donor: %d, recips %d, inerts: %d\n',...
                   nnz(assby_data.col_sum_diff>0),nnz(assby_data.col_sum_diff<0),nnz(assby_data.col_sum_diff==0));
              
            
            end
            
            
            
       else
           
           %insufficent pairs for swapping...
            if debug_level>=1
                  fprintf('\tinsufficient pairs for swapping...\n');
            end
         

       end
       
         
       num_attempt_swap_rounds = num_attempt_swap_rounds + 1;
          
        if debug_level>=1

           if pick_type ==1 
                fprintf('\n\t%d inert & %d recip remain after %d rounds...\n',...
                    inert_swap_data.num_inert_topick,inert_swap_data.num_recips_topick,num_attempt_swap_rounds);
                if debug_level>=2
                    recip_swap_list_topick
                end
           elseif pick_type == 2
                fprintf('\n\t%d inert & %d donors remain after %d rounds...\n',...
                    inert_swap_data.num_inert_topick,inert_swap_data.num_donors_topick,num_attempt_swap_rounds);
                if debug_level>=2
                    donor_swap_list_topick
                end
           end

        end    
        
       %if we still have target (recip or donor) cols to swap, and sufficient pool of inerts, flag 
       %to try another round...
       try_swap = 0;
       if inert_swap_data.num_inert_topick>0
           if (pick_type == 1) && (inert_swap_data.num_recips_topick>0)
               try_swap = 1;
           elseif (pick_type == 2) && (inert_swap_data.num_donors_topick>0)
               try_swap = 1;
%            else
%                try_swap = 0;
           end
%        else
%            try_swap = 0;
       end        
       
       
    end
    
    
    %update overall swap tracking
    inert_swap_data.tot_num_swaps = inert_swap_data.tot_num_swaps + inert_swap_data.num_swaps;
    
    if debug_level >= 1
        
       fprintf('\n\tpost-inert shuffle: (assby_data) num donor: %d, recips: %d, inerts: %d\n',assby_data.num_donors,assby_data.num_recips,assby_data.num_inerts);
       fprintf('\t (logic) num donor: %d, recips %d, inerts: %d\n',...
           nnz(assby_data.col_sum_diff>0),nnz(assby_data.col_sum_diff<0),nnz(assby_data.col_sum_diff==0));
        
    end

end

function [out_flag, out_msg, open_edge_data, assby_data, performance] = permA_if_open_edge(matrix_data,assby_data,auto_edge_data,open_edge_data,performance)
%internal function searches rows of a given donor node/column until find an 'open' edge for swapping
%with in corresponding recipient cols; 
%note: flag setting for search mode in open_edge_data.search_mode; default 'exhaustive' search <->
%1, 'quick' sub-matrix search <-> '2'

    global N debug_level
    
    out_flag = 0;
    out_msg = 'permA_if_open_edge';
    
    %flag returns either open edge found '1' or not '0'
    open_edge_data.flag = 0;    
    
    %update history tracker - shift previous effort back in history
    open_edge_data.search_history = circshift(open_edge_data.search_history,1);
    
    
    
    if debug_level >= 1
        fprintf('permA_if_open_edge called...\n');
        
        %if weight calculations requested, notify
        if open_edge_data.calc_weights
           fprintf('\n\tWeight calculations requested >>> forcing "quick" edge search...\n') 
        end
    end
    
    open_edge_data.donor_edge_rows = matrix_data.node_idx(logical(matrix_data.Ai(:,assby_data.donor_col)));
    open_edge_data.num_donor_edge_rows = length(open_edge_data.donor_edge_rows);
    
%     %test alt method 1 (for memory efficiency)
%     open_edge_data.donor_edge_rows_test_idx = logical(matrix_data.Ai(:,assby_data.donor_col));
%     open_edge_data.num_donor_edge_rows_test = nnz(open_edge_data.donor_edge_rows_test_idx);
%     open_edge_data.donor_edge_rows_test = matrix_data.node_idx(open_edge_data.donor_edge_rows_test_idx); 
    
%     %open_edge_data.num_donor_edge_rows = nnz(logical(matrix_data.Ai(:,assby_data.donor_col)));
%     open_edge_data.donor_edge_rows_idx = logical(matrix_data.Ai(:,assby_data.donor_col));
%     open_edge_data.num_donor_edge_rows = nnz(open_edge_data.donor_edge_rows_idx);
%     open_edge_data.donor_edge_rows(1:open_edge_data.num_donor_edge_rows) = ...
%         matrix_data.node_idx(open_edge_data.donor_edge_rows_idx);  
    
%     %testing alternative method
%     open_edge_data.num_donor_edge_rows_test = open_edge_data.num_donor_edge_rows;
%     open_edge_data.donor_edge_rows_test(1:open_edge_data.num_donor_edge_rows_test) = ...
%         open_edge_data.donor_edge_rows(1:open_edge_data.num_donor_edge_rows);
    
%     open_edge_data.num_donor_edge_rows = length(open_edge_data.donor_edge_rows);
%     
%     open_edge_data.num_donor_edge_rows_test = open_edge_data.num_donor_edge_rows;
    
    if debug_level >= 2
       fprintf('\n\trows with edges: \n');
       open_edge_data.donor_edge_rows
    end
    
    
    %user weight calculation check - if user requests weights returned, force 'quick' edge search
    %since only via the sub-matrices of A can we compute the number available 'rows' and 'cols'
    %sans exhaustive testing for available options..
    if open_edge_data.calc_weights
        open_edge_data.search_mode = 2;
        
    end
    
    %'exhaustive' / standard search
    if open_edge_data.search_mode == 1
        
           performance.open_search_std_start = tic;
        
           %sample rows from the current passed donor until we find a open edge to swap
           %open_edge = 0;
           max_row_samples = nnz(open_edge_data.donor_edge_rows);
           cur_row_sample = 1;
        
           if debug_level >= 2
               fprintf('\n\tsearching for open edges along cur_donor rows; max_row_samples: %d \n',max_row_samples);

           end

           while ~open_edge_data.flag && cur_row_sample <= max_row_samples            
            
               %below notes from permute_a_assby_t1p6p1
               %>>> when number of recipients falls this procedure may
               %exhaustively search all possible rows for open edges;
               %such a search may be unneccessary: there may be NO open edges
               %between the donor and recipient:
               %try testing for any available open edges from donor -> recip
               %and if not skip the exhaustive search...

               %>>> if the difference between donor and recipient columns has NO
               %positive entries, then no swap is possible - perhaps after
               %exhaustive searches start failing, then perform difference
               %calculations (maybe more expensive to do at first) and skip the
               %exhaustive search...

               %for now just test if multi-edge swapping works then optimise as suggested

               %select one row from the donor column
               %row_select = randperm(open_edge_data.num_donor_edge_rows,1);
               
               %test randi
               row_select = randi(open_edge_data.num_donor_edge_rows);

               %get this row's index in A
               open_edge_data.row_select_idx = open_edge_data.donor_edge_rows(row_select);

               if debug_level>=3
                   fprintf('\n\tselected row %d\n',open_edge_data.row_select_idx);
               end

              %this row selected may be in the recips list; if so, remove it to prevent
               %Auto-connects (depending on permit_auto_connects flag)

               assby_data.recips_trim_trimmed_flag = 0;
               if ~auto_edge_data.flag
%                     if debug_level>=2
%                         fprintf('\n\tsize assby_data.recips_trim(1:assby_data.num_recips-1):\n')
%                         size(assby_data.recips_trim(1:assby_data.num_recips-1))
%                         fprintf('\tsize assby_data.recips(assby_data.recips(1:assby_data.num_recips) ~= row_select_idx):\n')
%                         
%                         size(assby_data.recips(assby_data.recips(1:assby_data.num_recips) ~= row_select_idx))
%                     end
                    %load recips_trim with logic array, then load values (sized per nnz entries)
%                    if debug_level>=3
%                     
%                     fprintf('\nremoving element %d from assby_data.recips_trim (len: %d, num_recips: %d); before:\n',...
%                         open_edge_data.row_select_idx,assby_data.recips_trim_len,assby_data.num_recips)
%                         assby_data.recips_trim
%                    end     
                   
%                     assby_data.recips_trim_idx = assby_data.recips(1:assby_data.num_recips) ~= open_edge_data.row_select_idx;
%                     assby_data.recips_trim_len = nnz(assby_data.recips_trim_idx);
%                     %nnz(assby_data.recips_trim_idx)
%                     assby_data.recips_trim(1:assby_data.recips_trim_len) = assby_data.recips(assby_data.recips_trim_idx);
                    
%                    if debug_level>=3
%                         fprintf('after (len: %d):\n',assby_data.recips_trim_len)
%                         assby_data.recips_trim(1:assby_data.recips_trim_len)
%                         assby_data.recips_trim_idx
%                    end
                    
                    %testing alternative method
%                     assby_data.recips_trim(1:assby_data.num_recips) = assby_data.recips(1:assby_data.num_recips);
%                     assby_data.recips_trim(assby_data.recips_trim == open_edge_data.row_select_idx) = [];
%                     assby_data.recips_trim_len = assby_data.num_recips-1;
                    %remove current selected item from recips list
                    %size(assby_data.recips_trim_test)
                    
                    %tracking length of this array should be init / updated in functions since it
                    %may shrink during open edge search...
                    %assby_data.recips_trim_test_len = assby_data.num_recips;
                    
                    %assby_data.recips_trim_trimmed_flag = 0;
                    %if current selected row in recips list, exclude it
                    if any(assby_data.recips_trim(1:assby_data.recips_trim_len) == open_edge_data.row_select_idx)
                        
                       if debug_level>=3
                        
                            fprintf('removing element %d from recips_trim (len: %d); before:\n',...
                                open_edge_data.row_select_idx,assby_data.recips_trim_len)
                            %assby_data.recips_trim_test(1:assby_data.recips_trim_test_len)
                            assby_data.recips_trim
                       end
                       
                        %flag removal of this row_select_idx; in case we need insert it back...
                        assby_data.recips_trim_trimmed_flag = 1;
                        
                        assby_data.recips_trim(assby_data.recips_trim == ...
                            open_edge_data.row_select_idx) = [];
                        assby_data.recips_trim_len = assby_data.recips_trim_len - 1;
                        
                        
                       if debug_level>=3
                        
                            fprintf('after (len: %d):\n',assby_data.recips_trim_len)
                            %assby_data.recips_trim_test(1:assby_data.recips_trim_test_len)
                            assby_data.recips_trim
                       end                       
                        %size(assby_data.recips_trim_test)
                    end
                    
%                     if debug_level>=1
%                         %fprintf('\n\trecips_trim list comp: ');
%                         assby_data.num_recips;
%                         temp_trim = sort(assby_data.recips_trim(1:assby_data.recips_trim_len));
%                         temp_trim_test = sort(assby_data.recips_trim_test(1:assby_data.recips_trim_test_len));
%                         %size(temp_trim)
%                         %size(temp_trim_test)
%                         
%                         
%                         try
%                             trim_list_diff = nnz(temp_trim - temp_trim_test);
%                         catch 
%                             trim_list_diff = -1;
%                             warning('dimensions mismatch: temp_trim and temp_trim_test');
%                             size(temp_trim)
%                             size(temp_trim_test)
%                             open_edge_data.trim_test_warning = open_edge_data.trim_test_warning + 1;
%                             open_edge_data.trim_test_donor{open_edge_data.trim_test_warning} = assby_data.donor_col;
%                             open_edge_data.trim_test_row_select_idx{open_edge_data.trim_test_warning} = open_edge_data.row_select_idx;
%                             open_edge_data.temp_trim_debug{open_edge_data.trim_test_warning} = ...
%                                 assby_data.recips_trim;
%                             open_edge_data.recips_trim_len_debug{open_edge_data.trim_test_warning} = assby_data.recips_trim_len;
%                             open_edge_data.temp_trim_test_debug{open_edge_data.trim_test_warning} = ...
%                                 assby_data.recips_trim_test;
%                             open_edge_data.recips_trim_test_len_debug{open_edge_data.trim_test_warning} = assby_data.recips_trim_test_len;
%                             open_edge_data.cur_row_sample{open_edge_data.trim_test_warning} = cur_row_sample;
%                             
%                         end
%                         %fprintf(' %d\n',trim_list_diff)
% %                         if (trim_list_diff ~= 0)
% %                             temp_trim
% %                             temp_trim_test
% %                             error('recips_trim error')
% %                         end
%                     end
                    
                    if debug_level>=3
                        fprintf('\n\tAuto-connects disallowed; excluding %d from recip\n',open_edge_data.row_select_idx);
                        fprintf('\tassby_data.recips_trim(1:assby_data.recips_trim_len):\n')
                        assby_data.recips_trim(1:assby_data.recips_trim_len)
                    end

               else
                    assby_data.recips_trim = assby_data.recips;
                    assby_data.recips_trim_len = nnz(assby_data.recips_trim);
                    if debug_level>=3
                        fprintf('\n\tAuto-connects allowed; including %d in recip\n',open_edge_data.row_select_idx);
                        fprintf('\tassby_data.recips_trim(1:assby_data.recips_trim_len):\n')
                        assby_data.recips_trim(1:assby_data.recips_trim_len)
                    end

               end

              %from this row, find recipient cols with no edge (0-entry)

               %first, test if we have a 'full-row' i.e., no 0-entries except auto-connect

               %compute the _difference_ between donor entry (A(row_selected,donor_col)) and other
               %entries in the row. We seek 'open' edges such that recipients are 1-edge less than 
               %the donor entry; we do this for maintaining #of multi-edges while performing a
               %finer-grained swapping instead of shifting whole multi-edge blocks at once, e.g.
               %2+ donor edges -> 0 open recip edges appears too coarse; instead do: 2, 3, 4 edges -> 1, 2, 3,
               %etc...
               %current_row_diff = A(row_select_idx,donor_col) - A(row_select_idx,recips_trim);

               %the difference above is parsed depending on the number of edges in this donor selection:
               %e.g, if A(row_select_idx,donor_col) == 1, then we can only transfer this solo edge to an
               %empty open edge in the recipient list, e.g., only if current_row_diff = 1
               %OR, if A(row_select_idx,donor_col) > 1, then we can only transfer this MULTI edge to
               %another edge such that we don't lose the multi edge; e.g., only if current_row_diff == 
               %A(row_select_idx,donor_col), or simply the recipient is an open, zero, edge...
               
               if matrix_data.Ai(open_edge_data.row_select_idx,assby_data.donor_col) == 1  %have solo edge in donor; can only swap to open recips
                   if debug_level >= 2
                       fprintf('\n\tDonor edge "solo"; selecting only empty-edge recips...\n');
                   end
%                    current_row_diff = matrix_data.Ai(open_edge_data.row_select_idx,assby_data.donor_col) - ...
%                        matrix_data.Ai(open_edge_data.row_select_idx,assby_data.recips_trim(1:assby_data.recips_trim_len));
%                    open_edge_data.current_row_select = current_row_diff == 1;
                   
                   %instead, only test for '0' entries since we can only
                   %xfer from solo to open edges
%                    current_row_diff_test = ...
%                        matrix_data.Ai(open_edge_data.row_select_idx,assby_data.recips_trim(1:assby_data.recips_trim_len)) == 0;
%                    open_edge_data.current_row_select = ...
%                        matrix_data.Ai(open_edge_data.row_select_idx,assby_data.recips_trim(1:assby_data.recips_trim_len)) == 0;

                   open_edge_data.current_row_select = ...
                       matrix_data.Ai(open_edge_data.row_select_idx,assby_data.recips_trim(1:assby_data.recips_trim_len)) == 0;
                   
                   
%                    if nnz(open_edge_data.current_row_select - current_row_diff_test) >0
%                        error('current_row_diff_test deviated from current_row_diff...')
%                    end
                   
%                    %~logical way
%                    current_row_diff_test2 = ...
%                        ~logical(matrix_data.Ai(open_edge_data.row_select_idx,assby_data.recips_trim(1:assby_data.recips_trim_len)));
                   
%                    if nnz(open_edge_data.current_row_select - current_row_diff_test2) >0
%                        error('current_row_diff_test2 deviated from current_row_diff...')
%                    end
                   

               else %should be > 1; have multi-edge in donor; can only swap to multi recips
                   if debug_level >= 2
                       fprintf('\n\tDonor edge "multi"; selecting only filled-edge recips...\n');
                   end

                   open_edge_data.current_row_select = ...
                       matrix_data.Ai(open_edge_data.row_select_idx,assby_data.recips_trim(1:assby_data.recips_trim_len)) ~= 0;


               end
               
               if debug_level >= 2
                   %fprintf('\n\tcurrent row difference: \n')
                   %current_row_diff
                   fprintf('\n\tcurrent row diff ones: \n')
                   open_edge_data.current_row_select
               end
               
                %check current_row_diff for any values >= 1 => the donor edge has surplus edges
                %over the recipients in this row
                %if max(current_row_diff) <= 0
                if ~any(open_edge_data.current_row_select)  %'any' returns true if any entry is nonzero or 1
                   if debug_level>=1
                        %fprintf('\n\trow selected (%d) full of edges\n',row_select_idx);
                        fprintf('\n\trow selected (%d) without available recipients...\n',open_edge_data.row_select_idx);
                   end
                   %remove this current row from pool
                   open_edge_data.donor_edge_rows = open_edge_data.donor_edge_rows(...
                       open_edge_data.donor_edge_rows ~= open_edge_data.row_select_idx);
                   open_edge_data.num_donor_edge_rows = open_edge_data.num_donor_edge_rows - 1;
                   
%                    %test removing alt method
%                    open_edge_data.donor_edge_rows_test(open_edge_data.donor_edge_rows_test == open_edge_data.row_select_idx) = [];
%                    open_edge_data.num_donor_edge_rows_test = open_edge_data.num_donor_edge_rows_test - 1;
                   
%                    fprintf('open_edge_data.donor_edge_rows (len %d)\n',open_edge_data.num_donor_edge_rows)
%                    open_edge_data.donor_edge_rows
%                    
%                    fprintf('open_edge_data.donor_edge_rows_test (len %d)\n',open_edge_data.num_donor_edge_rows_test)
%                    open_edge_data.donor_edge_rows_test
%                    
%                    
%                    if nnz(open_edge_data.donor_edge_rows(1:open_edge_data.num_donor_edge_rows) - ...
%                            open_edge_data.donor_edge_rows_test(1:open_edge_data.num_donor_edge_rows_test))>0
%                        error('donor edge row deviation...')
%                    end
                   
                   cur_row_sample = cur_row_sample + 1;
                   
                   %for test recips_trim, put this selected row idx back into list (if removed due
                   %to no auto-edges)
                   if assby_data.recips_trim_trimmed_flag
%                        fprintf('\n\tInserting removed row item %d back into assby_data.recips_trim_test\n',open_edge_data.row_select_idx)
                       assby_data.recips_trim_len = assby_data.recips_trim_len + 1;
                       assby_data.recips_trim(assby_data.recips_trim_len) = open_edge_data.row_select_idx;
                       
                   end
                   
                else
    %                if (A(row_select_idx,recips_trim))
    %                    error('logical if on A(row...): row full of edges but nnz says not...\n')
    %                end

                  if debug_level>=1

                       %num_open_edges = length(recips_trim) - nnz(A(row_select_idx,recips_trim));

                       %we define 'open_edges' as those recipients with suitable number of edges 
                       %compared to the current donor; e.g., # of donor edges > recipient...
                       %num_open_edges = nnz(current_row_diff>0);
                       num_open_edges = nnz(open_edge_data.current_row_select);

                       fprintf('\n\trow selected (%d) with %d open edge(s)...continuing...\n',...
                           open_edge_data.row_select_idx,num_open_edges)
                  end
                   open_edge_data.flag = 1;
                   
                   %return list of open edges in recips trim:
                   open_edge_data.recips_trim_open_edges = assby_data.recips_trim(open_edge_data.current_row_select);
                   
                   
                end
           end
           
        performance.open_search_std_num = performance.open_search_std_num + 1;   
        performance.open_search_std_time(performance.open_search_std_num) = toc(performance.open_search_std_start);
        open_edge_data.std_search_flag(performance.open_search_std_num) = open_edge_data.flag;
        
        
        %update history of searches: we use ratio of samples/max to indicate whether to switch
        %methods
        %if cur_row_sample exceeded max, trim it to ratio of 1otherwise
        if cur_row_sample > max_row_samples
            open_edge_data.std_search_try_ratio(performance.open_search_std_num) = 1;
            open_edge_data.search_history(1) = 1;
        else
            open_edge_data.std_search_try_ratio(performance.open_search_std_num) = cur_row_sample/max_row_samples;
            open_edge_data.search_history(1) = cur_row_sample/max_row_samples;
        end
        
        %save for debug
        open_edge_data.std_search_history{performance.open_search_std_num} = open_edge_data.search_history;

    %"quick" search mode: extracts a sub-A and pulls out potential open edges from reduced set
    %this may or may not actually _be_ quick compared to above exhaustive, but in general as the
    %number of available open edges drops precipitously (i.e., towards the end of assembly) this
    %method can be dramatically faster
    else
        
        performance.open_search_alt_start = tic;
        
        if debug_level >= 1
            if open_edge_data.calc_weights
                %this mode is _forced_
                fprintf('permA_if_open_edge: "quick" search mode >>forced<< for calculating weights\n')
            else
                fprintf('permA_if_open_edge: "quick" search mode requested\n')
            end
        end
        
        %we assemble a sub-matrix from the full A filtered
        %through this current donor-node's filled edge rows, and available recipient columns
        %this method then tests all possible rows from current donor for open edges - either solo ->
        %empty, or multi->filled
        
       %for ease in later preventing auto-edges, split into solo and multi-edge components
       donor_edge_rows_solo = matrix_data.node_idx(matrix_data.Ai(:,assby_data.donor_col)==1);
       donor_edge_rows_multi = matrix_data.node_idx(matrix_data.Ai(:,assby_data.donor_col)>1);
       
       %sift out sub-matrix for both solo and multi components
       A_sub_solo = matrix_data.Ai(donor_edge_rows_solo,assby_data.recips(1:assby_data.num_recips));
       A_sub_multi = matrix_data.Ai(donor_edge_rows_multi,assby_data.recips(1:assby_data.num_recips));

       %if auto-edges disallowed, need to modify above sub-matrices:
       %solo edges in donor col: prevent auto-connects via filling (with non-zero '1')
       %multi edges in donor col: prevented via emptying with zero entries
       
       if ~auto_edge_data.flag
       
           %below sub-matrices simply map over auto-connect diagonals to the above sub-matrix indexing
           A_sub_solo_eye = logical(open_edge_data.auto_edge_fill(donor_edge_rows_solo,assby_data.recips(1:assby_data.num_recips)));
           A_sub_multi_eye = logical(open_edge_data.auto_edge_fill(donor_edge_rows_multi,assby_data.recips(1:assby_data.num_recips)));

           %modify the sub-matrices with entries <-> auto-fills as per above solo / multi reqs
           A_sub_solo(A_sub_solo_eye) = 1;
           A_sub_multi(A_sub_multi_eye) = 0;
       
       end
       
       %test if any edges in either solo or multi-sub arrays available:
       %solo: must be empty or '0'; multi: must be filled or at least '1'
       A_sub_solo_open = ~all(A_sub_solo,2); %'all' returns true or 1 if all elements are filled/nonzero
                                                %we negate to indicate 'open' elements: solo-edges can only swap
                                                %with open entries
       A_sub_multi_open = any(A_sub_multi,2); %'any' returns true or 1 is any elements are non-zero or filled; multi-edges
                                                %can only swap with filled entries
        
        %get index of these open-edge rows and the sub-A matrix
        A_sub_solo_open_idx = matrix_data.node_idx(A_sub_solo_open);
        A_sub_multi_open_idx = matrix_data.node_idx(A_sub_multi_open);
       
        %test if any open edges in either solo or multi idx
       open_edge_data.A_sub_flag = 0;
       if any(A_sub_solo_open) || any(A_sub_multi_open)
            open_edge_data.flag = 1;
            %open_edge_data.A_sub_flag = 1;
           
           %since have open edge, select one of these rows in either solo / multi for swapping
           %get number of rows with open edges 
           A_sub_solo_open_num = nnz(A_sub_solo_open);
           A_sub_num_multi_open = nnz(A_sub_multi_open);

           %select random row among both the solo and multi-open sub matrices
           A_sub_row_select_open = randi(A_sub_solo_open_num+A_sub_num_multi_open);
           
           %determine if a solo or multi selection
           A_sub_select_type = 0;  %'0' <-> solo, '1' <-> multi
           if (A_sub_row_select_open > A_sub_solo_open_num)
               A_sub_select_type = 1;
               
               %get index from the row select for open sub to the full sub-A; shifted back by number
               %of open solo's since randi gave us value in range of multi's
               A_sub_row_select_open_idx = A_sub_multi_open_idx(A_sub_row_select_open - A_sub_solo_open_num);
               
               %with this selected row, sift out 'open' edges; for multi-these need to be non-zero
               %note this should already be modified to avoid self-connects
               %A_sub_recips_trim_open_edges = assby_data.recips(logical(A_sub_multi(A_sub_row_select_open_idx,:)));
                open_edge_data.recips_trim_open_edges = assby_data.recips(logical(A_sub_multi(A_sub_row_select_open_idx,:)));
%                open_edge_data.A_sub_recips_trim_open_edges = assby_data.recips(logical(A_sub_multi(A_sub_row_select_open_idx,:)));
               
               %map this selected row back to full A index
               %A_sub_row_select_idx = donor_edge_rows_multi(A_sub_row_select_open_idx);
               open_edge_data.row_select_idx = donor_edge_rows_multi(A_sub_row_select_open_idx);
%               open_edge_data.A_sub_row_select_idx = donor_edge_rows_multi(A_sub_row_select_open_idx);
           
               
           else
               
               %get index from the row select for open sub to the full sub-A
               A_sub_row_select_open_idx = A_sub_solo_open_idx(A_sub_row_select_open);
               
               %with this selected row, sift out 'open' edges; for solo-these need to be zeros
               %note this should already be modified to avoid self-connects
               %A_sub_recips_trim_open_edges = assby_data.recips(~logical(A_sub_solo(A_sub_row_select_open_idx,:)));
                open_edge_data.recips_trim_open_edges = assby_data.recips(~logical(A_sub_solo(A_sub_row_select_open_idx,:)));
%                open_edge_data.A_sub_recips_trim_open_edges = assby_data.recips(~logical(A_sub_solo(A_sub_row_select_open_idx,:)));
               
               %map this selected row back to full A index
               %A_sub_row_select_idx = donor_edge_rows_solo(A_sub_row_select_open_idx);
                open_edge_data.row_select_idx = donor_edge_rows_solo(A_sub_row_select_open_idx);
%                open_edge_data.A_sub_row_select_idx = donor_edge_rows_solo(A_sub_row_select_open_idx);
               
               
               
           end
           
       end
       
        performance.open_search_alt_num = performance.open_search_alt_num + 1;   
        performance.open_search_alt_time(performance.open_search_alt_num) = toc(performance.open_search_alt_start);
        open_edge_data.alt_search_flag(performance.open_search_alt_num) = open_edge_data.A_sub_flag;
       
       
        %for debug_tracking, save the solo and multi open edge flag arrays
        if open_edge_data.debug_tracking
           
        %if open_edge_data.debug_tracking
        %end
        
%            if debug_level >=2 
%               fprintf('\n\tdebug_tracking: total number open edges this attempt: ')
%                
%            end
           
           %update counts and compute total possible open edges 
           open_edge_data.track_count = open_edge_data.track_count + 1;
           open_edge_data.track_donor_col(open_edge_data.track_count) = assby_data.donor_col;
           open_edge_data.track_tot_open_edge(open_edge_data.track_count) = ...
               nnz(1 - A_sub_solo);
           
%            if debug_level >=2 
%            
%                 fprintf(' %d edges found\n',open_edge_data.tot_open_edge{open_edge_data.track_count});
%            end           
           
           %as well as tot number open rows - that we selected one of (if
           %open edges)
           %and tot number open (recipient) cols
           if open_edge_data.flag
               open_edge_data.track_num_open_rows(open_edge_data.track_count) = A_sub_solo_open_num;
               open_edge_data.track_num_open_cols(open_edge_data.track_count) = length(open_edge_data.recips_trim_open_edges);
           else
               open_edge_data.track_num_open_rows(open_edge_data.track_count) = 0;
               open_edge_data.track_num_open_cols(open_edge_data.track_count) = 0;
           end
           
           if debug_level >=2 
              fprintf('\n\tdebug_tracking: tot open edges %d, selected row (%d) with %d\n',...
                  open_edge_data.track_tot_open_edge(open_edge_data.track_count),open_edge_data.row_select_idx,...
                  open_edge_data.track_num_open_rows(open_edge_data.track_count))

%                   open_edge_data.tot_open_edge(open_edge_data.track_count),open_edge_data.row_select_idx,...
%                   open_edge_data.track_num_open_rows(open_edge_data.track_count))
               
           end
           
           
           
           if debug_level >=2 
              fprintf('\n\tdebug_tracking: Saving A_sub_solo_open and A_sub_multi_open flag arrays...\n')
               
           end
           
           
           open_edge_data.tracking{assby_data.num_updates}.A_sub_solo = A_sub_solo;
           open_edge_data.tracking{assby_data.num_updates}.A_sub_multi = A_sub_multi;
           %open_edge_data.tracking{assby_data.num_updates}.A_sub_solo_open = A_sub_solo_open;
           %open_edge_data.tracking{assby_data.num_updates}.A_sub_multi_open = A_sub_multi_open;
           
           %save the size of potential rows for selection - across all open cols
           open_edge_data.tracking{assby_data.num_updates}.open_row_solo_sums = sum(~logical(A_sub_solo));
           open_edge_data.tracking{assby_data.num_updates}.open_row_multi_sums = sum(logical(A_sub_multi));
           %and potential cols for selection
           open_edge_data.tracking{assby_data.num_updates}.num_open_solo_cols = ~all(A_sub_solo);
           open_edge_data.tracking{assby_data.num_updates}.num_open_multi_cols = any(A_sub_multi);
            
            
        end
        
        %if user-requested weight calculation returns do so here
        if open_edge_data.calc_weights
            
            %bump counter
            open_edge_data.weight_count = open_edge_data.weight_count + 1;
            
            %append these row and col values - if we have an open edge, if not simply append '0'
            if open_edge_data.flag
                open_edge_data.weight_num_open_rows(open_edge_data.weight_count) = ...
                    A_sub_solo_open_num;
                open_edge_data.weight_num_open_cols(open_edge_data.weight_count) = ...
                    length(open_edge_data.recips_trim_open_edges);
                
            else
                
                open_edge_data.weight_num_open_rows(open_edge_data.weight_count) = 0;
                open_edge_data.weight_num_open_cols(open_edge_data.weight_count) = 0;
                
                
            end
            
            
        end
       
    end
    
    
    %update history tracker
    %open_edge_data.search_history(1) = open_edge_data.flag;
    %open_edge_data.search_history(1) = open_edge_data.flag;




end

function [out_flag, out_msg, assby_data, inert_swap_data, matrix_data] = permA_update_state_Ai(kout,matrix_data,assby_data,inert_swap_data,state)
%internal function calculating current column sum, or out-degrees, for interim Ai, and L2 norm of
%distance from target kout, depending on 'state':
%   state: 1 <-> 'initialise'; 2 <-> 'update col_sum_diff, recip list';  3 <-> 'update donor list';

    global N debug_level
    
    out_flag = 0;
    out_msg = 'permA_update_state_Ai';
    out_data = [];

    if debug_level >= 1
        fprintf('permA_update_state_Ai called...\n');
    end
    
    
    %if assby_data.num_updates < 1
    switch state
        
        case 1 %if state flag == 1; initialise
        
        if debug_level >= 1
            fprintf('\n\t\tInitial update/init A_col_sum and node type...\n')
            
        end

        %get current out-degrees for this Ai - via column sum
        assby_data.A_col_sum = sum(matrix_data.Ai);
        assby_data.col_sum_diff = assby_data.A_col_sum - kout';
        %assby_data.col_sum_diff_sum = sum(assby_data.col_sum_diff); 

        assby_data.num_updates = assby_data.num_updates + 1;

        %update L2 norm distance from target kout
        assby_data.norm_targ_dist(assby_data.num_updates) = norm(assby_data.col_sum_diff);

        %get indices to donors and recipients
        assby_data.donor_idx = assby_data.col_sum_diff > 0;
        assby_data.recips_idx = assby_data.col_sum_diff < 0;
        
        %and inerts
        assby_data.inert_idx = assby_data.col_sum_diff == 0;
    
        assby_data.num_donors = nnz(assby_data.donor_idx);
        assby_data.prev_num_donors = assby_data.num_donors;
        assby_data.num_recips = nnz(assby_data.recips_idx);
        
        assby_data.num_inerts = nnz(assby_data.inert_idx);
        
        assby_data.donors(1:assby_data.num_donors) = matrix_data.node_idx(assby_data.donor_idx);
        assby_data.recips(1:assby_data.num_recips) = matrix_data.node_idx(assby_data.recips_idx);
        
        assby_data.inerts(1:assby_data.num_inerts) = matrix_data.node_idx(assby_data.inert_idx);
        
        %initialise recips_trim to same recipient list and exclude trimmed out nodes during
        %open edge finding
        assby_data.recips_trim(1:assby_data.num_recips) = assby_data.recips(1:assby_data.num_recips);
        assby_data.recips_trim_len = assby_data.num_recips;
        
        %map over donor and recip lists for inert swap 
        inert_swap_data.inert_idx_topick(1:assby_data.num_inerts) = assby_data.inerts(1:assby_data.num_inerts);
        inert_swap_data.recip_idx_topick(1:assby_data.num_recips) = assby_data.recips(1:assby_data.num_recips);
        inert_swap_data.donor_idx_topick(1:assby_data.num_donors) = assby_data.donors(1:assby_data.num_donors);
        
        inert_swap_data.num_inert_topick = assby_data.num_inerts;
        inert_swap_data.num_donors_topick = assby_data.num_donors;
        inert_swap_data.num_recips_topick = assby_data.num_recips;
        
%         %testing arrays
%         %assby_data.donors_test = assby_data.donors;
%         %assby_data.donor_test_idx = assby_data.donor_idx;
         assby_data.num_recips_test = assby_data.num_recips;
        assby_data.recips_test = assby_data.recips;
        assby_data.recips_test_idx = assby_data.recips_idx;
        
%         assby_data.recips_trim_test(1:assby_data.num_recips) = assby_data.recips(1:assby_data.num_recips);
%         assby_data.recips_trim_test_len = assby_data.num_recips;
        
        %if debug tracking is on, save this current 'Ai'
        if matrix_data.debug_tracking
            if debug_level >= 2
                fprintf('\n\tInitialising matrix_data.tracking{%d}...\n\n',assby_data.num_updates)
            end
            %matrix_data.Ai
            
            matrix_data.tracking{assby_data.num_updates} = matrix_data.Ai;
            
            %and donor / recipient lists
            assby_data.tracking{assby_data.num_updates}.donor_idx = assby_data.donor_idx;
            assby_data.tracking{assby_data.num_updates}.recip_idx = assby_data.recips_idx;
            
        end

        
        
        
    %otherwise, update individual nodes for speed
    %else
        %update col_sum_diff, recip list
        case 2
        
        if debug_level >= 1
            fprintf('\n\t\tUpdating A_col_sum and recip lists...\n')
            
        end
        
        %if debug_level >= 3
            
%             assby_data.test_A_col_sum = sum(matrix_data.Ai);
%             assby_data.test_col_sum_diff = assby_data.A_col_sum - kout';
%             
%             if nnz(assby_data.test_col_sum_diff - assby_data.col_sum_diff) > 0
%                
%                 testdiff = assby_data.test_col_sum_diff - assby_data.col_sum_diff;
%                 testdiff_idx = testdiff ~=0
%                error('col sum diff deviation (test1)') 
%                
%                 
%             end
            
            
        %end
        
%         fprintf('A_col_sum pre-update\n')
%         assby_data.A_col_sum

        
        %increment / decrement the donor and target edge col tracking
        assby_data.A_col_sum(assby_data.donor_col) = assby_data.A_col_sum(assby_data.donor_col) - 1;
        assby_data.A_col_sum(assby_data.target_open_edge_col) = assby_data.A_col_sum(assby_data.target_open_edge_col) + 1;
        
%         fprintf('A_col_sum post-update\n')        
%         assby_data.A_col_sum
% 
%         assby_data.A_col_sum - kout'
        
%         fprintf('col_sum_diff pre-update\n')        
%         assby_data.col_sum_diff
        
        assby_data.col_sum_diff(assby_data.donor_col) = assby_data.col_sum_diff(assby_data.donor_col) - 1;
        assby_data.col_sum_diff(assby_data.target_open_edge_col) = assby_data.col_sum_diff(assby_data.target_open_edge_col) + 1;
        
%         fprintf('col_sum_diff post-update\n')                
%         assby_data.col_sum_diff        
        
        if debug_level >= 3
            fprintf('\n\tupdated col_sum_diff: \n');
            assby_data.col_sum_diff
        end
        
        %if debug_level >= 3
        
%             assby_data.test_A_col_sum = sum(matrix_data.Ai);
%             assby_data.test_col_sum_diff = assby_data.A_col_sum - kout';
%         
%         
%             if nnz(assby_data.test_col_sum_diff - assby_data.col_sum_diff) > 0
%                 testdiff = assby_data.test_col_sum_diff - assby_data.col_sum_diff;
%                 testdiff_idx = matrix_data.node_idx(testdiff ~=0)
%                error('col sum diff deviation (test2)') 
% 
%             end
            
        %end
        
        
        %update tracking vars
        assby_data.num_updates = assby_data.num_updates + 1;
        
        %if debug tracking is on, save this current 'Ai'
        if matrix_data.debug_tracking
            if debug_level >= 2
                fprintf('\n\tUpdating matrix_data.tracking{%d}...\n\n',assby_data.num_updates)
            end
            
            matrix_data.tracking{assby_data.num_updates} = matrix_data.Ai;

            %and donor / recipient lists
            assby_data.tracking{assby_data.num_updates}.donor_idx = assby_data.donor_idx;
            assby_data.tracking{assby_data.num_updates}.recip_idx = assby_data.recips_idx;
            
            
        end
        
        

        %update L2 norm distance from target kout
        assby_data.norm_targ_dist(assby_data.num_updates) = norm(assby_data.col_sum_diff);
        
        if debug_level >=2
           fprintf('\n\tpermA_update_state_Ai: distance: %2.3f...\n',assby_data.norm_targ_dist(assby_data.num_updates))
        end
        
        %update recipient list - this swap may fulfil current recipient
        %target kout
        %may update only on current recip; for now do all
        %first index:
%         assby_data.recips_idx = assby_data.col_sum_diff < 0;
%         assby_data.recips = matrix_data.node_idx(assby_data.recips_idx);
%         assby_data.num_recips = nnz(assby_data.recips_idx);
        
        %update recipient list based on current change to 'target_open_edge_col'; if
        %col_sum_diff(target_open_edge_col) == 0, then remove from recip list
        if assby_data.col_sum_diff(assby_data.target_open_edge_col) == 0
            
           if debug_level>=1
                fprintf('\n\tpermA_update_state_Ai: recipient col %d edge count filled...removing from list\n',assby_data.target_open_edge_col)
                size(assby_data.recips);
           end
           
            %test for num_recip error:
            pre_num_recip_col_sum = nnz(assby_data.col_sum_diff < 0);
            pre_num_donor_col_sum = nnz(assby_data.col_sum_diff > 0);
    %         fprintf('num recips (logic test): %d, num recips (assby_data): %d\n',num_recip_col_sum,assby_data.num_recips)
    %         fprintf('num donors (logic test): %d, num donors (assby_data): %d\n',num_donor_col_sum,assby_data.num_donors)

%             if (num_recip_col_sum ~= assby_data.num_recips)
%                  fprintf('pre update: num recips (logic test): %d, num recips (assby_data): %d\n',num_recip_col_sum,assby_data.num_recips)
%                  fprintf('pre update: num donors (logic test): %d, num donors (assby_data): %d\n',num_donor_col_sum,assby_data.num_donors)
% 
%                 %error('recip number deviation')
%             end
           
           
%            size(assby_data.recips(1:assby_data.num_recips))
%            size(assby_data.recips(assby_data.recips(1:assby_data.num_recips) ~= assby_data.target_open_edge_col))
           
%            assby_data.recips(1:assby_data.num_recips-1) = ...
%                assby_data.recips(assby_data.recips(1:assby_data.num_recips) ~= assby_data.target_open_edge_col);

           %
           %below statement changes size of recips array...this is not intentional... ToDo list:
           %change removal of this current filled recipient such that array size is fixed avoiding
           %memory resizing...
           %
           assby_data.recips(assby_data.recips == assby_data.target_open_edge_col) = [];


%            %try instead:
%            assby_data.recips_test(assby_data.recips_test == assby_data.target_open_edge_col) = [];
%            assby_data.num_recips_test = assby_data.num_recips_test - 1;
%            assby_data.recips_test_idx(assby_data.target_open_edge_col) = 0;
           
           % update / decrement the num_recips tracker
           assby_data.num_recips = assby_data.num_recips - 1;
           
           %further swap flag for this recip col to '0' removing from recips_idx list
           assby_data.recips_idx(assby_data.target_open_edge_col) = 0;
           
           %alternative method: update idx vector first, then recips list
           assby_data.recips_test_idx(assby_data.target_open_edge_col) = 0;
           assby_data.num_recips_test = nnz(assby_data.recips_test_idx);
           assby_data.recips_test(1:assby_data.num_recips_test) = matrix_data.node_idx(assby_data.recips_test_idx);
           
           %simply map over the recips list to recips_trim
           assby_data.recips_trim(1:assby_data.num_recips) = assby_data.recips(1:assby_data.num_recips);
           assby_data.recips_trim_len = assby_data.num_recips;
%            assby_data.recips_trim_test(1:assby_data.num_recips) = assby_data.recips(1:assby_data.num_recips);
           
%            if debug_level >= 2
%                
%               fprintf('recips_list test: num_recips / test: %d/%d\n',assby_data.num_recips,assby_data.num_recips_test)
%               recips_diff = assby_data.recips(1:assby_data.num_recips) - assby_data.recips_test(1:assby_data.num_recips_test);
%               fprintf('diff recips list nnz: %d\n',nnz(recips_diff)); 
%               
%            end

           
        else
           if debug_level>=1
                fprintf('\n\tpermA_update_state_Ai: recip col %d edge fill remains...retaining in list\n',assby_data.target_open_edge_col)
           end
        end
        
%         % map over the recips list to recips_trim; looks like the 'test' version must be updated
%         % each time...
%         %assby_data.recips_trim(1:assby_data.num_recips) = assby_data.recips(1:assby_data.num_recips);
%         assby_data.recips_trim_test(1:assby_data.num_recips) = assby_data.recips(1:assby_data.num_recips);
%         assby_data.recips_trim_test_len = assby_data.num_recips;

        
        
        %cannot do below donor update combined with recips update _inside_ the donor loop!
        
%         %update donor list here since we now have column utilised - easier to track!
%         if assby_data.col_sum_diff(assby_data.donor_col) == 0
%             %we've exhausted the surplus edges of this donor column - remove from donor list
%             if debug_level >= 1
%                fprintf('\n\tpermA_update_state_Ai: donor %d edge surplus exhausted; removing from donor list...\n',assby_data.donor_col)
%                 
%             end
%             
% %             assby_data.donors(1:assby_data.num_donors-1) = ...
% %                 assby_data.donors(assby_data.donors(1:assby_data.num_donors) ~= assby_data.donor_col);
%             %try instead:
%             assby_data.donors(assby_data.donors_test == assby_data.donor_col) = [];
%             
%             assby_data.prev_num_donors = assby_data.num_donors;
%             assby_data.num_donors = assby_data.num_donors - 1;
%             
%             %turn off flag for this donor in index 
%             assby_data.donor_idx(assby_data.donor_col) = 0;
% %             assby_data.donor_test_idx(assby_data.donor_col) = 0;
%             
%         end
        
        %test for num_recip error:
        %note this only seems to happen in the matlab gui and not batch
        %mode at command line...?!
        num_recip_col_sum = nnz(assby_data.col_sum_diff < 0);
        num_donor_col_sum = nnz(assby_data.col_sum_diff > 0);
%         fprintf('num recips (logic test): %d, num recips (assby_data): %d\n',num_recip_col_sum,assby_data.num_recips)
%         fprintf('num donors (logic test): %d, num donors (assby_data): %d\n',num_donor_col_sum,assby_data.num_donors)
        
        if (num_recip_col_sum ~= assby_data.num_recips)
             fprintf('pre update: num recips (logic col sum diff): %d\n',pre_num_recip_col_sum);
             fprintf('post update: num recips (logic test): %d, num recips (assby_data): %d\n',num_recip_col_sum,assby_data.num_recips)
             fprintf('post update: num donors (logic test): %d, num donors (assby_data): %d\n',num_donor_col_sum,assby_data.num_donors)
            
            error('recip number deviation -- try using command line MATLAB and NOT the GUI')
        end


        %update donor list
        case 3
        
            if debug_level >= 1
                fprintf('\n\t\tUpdating donor lists...\n')

            end
            assby_data.prev_num_donors = assby_data.num_donors;
            assby_data.donor_idx = assby_data.col_sum_diff > 0;
            assby_data.num_donors = nnz(assby_data.donor_idx);
            assby_data.donors(1:assby_data.num_donors) = matrix_data.node_idx(assby_data.donor_idx);
            
            if debug_level >= 1
                fprintf('\t\tprevious donor num: %d, post-loop donor num: %d\n',assby_data.prev_num_donors,assby_data.num_donors)
            end
            
            
        %update inert swap lists
        case 4
            
            if debug_level >= 1
                fprintf('\n\t\tUpdating inert swap lists...\n')

            end
            
            %update inert swapping lists (may invoke this only prior to calling inert_swap)
            inert_swap_data.num_donors_topick = assby_data.num_donors;
            inert_swap_data.donor_idx_topick(1:assby_data.num_donors) = assby_data.donors(1:assby_data.num_donors);
            
            inert_swap_data.num_recips_topick = assby_data.num_recips;
            inert_swap_data.recip_idx_topick(1:assby_data.num_recips) = assby_data.recips(1:assby_data.num_recips);
            
            assby_data.inert_idx = assby_data.col_sum_diff == 0;
            assby_data.num_inerts = nnz(assby_data.inert_idx);
            assby_data.inerts(1:assby_data.num_inerts) = matrix_data.node_idx(assby_data.inert_idx);
            
            inert_swap_data.inert_idx_topick(1:assby_data.num_inerts) = assby_data.inerts(1:assby_data.num_inerts);
            inert_swap_data.num_inert_topick = assby_data.num_inerts;
            
          
            
    end

end

function [out_flag, out_msg, matrix_data] = permA_permute_A1(kin,matrix_data,auto_edge_data)
%internal function taking A0 as input and permuting its rows into the method's A1 matrix - the first
%step towards the 'Af' satisfying kin _and_ kout

    global N debug_level
    
    out_flag = 0;
    out_msg = 'assemble_A0: ';
    out_data = [];
    
    if debug_level >= 1
        fprintf('permA_permute_A1 called...\n');
    end
    
    %pre-allocate vector for node-indexing if auto-connects/self-loops disallowed
    if ~~auto_edge_data.flag
        node_idx_cur = zeros(1,N-1);
    end
    
    %simply loop over all nodes (rows) of A0, permuting them into (uniform) random distribution of
    %its entries into the method's A1 matrix
    for cur_node = 1:N
        
       %if auto-connects / self-loops disallowed, exclude current node from permuting
       if ~auto_edge_data.flag
           node_idx_cur = matrix_data.node_idx(matrix_data.node_idx ~= cur_node);
           %note: if no auto-connects, then max(kin) = N-1...
           matrix_data.A1(cur_node,node_idx_cur) = matrix_data.A0(cur_node,randperm(N-1));
       else
           matrix_data.A1(cur_node,:) = matrix_data.A0(cur_node,randperm(N));
           
       end
        
        
    end
    
    %quick error check: above permutation should not change kin...
    if (max(sum(matrix_data.A1,2) - kin) > 0)
        error('permuted A1 kin error...\n');
    end
    

end

function [out_flag, out_msg, matrix_data, multi_edge_data] = permA_assemble_A0(kin,kout,matrix_data,multi_edge_data,auto_edge_data)
%internal function generating initial matrix 'A0'

    global N debug_level

    out_flag = 0;
    out_msg = 'assemble_A0: ';
    out_data = [];
    
    %define emtpy A0 given N
    %A0 = zeros(N,N);
    
    if debug_level >=1
        if multi_edge_data.flag
            fprintf('\n\t\tMulti-edges requested: %2.3f percent total\n',multi_edge_data.targ_prop);
        else
            fprintf('\n\t\tMulti-edges disallowed; Solo-connections only\n')
        end
        if auto_edge_data.flag
            fprintf('\n\t\tAuto-edges permitted\n');
        else
            fprintf('\n\t\tAuto-edges disallowed\n');
        end
    end
    
%     %assemble node edge list (solo-only for start)
%     for cur_node = 1:N
%         
%         matrix_data.node_edges{cur_node} = zeros(1,N);
%         matrix_data.node_edges{cur_node}(1:kin(cur_node)) = 1;
%         
%     end
    
    %given target percentage of multi-edges, compute how many each node should have
    multi_edge_data.kin_mult = round(kin*multi_edge_data.targ_prop);
    %get maximum number of multi-edges
    multi_edge_data.kin_mult_max = max(multi_edge_data.kin_mult);
    
    if debug_level >= 2
        fprintf('\n\tmulti_edge_data.kin_mult: Num multi-edge: %d\n',sum(multi_edge_data.kin_mult))
        multi_edge_data.kin_mult'
        
    end

    multi_edge_data.skip_multi = 0;
    %test for zero multi: if percentage low enough, may not have any and skip all this multi-edge stuff
    if nnz(multi_edge_data.kin_mult) < 1
        if multi_edge_data.flag
            fprintf('\n\tZero multi-edges result given target multi-edge percentage %2.3f...\n',multi_edge_data.targ_prop);
            fprintf('\n\t>>> Skipping multi-edge processing...\n')
        end
        multi_edge_data.skip_multi = 1;
        multsumright = 1;
    end
    
    %if sufficient multi-edges process them here
    if ~multi_edge_data.skip_multi
        
        if debug_level >= 1
            fprintf('\n\tprocessing multi-edges...\n')
        end

        
        %loop over all nodes assembling list of multi-edge values for each
        %rand selection from: [2, #mult-edge-1] such that Sum of list = 
        %#mult-edge

        num_mult_pick_fails = 0;
        most_pick_tries = -1;
        
        %in-degree error count
        kin_err = 0;
        kin_err_num = 0;
        kin_err_node = zeros(1,N);
        
        %preallocate space for multi-edge list
        mult_list2 = zeros(1,multi_edge_data.kin_mult_max);
        
         for cur_node = 1:N
            if debug_level >= 1
                fprintf('\n\tprocessing node %d with %d multi-edges...\n',cur_node,multi_edge_data.kin_mult(cur_node))
            end
             
             
            %for this node's mult number, assemble vector of possible mult-edge values
            cur_node_mult_edge_val = [1:multi_edge_data.kin_mult(cur_node)];
      
            %compute a distribution for these potential multi-edge vals
            %given exponential parm 'lambda'
            cur_node_mult_edge_val_distn = ...
                (1/multi_edge_data.lambda)*exp(-(1/multi_edge_data.lambda)*cur_node_mult_edge_val);
            
            if debug_level >= 2
                fprintf('\n\tCur nodes multi-edge vals: \n')
                cur_node_mult_edge_val
                fprintf('\n\tCur nodes multi-edge distn: \n')
                cur_node_mult_edge_val_distn
            end
            
            
           %randomly pick one of these values, add to list, test sum, keep going if not right
           multsumright = 0;
           num_cur_mult = 0;
           %mult_list = [];
           mult_list_sum = 0;
           max_tries = multi_edge_data.max_pick_tries;
           num_tries = 0;
           
           %generate pool of randoms shaped to distn prior to looping
           max_cur_mult_edge_val = max(cur_node_mult_edge_val_distn);
           rand_picks = max_cur_mult_edge_val*rand(1,max_tries+1);
           
           while ~multsumright && num_tries <= max_tries
               %pick from list using the 'exponential' distn 
               %temprand = rand(1,1);
               tempmult_idx = nnz(rand_picks(num_tries+1) < cur_node_mult_edge_val_distn);
               %using find instead; this appears consistently slower than above nnz() method
%               tempmult_idx = find(rand_picks(num_tries+1)<cur_node_mult_edge_val_distn,1,'last');
              
               tempmult = cur_node_mult_edge_val(tempmult_idx);
               
                if debug_level >= 2
                    %fprintf('\n\ttempmult_idx way 1: %d, way 2: %d...\n',tempmult_idx,tempmult_idx2);
                    fprintf('\n\tcurrent selection "tempmult" %d\n',tempmult)
%                     if (tempmult_idx - tempmult_idx2) ~= 0
%                         error('tempmult_idx different...')
%                     end
                end
               

               %tempmultsum = sum([mult_list tempmult]);
               tempmultsum = mult_list_sum + tempmult;
               if tempmultsum == multi_edge_data.kin_mult(cur_node)
                   %hit the right sum!
                   %add this mult to the list
                   %mult_list = [mult_list tempmult];
                   
                   num_cur_mult = num_cur_mult + 1;
                   mult_list2(num_cur_mult) = tempmult;

                   mult_list_sum = mult_list_sum + tempmult;
                   multsumright = 1;
                   
                    if debug_level >= 2
                        fprintf('\n\ttempmultsum (%d) correct!\n',tempmultsum)
                    end
                   
               elseif tempmultsum > multi_edge_data.kin_mult(cur_node)
                   %sum too large - reject this pick
                   num_tries = num_tries + 1;
                   
                    if debug_level >= 2
                        fprintf('\n\ttempmultsum (%d) too large, rejecting...\n',tempmultsum)
                    end
                   
               else
                   %sum too small - add this pick
                   %mult_list = [mult_list tempmult];
                   
                   num_cur_mult = num_cur_mult + 1;
                   mult_list2(num_cur_mult) = tempmult;
                   
                   %update running sum of mult_list with this pick
                   mult_list_sum = mult_list_sum + tempmult;
                   num_tries = num_tries + 1;

                   if debug_level >= 2
                       fprintf('\n\ttempmultsum (%d) too small, appending and continuing...\n',tempmultsum)
                   end
                   

               end
               
               
           end
           
           %if got the sum right, add list to the node_edge data
           if multsumright
               
                if debug_level >= 1
                    fprintf('\n\tMultisum right; final multi-edge list:\n')
                    %mult_list;
                    mult_list2(1:num_cur_mult)
                    
%                     if nnz(mult_list - mult_list2(1:num_cur_mult)) > 0
%                         error('mult_list & mult_list2 mismatch')
%                     end
                end
               
               
               %compute number of solo edges given this mult_list
               %num_solo = kin(cur_node) - (multi_edge_data.kin_mult(cur_node) + length(mult_list));
               
               num_solo2 = kin(cur_node) - (multi_edge_data.kin_mult(cur_node) + num_cur_mult);
               
               %num_solo_diff = num_solo - num_solo2;
%                if debug_level >=1
%                    fprintf('num_solo: %d, num_solo2: %d\n',num_solo,num_solo2);
%                    if num_solo_diff ~= 0
% 
%                       error('num_solo diff error.') 
%                    end
%                end
               %if above num_solo is negative, we have too many edges for this
               %node; i.e., we are violating kin_c
               if num_solo2 < 0
                   
                  if debug_level >= 1 
                      fprintf('\n\tCurrent node %d with too many edges (num_solo = %d); correcting...\n',...
                          cur_node,num_solo2)
                  end

                  %we remove edges here - number removed == num_solo...
                  %loop back over mult_list removing edges until num_solo = 0;
                  %this won't work if we have insufficient multi-edges - must be at
                  %least one!
                  %len_multi_list = length(mult_list);
                  %if len_multi_list < 1
                  if num_cur_mult < 1
                     %fprintf('\n\tInsufficient multi-edges to repair edge-count...\n');
                     mult_edge_repair_failed = 1;
                     
                  else
                     mult_edge_repair_failed = 0;
                     
                     %get sorted list of the multi-edges; ascending
                     %mult_list_sorted = sort(mult_list);
                     
                     mult_list2(1:num_cur_mult) = sort(mult_list2(1:num_cur_mult));
                     
%                     if debug_level >= 1
%                         if nnz(mult_list_sorted - mult_list2(1:num_cur_mult)) > 0
%                             error('mult_list_sorted & mult_list2 mismatch')
%                         end
%                     end
                     

                     num_mult_edge_remove = 0;
                     %removed_mult_edge = [];
                     removed_mult_edge2 = zeros(1,num_cur_mult);
                     %cur_multi_fix = len_multi_list;
                     %multi_edge_num_over = num_solo;
                     
                     multi_edge_num_over2 = num_solo2;
                     
                     %while (num_mult_edge_remove <= len_multi_list) && (multi_edge_num_over <0)
                     while (num_mult_edge_remove <= num_cur_mult) && (multi_edge_num_over2 <0)
                     
                        num_mult_edge_remove = num_mult_edge_remove + 1; 
                        %cur_multi_edge = mult_list_sorted(num_mult_edge_remove);
                        cur_multi_edge2 = mult_list2(num_mult_edge_remove);
                        %fprintf('\n\tcurrent multi-edge: %d\n',cur_multi_edge)
                        %add this multi-edge targeted for removal to removed list

                        %removed_mult_edge(num_mult_edge_remove) = cur_multi_edge;
                        removed_mult_edge2(num_mult_edge_remove) = cur_multi_edge2;

                        %get sum of total edges targeted for removal
                        %sum_multi_edges_removal = sum(1+removed_mult_edge);
                        sum_multi_edges_removal2 = sum(1+removed_mult_edge2(1:num_mult_edge_remove));

                        %and number of multi-edges 'too much'
                        %multi_edge_num_over = num_solo + sum_multi_edges_removal;
                        multi_edge_num_over2 = num_solo2 + sum_multi_edges_removal2;
                         
                     
                     end
                     
%                      if debug_level >=1 
%                          fprintf('sum_multi_edges_removal: %d, sum_multi_edges_removal2: %d\n',...
%                             sum_multi_edges_removal,sum_multi_edges_removal2)
% 
%                      end
                     
                     %swap over these lists and values
%                      old_mult_list = mult_list;
%                      old_num_solo = num_solo;

                     %mult_list = mult_list_sorted((num_mult_edge_remove+1):end);
                     mult_list2 = mult_list2((num_mult_edge_remove+1):num_cur_mult);
                     
                     num_cur_mult = num_cur_mult - num_mult_edge_remove;
                     
                     %num_solo = sum_multi_edges_removal + num_solo;
                     num_solo2 = sum_multi_edges_removal2 + num_solo2;
                     
                    if debug_level >= 2
                        fprintf('\n\tCorrected multi-edge list:\n')
                        mult_list2
                        
                    elseif debug_level >= 1
                        fprintf('\n\tCorrected edge counts: multi %d, solo %d, kin: %d\n',...
                            sum(1+mult_list2(1:num_cur_mult)),num_solo2,kin(cur_node))
                        
%                        diff_mult_list = nnz(mult_list - mult_list2(1:num_cur_mult));
%                         if diff_mult_list > 0
%                             error('mult list difference > 0')
%                         end
%                         fprintf('sum_multi_edges_removal: %d, sum_multi_edges_removal2: %d\n',...
%                             sum_multi_edges_removal,sum_multi_edges_removal2)
                    end
                         
                  end

               end
               
               most_pick_tries = max(most_pick_tries,num_tries);

           else
               %mult_list = -1;
               num_mult_pick_fails = num_mult_pick_fails + 1;
               most_pick_tries = multi_edge_data.max_pick_tries;
               
           
           end
           
           %error check: number of edges should be correct for kin still
           %edge_num = sum(1+mult_list) + num_solo;
           edge_num2 = sum(1+mult_list2(1:num_cur_mult)) + num_solo2;
           
%            if debug_level>=1
%               fprintf('num_solo: %d, num_solo2: %d\n',num_solo,num_solo2)
%               edge_num_diff = edge_num - edge_num2;
%               fprintf('edgenum: %d, edgenum2: %d\n',edge_num,edge_num2)
%               if edge_num_diff ~= 0
%                   error('edge num diff error')
%               end
%               
%               %test difference in multi_lists
%               mult_list_diff = mult_list - mult_list2(1:num_cur_mult);
%               fprintf('mult_list_diff nnz: %d\n',nnz(mult_list_diff));
%               mult_list
%               mult_list2(1:num_cur_mult)
%                
%            end

            %kin_err = 0;
            %if above not == kin_c then have problem...
            if kin(cur_node) - edge_num2 ~= 0
%                 fprintf('\n\tError: node %d: edge num (%d) ~= kin (%d)...\n',cur_node,...
%                     edge_num2,kin(cur_node));
%             if kin(cur_node) - edge_num2 ~= 0
%                 fprintf('\n\tError: node %d: edge num (%d) ~= kin (%d)...\n',cur_node,...
%                     edge_num2,kin(cur_node));

                kin_err = 1;
                kin_err_num = kin_err_num + 1;
                kin_err_node(kin_err_num) = cur_node;
            end
            
            %load this nodes edges into A0
            
           %first load up multi-edge list
%            if debug_level >=1
%               fprintf('loading mult lists into A0:\n');
%               1+mult_list
%               1+mult_list2(1:num_cur_mult)
% %               fprintf('pre-load A0 values:\n')
% %               matrix_data.A0(cur_node,:)
% %               matrix_data.A02(cur_node,:)
%                
%            end
           %matrix_data.A0(cur_node,1:length(mult_list)) =  1+mult_list;
           matrix_data.A0(cur_node,1:num_cur_mult) =  1+mult_list2(1:num_cur_mult);

           %then pad remaining edges with solo's
%            if debug_level >=1
%               fprintf('loading solo lists into A0:\n');
%               index = [length(mult_list)+1:(length(mult_list)+num_solo)]
%               index2 = [num_cur_mult+1:(num_cur_mult+num_solo2)]
%                
%            end
           
%            matrix_data.A0(cur_node,(length(mult_list)+1):...
%                (length(mult_list)+num_solo)) =  1;
           matrix_data.A0(cur_node,(num_cur_mult+1):...
               (num_cur_mult+num_solo2)) =  1;
           
            
           
           
           
         end
         
    else
        
        %in-degree error count
        kin_err = 0;
        kin_err_num = 0;
        kin_err_node = zeros(1,N);
        
        %solo edge only assembly...
        for cur_node = 1:N
            matrix_data.A0(cur_node,1:kin(cur_node)) = 1;
            
            %error check - edge count this row / node should == kin(cur_node)
            if kin(cur_node) ~= sum(matrix_data.A0(cur_node,:))
                kin_err = 1;
                kin_err_num = kin_err_num + 1;
                kin_err_node(kin_err_num) = cur_node;
            end
            
        end
        
    end
    
    %if kin_err's flagged, display error message
    if kin_err
       fprintf('\n\t\tWARNING: Edge-count error: kin mismatch for %d nodes...',kin_err_num) 
       fprintf('\n\t\t\tRejecting A0...\n')
        
        
    end

    
    %post-assembly tests:
    %correct in-degree
    matrix_data.A0_kin = sum(matrix_data.A0,2);
    matrix_data.A0_kin_diff = (matrix_data.A0_kin - kin) ~= 0;
    matrix_data.A0_kin_errnum = nnz(matrix_data.A0_kin_diff);
    
    %disregard kout due to A0 design (it should not satisfy kout unless
    %we're quite lucky!)
    
    %calculate number multi-edges, if any
    matrix_data.A0_num_multi = sum(sum(matrix_data.A0(matrix_data.A0>1)-1));
    
    if debug_level >= 1
        fprintf('\n\tNum multi in A0: %d; num multi in kin_multi: %d\n',matrix_data.A0_num_multi,sum(multi_edge_data.kin_mult))
        
    end
    
    
    %resulting percentage of edges
     multi_edge_data.A0_pct_multi =   matrix_data.A0_num_multi/matrix_data.sum_kin;

    
    
    
end

function [out_flag, out_msg] = permA_test_kin_kout(kin,kout,matrix_data)
%internal function performing tests on given in & out degree sequences

    out_flag = 0;
    out_msg = 'graphicality test for kin & kout: ';
    
    global N
    
    %test (1) sum(kin) = sum(kout) 
%     sum_kin = sum(kin);
%     sum_kout = sum(kout);
    
    sum_test = 0;
    if (matrix_data.sum_kin == matrix_data.sum_kout)
        sum_test = 1;
        out_msg_str1 = ' ';
    else
        out_msg_str1 = 'failed sum test ';
    end
    
    %test (2) aka the Gale-Ryser Inequalities
    %
    %       Sum_i=1...n [min(kout(i),j]  >=   Sum_i=1...j [kin(i)]
    %
    %       for all j in [1,...,n-1]
    %loop over all columns
    for cur_col = 1:N-1

        %get Sum_i=1...j [kin(i)]
        kin_sum_j(cur_col) = sum(kin(1:cur_col));

        %get Sum_i=1...n [min(kout(i),j]
        kout_sum_min_i(cur_col) = 0;
        for cur_row = 1:N
           cur_term = min(kout(cur_row),cur_col);
           kout_sum_min_i(cur_col) = kout_sum_min_i(cur_col) + cur_term;

        end

    end

    %test for satisfying inequality
    ineq_kin_kout = nnz(kout_sum_min_i >= kin_sum_j); 
    
    ineq_test = 0;
    if ineq_kin_kout == length(kout_sum_min_i)
        ineq_test = 1;
        out_msg_str2 = ' ';
    else
        out_msg_str2 = 'failed inequality test ';
    end
    
    %final result
    if sum_test && ineq_test
        out_flag = 1;
        out_msg = [out_msg 'success'];
    else
        out_flag = 0;
        out_msg = [out_msg out_msg_str1 out_msg_str2];
    end
    



end

function [out_flag, out_msg, passed_parms] = ...
    permA_parse_args(kin,kout,nargin,varargin_local,nargout)
%internal function for parsing out input arguments into internal structures

    %fprintf('permA_parse_args called...\n') 

    out_flag = [];
    out_msg = [];
    passed_parms = [];
    
    %set option indices here 
    %we define the indices thus:
    %   1 <-> targ_prop
    %   2 <-> max_pct_attempts
    %   3 <-> auto_connects included or no
    %   4 <-> user-supplied A1
    %   5 <-> calculate weights
    %   6 <-> undefined
    OPT_TARG_PROP_IDX = 1;
    OPT_MAX_ATT_IDX = 2;
    OPT_AUTO_CONNECT_IDX = 3;
    OPT_USER_A1_IDX = 4;
    OPT_CALC_WEIGHTS_IDX = 5;
    OPT_UNDEF_IDX = 6;
    
        
    %below 'std' full-size strings
    passed_parms.multi_edge.option_str{OPT_TARG_PROP_IDX} = 'multi_edge_target_proportion';
    passed_parms.multi_edge.option_str{OPT_MAX_ATT_IDX} = 'max_targ_attempts';
    passed_parms.multi_edge.option_str{OPT_AUTO_CONNECT_IDX} = 'auto_connects_permitted';
    %note below string 'lower case' for parsing...
    passed_parms.matrix_data.user_A1_str{OPT_USER_A1_IDX} = 'user_a1';
    passed_parms.open_edge_data.option_return_weights_str = 'calculate_sampling_weights';
    
    %below 'abbreviated' strings
    passed_parms.multi_edge.option_str_abb{OPT_TARG_PROP_IDX} = 'targ_prop';
    passed_parms.multi_edge.option_str_abb{OPT_MAX_ATT_IDX} = 'max_att';
    passed_parms.multi_edge.option_str_abb{OPT_AUTO_CONNECT_IDX} = 'self_loops';
    %note below string 'lower case' for parsing...
    passed_parms.matrix_data.user_A1_str_abb{OPT_USER_A1_IDX} = 'a1';
    passed_parms.open_edge_data.option_return_weights_str_abb = 'weights';

    %number of acceptable strings 
    num_acceptable_strs = 2*length(passed_parms.multi_edge.option_str); 
    
    passed_parms.multi_edge.targ_prop_flag = 0;
    passed_parms.multi_edge.max_pct_attempts_flag = 0;
    passed_parms.auto_edge.user_flag = 0;
    passed_parms.matrix_data.user_A1_flag = 0;
    passed_parms.open_edge_data.calc_weights_flag = 0;
    passed_parms.unrecognised_flag = 0;
    
    %below array for tracking values passed
    %   1 <-> targ_prop
    %   2 <-> max_pct_attempts
    %   3 <-> auto_connects
    %   4 <-> user A1
    %   5 <-> calculate weights
    passed_parms.val{OPT_TARG_PROP_IDX} = -1; %init to '-1'; later set to default if not passed by user
    passed_parms.val{OPT_MAX_ATT_IDX} = -1; 
    passed_parms.val{OPT_AUTO_CONNECT_IDX} = -1;
    passed_parms.val{OPT_CALC_WEIGHTS_IDX} = -1;
    passed_parms.val{OPT_USER_A1_IDX} = -1;
    
    
    %assemble array of string identifiers for these options (internal refs)
    passed_parms.type_str{OPT_TARG_PROP_IDX} = 'targ_prop';
    passed_parms.type_str{OPT_MAX_ATT_IDX} = 'max_att';
    passed_parms.type_str{OPT_AUTO_CONNECT_IDX} = 'auto_edge';
    passed_parms.type_str{OPT_USER_A1_IDX} = 'user_A1';
    passed_parms.type_str{OPT_CALC_WEIGHTS_IDX} = 'weights';
    passed_parms.type_str{OPT_UNDEF_IDX} = 'undef';
    
    passed_parms.num_accepted_opts = length(passed_parms.type_str);
    
    %vector of flags for options '0' not set, '1' passed...
    passed_parms.opts_flag = zeros(1,passed_parms.num_accepted_opts);
    
    %unrecognised or incorrect passed options index
    passed_parms_unrecognised_idx = [];

    %fprintf('Total number of inputs = %d\n',nargin);
    
    
    %test number of var-args...
    nVarargs = length(varargin_local);
    %require pairs of arguments: string,value
    %if length is not even, reject
    if mod(nVarargs,2) ~= 0
        error('Require pairs of arguments: string,value...')
    end
    
    num_options = nVarargs/2;
    tempstr_idx = zeros(1,num_options);
    num_string = 0;
    %using tempvals as a cell array now to handle different data types instead of vector
    %tempvals = zeros(1,num_options);
    tempval_idx = zeros(1,num_options);
    cur_numeric_value = 0;
    num_valid = 0;
    num_unrecognised = 0;
    num_invalid = 0;
    temp_invalid = [];
    temp_invalid_idx = [];
    
    %vect for mapping which option passed from command line to our index
    passed_parm_map = zeros(1,num_options);
    
    %fprintf('Inputs in varargin_local(%d):\n',nVarargs)
    
    %test and flag string arguments against acceptable inputs (set above)
    for k = 1:nVarargs
      
      %if passed a string, test for proper input
      if ischar(varargin_local{k})
          
          %varargin_local{k}
      
          %build list of passed parms according to type
          switch lower(varargin_local{k})
              case {passed_parms.multi_edge.option_str{OPT_TARG_PROP_IDX},...
                      passed_parms.multi_edge.option_str_abb{OPT_TARG_PROP_IDX}}
                  %fprintf('\n\t"%s" passed...flagging\n',passed_parms.multi_edge.option_str{1})
                  passed_parms.multi_edge.targ_prop_flag = 1;
                  %passed_parms.multi_edge.targ_prop_idx = k;
                  num_string = num_string + 1;
                  tempstr_idx(num_string) = k;
                  passed_parms.varargin_idx(OPT_TARG_PROP_IDX) = k;
                  passed_parms.opts_flag(OPT_TARG_PROP_IDX) = 1;
                  num_valid = num_valid + 1;
                  passed_parm_map(num_valid) = OPT_TARG_PROP_IDX;
              case {passed_parms.multi_edge.option_str{OPT_MAX_ATT_IDX},...
                      passed_parms.multi_edge.option_str_abb{OPT_MAX_ATT_IDX}}
                  %fprintf('\n\t"%s" passed...flagging\n',passed_parms.multi_edge.option_str{2})
                  passed_parms.multi_edge.max_pct_attempts_flag = 1;
                  %passed_parms.multi_edge.max_pct_attempts_idx = k;
                  num_string = num_string + 1;
                  tempstr_idx(num_string) = k;
                  passed_parms.varargin_idx(OPT_MAX_ATT_IDX) = k;
                  passed_parms.opts_flag(OPT_MAX_ATT_IDX) = 1;
                  num_valid = num_valid + 1;
                  passed_parm_map(num_valid) = OPT_MAX_ATT_IDX;
              case {passed_parms.multi_edge.option_str{OPT_AUTO_CONNECT_IDX},...
                      passed_parms.multi_edge.option_str_abb{OPT_AUTO_CONNECT_IDX}}
                  %fprintf('\n\t"%s" passed...flagging\n',passed_parms.type_str{3})
                  passed_parms.auto_edge.user_flag = 1;
                  num_string = num_string + 1;
                  tempstr_idx(num_string) = k;
                  passed_parms.varargin_idx(OPT_AUTO_CONNECT_IDX) = k;
                  passed_parms.opts_flag(OPT_AUTO_CONNECT_IDX) = 1;
                  num_valid = num_valid + 1;
                  passed_parm_map(num_valid) = OPT_AUTO_CONNECT_IDX;     
              case {passed_parms.matrix_data.user_A1_str{OPT_USER_A1_IDX},...
                      passed_parms.matrix_data.user_A1_str_abb{OPT_USER_A1_IDX}}
                  %fprintf('\n\t"%s" passed...flagging\n',passed_parms.type_str{3})
                  passed_parms.matrix_data.user_A1_flag = 1;
                  num_string = num_string + 1;
                  tempstr_idx(num_string) = k;
                  passed_parms.varargin_idx(OPT_USER_A1_IDX) = k;
                  passed_parms.opts_flag(OPT_USER_A1_IDX) = 1;
                  num_valid = num_valid + 1;
                  passed_parm_map(num_valid) = OPT_USER_A1_IDX;                  
                  
              case {passed_parms.open_edge_data.option_return_weights_str,...
                      passed_parms.open_edge_data.option_return_weights_str_abb}
%                   fprintf('\n\t"%s" passed...flagging\n',...
%                       passed_parms.type_str{OPT_CALC_WEIGHTS_IDX})
                  passed_parms.open_edge_data.calc_weights_flag = 1;
                  num_string = num_string + 1;
                  tempstr_idx(num_string) = k;
                  passed_parms.varargin_idx(OPT_CALC_WEIGHTS_IDX) = k;
                  passed_parms.opts_flag(OPT_CALC_WEIGHTS_IDX) = 1;
                  num_valid = num_valid + 1;
                  passed_parm_map(num_valid) = OPT_CALC_WEIGHTS_IDX;                  
                  
                  
                  
              otherwise
                  %fprintf('\n\t"%s" unrecognised; rejecting...\n',varargin_local{k});
                  passed_parms.unrecognised_flag = 1;
                  num_unrecognised = num_unrecognised + 1;
                  passed_parms.unrecognised_str{num_unrecognised} = varargin_local{k};
                  %passed_parms.multi_edge.unrecognised_idx = k;
                  
                  passed_parms_unrecognised_idx(num_unrecognised) = k;

          end
      
      elseif isnumeric(varargin_local{k})
          %load up values into temporary array for later matching with
          %option string
          
          cur_numeric_value = cur_numeric_value + 1;
          tempvals{cur_numeric_value} = varargin_local{k};
          tempval_idx(cur_numeric_value) = k;
          
      else 
          %unacceptable inputs - not string or numeric...
          num_invalid = num_invalid + 1;
          temp_invalid{num_invalid} = varargin_local{k};
          temp_invalid_idx(num_invalid) = k;
          
      end
      
    end
    
    %if unrecognised string, print set of recognised/acceptable strings and exit
    if num_unrecognised>0
       fprintf('\n\tUnrecognised option "%s" passed\n\tacceptable options are:\n\n',passed_parms.unrecognised_str{1})
       for cur_opt = 1:0.5*num_acceptable_strs
           
           fprintf('\t\t"%s"\n',passed_parms.multi_edge.option_str{cur_opt});
           fprintf('\t\t\talternate: "%s"\n',passed_parms.multi_edge.option_str_abb{cur_opt});
           
       end
       fprintf('\n')
       error(' ');
        
    end
    
%     fprintf('\n\tNumber valid option-pairs passed: %d\n\tnumber unrecognised: %d\n\tnumber invalid: %d\n',...
%         num_valid,num_unrecognised,num_invalid)
%     
%     passed_parms.varargin_idx
%     
%     tempstr_idx
%     tempval_idx
%     temp_invalid_idx
    
    %quick error check: for each valid string argument passed, we require a
    %valid numerical argument; ensure we have suitable pairs:
    if nnz(tempstr_idx) ~= nnz(tempval_idx)
        error('String / numeric argument pair mismatch; %d string, %d numeric...',...
            nnz(tempstr_idx),nnz(tempval_idx));
    end
    %and ensure order of pairs: index for strings should be 1 less than
    %numeric paired value
    diff_idx = nnz((tempstr_idx+1) - tempval_idx);
    if diff_idx > 0
        error('String / numeric argument pairs must be sequence: "string1", "numeric1", ..., "stringn", "numericn",...');
    end
    
    %passed_parm_map
    
    %loop over number of pairs; for each flagged item load numeric value
    for cur_pair = 1:num_options
        
        %map over from passed option order of command line to our internal index
        cur_pair_internal = passed_parm_map(cur_pair);
        
        %key on whether indices for passed_parms are set, e.g., non-zero
        %values
        switch cur_pair_internal
            
            %we define the indices thus:
            %   1 <-> targ_prop
            %   2 <-> max_pct_attempts
            %   3 <-> auto_edges
            %   4 <-> user passed A1
            %   5 <-> calculate weights
            
            
            case 0
                %fprintf('Parameter for item# %d not set...\n',cur_pair)
        
            case 1  %
                 %fprintf('Parameter for item# %d "%s" set...\n',cur_pair_internal,passed_parms.type_str{cur_pair_internal})
                 
                 %test value passed: should be between 0, 1.0 for proportion of multi-edges
                 %fprintf('\tvalue passed: %2.3f...',tempvals(cur_pair_internal));
                 %fprintf('valid\n');
                 if (tempvals{cur_pair} >= 0) && (tempvals{cur_pair} <= 1.0)
                     passed_parms.val{1} = tempvals{cur_pair};
                     passed_parms.multi_edge.targ_prop = tempvals{cur_pair};
                     
                 else
                     error('Value %2.3f for option "%s" invalid; must be between 0 and 1.0',tempvals{cur_pair},passed_parms.type_str{cur_pair_internal})
                 end
                
            case 2
                 %fprintf('Parameter for item# %d "%s" set...\n',cur_pair_internal,passed_parms.type_str{cur_pair_internal})
                 
                 %test value passed: should be positive integer
                 %fprintf('\tvalue passed: %2.3f...',tempvals(cur_pair_internal));
                 %test for integer
                 diffval = floor(tempvals{cur_pair}) - tempvals{cur_pair};
                 if (tempvals{cur_pair} >= 1) && (diffval == 0)
                     %fprintf('valid\n');
                     passed_parms.val{2} = tempvals{cur_pair};
                     passed_parms.multi_edge.max_pct_attempts = tempvals{cur_pair};
                     
                 else
                     error('Value %2.3f for option "%s" invalid; must be positive integer',tempvals{cur_pair},passed_parms.type_str{cur_pair_internal})
                 end
                 
            case 3
                  %fprintf('Parameter for item# %d "%s" set...\n',cur_pair_internal,passed_parms.type_str{cur_pair_internal})
                  
                  %test value passed: should simply be binary '1' on/allowed, '0' off/disallowed
                  %fprintf('\tvalue passed: %2.3f...',tempvals{cur_pair});
                  if (tempvals{cur_pair} == 0) || (tempvals{cur_pair} == 1)
                     %fprintf('valid\n');
                     passed_parms.val{3} = tempvals{cur_pair};
                     passed_parms.auto_edge.permit_flag = tempvals{cur_pair};
                  else
                     error('Value %2.3f for option "%s" invalid; must be binary "0" (none) / "1" (permit)',...
                         tempvals{cur_pair},passed_parms.type_str{cur_pair_internal})
                  
                  end
                  
            case 4
                %fprintf('Parameter(s) for item # %d "%s" set...\n',cur_pair_internal,passed_parms.type_str{cur_pair_internal})
                
                %can only pre-test this item; don't know 'N'; but can test
                %for consistency with kin...
                if nnz(sum(tempvals{cur_pair},2) - kin) == 0
                    %fprintf('valid so far\n')
                    passed_parms.val{4} = tempvals{cur_pair};
                    passed_parms.matrix_data.user_A1_flag = 1;
                    passed_parms.matrix_data.A1 = tempvals{cur_pair};
                else
                    error('User-supplied "A1" permuted matrix must satisfy sum(A1,2) == kin')
                end
                
                %further check for negative values
                if any(any(tempvals{cur_pair} < 0))
                    error('User-supplied "A1" permuted matrix entries must be zero or positive integers')
                end
                
            case 5
%                 fprintf('Parameter for item # %d "%s" set...\n',...
%                     cur_pair_internal,passed_parms.type_str{cur_pair_internal})
                
                  %test value passed: should simply be binary '1' on/allowed, '0' off/disallowed
                  %fprintf('\tvalue passed: %2.3f...',tempvals{cur_pair});
                  if (tempvals{cur_pair} == 0) || (tempvals{cur_pair} == 1)
                     %fprintf('valid\n');
                     passed_parms.val{3} = tempvals{cur_pair};
                     passed_parms.open_edge_data.calc_weights_flag = tempvals{cur_pair};
                  else
                     error('Value %2.3f for option "%s" invalid; must be binary "0" (off) / "1" (calc & return)',...
                         tempvals{cur_pair},passed_parms.type_str{cur_pair_internal})
                  
                  end
                  
                  %check output variable containers: should be an additional output variable to
                  %accept the user-requested weights...
                  if nargout < 3
                      fprintf('\n\tUser requested weight calculations, but no output variable to accept data...\n')
                      fprintf('\t\tProvide additional output variable to continue\n\n');
                      error(' ')
                  end
                
                
            otherwise
                fprintf('I dont know what to do here...\n')
                
        end
        
    end

end

function [out_flag, out_msg, matrix_data, multi_edge_data, auto_edge_data, assby_data, open_edge_data, inert_swap_data, performance] = ...
    permA_init_structs(kin,kout,passed_parms, performance)
%internal function for parsing out input arguments into internal structures

    global debug_level N

    if debug_level >=1 
        fprintf('permA_init_structs called...\n') 
    end
    
    out_flag = [];
    out_msg = [];
    
    %
    %Initialise matrix_data fields
    %
    matrix_data.Nkin = length(kin);
    matrix_data.Nkout = length(kout);
    matrix_data.sum_kin = sum(kin);
    matrix_data.sum_kout = sum(kout);
    matrix_data.num_multi = 0;
    matrix_data.num_auto = 0;
    
    %quick test length kin & kout, set N if pass
    if matrix_data.Nkin == matrix_data.Nkout
        N = matrix_data.Nkin;
    else
        error('Length kin must equal length kout')
    end
    
    %set defaults for parms of interest
    
    DEFAULT_VERBOSITY = 1.0;            %level of output; 1.0 <-> normal; 0.0 <-> nothing...
    
    %A0 Assembly parameters
    DEFAULT_MULT_EDGE_PROP = 0.0;       %percentage multi-edges
    DEFAULT_MAX_TARG_ATTEMPTS = 4;      %A0 assby: attempts to meet above percentage via increasing lambda
    DEFAULT_MAX_PICK_TRIES = 1000;      %A0 assby: number picks from multi-edge distn to satisfy target mult-edge for node
    DEFAULT_DISTN_LAMBDA = 16;          %A0 assby: parm for exponential distn multi-edges
    DEFAULT_MULT_PCT_EPSILON = 0.01;    %A0 assby: acceptable deviation from target pct
    
    %Af assembly parameters
    DEFAULT_MAX_DONOR_LOOPS = 2*N;    %maximum number loops over all 'donor' nodes during Af assembly
    %DEFAULT_MAX_DONOR_LOOPS = 50;    %debug
    DEFAULT_INERT_SWAP_PICK_THRESH = 5; %threshold for switching pick type given ineffective inert swapping
    DEFAULT_INERT_SWAP_PICK_TYPE = 1;   %default to '1' <-> recipient, '2' <-> donor
    
    DEFAULT_INERT_SWAP_ACTIVATE_RATIO = 0.9;    %ratio of donor->recip swaps to total #donors when inert shuffling activated
                                                %set this higher (up to 1) for more frequent inert
                                                %shuffling
    
    DEFAULT_OPEN_EDGE_SEARCH_SWITCH_THRESH = 0.8;  %threshold of ratio (number row samples / max possible) search for switching methods
                                                   %set this lower (1 -> 0) for faster swap from
                                                   %exhaustive to 'quick' search method; note, the
                                                   %'quick' method is far less efficient when the
                                                   %ratio is low, but if this threshold is too
                                                   %high, and the 'quick' method doesn't activate
                                                   %the exhaustive method is substantially
                                                   %slower...your mileage may vary
    DEFAULT_OPEN_EDGE_SUCCESS_MEMORY = 10;       %history length for tracking open-edge search attempts (used for above switch)
    DEFAULT_OPEN_EDGE_SUCCESS_NUM_TEST = 3;     %number of history to test for failure
   
    %debug tracking: this flag determines whether to save _all_ iterations of Ai - very memory
    %intensive for large N!
    DEFAULT_DEBUG_TRACKING_FLAG = 0;
    
    %matrix debug tracking
    matrix_data.debug_tracking = DEFAULT_DEBUG_TRACKING_FLAG;
    
    %
    %Multi-edge initialisations
    %
    multi_edge_data.flag = 0;                                  %default: no multi-edges
    multi_edge_data.targ_prop = DEFAULT_MULT_EDGE_PROP;        %default: zero proportion multi-edges
    
    %test if user requested multi-edge proportion and proportion non-zero
    if passed_parms.multi_edge.targ_prop_flag &&...
        passed_parms.multi_edge.targ_prop > 0 
        multi_edge_data.flag = 1;     
        multi_edge_data.targ_prop = passed_parms.multi_edge.targ_prop;
    end

    %A0 assembly filds
    multi_edge_data.max_pick_tries = DEFAULT_MAX_PICK_TRIES;  %maximum attempts assembling multi-edge lists
    
    multi_edge_data.max_pct_attempts = DEFAULT_MAX_TARG_ATTEMPTS;  %default maximum attempts to meet target percentage
    if passed_parms.multi_edge.max_pct_attempts_flag %user passed parameter
        multi_edge_data.max_pct_attempts = passed_parms.multi_edge.max_pct_attempts;
    end
    
    multi_edge_data.distn_type = 1;                 %assume exponential distn
    multi_edge_data.lambda = DEFAULT_DISTN_LAMBDA;  %exponential distn parameter
    multi_edge_data.epsilon = DEFAULT_MULT_PCT_EPSILON;  %acceptable deviation actual pct multi from target
    
    %A1 flag - either user-supplied or no
    matrix_data.user_A1_flag = 0;
    if passed_parms.matrix_data.user_A1_flag
        matrix_data.user_A1_flag = 1;
        matrix_data.A1 = passed_parms.matrix_data.A1;
    end
    
    
    %
    %initialise Af assembly fields
    %
    
    %debug tracking
    assby_data.debug_tracking = DEFAULT_DEBUG_TRACKING_FLAG;
    
    %flags for success type: either A1 happens to be an Af, or went through full assembly
    %only really relevant for computing weights (see open_edge_data)
    %success types: 0 <-> no success, 1 <-> full assembly (A0 -> A1 -> Ai -> Af)
    %2 <-> A0 -> A1 == Af
    assby_data.success_type = -1;
    
    %col_sum_diff tracking
    assby_data.A_col_sum = zeros(1,N);
    assby_data.col_sum_diff = zeros(1,N);
    assby_data.col_sum_diff_sum = 0;
    
    assby_data.num_updates = 0;
    
    assby_data.num_donor_skips = 0;
    
    assby_data.tot_swaps = 0;
    
    assby_data.max_donor_loops = DEFAULT_MAX_DONOR_LOOPS;
    
    %target 'distance' tracking; e.g., deviation from target kout
    assby_data.norm_targ_dist = zeros(1,DEFAULT_MAX_DONOR_LOOPS);
    
    %donor loop data tracking
    assby_data.donor_loop_num_swaps = zeros(1,DEFAULT_MAX_DONOR_LOOPS);
    assby_data.donor_loop_norm_tracking = zeros(1,DEFAULT_MAX_DONOR_LOOPS);
    assby_data.donor_loop_inert_swaps = zeros(1,DEFAULT_MAX_DONOR_LOOPS);
    assby_data.skip_donor_tracking = zeros(1,DEFAULT_MAX_DONOR_LOOPS);
    
    %node type indexing
    assby_data.donor_idx = zeros(1,N);
    %assby_data.donor_test_idx = zeros(1,N);
    
    assby_data.recips_idx = zeros(1,N);  %this array with flags (0/1) for columns of recipient nodes
    
    assby_data.inert_idx = zeros(1,N);
    
    %node lists by type: stores current IDs each type out to number of type, vectors padded with
    %zeros after
    assby_data.donors = zeros(1,N);
    %assby_data.donors_test = zeros(1,N);
    
    assby_data.recips = zeros(1,N);     %this array for storing column index numbers of recipient nodes; padded zeros after length(recips)
    
    assby_data.inerts = zeros(1,N);
    
    assby_data.recips_trim = zeros(1,N);    %this array for storing col idx of trimmed recipient lists sans auto-edge nodes; padded zeros after length
    assby_data.recips_trim_idx = zeros(1,N);     %thhis array for flags (0/1) for columns of trimmed recipient list
    
%     %debug testing
%     assby_data.recips_test = zeros(1,N);
%     assby_data.recips_test_idx = zeros(1,N);
%     assby_data.recips_trim_test = zeros(1,N); 
%     assby_data.recips_trim_test_idx = zeros(1,N);
    
    %inert swapping parms & lists - inert swapping uses own local copy of donors & recips lists
%     assby_data.inert_swap_thresh = DEFAULT_INERT_SWAP_PICK_THRESH;
%     assby_data.inert_pick_type = DEFAULT_INERT_SWAP_PICK_TYPE;

    %debug tracking
    inert_swap_data.debug_tracking = DEFAULT_DEBUG_TRACKING_FLAG;
    
    inert_swap_data.inert_swap_thresh = DEFAULT_INERT_SWAP_PICK_THRESH;
    inert_swap_data.inert_pick_type = DEFAULT_INERT_SWAP_PICK_TYPE; %deprecated?
    inert_swap_data.activate_ratio = DEFAULT_INERT_SWAP_ACTIVATE_RATIO;
    inert_swap_data.inert_idx_topick = zeros(1,N);
    inert_swap_data.recip_idx_topick = zeros(1,N);
    inert_swap_data.donor_idx_topick = zeros(1,N);
    
    inert_swap_data.tot_num_swaps = 0;
    
    inert_swap_data.auto_edge_err_num = 0;
    
    %initialise auto-edge data fields
    auto_edge_data.flag = 0;            %auto-edges disallowed by default
    
    if passed_parms.auto_edge.user_flag   %user passed auto-edge option
        auto_edge_data.flag = passed_parms.auto_edge.permit_flag;
        
    end
    
    %debug tracking
    open_edge_data.debug_tracking = DEFAULT_DEBUG_TRACKING_FLAG;
    %for use with above
    open_edge_data.track_count = 0;
    
    open_edge_data.calc_weights = 0;
    %check if user requested weight calculation
    if passed_parms.open_edge_data.calc_weights_flag
        open_edge_data.calc_weights = 1;
        %fprintf('\n\tUser requested weight calculation...flagging\n')
        
        open_edge_data.weight_count = 0;
        
        %init the final weight value for return
        open_edge_data.weights = -1;
        
    end
    
    
    %for disallowed auto-edges, init auto-edge-fill matrix for their exclusion 
    open_edge_data.auto_edge_fill = eye(N);

    %open_edge flags, indexes, etc
    open_edge_data.search_mode = 1; %default exhaustive; '2' <-> sub-matrix search
    open_edge_data.donor_edge_rows = zeros(1,N);  %may init to zero-vector sized N...
    %open_edge_data.donor_edge_rows_idx = zeros(1,N);
    open_edge_data.num_donor_edge_rows = 0;
    open_edge_data.flag = 0;
    open_edge_data.current_row_select = 0;
    open_edge_data.row_select_idx = 0;
    
    open_edge_data.search_history_length = DEFAULT_OPEN_EDGE_SUCCESS_MEMORY;
    open_edge_data.search_history_test_length = DEFAULT_OPEN_EDGE_SUCCESS_NUM_TEST;
    open_edge_data.search_switch_thresh = DEFAULT_OPEN_EDGE_SEARCH_SWITCH_THRESH;
    
    %set below to zero since tracking ratio attempted/max searches
    open_edge_data.search_history = zeros(1,open_edge_data.search_history_length);
    
    %below for testing alternative trim method
    %open_edge_data.trim_test_warning = 0;
    
%     open_edge_data.donor_edge_rows_test = zeros(1,N);
%     open_edge_data.donor_edge_rows_test_idx = zeros(1,N);
%     open_edge_data.num_donor_edge_rows_test = 0;
    

    %performance tracking - 
    performance.open_search_std_num = 0;
    performance.open_search_alt_num = 0;
    
    

end

function [out_flag, out_msg] = permA_errcheck_Ai(kin,matrix_data,assby_data)
%internal function for testing current Ai build for errors: e.g., violating kin, mismatch sum of
%donors & recipients (sum(col_sum_A) must = 0) etc.

    global N debug_level
    
    out_flag = 0;
    out_msg = 'permA_errcheck_Ai';
    out_data = [];

    if debug_level >= 1
        fprintf('permA_errcheck_Ai called...\n');
    end
    
    %compare current in-degrees for Ai with prescribed kin
    matrix_data.Ai_kin_diff = kin - sum(matrix_data.Ai,2);
    matrix_data.Ai_kin_err = 0;
    if nnz(matrix_data.Ai_kin_diff) > 0
        matrix_data.Ai_kin_err = 1;
        out_flag = -1;
        out_msg = [out_msg 'kin error'];
    end
    
    %number of donors must == number of recipients:
    assby_data.Ai_col_sum_err = 0;
    assby_data.col_sum_diff_sum = sum(assby_data.col_sum_diff); 
    if assby_data.col_sum_diff_sum ~= 0
        assby_data.Ai_col_sum_err = 1;
        out_flag = -1;
        out_msg = [out_msg 'Ai col sum error'];
        
    end

end