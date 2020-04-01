%
%Control script demonstrating generation of sample weights and estimation of spectral radius for
%space of potential adjacency matrices for sequences kin = [1 2 1] = kout
%as shown in Figure 5 (Distribution 3x3 Matrices) and Figure 7 (Comparison spectral radius)
%from "A Permutation Method for Network Assembly"
%
%2020 - Shawn Means
%School of Natural and Computational Sciences
%Massey University, New Zealand
%s.means@massey.ac.nz
%
%this .m file is part of a .zip archive providing access to functional Matlab implementation of the
%'Permutation Method' as described in "A Permutation Method for Network Assembly" (in submission,
%2020). If you have any problems with these items feel free to contact at email above.

%It is recommended that the variable 'verbosity_level' is set to '0' or off in the permuteA_beta3.m function
%for this suite of tests

clear;close all;

%set number of tests to perform - at least about 50,000 tests are needed to converge to expected
%distribution of final outcome Af's
num_tests = 10000;

%Given kin = [1 2 1] = kout, the permutation method assembles an A0 = [1 0 0;1 1 0;1 0 0] and
%permutations of this A0 -> A1 gives 27 possible variants; all these forms are in the following
%datafile
load ./permute_assby_3x3_A1_types_wselfloop.mat;

num_A1_types = length(A1_types);
A1_type_count = zeros(1,num_A1_types);

%There are only 5 possible final A forms (Af) satisfying the sequence kin = [1 2 1] = kout; we set
%these types here
Af_type{1} = [1 0 0;0 1 1;0 1 0];
Af_type{2} = [0 1 0;1 1 0;0 0 1];
Af_type{3} = [0 1 0;1 0 1;0 1 0];
Af_type{4} = [0 1 0;0 1 1;1 0 0];
Af_type{5} = [0 0 1;1 1 0;0 1 0];

num_A_types = length(Af_type);
A_type_count = zeros(1,num_A_types);

%Analysis of all possible paths from 27 A1 variants to the 5 Af forms, the following expected
%proportions of Af are calculated:
Af_type_expected = [0.2099	0.2006	0.1790	0.2099	0.2006];

%setting kin and kout and N
kin = [1 2 1];
kout = [1 2 1];

N = length(kin);

%set spectral radius values for the Af_types here - as calculated via max(abs(eig(Af_type{n})))
Af_spectral_r = [1.6180	1.6180	1.4142	1.4656	1.4656];

%calculate actual 'theoretical' or arithmetic mean of spectral radius for all Af_types
spectral_r_mean_t = mean(Af_spectral_r);

%allocate for A1 type, Af type, weights and weighted Q values tracking
A1_type_tracking = zeros(1,num_tests);
Af_type_tracking = zeros(1,num_tests);
weight_tracking = zeros(1,num_tests);
weight_tracking_returned = zeros(1,num_tests);
weight_mean_num = zeros(1,num_tests);
weight_mean_den = zeros(1,num_tests);
weight_mean = zeros(1,num_tests);

%printout frequency for progress header
header_frequency = num_tests/10;

processing_start = tic;

fprintf('\n\n\t\t%s\n\n',datestr(now));


    process_A1_start = tic;
    
    %loop over A assembly saving each realisation
    num_success_tests = 0;
    num_failure_tests = 0;
    
    %track successful A1 == Af runs
    num_A1_is_Af = 0;
    
    %preallocate success run and A_type tracking
    success_test_ID = zeros(1,num_tests);

    A_type_ID = zeros(1,num_tests);
    
    %and product of open edges for our weights
    open_edge_product = zeros(1,num_tests);

    %open_edge_product2 = zeros(1,num_tests);
    

    
    
    for cur_test = 1:num_tests
        
        
        if mod(cur_test,header_frequency) == 0
            fprintf('\n\tPerforming test %d/%d\n',cur_test,num_tests)
        end

        
        [testA{cur_test},success{cur_test},weights{cur_test},testdata{cur_test}] = ...
            permuteA_beta3(kin',kout','targ_prop',0,'self_loops',1,'weights',1);

        
        %if a successful A returned;
        if success{cur_test}
        
            num_success_tests = num_success_tests + 1;

            %save this run ID 
            success_test_IDs(num_success_tests) = cur_test;
            
            %search for matching Af_type from this current returned A
            Af_matched = 0;
            cur_type_test = 1;
            while ~Af_matched && (cur_type_test <= num_A_types)

               if nnz(testA{cur_test} - Af_type{cur_type_test}) == 0
                  %fprintf('\n\t\tMatch found with canonical A #%d...\n',cur_type_test)
                  A_type_count(cur_type_test) = A_type_count(cur_type_test) + 1;
                  testdata{cur_test}.A_type = cur_type_test;

                  %A_type_ID(num_success_tests) = cur_type_test;
                  Af_type_tracking(num_success_tests) = cur_type_test;

                  Af_matched = 1;

               else
                  %not a match bump cur_type_test continue...
                  cur_type_test = cur_type_test + 1;
               end
            end  
            
            %same thing but search for A1_type
            A1_matched = 0;
            cur_type_test = 1;
            while ~A1_matched && (cur_type_test <= num_A1_types)

               if nnz(testdata{cur_test}.matrix_data.A1 - A1_types{cur_type_test}) == 0
                  %fprintf('\n\t\tMatch found with canonical A #%d...\n',cur_type_test)
                  A1_type_count(cur_type_test) = A1_type_count(cur_type_test) + 1;
                  %testdata{cur_test}.A_type = cur_type_test;

                  %A_type_ID(num_success_tests) = cur_type_test;
                  A1_type_tracking(num_success_tests) = cur_type_test;

                  A1_matched = 1;

               else
                  %not a match bump cur_type_test continue...
                  cur_type_test = cur_type_test + 1;
               end
            end
            
            %if the provided A1 _is_ an Af_type, then permuteA does no work
            %and below data not returned; assby_data should return in this
            %case with success_type == 1
            if testdata{cur_test}.assby_data.success_type == 1
            
                %save built-in weight calc
                weight_tracking_returned(num_success_tests) = testdata{cur_test}.open_edge_data.weights;
                        
            %otherwise, success type without such luck    
            elseif testdata{cur_test}.assby_data.success_type == 2

                %save built-in weight calc
                weight_tracking_returned(num_success_tests) = testdata{cur_test}.open_edge_data.weights;

                %track this run - which A1->Af and number of
                num_A1_is_Af = num_A1_is_Af + 1;
                A1_is_Af_tracking(1,num_A1_is_Af) = cur_type_test;
                
            end
            
                %we keep running sums over numerator and denominator since weighted mean =
                %sum(...)/sum(...)
                if num_success_tests == 1
                    weight_mean_num(num_success_tests) = ...
                        weight_tracking_returned(num_success_tests)*Af_spectral_r(Af_type_tracking(num_success_tests));

                    weight_mean_den(num_success_tests) = weight_tracking_returned(num_success_tests);

                else

                    weight_mean_num(num_success_tests) = weight_mean_num(num_success_tests-1) + ...
                        weight_tracking_returned(num_success_tests)*Af_spectral_r(Af_type_tracking(num_success_tests));

                    weight_mean_den(num_success_tests) = weight_mean_den(num_success_tests-1) + ...
                        weight_tracking_returned(num_success_tests);

                end
                    %compute mean given above running sums
                    weight_mean(num_success_tests) = weight_mean_num(num_success_tests)/weight_mean_den(num_success_tests);
                
            
                
        else
            num_failure_tests = num_failure_tests + 1;
            
            
        end
        
    end
    
    fprintf('num successes: %d\n',num_success_tests)

%compute actual proportions of Af_types for this run
Af_type_actual = A_type_count./num_tests;

%pull out strings for labels from these actual vals
xtick_prop_str{1} = num2str(Af_type_actual(1));
xtick_prop_str{2} = num2str(Af_type_actual(2));
xtick_prop_str{3} = num2str(Af_type_actual(3));
xtick_prop_str{4} = num2str(Af_type_actual(4));
xtick_prop_str{5} = num2str(Af_type_actual(5));
   
    
process_A1_finish = toc(process_A1_start);

%plot progression of mean spectral radius - arithmetic and sampled
figure
set(gcf,'position',[298   536   560   420]);

subplot(2,1,1)

semilogx([1:num_success_tests],weight_mean,'b-','linewidth',2);
hold on; grid on
semilogx([1:num_success_tests],spectral_r_mean_t*ones(1,num_success_tests),'r-','linewidth',2)

legend('Weighted','Arithmetic','Location','best')

%xlabel('Number of Trials');
ylabel('Mean Spectral Radius');

title('Mean Spectral Radius Comparison')

set(gca,'fontsize',14)
set(gca,'fontweight','bold')

subplot(2,1,2)

%plot relative error from arithmetic mean
loglog([1:num_success_tests],abs(weight_mean - spectral_r_mean_t)./spectral_r_mean_t,'b-','linewidth',2);
hold on;grid on

xlabel('Number of Trials');
ylabel('Relative Error');

set(gca,'fontsize',14)
set(gca,'fontweight','bold')

%plot histogram of distributed Af results against expected
figure
set(gcf,'position',[871   536   560   420]);

bar(A_type_count);

titlestr = sprintf('%d tests (%d/%d fail/success)',num_tests,num_failure_tests,num_success_tests);

title(titlestr,'fontsize',14,'fontweight','bold');

%compute expected # each type given num success - for perfect uniform distn
expected_num = num_success_tests/num_A_types;

%plot comparison line
cur_xlim = get(gca,'xlim');
hold on;grid on;
xvect = [cur_xlim(1):cur_xlim(end)];
plot(xvect,expected_num*ones(1,length(xvect)),'r-');

%append xlabel tick marks with actual proportions
set(gca,'xticklabel',xtick_prop_str)

%plot expected values for reference
plot([1:5],num_tests*Af_type_expected,'r+')
