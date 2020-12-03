%creates cell with objects needed for estimation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of draws after burn-in 
setup.number_of_draws=30000;%500000



%number of draws for choosing the scaling matrix in case of standard RW
%proposal
setup.scaling_draws=10000;%10000
%scaling is recomputed after check_scaling draws
setup.check_scaling=500;%200
%initial scaling for standard RW proposal
setup.initial_scaling=[5 .001 .001 .1 .0001 .001]';
%proposal=1 ->standard RW
%log likelihood computation
%likelihood=1 -> SS KF (only uses the SS covariance matrix, initial mean
%still has to be provided)
%likelihood=2 ->  KF

setup.likelihood=2;
%name of function that evaluates log prior
setup.prior_function='prior'; %right now prior function has to be called prior.m!
%initial value for the state and covariance of the state
setup.skip_opt=1; %skip optimization and go directly to MCMC (previously saved values are used then)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


setup.TVKF=0; %does the state space model feature time-varying matrices?


%should additional matrices be stored
setup.add_matrices=1;
%dimension of those variables
setup.dim_add_matrices=[5 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%name of dataset - see main_example.m for construction - different models
%can have different datasets
setup.data=['data.mat'];



%names of files that map parameters into state space matrices (i.e. solve models)
setup.wrapper{1}='benchmark'; 
setup.wrapper{2}='exp_utility'; 
setup.wrapper{3}='rot'; 
setup.wrapper{4}='rbc'; %large NK model
setup.wrapper{5}='liquid'; %large NK model
%=1 if initial values (means and covariance matrix) for Kalman filter are provided
setup.initial_provided=1;

%initial values




setup.state_initial{1}=[zeros(7,1);1]; 
setup.state_initial{2}=[zeros(8,1);1]; 

setup.state_initial{3}=[zeros(10,1);1];
setup.state_initial{4}=[zeros(13,1);1];
setup.state_initial{5}=[zeros(25,1);1];




for jj=1:5
    temp_cov=zeros(length(setup.state_initial{jj}),length(setup.state_initial{jj}));
temp_cov(1:length(setup.state_initial{jj})-1,1:length(setup.state_initial{jj})-1)=eye(length(setup.state_initial{jj})-1);
setup.cov_initial{jj}=temp_cov; 
end

%storage and display options

%display every disp iter draw
setup.disp_iter=500;
% keep every keep_draw th draw
setup.keep_draw=1;




%using uniform priors throughout
setup.index_normal=[1:42];
setup.normal_prior_means=[.7 .02 1.5 1 1 .01 .02 .02 1.5 1 1  .02 .02 .7 1 1 .02 1.5 1 1 .01 .02 .02 .3 1.5 1.3 .01 .02 2 1 2 .02 .02 .02  .5 1 1 1 1 1 1 1]'; % index for parameters that have normal priors
setup.normal_prior_std=[.1 .01 .3 2 2 .005 .01 .01 .3 2 2  .01 .01 .1 1 1 .01 .3 2 2 .005 .01 .01 .3 2 2 .005 .01 1 1 1 .01 .01 .01 .2  2 2 2 2 2 2 2]';
setup.index_gamma=[];
setup.initial_parameter=setup.normal_prior_means;
%parameter transformations for proposal (we want proposal to be unrestricted)
setup.index_logit_general=[2:31 36:42];
 setup.logit_general_lb=zeros(37,1);
 setup.logit_general_lb([2 8 17])=1; %a_bar>1 to get determinacy (note that a_bar here is the total coefficient on a_t-1 in the budget constraint
setup.logit_general_ub=10*ones(37,1);
setup.logit_general_ub([2 8 17])=100000;
setup.length_logit_general=37;


setup.index_log=[32 33 34];
setup.length_log=length(setup.index_log);

setup.index_logit=[35]; 
setup.length_logit=length(setup.index_logit);









%set parameter blocks


setup.number_blocks=7;
setup.index_block{1}=[1]';
setup.index_block{2}=[2:7]';
setup.index_block{3}=[8:16]';
setup.index_block{4}=[17:22]';
setup.index_block{5}=[23:31]';
setup.index_block{6}=[32:42]';
setup.index_block{7}=[43:47]';

%information on the models 
setup.number_models=5; %number of models
setup.index_model{1}=[1]';
setup.index_model{2}=[2:7];
setup.index_model{3}=[8:16];
setup.index_model{4}=[17:22];
setup.index_model{5}=[23:31]';
setup.index_model{6}=[32:42]';
setup.index_model{7}=[43:47]';

setup.weight_index=[43;44;45;46;47];
setup.dirichlet_scaling=50; %scaling for Dirichlet proposal

temp_prior_weights=40*[1/5;1/5;1/5;1/5;1/5];
setup.dirichlet_prior_parameters=temp_prior_weights;
%load starting_value max of posterior from previous run

load max_draw_fixed
setup.initial_parameter=max_draw_fixed;
load max_draw_model7;

setup.initial_parameter=[setup.initial_parameter;max_draw_model7(2:end);temp_prior_weights/sum(temp_prior_weights)];




setup.length_param_vector=length(setup.initial_parameter);


%legacy options - do not change
setup.proposal=1;
