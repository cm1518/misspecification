%creates cell with objects needed for estimation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of draws after burn-in 
setup.number_of_draws=50000;%500000
%adaptive MH not implemented in this code
setup.proposal=1;

%number of draws for choosing the scaling matrix in case of standard RW
%proposal
setup.scaling_draws=30000;%10000
%scaling is recomputed after check_scaling draws
setup.check_scaling=500;%200
%initial scaling for standard RW proposal
setup.initial_scaling=[15 5 .1]'; %note that right now there is no scaling parameter for the dirichlet params
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


%name of dataset - see main_example.m for construction - different models
%can have different datasets
setup.data=['data.mat'];



%names of files that map parameters into state space matrices (i.e. solve models)
setup.wrapper{1}='HS'; %NK model without wages
setup.wrapper{2}='modelJPT'; %large NK model
setup.initial_provided=1;

%initial values

setup.state_initial{1}=[zeros(8,1);1]; 


setup.state_initial{2}=[zeros(56,1);1]; 
temp_cov=zeros(9,9);
temp_cov(1:8,1:8)=eye(8);

setup.cov_initial{1}=temp_cov; 

temp_cov=zeros(57,57);
temp_cov(1:56,1:56)=eye(56);
setup.cov_initial{2}=temp_cov; 

%storage and display options

%display every disp iter draw
setup.disp_iter=500;
% keep every keep_draw th draw
setup.keep_draw=1;





load initial_values_and_priors_HS
setup.initial_parameter=starting_value;


%priors - I use normal priors in this example. Right now the code only
%allows for normal and gamma priors, but than can easily be modified
setup.index_normal=1:length(starting_value); % index for parameters that have normal priors
setup.normal_prior_means=normal_means; %prior means
setup.normal_prior_std=normal_std;
 
setup.index_gamma=[];

%parameter transformations for proposal (we want proposal to be unrestricted)







setup.index_logit_general=logit_general_total;
 setup.logit_general_lb=lb;
setup.logit_general_ub=ub;
setup.length_logit_general=length(setup.index_logit_general);


setup.index_log=log_total;
setup.length_log=length(setup.index_log);

setup.index_logit=logit_total; 
setup.length_logit=length(setup.index_logit);











setup.number_blocks=4;
setup.index_block{1}=[1:2]';
setup.index_block{2}=[3:13]';
setup.index_block{3}=[14:46]';
setup.index_block{4}=[47:48]';
%information on the models
setup.number_models=2; %number of models
setup.index_model{1}=[1:2]';
setup.index_model{2}=[3:13];
setup.index_model{3}=[14:46]';
setup.index_model{4}=[47:48]';
%setup.model_weights=[.25 .25 .25 .25]; %model weights
setup.weight_index=[47;48];
setup.dirichlet_scaling=800;

temp_prior_weights=800*[1/2;1/2];
setup.dirichlet_prior_parameters=temp_prior_weights;
%load starting_value max of posterior from previous run
setup.initial_parameter=[setup.initial_parameter;temp_prior_weights/sum(temp_prior_weights)];




setup.length_param_vector=length(setup.initial_parameter);

%the next 2 values should be kept fixed 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%should additional matrices be stored
setup.add_matrices=1;
%dimension of those variables
setup.dim_add_matrices=[2 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


