%creates cell with objects needed for estimation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of draws after burn-in 
setup.number_of_draws=5000;%500000

%number of draws for choosing the scaling matrix in case of standard RW proposal
setup.scaling_draws=1000;%10000
%scaling is recomputed after check_scaling draws
setup.check_scaling=100;%200
%initial scaling for standard RW proposal
setup.initial_scaling=[1]';
%proposal=1 ->standard RW
%adaptive MH not implemented in this code
setup.proposal=1;
%log likelihood computation
%likelihood=1 -> SS KF (only uses the SS covariance matrix, initial mean
%still has to be provided)
%likelihood=2 ->  KF

setup.likelihood=2;
%name of function that evaluates log prior
setup.prior_function='prior'; %right now prior function has to be called prior.m!
%initial value for the state and covariance of the state
setup.skip_opt=0; %skip optimization and go directly to MCMC (previously saved values are used then)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


setup.TVKF=0; %does the state space model feature time-varying matrices?

%the next 4 values should be kept fixed - ultimately I will modify the code
%so that the filtered states for each model are saved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setup.state_size=0;
setup.sample_size=0;
%should additional matrices be stored
setup.add_matrices=1;
%dimension of those variables
setup.dim_add_matrices=[2 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%name of dataset - see main_example.m for construction - different models
%can have different datasets
setup.data=['data.mat'];



%names of files that map parameters into state space matrices (i.e. solve models)
setup.wrapper{1}='AR_1';
setup.wrapper{2}='MA_1';
%setup.wrapper{3}='ARMA_11';

%=1 if initial values (means and covariance matrix) for Kalman filter are provided
setup.initial_provided=1;

%initial values
setup.state_initial{1}=0; %AR(1)
setup.state_initial{2}=[0;0]; %MA(1)
%setup.state_initial{3}=[0;0;0]; %ARMA(1)

setup.cov_initial{1}=1; %AR(1)
setup.cov_initial{2}=eye(2); %MA(1)
%setup.cov_initial{2}=eye(3); %ARMA(1)


%storage and display options

%display every disp iter draw
setup.disp_iter=5000;
% keep every keep_draw th draw
setup.keep_draw=1;





%priors - I use normal priors in this example. Right now the code only
%allows for normal and gamma priors, but than can easily be modified

setup.index_normal=[1;2;3]; % index for parameters that have normal priors
setup.normal_prior_means=[0;0;0]; %prior means
setup.normal_prior_std=[1;.5;1]; % prior std


setup.index_gamma=[];

%setup.index_uniform=[1;2;3]; % index for parameters that have normal priors
%setup.uniform_prior_means=[0;0;0]; %prior means
%setup.uniform_prior_std=[1;.5;1]; % prior std

%parameter transformations for proposal (we want proposal to be unrestricted)


setup.length_logit_general=0;
setup.index_logit_general=[];
 setup.logit_general_lb=[]';
setup.logit_general_ub=[]';


setup.length_log=2;
setup.index_log=[1 3];

setup.length_logit=1;
setup.index_logit=[2];

%the code can handle multiple block MH, but right now that would be very
%inefficient since each model and the associated likelihoods are recomputed
%in every block. So stick to 1 block for now.

setup.number_blocks=1;
setup.index_block{1}=[1:3]';

%information on the models
setup.number_models=2; %number of models
setup.index_model{1}=[1]; %parameter indices for different models - there is one more set of indices than there are 
%models - the first set of indices refers to the joint parameters
setup.index_model{2}=[2];
setup.index_model{3}=[3];
setup.model_weights=[.5 .5]; %model weights

setup.initial_parameter=[2;.5;.5];
setup.length_param_vector=length(setup.initial_parameter);
