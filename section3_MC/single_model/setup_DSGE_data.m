%creates cell with objects needed for estimation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of draws after burn-in 
setup.number_of_draws=10000;%500000

%number of draws for choosing the scaling matrix in case of standard RW
%proposal
setup.scaling_draws=1000;%10000
%scaling is recomputed after check_scaling draws
setup.check_scaling=500;%200
%initial scaling for standard RW proposal
setup.initial_scaling=[2.4327e-04]';
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


setup.state_size=57;
setup.sample_size=230;
%should additional matrices be stored
setup.add_matrices=0;
%dimension of those variables
setup.dim_add_matrices=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%name of dataset - see main_example.m for construction - different models
%can have different datasets
setup.data=['data.mat'];



%names of files that map parameters into state space matrices (i.e. solve models)
setup.wrapper='modelJPT'; 



%=1 if initial values (means and covariance matrix) for Kalman filter are provided
setup.initial_provided=1;

%initial values

setup.state_initial=[zeros(56,1);1]; %nk1 - 1 is for the intercept


temp_cov=zeros(57,57);
temp_cov(1:56,1:56)=eye(56);

setup.cov_initial=temp_cov; %nk1



%storage and display options

%display every disp iter draw
setup.disp_iter=500;
% keep every keep_draw th draw
setup.keep_draw=1;





%priors - I use normal priors in this example. Right now the code only
%allows for normal and gamma priors, but than can easily be modified

setup.index_normal=[1:11 12 13:35]; % index for parameters that have normal priors



setup.normal_prior_means=[.3 .2 .5 .5 .5 .15 .15 0 .5 .25 2 .2 .66 5 4 1.7 0.13 0.13 0.5 .5 .5 .5 .5 .5 .5 .5...
 .5 .5 .1 .5 .5 .5 .1 .1 .1]'; %prior means
setup.normal_prior_std=.25*setup.normal_prior_means; % prior std
setup.normal_prior_std(8)=.5;
setup.normal_prior_std(12)=.5;
setup.index_gamma=[];



%parameter transformations for proposal (we want proposal to be unrestricted)


setup.length_logit_general=2;
setup.index_logit_general=[2 16];
 setup.logit_general_lb=[0 1]';
setup.logit_general_ub=[.5 10]';


setup.length_log=11;
setup.index_log=[1 10 11 12 17 29 30 31 32 33 34 35];


setup.index_logit=[3 5 13 19:28];
setup.length_logit=length(setup.index_logit);


setup.number_blocks=1;
setup.index_block{1}=[1:35]';


load draw_max;

%setup.initial_parameter=[.3 .2 .5 .5 .5 .15 .15 .4 .5 .25 2 .4  .66 5 4 1.7 0.13 0.13 0.6 .6 .6 .6 .6 .6 .6 .6...
% .5 .5 .1 .5 .5 2 .1 .1 .1]';
setup.initial_parameter=draw_max;
setup.initial_parameter(17)=0.01;
%setup.initial_parameter(21)=.8;
setup.length_param_vector=length(setup.initial_parameter);
