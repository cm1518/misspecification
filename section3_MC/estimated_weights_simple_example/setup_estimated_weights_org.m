%creates cell with objects needed for estimation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of draws after burn-in 
setup.number_of_draws=1000;%500000
%adaptive MH not implemented in this code
setup.proposal=1;

%number of draws for choosing the scaling matrix in case of standard RW
%proposal
setup.scaling_draws=1000;%10000
%scaling is recomputed after check_scaling draws
setup.check_scaling=250;%200
%initial scaling for standard RW proposal
setup.initial_scaling=[5 .0057 .01 .0001 0.012]'; %note that right now there is no scaling parameter for the dirichlet params
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

%the next 4 values should be kept fixed - ultimately I will modify the code
%so that the filtered states for each model are saved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setup.state_size=0;
setup.sample_size=0;
%should additional matrices be stored
setup.add_matrices=1;
%dimension of those variables
setup.dim_add_matrices=[4 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%name of dataset - see main_example.m for construction - different models
%can have different datasets
setup.data=['data.mat'];



%names of files that map parameters into state space matrices (i.e. solve models)
setup.wrapper{1}='RRR'; %NK model with wages
setup.wrapper{2}='RRR_small'; %NK model without wages
setup.wrapper{3}='modelJPT'; %large NK model
setup.wrapper{4}='CK3'; %search & matching model
%=1 if initial values (means and covariance matrix) for Kalman filter are provided
setup.initial_provided=1;

%initial values
for jj=1:2
setup.state_initial{jj}=[zeros(16,1);1]; 
end

setup.state_initial{3}=[zeros(56,1);1]; 
setup.state_initial{4}=[zeros(31,1);1];
temp_cov=zeros(17,17);
temp_cov(1:16,1:16)=eye(16);
for jj=1:2
setup.cov_initial{jj}=temp_cov; 
end
temp_cov=zeros(57,57);
temp_cov(1:56,1:56)=eye(56);
setup.cov_initial{3}=temp_cov; 

temp_cov=zeros(32,32);
temp_cov(1:31,1:31)=eye(31);

setup.cov_initial{4}=temp_cov;

%storage and display options

%display every disp iter draw
setup.disp_iter=500;
% keep every keep_draw th draw
setup.keep_draw=1;





%priors - I use normal priors in this example. Right now the code only
%allows for normal and gamma priors, but than can easily be modified

setup.index_normal=[1:26 27:(27+33) 61:79]; % index for parameters that have normal priors
setup.normal_prior_means=[.2 ;1;0.5;1.5;0.125;.5;.5/100 ;.5;.5/100;.5/100;.5/100;.2;.2; .4;1;0.5;1.5;0.125;.5;.5/100;.5;.5/100;.5/100;.5/100;.2; .4;...
     .2 ;.5 ;.5 ;.5 ;.15; .15; 0 ;.5 ;.25; 2 ;.3 ;.66; 5 ;4; 1.7; 0.13; 0.13; 0.5; .5 ;.5; .5; .5; .5; .5; .5;.5 ;.5 ;.1 ;.5 ;.5; .5 ;.1 ;.1 ;.1;...
	 .2;11;0;1.5;2;0.5;0.2;0.03*3;0.85;1.5;0.5;0.9;0.2;0.7912;0.6712;0.3640/100;0.25/100;0.8716/100;0.6849/100]; %prior means
setup.normal_prior_std=zeros(size(setup.normal_prior_means));
 setup.normal_prior_std(1:60)=.25*setup.normal_prior_means(1:60); % prior std
setup.normal_prior_std(61:end)=.2*setup.normal_prior_means(61:end); % prior std
setup.normal_prior_std(find(setup.normal_prior_means==0))=.5;
setup.normal_prior_std(1)=.5;
setup.normal_prior_std([7 9 10 11 20 22 23 24])=.5/100;

setup.normal_prior_std(27+6)=.5;
setup.normal_prior_std(63)=.2;
setup.normal_prior_std(72:75)=.05;
setup.index_gamma=[];

%parameter transformations for proposal (we want proposal to be unrestricted)


setup.index_logit_general=[4 12 17 25 27 27+14];
 setup.logit_general_lb=[1 0 1 0 0 1]';
setup.logit_general_ub=[10 0.5 10 0.5 .5 10]';
setup.length_logit_general=length(setup.index_logit_general);


setup.index_log=[1;2;5;9;10;11;14;15;18;22;23;24;26; 10+25; 11+25; 12+25; 29+25; 30+25; 31+25; 32+25; 33+25; 34+25; 35+25;61;62;64;65;70;71;76;77;78;79];
setup.length_log=length(setup.index_log);

setup.index_logit=[3;6;7;8;13;16;19;20;21; 3+25; 5+25; 13+25; ((19:28)+25)';63;66;67;68;69;72;73;74;75]; 
setup.length_logit=length(setup.index_logit);









setup.number_blocks=6;
setup.index_block{1}=[1]';
setup.index_block{2}=[2:14]';
setup.index_block{3}=[15:26]';
setup.index_block{4}=[27:60]';
setup.index_block{5}=[61:79]';
setup.index_block{6}=[80:83]';
%information on the models
setup.number_models=4; %number of models
setup.index_model{1}=[1]';
setup.index_model{2}=[2:14];
setup.index_model{3}=[15:26];
setup.index_model{4}=[27:60];
setup.index_model{5}=[61:79];
setup.index_model{6}=[80:83];
%setup.model_weights=[.25 .25 .25 .25]; %model weights
setup.weight_index=[80;81;82;83];

setup.dirichlet_proposal_parameters=[1;1;1;1];
setup.dirichlet_prior_parameters=[2;2;4;2];
%load starting_value max of posterior from previous run
   load draw_max_fixed;
setup.initial_parameter=[draw_max_fixed;[.2;.2;.4;.2]];



setup.length_param_vector=length(setup.initial_parameter);
