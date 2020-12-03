function [ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH( setup )
%main function for the MH algorithm
data=load(setup.data);
data=data.data;





warning off;
options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-4,'MaxIter',500, 'MaxFunEvals',50000);

if setup.proposal==1
   temp=setup.scaling_draws/setup.check_scaling;
 
   
   
      
       post=@(params,setup,data) (posterior( params,setup,data ));
       postopt=@(params,setup,data) -(posterior( params,setup,data ));
       postopt2=@(params) -(posterior( params,setup,data ));
      if setup.skip_opt==0
   %initial minimization
   [ x0 ] = transform( setup.initial_parameter,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub  );
 xestimate1=x0;
   
      options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-4,'MaxIter',5000, 'MaxFunEvals',50000);
 
 %setting up constraints for fmincon (constraints on model weights)
 A=-eye(setup.number_models);
 A=[zeros(setup.number_models,setup.length_param_vector-setup.number_models) A];
 B=zeros(setup.number_models,1);
 Aeq=zeros(1,setup.length_param_vector);
 Aeq(end-setup.number_models+1:end)=1;
 Beq=1;
 
[xestimate1,functionvalue1,exitflag,output,lambda,grad,hessian]=fmincon(postopt2,xestimate1,A,B,Aeq,Beq,[],[],[],options);
 
       xestimate=xestimate1;
       functionvalue=functionvalue1;
 save temp_max xestimate1
    elseif setup.skip_opt==1
[ x0 ] = transform( setup.initial_parameter,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub  );
 xestimate=x0;
functionvalue=postopt2(xestimate);
end
  
  

 for bb=1:setup.number_blocks-1  
     try
        
        Hess_temp = hessian3(func2str(@postblock),xestimate(setup.index_block{bb}),setup,data,bb,xestimate);

   Hess_temp=reshape(-Hess_temp,length(setup.index_block{bb}),length(setup.index_block{bb}));
  Hess_temp=eye(size(Hess_temp,1))/Hess_temp; 
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %making sure the Hessian is positive definite
  [V,D]=eig(Hess_temp);

   d=diag(D);
   d(d<=0)=eps;

   Hess_temp= V*diag(d)*V';
   Hess{bb}=.75*Hess_temp+.25*max(diag(Hess_temp))*eye(length(setup.index_block{bb}));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     catch ME
         disp('Hessian computation in the following block did not work:')
         bb
Hess{bb}=.5*eye(length(setup.index_block{bb}));
     end
     
end  %end

 scaling=ones(setup.number_blocks-1,1).*setup.initial_scaling; %because the scaling only refers to the RW blocks

%set values for Hess and scaling for Dirichlet block so that code won't
%return an error
Hess{setup.number_blocks}=[];
scaling(setup.number_blocks)=NaN;

  
   old_posterior=-functionvalue;
   posterior_draw=xestimate;
   acc_rate=zeros(setup.number_blocks,1);
   for jj=1:temp %loop to adjust scaling matrix
       disp('iteration for adjustment of proposal scaling factor:')
       jj
    for bb=1:setup.number_blocks
       [ acc_rate(bb) ,posterior_draw,old_posterior] = acceptance_rate( posterior_draw,Hess{bb},scaling(bb),setup.check_scaling,old_posterior,setup,data,bb  );  
    
       disp('block')
       bb
       disp('acceptance rate during that iteration:')
       acc_rate(bb) 
       
     if acc_rate(bb)<.2 && bb<setup.number_blocks 
       scaling(bb)=.75*scaling(bb);
       elseif acc_rate(bb)>.5 && bb<setup.number_blocks 
           scaling(bb)=1.5*scaling(bb);
     end
	 
	 if acc_rate(bb)<.02 && bb==setup.number_blocks 
       setup.dirichlet_scaling=2*setup.dirichlet_scaling;
       elseif acc_rate(bb)>.95 && bb==setup.number_blocks 
           setup.dirichlet_scaling=.5*setup.dirichlet_scaling;
     end
	 
	 
    end 
    end
    disp('scaling:')
    scaling
    disp('acceptance rate:')
    acc_rate
    
    [ acc_rate ,draws,log_posteriors, statedraws, add_matrices] = draws_standard_RW( posterior_draw,Hess,scaling,old_posterior,setup,data  );
    
 else %
     
     disp('this code is only implemented for the standard MH algorithm for now')

end



end

