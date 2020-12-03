function [ covariance, mean ] = recursive_moments(index, data,covariance, mean,setup )
%recursive calulation of mean and covariance matrix including correction
meanold=mean;
mean=1/index*data+(index-1)/(index)*meanold;
covariance=(index-2)/(index-1)*covariance+setup.scaling_adaptive/(index-1)*((index-1)*(meanold*meanold')-(index)*(mean*mean')+data*data'+setup.eps_adaptive*eye(setup.length_param_vector));
end

