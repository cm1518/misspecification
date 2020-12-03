
test_vec=[0:.1:5];

for jj=1:length(test_vec)

[A B C D add_matrices]=wrapper_func([params(setup.index_model{1});test_vec(jj)],setup,data{mm});
[ llk(jj) xest] = KF(A,B,C,D,setup.state_initial{mm},setup.cov_initial{mm},data{mm});

end