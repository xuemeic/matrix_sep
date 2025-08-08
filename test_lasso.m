% test lasso
m = 1000;

p = m*2;
k = 2;

first_row = [-1 1 zeros(1, m - 2)]; % First row of the circulant matrix H
%A = toeplitz([first_row(1), fliplr(first_row(2:end))], first_row); 

A = randn(p, m);
% Make x0
s = 15;
x0 = zeros(m, k);
nonzero_idx = randsample(m*k, s*k);
x0(nonzero_idx) = rand(s*k, 1) + 1;

b = A*x0;
lam = 1e-2;

clear para
para.max_iter = 1000;
para.tol = 1e-7;

para.method = 'ADMM';
para.rho = 1;
fprintf("*********** ADMM ***********\n")
tic
[x, numi] = my_lasso(A, b, lam, para);
toc

fprintf("Ran %g many outer loops.\n", numi);
fprintf("relative error for ADMM is %e.\n", norm(x - x0)/norm(x0));



para.method = 'ISTA';
fprintf("*********** ISTA ***********\n")
tic
[xi, num2] = my_lasso(A, b, lam, para);
toc

fprintf("Ran %g many outer loops.\n", num2);
fprintf("relative error for ISTA is %e.\n", norm(xi - x0)/norm(x0));


para.method = 'FISTA';
fprintf("*********** FISTA ***********\n")
tic
[xf, nf] = my_lasso(A, b, lam, para);
toc
fprintf("Ran %g many outer loops.\n", nf);
fprintf("relative error for FISTA is %e.\n\n", norm(xf - x0)/norm(x0));