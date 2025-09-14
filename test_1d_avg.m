% test on recover H, L from M0 = L0 + H*S0
% Make H: averaging filter

n = 300;
m = n - 1;
p = m;

first_row = [-1 1 zeros(1, m - 2)]; % First row of the circulant matrix H
H = toeplitz([first_row(1), fliplr(first_row(2:end))], first_row); 

rng(2);

% Make S0
s = 15;
S0 = zeros(p, n);
nonzero_idx = randsample(p*n, s*n);
S0(nonzero_idx) = rand(s*n, 1) + 1;

% Make L
rank_r = floor(min(m, n)*0.05);
L1 = randn([m, rank_r]);  % m x r matrix
L1 = orth(L1);
L2 = randn([n, rank_r]);  % r x n matrix
L2 = orth(L2);
% Multiply to get the low-rank matrix (m x n)
L0 = L1*L2';

% recover
M0 = L0 + H*S0;
clear para
para.rho_outer = 1; 
para.lasso_rho = 1; 
para.max_iter = 500;
para.lasso_max_iter = 20;
para.tol_outer = 1e-7;
para.lasso_tol = 1e-5;
para.lasso_method = 'ADMM';
para.lasso_decomp = "svd"; 
para.preconditioned = false;
lam =  1/sqrt(n);

tic
output = gen_matrix_sep(M0, H, lam, para);
t = toc;
if output.para.preconditioned
    sp = "with preconditioning";
else
    sp = "no preconditioning";
end
fprintf("******* %s, lasso by %s ******\n", sp, output.para.lasso_method)
fprintf("Number of iterations: %g.\n", output.count_outer);
fprintf("Duration: %.3f seconds.\n", t)
rel_L_c = norm(output.L - L0, 'fro')/norm(L0, 'fro');
rel_S_c = norm(output.S - S0, 'fro')/norm(S0, 'fro');
fprintf("Relative error of recovering L0: %e\n", rel_L_c);
fprintf("Relative error of recovering S0: %e\n\n", rel_S_c);

tic
outputC = gen_matrix_sep_con(M0, H, lam, para);
t = toc;
if outputC.para.preconditioned
    sp = "with preconditioning";
else
    sp = "no preconditioning";
end
fprintf("******* %s, lasso by %s ******\n", sp, outputC.para.lasso_method)
fprintf("Number of iterations: %g.\n", outputC.count_outer);
fprintf("Duration: %.3f seconds.\n", t)
rel_L_c = norm(outputC.L - L0, 'fro')/norm(L0, 'fro');
rel_S_c = norm(outputC.S - S0, 'fro')/norm(S0, 'fro');
fprintf("Relative error of recovering L0: %e\n", rel_L_c);
fprintf("Relative error of recovering S0: %e\n\n", rel_S_c);



para.lasso_method = 'FISTA';
tic
outputf = gen_matrix_sep(M0, H, lam, para);
t = toc;
if outputf.para.preconditioned
    sp = "with preconditioning";
else
    sp = "no preconditioning";
end
fprintf("******* %s, lasso by %s ******\n", sp, outputf.para.lasso_method)
fprintf("Number of iterations: %g.\n", outputf.count_outer);
fprintf("Duration: %.3f seconds.\n", t)
rel_L_c = norm(outputf.L - L0, 'fro')/norm(L0, 'fro');
rel_S_c = norm(outputf.S - S0, 'fro')/norm(S0, 'fro');
fprintf("Relative error of recovering L0: %e\n", rel_L_c);
fprintf("Relative error of recovering S0: %e\n\n", rel_S_c);

tic
outputfC = gen_matrix_sep_con(M0, H, lam, para);
t = toc;
if outputfC.para.preconditioned
    sp = "with preconditioning";
else
    sp = "no preconditioning";
end
fprintf("******* %s, lasso by %s ******\n", sp, outputfC.para.lasso_method)
fprintf("Number of iterations: %g.\n", outputfC.count_outer);
fprintf("Duration: %.3f seconds.\n", t)
rel_L_c = norm(outputfC.L - L0, 'fro')/norm(L0, 'fro');
rel_S_c = norm(outputfC.S - S0, 'fro')/norm(S0, 'fro');
fprintf("Relative error of recovering L0: %e\n", rel_L_c);
fprintf("Relative error of recovering S0: %e\n\n", rel_S_c);


