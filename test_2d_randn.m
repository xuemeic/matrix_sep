m1 = 10;
p1 = 9;
m2 = 9;
p2 = 8;
K = 50;

rng(3);
G1 = randn(m1, p1);
G2 = randn(m2, p2);

% make S0
num_elements = p1 * p2;
percent_nonzero = 4;
num_nonzeros = K*round((percent_nonzero / 100) * num_elements);
nonzero_idx = randsample(p1*p2*K, num_nonzeros);
S0 = zeros(p1, p2, K);
S0(nonzero_idx) = rand(num_nonzeros, 1) + 1;

% make L0
rank_r = 1;
L1 = randn([m1*m2, rank_r]);  % m x r matrix
L1 = orth(L1);
L2 = randn([K, rank_r]);  % r x n matrix
L2 = orth(L2);
% Multiply to get the low-rank matrix (m x n)
L0 = L1*L2';
L0 = reshape(L0, [m1, m2, K]);

% make M0
M0 = L0 + pagemtimes(pagemtimes(G1, S0), G2');

% run
clear para
para.rho_outer = 1; 
para.rho_inner = 1; 
para.max_iter = 200;
para.lasso_max_iter = 50;
para.tol_outer = 1e-8;
para.lasso_tol = 1e-6;
para.is_E = false;
para.lasso_method = 'ADMM';
para.preconditioned = false;
lam =  1/sqrt(max(m1*m2, K));


tic
output = gen_matrix_sep_2d_con(M0, G1, G2, lam, para);
t = toc;

desired_print(output, t, L0, S0);




para.lasso_method = 'FISTA';
tic
outputC = gen_matrix_sep_2d_con(M0, G1, G2, lam, para);
t2 = toc;
desired_print(outputC, t2, L0, S0);

function desired_print(o, t, L0, S0)
if o.para.preconditioned
    sp = "with preconditioning";
else
    sp = "no preconditioning";
end
fprintf("******* Gi Gaussian, %s, lasso by %s ******\n", sp, o.para.lasso_method)
fprintf("Number of iterations: %g.\n", o.count_outer);
fprintf("Duration: %.3f seconds.\n", t)
rel_L = norm(o.L - L0, 'fro')/norm(L0, 'fro');
rel_S = norm(o.S - S0, 'fro')/norm(S0, 'fro');
fprintf("Relative error of recovering L0: %e\n", rel_L);
fprintf("Relative error of recovering S0: %e\n\n", rel_S);
end
