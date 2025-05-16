% test on recover H, L from M0 = L0 + H*S0
% no preconditioning is applied
% Make H: random filter

n = 50;
m = 45;
p = 48;

rng(30);
H = randn(m, p);

% Make S0
s = 2;
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
para.rho_outer = 0.5; 
para.rho_inner = 1; 
para.N_outer = 500;
para.N_inner = 30;
para.tol_outer = 1e-7;
para.tol_inner = 1e-5;
para.decomp = "svd"; 
lam =  1/sqrt(m);
tic
output = gen_matrix_sep(M0, H, lam, para);
toc

rel_L = norm(output.L - L0, 'fro')/norm(L0, 'fro');
rel_S = norm(output.S - S0, 'fro')/norm(S0, 'fro');
fprintf("Relative error of recovering L0: %e\n", rel_L);
fprintf("Relative error of recovering S0: %e\n", rel_S);
fprintf("Ran %g many outer loops.\n", output.count_outer);
fprintf("Input H was circulant? Answer: %g \n", output.isCirc)
fprintf("****** with preconditioning ******\n")
[CH, C] = precondition(H);
tic
outputC = gen_matrix_sep(C*M0, CH, lam, para);
toc
Lhat = M0 - H*outputC.S;
rel_L = norm(Lhat - L0, 'fro')/norm(L0, 'fro');
rel_S = norm(outputC.S - S0, 'fro')/norm(S0, 'fro');
fprintf("Relative error of recovering L0: %e\n", rel_L);
fprintf("Relative error of recovering S0: %e\n", rel_S);
fprintf("Ran %g many outer loops.\n", outputC.count_outer);
fprintf("Input H was circulant? Answer: %g \n", outputC.isCirc)
