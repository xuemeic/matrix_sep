% This is when both Gi are circulant
m1 = 20;
p1 = m1;
m2 = 18;
p2 = m2;
K = 50;

rng(3);
G1r1 = randn(1, m1);  % First row of the circulant matrix 
G1 = toeplitz([G1r1(1), fliplr(G1r1(2:end))], G1r1); 

G2r1 = rand(1, m2);  % First row of the circulant matrix 
G2 = toeplitz([G2r1(1), fliplr(G2r1(2:end))], G2r1); 

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
para.rho_outer = 1; % smaller rho results fewer iterations
para.rho_inner = 1; % 10 runs longer
para.N_outer = 200;
para.N_inner = 50;
para.tol_outer = 1e-8;
para.tol_inner = 1e-6;
para.is_E = false;
lam =  1/sqrt(min(m1*m2, K));


tic
output = gen_matrix_sep_2d(M0, G1, G2, lam, para);
toc

rel_L = norm(output.L - L0, 'fro')/norm(L0, 'fro');
rel_S = norm(output.S - S0, 'fro')/norm(S0, 'fro');
fprintf("Relative error of recovering L0: %e\n", rel_L);
fprintf("Relative error of recovering S0: %e\n", rel_S);
fprintf("Ran %g many outer loops.\n", output.count_outer);
fprintf("Input H was circulant? Answer: %g \n", output.isCirc)

fprintf("****** with preconditioning ******\n")
tic
outputC = gen_matrix_sep_2d_con(M0, G1, G2, lam, para);
toc

rel_Lc = norm(outputC.L - L0, 'fro')/norm(L0, 'fro');
rel_Sc = norm(outputC.S - S0, 'fro')/norm(S0, 'fro');
fprintf("Relative error of recovering L0: %e\n", rel_Lc);
fprintf("Relative error of recovering S0: %e\n", rel_Sc);
fprintf("Ran %g many outer loops.\n", outputC.count_outer);
fprintf("Input H was circulant? Answer: %g \n\n", outputC.isCirc)


