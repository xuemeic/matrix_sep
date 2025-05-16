n1 = 2;
rng(3);
E1 = rand(n1, n1);
E1 = E1./sum(E1, 2);
% E1 will have row sum to be 1

n2 = 3;
E2 = rand(n2, n2);
E2 = E2./sum(E2, 1);
% E2 will have column sum to be 1

m1 = 20;
m2 = 30;
p1 = m1;
p2 = m2;
K = 100;

G1 = kron(eye(m1/n1), E1);
G2 = kron(eye(m2/n2), E2);

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
para.E1 = E1;
para.E2 = E2;
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

[CG1, C1] = precondition(G1);
[CG2, C2] = precondition(G2);
CM0 = pagemtimes(pagemtimes(C1, M0), C2');

tic
outputC = gen_matrix_sep_2d(CM0, CG1, CG2, lam, para);
toc
HS_output = pagemtimes(pagemtimes(G1, outputC.S), G2');

L_output = M0 - HS_output;
rel_Lc = norm(L_output - L0, 'fro')/norm(L0, 'fro');
rel_Sc = norm(outputC.S - S0, 'fro')/norm(S0, 'fro');
fprintf("Relative error of recovering L0: %e\n", rel_Lc);
fprintf("Relative error of recovering S0: %e\n", rel_Sc);
fprintf("Ran %g many outer loops.\n", outputC.count_outer);
fprintf("Input H was circulant? Answer: %g \n", outputC.isCirc)


