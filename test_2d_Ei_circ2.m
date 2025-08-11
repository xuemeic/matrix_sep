% Gi = kron(eye(mi/ni), Ei);
% E1, E2 will be circulant
% but n1 and n2 are relatively big

m1 = 20;
m2 = 30;
p1 = m1;
p2 = m2;
K = 100;

n1 = 10;
rng(3);
E1r1 = rand(1, n1); 
E1r1 = E1r1/sum(E1r1); % First row of the circulant matrix 
E1 = toeplitz([E1r1(1), fliplr(E1r1(2:end))], E1r1); 
% E1 will have row sum to be 1

n2 = 10;
E2r1 = rand(1, n2); 
E2r1 = E2r1/sum(E2r1); % First row of the circulant matrix 
E2 = toeplitz([E2r1(1), fliplr(E2r1(2:end))], E2r1); 
E2 = E2';
% E2 will have column sum to be 1

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
clear para
para.rho_outer = 1; 
para.rho_inner = 1; 
para.max_iter = 200;
para.lasso_max_iter = 20;
para.tol_outer = 1e-8;
para.lasso_tol = 1e-6;
para.E1 = E1;
para.E2 = E2;
para.lasso_method = 'ADMM';
para.preconditioned = false;
lam =  1/sqrt(max(m1*m2, K));

para.is_E = true;
tic
output = gen_matrix_sep_2d(M0, G1, G2, lam, para);
t = toc;
desired_print(output, t, L0, S0);


tic
output = gen_matrix_sep_2d_con(M0, G1, G2, lam, para);
t = toc;
desired_print(output, t, L0, S0);

para.lasso_method = 'FISTA';
tic
outputC = gen_matrix_sep_2d(M0, G1, G2, lam, para);
t = toc;
desired_print(outputC, t, L0, S0);

tic
outputC = gen_matrix_sep_2d_con(M0, G1, G2, lam, para);
t = toc;
desired_print(outputC, t, L0, S0);



function desired_print(o, t, L0, S0)
if o.para.preconditioned
    sp = "with preconditioning";
else
    sp = "no preconditioning";
end
fprintf("******* Gi block circ, %s, lasso by %s ******\n", sp, o.para.lasso_method)
fprintf("Number of iterations: %g.\n", o.count_outer);
fprintf("Duration: %.3f seconds.\n", t)
rel_L = norm(o.L - L0, 'fro')/norm(L0, 'fro');
rel_S = norm(o.S - S0, 'fro')/norm(S0, 'fro');
fprintf("Relative error of recovering L0: %e\n", rel_L);
fprintf("Relative error of recovering S0: %e\n\n", rel_S);
end
