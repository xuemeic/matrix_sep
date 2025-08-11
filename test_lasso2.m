% This is when both Gi are circulant
% not sure why it does not work well with FISTA
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

%G1 = eye(m1) + 0.1;

%G2 = eye(m2) + 0.15;
% make S0
num_elements = p1 * p2;
percent_nonzero = 4;
num_nonzeros = K*round((percent_nonzero / 100) * num_elements);
nonzero_idx = randsample(p1*p2*K, num_nonzeros);
S0 = zeros(p1, p2, K);
S0(nonzero_idx) = rand(num_nonzeros, 1) + 1;

B = pagemtimes(pagemtimes(G1, S0), G2');

% run
clear para
para.method = 'ADMM';
para.rho = 1; 
para.max_iter= 500;
para.tol = 1e-5;
%para.is_E = false;
lam =  1e-4;

tic
[X, num_iter, para2] = my_lasso2(G1, G2, B, lam, para);
toc
norm(X - S0, 'fro')


