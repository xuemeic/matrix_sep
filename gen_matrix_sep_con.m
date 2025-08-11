function output = gen_matrix_sep_con(M, H, lam, para)
% preconditioned generalized matrix separation
% min||L||_* + lam||S||_1, subject to, L + HS = M
% input: M: given matrix: m x n
% input: H: given matrix: m x p, not necessarily circulant
% input: lam
% input: para has fields: 
%        - rho_outer 
%        - rho_inner
%        - N_outer
%        - N_inner
%        - tol_outer
%        - tol_inner
%        - decomp: 'svd' or 'chol'

% output has fields
% - L: m x n
% - S: p x n
% - count_outer
% - para
% created on 7/18/2025, Xuemei Chen
% updated on 8/5/2025


[u, s, v] = svd(H);
k = rank(H);
uk = u(:,1:k); % m by k
vk = v(:,1:k); % p by k

s = diag(s); % a vector
s = s(1:k);
s_inv = 1./s;

C = uk*diag(s_inv)*uk';

CH = uk*vk';
para.preconditioned = true;
para.V_H = v;

[~, p] = size(H);
ss = zeros(p, 1);
ss(1:k) = 1;
para.S2_H = ss;
outputC = gen_matrix_sep(C*M, CH, lam, para);


output.S = outputC.S;
output.L = M - H*outputC.S;
output.count_outer = outputC.count_outer;
%output.isCirc = outputC.isCirc;
output.para = outputC.para;