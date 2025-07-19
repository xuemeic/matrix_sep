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
% - isCirc: TRUE or FALSE. TRUE if H is circulant.
% created on 7/18/2025

[CH, C] = precondition(H);
outputC = gen_matrix_sep(C*M, CH, lam, para);


output.S = outputC.S;
output.L = M - H*outputC.S;
output.count_outer = outputC.count_outer;
output.isCirc = outputC.isCirc;