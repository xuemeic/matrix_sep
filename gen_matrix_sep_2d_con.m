function output = gen_matrix_sep_2d_con(M, G1, G2, lam, para)
% preconditioned generalized matrix separation 2d
% min||L||_* + lam||S||_1, subject to, L + HS = M, 
% H = G2 otimes G1 = kron(G2, G1)
% H is theoretically m1m2 x p1p2 
% input: 
%      - M: given 3D matrix: m1 x m2 x K
%      - Gi: given matrix: mi x pi
%      - lam: positive scalar, as described above
%      - para has fields: 
%        - rho_outer
%        - rho_inner
%        - N_outer
%        - N_inner
%        - tol_outer
%        - tol_inner
%        - is_E: if true, then E1, E2 need to be provided
%        - E1
%        - E2


% output has fields:
% - L: m1 x m2 x K
% - S: p1 x p2 x K
% - count_outer
% - isCirc

% created on 7/19/2025

if para.is_E
    E1 = para.E1;
    E2 = para.E2;
    [n1, ~] = size(E1);
    [n2, ~] = size(E2);
    [m1, m2, ~] = size(M);
    [CE1, C_e1] = precondition(E1);
    
    CG1 = kron(eye(m1/n1), CE1);
    C1 = kron(eye(m1/n1), C_e1);
    [CE2, C_e2] = precondition(E2);
    CG2 = kron(eye(m2/n2), CE2);
    C2 = kron(eye(m2/n2), C_e2);

    % E1, E2 needs to be preconditioned
    para.E1 = CE1;
    para.E2 = CE2;
else
    [CG1, C1] = precondition(G1);
    [CG2, C2] = precondition(G2);
    
end

CM = pagemtimes(pagemtimes(C1, M), C2');
outputC = gen_matrix_sep_2d(CM, CG1, CG2, lam, para);
output.S = outputC.S;
HS_output = pagemtimes(pagemtimes(G1, outputC.S), G2');
output.L = M - HS_output;
output.count_outer = outputC.count_outer;
output.isCirc = outputC.isCirc;

