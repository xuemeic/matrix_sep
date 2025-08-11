function output = gen_matrix_sep_2d_con(M, G1, G2, lam, para)
% preconditioned generalized matrix separation 2d, solved by ADMM
% min||L||_* + lam||S||_1, subject to, L + HS = M, 
% H = G2 otimes G1 = kron(G2, G1)
% H is theoretically m1m2 x p1p2 
% input: 
%      - M: given 3D matrix: m1 x m2 x K
%      - Gi: given matrix: mi x pi
%      - lam: positive scalar, as described above
%      - para has fields: 
%       .rho_outer: step size for ADMM, provided by user
%       .max_iter: max number of iterations for ADMM
%       .tol_outer: convergence tolerance for ADMM
%       .is_E
%       .lasso_method: {'ADMM', 'FISTA'}
%       .lasso_rho: step size for lasso if lasso_method is ADMM
%       .lasso_max_iter: max num of iterations for lasso
%       .lasso_tol: convergence tolerance for lasso

%%%%%%%%%%%%% notes %%%%%%%%%%%%%%%%
% para.is_E needs to be provided
%%%%%%%%%%%%

% output has fields:
% - L: m1 x m2 x K
% - S: p1 x p2 x K
% - count_outer
% - para

% created on 7/19/2025 by XC
% updated on 8/10/2025

para.preconditioned = true;
[~, p1] = size(G1);
[~, p2] = size(G2);


if para.is_E
    E1 = para.E1;
    E2 = para.E2;
    [n1, ~] = size(E1);
    [n2, ~] = size(E2);
    [m1, m2, ~] = size(M);
    %[CE1, C_e1] = precondition(E1);
    [C_e1, CE1, s1, v1, k1] = pre_con(E1); % v1 is n1 by n1
    %{
    [u1, s1, v1] = svd(E1);
    k1 = rank(E1);
    uk1 = u1(:,1:k1); % m by k
    vk1 = v1(:,1:k1); % p by k
    s1 = diag(s1); % a vector
    s1 = s1(1:k1);
    s1_inv = 1./s1;
    CE1 = uk1*vk1';
    C_e1 = uk1*diag(s1_inv)*uk1';
    %}
    [C_e2, CE2, s2, v2, k2] = pre_con(E2);
    
    CG1 = kron(eye(m1/n1), CE1);
    C1 = kron(eye(m1/n1), C_e1);
    
    CG2 = kron(eye(m2/n2), CE2);
    C2 = kron(eye(m2/n2), C_e2);

    % E1, E2 needs to be preconditioned
    para.E1 = CE1;
    para.E2 = CE2;

    %%%%%% dont know why the following two lines don't work
    para.V_G1 = kron(eye(m1/n1), v1);
    para.V_G2 = kron(eye(m2/n2), v2);

    [C2, CG2, ~, v2, k2] = pre_con(G2);
    [C1, CG1, ~, v1, k1] = pre_con(G1);
    para.V_G1 = v1;
    para.V_G2 = v2;
    
    s1 = zeros(p1, 1);
    s1(1:k1) = 1;
    
    s2 = zeros(p2, 1);
    s2(1:k2) = 1;
    para.Sig = s1*s2';


   

else
    [C2, CG2, ~, v2, k2] = pre_con(G2);
    [C1, CG1, ~, v1, k1] = pre_con(G1);
    para.V_G1 = v1;
    para.V_G2 = v2;

    s1 = zeros(p1, 1);
    s1(1:k1) = 1;
    
    s2 = zeros(p2, 1);
    s2(1:k2) = 1;
    para.Sig = s1*s2';

    
end



CM = pagemtimes(pagemtimes(C1, M), C2');
outputC = gen_matrix_sep_2d(CM, CG1, CG2, lam, para);
output.S = outputC.S;
HS_output = pagemtimes(pagemtimes(G1, outputC.S), G2');
output.L = M - HS_output;
output.count_outer = outputC.count_outer;
output.para = outputC.para;


end


function [C, CA, s, v, k] = pre_con(A)
    % A is m by p
    [u, s, v] = svd(A);
    k = rank(A);
    uk = u(:,1:k); % m by k
    vk = v(:,1:k); % p by k
    s = diag(s); % a vector
    sk = s(1:k);
    s_inv = 1./sk;
    CA = uk*vk';
    C = uk*diag(s_inv)*uk';
end

