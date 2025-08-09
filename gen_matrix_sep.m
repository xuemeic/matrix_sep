function output = gen_matrix_sep(M, H, lam, para)
% generalized matrix separation
% use ADMM to solve the following
% min||L||_* + lam||S||_1, subject to, L + HS = M
% input: 
%   M: given matrix: m x n
%   H: given matrix: m x p, not necessarily circulant
%   lam: positive scalar
%   para has fields: 
%       .rho_outer: step size for ADMM, provided by user
%       .max_iter: max number of iterations for ADMM
%       .tol_outer: convergence tolerance for ADMM
%       .preconditioned: true or false. If true, H has been preconditioned
%       with all singular values to be 1
%       .lasso_decomp: {'svd', 'chol'}
%       .lasso_method: {'ADMM', 'FISTA'}
%       .lasso_rho: step size for lasso if lasso_method is ADMM
%       .lasso_max_iter: max num of iterations for lasso
%       .lasso_tol: convergence tolerance for lasso

% para.preconditioned needs to be provided by user!

% output has fields
% - L: m x n
% - S: p x n
% - count_outer
% - isCirc: TRUE or FALSE. TRUE if H is circulant.
% updated on 4/9/2025
% updated on 4/21/2025
% updated on 7/5/2025: precompute A^Tb outside of loop
% updated on 8/5/2025: use my_lasso() and take advantage of V if
% preconditioning H

[m, n] = size(M);
[~, p] = size(H);
%% default parameter values
rho_outer_default = 1;
max_iter_default = 200;
tol_outer_default = 1e-7;
lasso_method_default = 'FISTA';
lasso_decomp_default = 'svd';
%% pass the parameters
if ~isfield(para, 'preconditioned')
    error('Please specify "para.preconditioned" value.')
end

if ~isfield(para,'rho_outer')
    para.rho_outer = rho_outer_default;        
end

if ~isfield(para, 'max_iter')
    para.max_iter = max_iter_default;    
end

if ~isfield(para, 'tol_outer')
    para.tol_outer = tol_outer_default;   
end

if ~ isfield(para, 'lasso_method')
    para.lasso_method = lasso_method_default;   
end

if ~ isfield(para, 'lasso_decomp')
    para.lasso_decomp = lasso_decomp_default; 
end

rho_outer = para.rho_outer;
N_outer = para.max_iter;
tol_outer = para.tol_outer;
lasso_method = para.lasso_method;
lasso_decomp = para.lasso_decomp;
%% parameters for my_lasso()
para_lasso.max_iter = para.lasso_max_iter;
para_lasso.rho = para.lasso_rho; % if method is ADMM
para_lasso.tol = para.lasso_tol;
para_lasso.method = lasso_method;
para_lasso.decomp = lasso_decomp; % default value

% para_lasso.isCirculant is determined later

if strcmp(lasso_method, 'FISTA')
    % FISTA
    % H circulant does not matter
    %isCirculant = [];
    if para.preconditioned
        para_lasso.L = 1; % this is biggest singular value of H'H
    else
        s = svd(H'*H);
        s = diag(s);
        para_lasso.L = s(1);
    end
elseif strcmp(lasso_method, 'ADMM')
    % ADMM
    if is_circulant(H)
        isCirculant = true; 
        para_lasso.isCirculant = isCirculant; 
        d = abs(fft(H(:,1)));
        para_lasso.coef = d.^2;
    else
        % more common case for ADMM
        isCirculant = false;
        para_lasso.isCirculant = isCirculant;
        if (para_lasso.decomp == "svd") && para.preconditioned
            % singular values of H'H
            para_lasso.S2_A = para.S2_H; % p by 1
            % singular vectors of H or H'H
            para_lasso.V_A = para.V_H; % p by p
        elseif (para_lasso.decomp == "svd") && (~para.preconditioned)
            [para_lasso.V_A, para_lasso.S2_A] = pref(H);
        elseif para_lasso.decomp == "chol"
            L_H = chol(H'*H + para.lasso_rho*eye(p), 'lower');
            para_lasso.L = L_H;
        end

            
    end
end


%% initialization and main loop
L = zeros(m, n);
S = zeros(p, n);
U = zeros(m, n);

RelChg = 1;
count_outer = 0;
while RelChg > tol_outer && count_outer < N_outer
    Slast = S;
    Llast = L;
    L = SVT(M - H*S - U, 1/rho_outer);
    [S, ~] = my_lasso(H, M - U - L, lam/rho_outer, para_lasso);  
    
    U = U + L + H*S - M;
    count_outer = count_outer + 1;

    % Check convergence
    % finding the relative error
    Ldn = norm(L - Llast, 'fro');
    Sdn = norm(S - Slast, 'fro');
    Ln = norm(Llast, 'fro');
    Sn = norm(Slast, 'fro');

    % updating stopping critera
    RelChg = (Ldn^2 + Sdn^2)^0.5 / ((Ln^2 + Sn^2)^0.5 + 1);
end

output.L = L;
output.S = S;
output.count_outer = count_outer;
output.para = para;
end

%% helper functions
function flag = is_circulant(A)
    [m, n] = size(A);
    
    % we need A to be square if circulant
    if m ~= n
        flag = false;
        return;
    end

    first_row = A(1, :);
    
    % check if each row is a shift of the first row
    for i = 2:m
        expected_row = circshift(first_row, [0, i-1]);
        if any(A(i, :) ~= expected_row)
            flag = false;
            return;
        end
    end
    flag = true;
end

function [V, sig] = pref(A)
% prefactorization of A
% A is m by p
% [~, sig, V] = svd(A'*A)
% V is square p by p
% sig is p by 1
% this won't be used if H is preconditioned

[m, p] = size(A);

if p < m % tall
   [~, sig, V] = svd(A'*A); % V is p by p       
   sig = diag(sig); % p by 1
            
else % fat, m smaller
   [~, s, V] = svd(A); % V is p by p 
   s = diag(s); % m by 1
   sig = zeros(p, 1);
   sig(1:length(s)) = (abs(s)).^2;
end
end



