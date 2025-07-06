function output = gen_matrix_sep(M, H, lam, para)
% generalized matrix separation
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
% updated on 4/9/2025
% updated on 4/21/2025
% updated on 7/5/2025: precompute A^Tb outside of loop

%%%%% pass the parameters
rho_outer = para.rho_outer;
tol_outer = para.tol_outer;
N_outer = para.N_outer;

%%%%% parameters for lasso1()
para_lasso.max_iter = para.N_inner;
para_lasso.rho = para.rho_inner;
para_lasso.tol = para.tol_inner;

%%%%% initialization
[m, n] = size(M);
[~, p] = size(H);

L = zeros(m, n);
S = zeros(p, n);
U = zeros(m, n);

RelChg = 1;
count_outer = 0;
if is_circulant(H)
    isCirculant = true;
    
    para_lasso.isCirculant = isCirculant;
    
    d = abs(fft(H(:,1)));
    para_lasso.coef = d.^2 + para.rho_inner;
    
else
    isCirculant = false;
    para_lasso.isCirculant = isCirculant;
    if para.decomp == "chol"
        L_H = chol(H'*H + rho_inner*eye(m), 'lower');
        para_lasso.decomp = "chol";
        para_lasso.L = L_H;
    elseif para.decomp == "svd"
        [~, sig, V_H] = svd(H'*H);
        para_lasso.decomp = "svd";
        para_lasso.S = sig;
        para_lasso.A2_V = V_H;
    end
end


while RelChg > tol_outer && count_outer < N_outer
    Slast = S;
    Llast = L;
    L = SVT(M - H*S - U, 1/rho_outer);

    [~, S] = lasso1(H, M - U - L, lam/rho_outer, para_lasso);
    
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
output.isCirc = isCirculant;
end


%%%%%%%%%%%%% lasso %%%%%%%%%%%%%%%%%
function [count, s] = lasso1(A, b, a, para_lasso)
% A is H:  p x m
% b can have more than 1 columns; b is p by k
% para_lasso has fields
%       - isCirculant
%       - rho
%       - tol
%       - max_iter
%       - coef: if circulant
%       - A2_V: The V matrix from "A^T*A = V*S*V^T" (pre-compute SVD if not circulant and using svd)
%       - S: THe S matrix from "A^T*A = V*S*V^T".
%       - L: pre-compute cholesky if not circulant and using chol
% output s: m x k
% s = argmin_x {a||x||_1 + 0.5||Ax - b||_2^2}
% update x: solve {A'*A + rho*eye(n)}x = A'*b + rho*(z - u)

max_iter = para_lasso.max_iter;
rho = para_lasso.rho;

[~, k] = size(b);
[~, m] = size(A);
x = zeros(m, k);
z = zeros(m, k);
u = zeros(m, k);
count = 0;
RelChg = 1;
eps = para_lasso.tol;
Atb = A'*b;

while count < max_iter && RelChg > eps
    zlast = z;
    xlast = x;

    rhs = Atb + rho*(z - u);
    % update x
    if para_lasso.isCirculant
    
    coef = para_lasso.coef;
    % [m,n]./[m,1] will divide each col
    x = ifft(fft(rhs)./coef);
    x = real(x);
    elseif (~para_lasso.isCirculant) && (para_lasso.decomp == "chol")
        L = para_lasso.L;
        y = L \ rhs;                    % y is lower triangular: n x p

        
        x = L' \ y;  
    elseif (~para_lasso.isCirculant) && (para_lasso.decomp == "svd")
        V = para_lasso.A2_V;
        sig = diag(para_lasso.S); % column vector
        x = (V'*rhs)./(sig + rho);
        x = V*x;
    end


    % update z
    z = SoftThresh(x + u, a/rho);

    % update u
    u = u + x - z;

    % Check convergence
    % finding the relative error
    xdn = norm(x - xlast, 'fro');
    zdn = norm(z - zlast, 'fro');
    xn = norm(xlast, 'fro');
    zn = norm(zlast, 'fro');

    % updating stopping critera
    RelChg = (xdn^2 + zdn^2)^0.5 / ((xn^2 + zn^2)^0.5 + 1);

    count = count + 1;
end

s = x;


end

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

