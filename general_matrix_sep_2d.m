function output = general_matrix_sep_2d(M, G1, G2, lam, para)
% generalized matrix separation
% min||L||_* + lam||S||_1, subject to, L + HS = M, H = G2 otimes G1
% input: M: given matrix: m1 x m2 x K
% input: Gi: given matrix: pi x mi
% input: lam
% input: para has fields: 
% - rho_outer, rho_inner, N_outer, N_inner, tol, decomp

% output has fields
% - L: m1 x m2 x K
% - S: p1 x p2 x K
% - count_outer

% written by Owen Deen, 11/23/2024

rho_outer = para.rho_outer;

para_lasso.rho = para.rho_inner;

% if there is E this is only in block circulant case
E = para.E;

N_outer = para.N_outer;

para_lasso.max_iter = para.N_inner;
%para_lasso.tol = para.tol;
para_lasso.tol = para.inner_tol;

tol = para.tol;

[m1, m2, K] = size(M);
[p1, m1] = size(G1);
[p2, m2] = size(G2);

[n1, ~] = size(E);


L = zeros(m1, m2, K);
S = zeros(p1, p2, K);
U = zeros(m1, m2, K);

RelChg = 1;
count_outer = 0;


if is_circulant(G1) && is_circulant(G2)
    disp("at 1st block");
    isCirculant = true;
    
    para_lasso.isCirculant = isCirculant;
    
    % Since we are computing F G^T G F^* = D

    D1 = abs(fft(G1(:,1))).^2;

    D2 = abs(fft(G2(:,1))).^2;

    para_lasso.D = D1 * D2' + para.rho_inner;

elseif is_circulant(E)
    disp("at 2nd block");
    isCirculant = true;
    para_lasso.isCirculant = isCirculant;
    

    a = abs(fft(E(:,1))).^2;

    b = abs(fft(E(:,1))).^2;

    Lc = kron(ones(m1/n1, m2/n1), a * b');

    para_lasso.D = Lc + para.rho_inner;
    





else
    disp("at 3rd block");
    isCirculant = false;
    para_lasso.isCirculant = isCirculant;

    % calculating svd of G1 and G2
    [~, S1, para_lasso.V1] = svd(G1' * G1);
    S1 = diag(S1); % m1
 

    [~, S2, para_lasso.V2] = svd(G2' * G2);
    S2 = diag(S2); % m2

    para_lasso.Sig = S1 * S2'; %m1 by m2


    para_lasso.decomp = "svd";
   
end



while RelChg > tol && count_outer < N_outer
    disp(count_outer)
%while (r_primal > tol || s_dual > tol) && count_outer < N_outer
    Slast = S;
    Llast = L;
    
    HS = framewise_mult(S, G1, G2');
    L = SVT(reshape(M - HS - U, [m1*m2, K]), 1/rho_outer);
    L = reshape(L, [m1, m2, K]);
    
    [count_inner,RelChg_lasso, S] = lasso2(G1, G2, M - U - L, lam/rho_outer, para_lasso);
    
    HS = framewise_mult(S, G1, G2');
    U = U + L + HS - M;
    count_outer = count_outer + 1;

    % 
    % 
    % % Check convergence
    % % finding the relative error

    
    Ldn = norm(L - Llast, 'fro');
    Sdn = norm(S - Slast, 'fro');
    Ln = norm(Llast, 'fro');
    Sn = norm(Slast, 'fro');
    % 
    % % updating stopping critera
    top = (Ldn^2 + Sdn^2)^0.5;
    bottom = ((Ln^2 + Sn^2)^0.5 + 1);
    RelChg =  top / bottom;
    % 
    % 

end

output.L = L;
output.S = S;
output.count_outer = count_outer;
output.count_inner = count_inner;
output.relchg = RelChg;
output.relchg_top = top;
output.relchg_bottom = bottom;
output.Ldn = Ldn;
output.Sdn = Sdn;
output.Ln = Ln;
output.Sn = Sn;
output.RelChg_lasso = RelChg_lasso;
end

function [count,RelChg, S] = lasso2(G1, G2, Video, lambda, para_lasso)
% A is H: 
% G1, G2 are from H
% B can have more than 1 columns; b is m1 by m2 by k
% output S will have the same size as B

% s = argmin_x {a||x||_1 + 0.5||HS - L - U + M||}
% update x: solve {A'*A + rho*eye(n)}x = A'*b + rho*(z - u)

max_iter = para_lasso.max_iter;
rho = para_lasso.rho;

[p1, p2, k] = size(Video);
[~, m1] = size(G1);
[~, m2] = size(G2);
X = zeros(m1, m2, k);
Z = X;
U = X;
count = 0;
RelChg = 1;
eps = para_lasso.tol;
%coef = A'*A + rho*eye(n);

while count < max_iter && RelChg > eps
    Zlast = Z;
    Xlast = X;

    rhs = framewise_mult(Video, G1', G2) + rho*(Z - U);
 

    % update X
    if para_lasso.isCirculant
    
        D = para_lasso.D;
         
        % % 
        temp = fft2(rhs);
        % % 
        % % 
        temp2 =  temp ./ D;
        % % 
        % % 
        X = real(ifft2(temp2));
        


      

    else 
        % this is for the 1d case
        % V = para_lasso.V;
        % sig = diag(para_lasso.sig); % column vector
        % x = (V'*rhs)./(sig + rho);
        % x = V*x;
        V1 = para_lasso.V1;
        S = para_lasso.Sig; % m1 by m2
        V2 = para_lasso.V2;

        %U1 = para_lasso.U1;
        %U2 = para_lasso.U2;
        

        % needs to be 3d

        % Sig is 10 * 20

        %Sig = abs(S1) * abs(S2');
        
        temp = framewise_mult(rhs, V1', V2);

        temp2 = temp ./ (S + rho);

        
        X = framewise_mult(temp2, V1, V2');

      

    end


    % solve coef*x = rhs
    %x = coef\rhs;

    % update z
    Z = SoftThresh(X + U, lambda/rho);

    % update u
    U = U + X - Z;

    % Check convergence
    % finding the relative error
    xdn = norm(X - Xlast, 'fro');
    %zdn = norm(Z - Zlast, 'fro');
    xn = norm(Xlast, 'fro');
    %zn = norm(Zlast, 'fro');

    % updating stopping critera
    %RelChg = (xdn^2 + zdn^2)^0.5 / ((xn^2 + zn^2) + 1)^0.5;

    RelChg = (xdn) / ((xn) + 1);

    count = count + 1;
end

S = X;


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


function v = framewise_mult(V, left, right)
    % performs the frame-wise multiplication on a tensor V.
    %
    % Inputs:
    %   V     - an m1 x m2 x k tensor 
    %   left  - a p x m1 matrix to be applied on the left of each frame
    %   right - a m2 x q matrix to be applied on the right of each frame
    %
    % Output:
    %   v     - a p x q x k tensor where each frame is transformed as left * V(:,:,i) * right

    % validate input dimensions
    [m1_V, m2_V, ~] = size(V);
    [~, m1_left] = size(left);
    [m2_right, ~] = size(right);

    if m1_left ~= m1_V
        error('Number of columns of left (%d) must match first dimension of V (%d).', m1_left, m1_V);
    end
    if m2_right ~= m2_V
        error('Number of rows of right (%d) must match second dimension of V (%d).', m2_right, m2_V);
    end

    LV = pagemtimes(left, V);        % Size: [p, m2_V, k]

    v = pagemtimes(LV, right);       % Size: [p, q, k]
end




















