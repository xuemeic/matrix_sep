function [Hc, C] = precondition(H)
% H is arbitrary matrix: m by n
% ker(H) = ker(Hc)
% cond(Hc) = 1
% C*H = Hc
[u, s, v] = svd(H);
k = rank(H);
u = u(:,1:k); % m by k
v = v(:,1:k); % n by k
Hc = u*v'; % m by n
s = diag(s); % a vector
s_inv = 1./s(1:k);

C = u*diag(s_inv)*u';