% based on test_1d_avg
% this is to test which rho's to pick
% updated on 5/27/2026

n = 100;
m = n - 1;
p = m;

first_row = [-1 1 zeros(1, m - 2)]; % First row of the circulant matrix H
H = toeplitz([first_row(1), fliplr(first_row(2:end))], first_row); 

rng(2);
% Make S0
s = 5;
S0 = zeros(p, n);
nonzero_idx = randsample(p*n, s*n);
S0(nonzero_idx) = rand(s*n, 1) + 1;

% Make L
rank_r = floor(min(m, n)*0.05);
L1 = randn([m, rank_r]);  % m x r matrix
L1 = orth(L1);
L2 = randn([n, rank_r]);  % r x n matrix
L2 = orth(L2);
% Multiply to get the low-rank matrix (m x n)
L0 = L1*L2';

rho_outer_list = 10.^linspace(-2, 2, 5);
rho_outer_list = reshape([1;5]*rho_outer_list, [1, 10]);

rho_inner_list = rho_outer_list(2:9);

% recover
M0 = L0 + H*S0;

para.N_outer = 100;
para.N_inner = 20;
para.tol_outer = 1e-7;
para.tol_inner = 1e-5;
para.decomp = "svd"; 
para.preconditioned = false;
para.lasso_max_iter = 20;
para.max_iter = 100;
para.lasso_method = 'ADMM';
lam =  1/sqrt(n);



ol = length(rho_outer_list);
il = length(rho_inner_list);

record = zeros(ol, il, 3);

for i = 1:ol
    for j = 1:il
        para.rho_outer = rho_outer_list(i); 
        para.lasso_rho = rho_inner_list(j); 
        tic
        output = gen_matrix_sep(M0, H, lam, para);
        tim = toc;
        rel_L = norm(output.L - L0, 'fro')/norm(L0, 'fro');
        rel_S = norm(output.S - S0, 'fro')/norm(S0, 'fro');
        record(i, j, :) = [rel_L, rel_S, tim];
    end
end


record_c = zeros(ol, il, 3);

for i = 1:ol
    for j = 1:il
        para.rho_outer = rho_outer_list(i); 
        para.lasso_rho = rho_inner_list(j); 
        tic
        output = gen_matrix_sep_con(M0, H, lam, para);
        tim = toc;
        rel_L = norm(output.L - L0, 'fro')/norm(L0, 'fro');
        rel_S = norm(output.S - S0, 'fro')/norm(S0, 'fro');
        record_c(i, j, :) = [rel_L, rel_S, tim];
    end
end

%% heatmap
figure(1)

subplot(1,2,1)
imshow(record(:,:,1) + record(:,:,2))
yticklabels(rho_outer_list) % i are for rows, y label
xticklabels(rho_inner_list)
colormap("parula");
set(gca,'YDir','normal')
axis on
colorbar;
clim([1e-5,.5]);
xl = xlabel('$\rho_I$', 'interpreter', 'latex');
yl = ylabel('$\rho_O$', 'interpreter', 'latex');
xl.FontSize = 15;
yl.FontSize = 15;
t1=title("Relative Error of Recovering $S_0, L_0$: no precprocessing", "Interpreter","latex");
t1.FontSize = 16;

subplot(1,2,2)
imshow(record_c(:,:,1) + record_c(:,:,2))
yticklabels(rho_outer_list) % i are for rows, y label
xticklabels(rho_inner_list)
colormap("parula");
set(gca,'YDir','normal')
axis on
colorbar;
clim([1e-5,.5]);
xl = xlabel('$\rho_I$', 'interpreter', 'latex');
yl = ylabel('$\rho_O$', 'interpreter', 'latex');
xl.FontSize = 15;
yl.FontSize = 15;
t1=title("Relative Error of Recovering $S_0, L_0$: preprocessing", "Interpreter","latex");
t1.FontSize = 16;

%%
%print("-f1", 'figs/heat_rho', '-djpeg', '-r500')



