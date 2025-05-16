load('data/video.mat'); % get V: 240 by 320 by 300
V = uint8(V); % values are integers from 0 to 255
V = im2double(V); % ranges from 0 to 1

[h, w, numFrames] = size(V);


% Calculate center crop indices
half_crop = 24;  
center_row = floor(h/2);
center_col = floor(w/2);

row_range = (center_row - half_crop + 1):(center_row + half_crop);
col_range = (center_col - half_crop + 1):(center_col + half_crop);

% Initialize cropped video array
V_cropped = V(row_range, col_range, :);

[m1, m2, K] = size(V);

n1 = 2;
rng(3);

E1r1 = rand(1, n1); 
E1r1 = E1r1/sum(E1r1); % First row of the circulant matrix 
E1 = toeplitz([E1r1(1), fliplr(E1r1(2:end))], E1r1); 
% E1 will have row sum to be 1

n2 = 4;
E2r1 = rand(1, n2); 
E2r1 = E2r1/sum(E2r1); % First row of the circulant matrix 
E2 = toeplitz([E2r1(1), fliplr(E2r1(2:end))], E2r1); 
E2 = E2';
% E2 will have column sum to be 1

G1 = kron(eye(m1/n1), E1);
G2 = kron(eye(m2/n2), E2);

% preconditioning is appliedm
[CG1, C1] = precondition(G1);
[CG2, C2] = precondition(G2);

M0 = pagemtimes(pagemtimes(G1, V), G2');
% M0 is very blurred

M = -M0 + pagemtimes(pagemtimes(G1, ones(size(M0))), G2');
CM = pagemtimes(pagemtimes(C1, M), C2');

para.rho_outer = 1;
para.rho_inner = 1;
% 48 seconds for 10, 3
% 108s for 10, 10
% 557s for 50, 10
% 1188s = 20m for 100, 10
para.N_outer = 100;
para.N_inner = 10;
para.tol_outer = 1e-8;
para.tol_inner = 1e-6;
para.is_E = false;
lam =  1/sqrt(min(m1*m2, K))*0.2;

tic
outputC = gen_matrix_sep_2d(CM, CG1, CG2, lam, para);
toc
HS_output = pagemtimes(pagemtimes(G1, outputC.S), G2');
L_output = M - HS_output;

%%
figure(1)
k = 10;
subplot(1, 4, 1)
imshow(V(:,:, k))
title('original')
subplot(1, 4, 2)
imshow(M0(:,:, k))
title('Blurred')
subplot(1, 4, 3)
imshow(1-outputC.S(:,:, k), [])
title('deblurred sparse component')
subplot(1, 4, 4)
imshow(-L_output(:,:, k), [])
title('background')
%%
%save('results/full_vid.mat', "outputC", "V", "M0", "L_output")
