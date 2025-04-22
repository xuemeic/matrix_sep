%% Make G1, G2
K = 700;
n = 320;
m = 240;

m1 = m; % Dimension of G1
p1 = m; % Dimension of G1
m2 = n;  % Dimension of G2
p2 = n;  % Dimension of G2

%p = p1 * p2;
%m = m1 * m2;

% Try on video with E = randn(2,2);
% E=[0.1,0.2;0.3,0.4]


method = 'video';

switch method
    case 'random'
        E = randn(2,2);
        % Random Gaussian matrices
        G1 = randn(m1, p1);
        G2 = randn(m2, p2);
    case 'orthogonal'
        E = randn(2,2);
        % Orthogonalized random matrices
        [Q1, ~] = qr(randn(m1, p1), 0);
        [Q2, ~] = qr(randn(m2, p2), 0);
        G1 = Q1;
        G2 = Q2;
    case 'circulant'
    E = randn(2,2);
    first_row = [-1 1 zeros(1, m1-2)]; % First row of the circulant matrix H
    G1 = toeplitz([first_row(1), fliplr(first_row(2:end))], first_row); 

    first_row = [-1 1 zeros(1, m2-2)]; % First row of the circulant matrix H
    G2 = toeplitz([first_row(1), fliplr(first_row(2:end))], first_row); 

    case 'video'
        %E = 0.5 * ones(2,2);
        E = [0.4,0.6;0.6,0.4];

        %G1 = kron(eye(m/2),E); %multiply on the left
        G1 = kron(randn(m/2, m/2), E);

        %G2 = kron(eye(n/2),E); % multiply on the right

        G2 = kron(randn(n/2, n/2), E);

end

% Normalize columns
% for j = 1:p1
%     G1(:, j) = G1(:, j) / norm(G1(:, j));
% end
% for j = 1:p2
%     G2(:, j) = G2(:, j) / norm(G2(:, j));
% end

% Check condition number of H
%H = kron(G2, G1);

%cond_num = cond(H);
%disp(['Condition number of H (' method '): ' num2str(cond_num)]);

%%

v = VideoReader('Converted.avi');  % or use 'Converted.mp4' if you re-encoded to mp4

% Preallocate cell array to store grayscale frames
gray_frames = cell(1, floor(v.Duration * v.FrameRate));
frameIdx = 1;

while hasFrame(v)
    frame = readFrame(v);  % Read RGB frame

    % NTSC grayscale conversion and cast to double
    gray = ...
        0.2989 * double(frame(:,:,1)) + ...
        0.5870 * double(frame(:,:,2)) + ...
        0.1140 * double(frame(:,:,3));

    gray_frames{frameIdx} = gray;
    frameIdx = frameIdx + 1;
end

% Stack into a 3D double array
gray_video = cat(3, gray_frames{:});

V = gray_video/255;
%V = gray_video;

%%
% Get original frame size
[h, w, numFrames] = size(V);
%numFrames = 400;

% Calculate center crop indices
half_crop = 20;  
center_row = floor(h/2);
center_col = floor(w/2);

row_range = (center_row - half_crop + 1):(center_row + half_crop);
col_range = (center_col - half_crop + 1):(center_col + half_crop);

% Initialize cropped video array
V_cropped = zeros(40, 40, numFrames);

% Crop each frame
for i = 1:numFrames
    frame = V(:, :, i);
    V_cropped(:, :, i) = frame(row_range, col_range);
end

%% Make G1, G2
m = 40; n = 40; K = numFrames;

%m = h; n = w;


E = 0.5 * ones(2,2);

%E = randn(2,2);
%E1 = [0.8,0.2;0.6,0.4];
%E = [0.4,0.6;0.6,0.4];

G1 = kron(eye(m/2), E); %multiply on the left
G2 = kron(eye(n/2), E); % multiply on the right



%%
i = 100;
subplot(1,2,1)
imshow(M(:,:,i),[])
subplot(1,2,2)
imshow(V_cropped(:,:,i),[])

%% Make S, L to make M

% Make S

num_elements = m * n * K;
percent_nonzero = 5;
num_nonzeros = round((percent_nonzero / 100) * num_elements);

row_indices = randi(m * n, num_nonzeros, 1);  % Random row indices
col_indices = randi(K, num_nonzeros, 1);  % Random column indices

values = rand(num_nonzeros, 1) + 1;

Strue = sparse(row_indices, col_indices, values, m * n, K);
Strue = full(Strue);

Strue = reshape(Strue, [m, n, K]);




%% Make L

rank_r = ceil(0.1 * (m * n));
L1 = rand([m*n, rank_r]);  % mn x r matrix
L2 = rand([rank_r, K]);  % r x K matrix

% Multiply to get the low-rank matrix (mn x K)
Ltrue = L1 * L2;
Ltrue = Ltrue * 10;

Ltrue = reshape(Ltrue, [m, n, K]);

%% Precondition 

% Pre-conditioning for G1
[U1, S1, V1] = svd(G1, 'econ');
k1 = rank(G1);
U1 = U1(:, 1:k1);
V1 = V1(:, 1:k1);
S2_1 = diag(1 ./ diag(S1(1:k1, 1:k1)));
G1_sharp = V1 * S2_1 * U1';
C1 = U1 * S2_1 * U1';
CG1 = U1 * V1';

% Pre-conditioning for G2
[U2, S2, V2] = svd(G2, 'econ');
k2 = rank(G2);
U2 = U2(:, 1:k2);
V2 = V2(:, 1:k2);
S2_2 = diag(1 ./ diag(S2(1:k2, 1:k2)));
G2_sharp = V2 * S2_2 * U2';
C2 = U2 * S2_2 * U2';
CG2 = U2 * V2';


%% Make M

%M0 = Ltrue + Strue;





% H_video = pagemtimes(pagemtimes(G1, ones(size(V_cropped))), G2');
% blurr_M = pagemtimes(pagemtimes(G1, V_cropped), G2');
% input_M = H_video - blurr_M;


%%%%%%%%%%%%%%%%%%%%%%%%%%

M0 = pagemtimes(pagemtimes(G1, V_cropped), G2');

%M0 = pagemtimes(pagemtimes(G1, V), G2');
M = -M0 + pagemtimes(pagemtimes(G1, ones(size(M0))), G2');



%M = pagemtimes(pagemtimes(G1, M0), G2');

%HS = pagemtimes(pagemtimes(G1, Strue), G2');
%M = Ltrue + HS;

%% Add the C
% 
% CM = pagemtimes(pagemtimes(C1, input_M), C2');
% 
% CM = pagemtimes(pagemtimes(C1, blurr_M), C2');

CM = pagemtimes(pagemtimes(C1, M), C2');


%% Testing the gen matrix 2d
para.rho_outer = 1;
para.rho_inner = 100;

para.E = E; % testing svd so this is not circulant
%para.E = randn(2,2);



para.N_inner = 170; % 170 works good for 40 x 40 x 400
para.N_outer = 7; % 7

para.inner_tol = 1e-6;
para.tol = 1e-6;

para.mode = false; % true is cholesky, false is svd

lam =  1/sqrt(min(m*n, K)); % since matrix is mn by K

%lam = 0.1/sqrt(min(m*n, K));
%lam = 1/sqrt(m);

% % trying normalization
%M = M / norm(M, 'fro');
% G1 = G1 / norm(G1, 'fro');
% G2 = G2 / norm(G2, 'fro');

% testing code

% -CM

%output = general_matrix_sep_2d(M, G1, G2, lam, para)

output = general_matrix_sep_2d(CM, CG1, CG2, lam, para)
HS_output = pagemtimes(pagemtimes(G1, output.S), G2');

L_output = M - HS_output;

%L_output = M0 - HS_output;

%L_output = output.L;

% S_output = output['S']
S_output = output.S; 

%% % Compute relative errors


rel_L = norm(L_output - Ltrue, 'fro') / norm(Ltrue, 'fro')
rel_S = norm(S_output - Strue, 'fro') / norm(Strue, 'fro')

% - M = -L + H(-S)

%% Save the variables
save('full-video-matrixsep-results2.mat', 'M', 'L_output', 'S_output', '-v7.3');




%% View the video frame 
i = 1;
subplot(1,3,1)
imshow(V_cropped(:,:,i),[])
subplot(1,3,2)
imshow(S_output(:,:,i),[])
subplot(1,3,3)
imshow(L_output(:,:,i),[])

%% Load the variables to view

%load('1234-E-matrixsep-results.mat');

load('full-video-matrixsep-results2.mat')


%%

numFrames = size(V, 3);

i = 1;

figure;

% -M=-L+H(-S)
%

while i <= numFrames

    subplot(1,4,1)
    imshow(V_cropped(:,:,i), []);
    title(['Frame ', num2str(i)]);


    subplot(1,4,2)
    imshow(M(:,:,i), []);
    title(['Blurred - Frame ', num2str(i)]);

    subplot(1,4,3)
    imshow(1 - S_output(:,:,i), []);
    title('1 - S');

    subplot(1,4,4)
    imshow(L_output(:,:,i), []);
    title('L');

    k = waitforbuttonpress;
    
    
    if k == 1 
        key = get(gcf, 'CurrentCharacter');
        if strcmp(key, char(27))  
            break;
        end
    end

    i = i + 1;
end


%% View the video frame 
i = 70;
subplot(1,3,1)
imshow(BV(:,:,i),[])
subplot(1,3,2)
imshow(1- S_output(:,:,i),[])
subplot(1,3,3)
imshow(-1-L_output(:,:,i),[])


%% Comparing with the 1d function

para.rho_outer = 1;
para.rho_inner = 1;

para.N_outer = 400;
para.N_inner = 100;
para.mode = false;

lam =  lam;
L0 = reshape(Ltrue, [m*n, K]);
S0 = reshape(Strue, [m*n, K]);
M0 = reshape(M, [m*n, K]);


output = gen_matrix_sep_chen(M0, H, lam, para)
rel_L = norm(output.L - L0, 'fro')/norm(L0, 'fro')
rel_S = norm(output.S - S0, 'fro')/norm(S0, 'fro')

%%
m = 5;
n = 10;
K = 3;


[G1, G2, Ltrue, Strue] = create_simple_data(m, n, K);
size(G1)      % 5x5
size(G2)      % 10x10
size(Ltrue)   % 5x10x3
size(Strue)   % 5x10x3












