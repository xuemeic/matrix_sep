% size of matrix
m = 100; % rows of L, M, H
n = 100; % columns of L, M, S
p = 100; % rows of S, column of H

max_rank = 30; % Maximum rank
max_sparsity = 30; % Maximum sparsity

%epsilon = 10e-2; % Threshold for relative error, successful trial
num_trials = 10; % Number of trials for each sparsity-rank pair

% choose the method to generate H
method = 'circulant'; 

% choose the method to generate S
sparse_method = 'uniform';


switch method
    case 'random'
        % Random Gaussian matrix

        % Get seed so that H is the same for each experiment
        seed = rng;

        % Set seed for H
        rng(100);

        H = randn(m, p);

        % Restore previous RNG state
        rng(seed);

    case 'circulant'
        first_row = [1 1 zeros(1, m - 2)]; % First row of the circulant matrix H
        H = toeplitz([first_row(1), fliplr(first_row(2:end))], first_row); 

    otherwise
        error('Unknown method for generating H.');
end







% Matrix Sep Function Parameters

para.rho_outer = 0.5; 
para.rho_inner = 1; 
para.N_outer = 50;
para.N_inner = 50;
para.tol_outer = 1e-7;
para.tol_inner = 1e-7;
para.decomp = "svd"; 
lam =  1/sqrt(m);


% Error matrix to store results
% error_matrix = zeros(max_sparsity, max_rank);
relative_error_matrix = cell(max_sparsity, max_rank, num_trials);


for s = 1:max_sparsity
    disp(s); % To show the iterations remaining
    for r = 1:max_rank
        errors = zeros(1, num_trials);
        
        for trial = 1:num_trials

            % Make L0
            % rank_r = floor(min(m, n)* r/100);
            % L1 = randn([m, rank_r]);  % m x r matrix
            % L1 = orth(L1);
            % L2 = randn([n, rank_r]);  % r x n matrix
            % L2 = orth(L2);
            % % Multiply to get the low-rank matrix (m x n)
            % L0 = L1*L2';
            rank_r = round(min(m, n) * r / 100);  

            U = randn(m, rank_r);                 % Left factor
            V = randn(n, rank_r);                 % Right factor

            L0 = U * V';                          % Resulting m x n matrix of rank r
            

            % Make S0
            mgL = max(abs(L0(:)));          % max absolute value of low-rank matrix
            S0 = zeros(p, n);
            spr = s/100;
            nonzero_idx = randsample(p*n, round(spr * n * p));

            
            switch sparse_method
                case 'gaussian'

                    %S0(nonzero_idx) = rand(s*n, 1) + 1;

                    % Gaussian
                    S0(nonzero_idx) = randn(length(nonzero_idx), 1);  % N(0,1)
                case 'uniform'

                    % Centered Uniform

                    a = mgL;  % or some fraction like 0.5 * mgL
            
                    % Centered uniform sparse values in [-a, a]
                    S0(nonzero_idx) = 2 * a * rand(length(nonzero_idx), 1) - a;
                
                case 'impulsive'

                    % Implusive

                    S0(nonzero_idx) = mgL * sign(randn(length(nonzero_idx), 1));  % Spiky noise

                
            
            end
            %%%%%%%%%%%%%%%%%%%%%

                
           

            
            % Create M0
            M0 = L0 + H * S0;
            
            % Using Precondition as default
            outputC = gen_matrix_sep_con(M0, H, lam, para);

            Lhat = M0 - H*outputC.S;
            rel_L = norm(Lhat - L0, 'fro')/norm(L0, 'fro');
            rel_S = norm(outputC.S - S0, 'fro')/norm(S0, 'fro');
            
            relative_error_matrix{s, r, trial} = [rel_L, rel_S];
            
        
            % % if rel_L < epsilon && rel_S < epsilon
            % %     errors(trial) = 1; % Successful
            % % else
            % %     errors(trial) = 0; % Failed
            % % end
        end
        
      
        % error_matrix(s, r) = mean(errors);
    end
end

% Save error matrix 

file_name = strcat('Relative Error Matrix ',sparse_method, datestr(now));
save(file_name, 'relative_error_matrix');


