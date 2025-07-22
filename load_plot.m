
% Change the file name to be the saved error matrix from generate_heatmap.m
file_name = 'Relative Error Matrixuniform22-Jul-2025 18:07:22.mat';
load(file_name, 'relative_error_matrix');

% Change new_epsilon to make plots at different thresholds
new_epsilon = 10e-3;

% These should be the same from generate_heatmap.m
max_rank = 30; % Maximum rank
max_sparsity = 30; % Maximum sparsity
num_trials = 10; % Number of trials for each sparsity-rank pair


new_error_matrix = zeros(max_sparsity, max_rank);

for s = 1:max_sparsity
    for r = 1:max_rank
        success_count = 0;
        for trial = 1:num_trials
            rel_errors = relative_error_matrix{s, r, trial};
            rel_L = rel_errors(1);
            rel_S = rel_errors(2);
            if rel_L <= new_epsilon && rel_S <= new_epsilon
                success_count = success_count + 1;
            end
        end
        new_error_matrix(s, r) = success_count / num_trials;
    end
end

% Save the new matrix if testing on different threshold
save(['Error Matrix eps=' num2str(new_epsilon) '.mat'], 'new_error_matrix');

%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
figure;
imagesc(new_error_matrix);
colorbar;
set(gca, 'YDir', 'normal');

% Axis labels and title with larger font
% Change title latex for different new_epsilon values
xlabel('Rank', 'FontSize', 64);
ylabel('Sparsity Ratio (%)', 'FontSize', 64);
title('Recoverability Success Rate for ${\mathrm{RelErr}} \leq 10^{-3}$', ...
    'Interpreter', 'latex', 'FontSize', 70);

% Tick marks
% Change title latex for different new_epsilon values
yticks = 5:5:max_sparsity;
xticks = 5:5:max_rank;
set(gca, 'YTick', yticks, 'YTickLabel', yticks, ...
         'XTick', xticks, 'XTickLabel', xticks, ...
         'FontSize', 14); % Tick label font size

% Export high-resolution image
%exportgraphics(gcf, 'figures/1D Experiments Circulant H/gaussian-4.png', 'Resolution', 1200);
exportgraphics(gcf, 'figures/1D Experiments Random H/impulsive-3.pdf');

