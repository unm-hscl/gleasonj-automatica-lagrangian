clear;
clc;
close all;
cov_dist = 1.6;
input_max = 0.5;
input_space = input_max * Polyhedron('lb', -ones(2, 1), 'ub', ones(2, 1));
gaussian_dist = RandomVector('Gaussian', zeros(2,1), cov_dist * eye(2));
sys = LtiSystem('StateMatrix', 0.8 * eye(2), ...
                'InputMatrix', eye(2), ...                
                'InputSpace', input_space, ...
                'Disturbance', gaussian_dist, ...
                'DisturbanceMatrix', eye(2));
            
prob_thresh = 0.8;            
time_horizon = 2;
safe_set = 2 * Polyhedron('lb', [-1,-1], 'ub', [1,1]);
safety_tube = Tube('viability', safe_set, time_horizon);

% %% Dynamic programming computation
% % Grid-based dynamic programming
% disp('Grid-based dynamic programming');
% dyn_prog_xinc = 0.05;
% dyn_prog_uinc = 0.15;
% timer_DP=tic;
% [prob_x, cell_of_xvecs] = SReachDynProg('term', sys, dyn_prog_xinc, ...
%     dyn_prog_uinc, safety_tube);
% elapsed_time_DP_recursion = toc(timer_DP);
% x = cell_of_xvecs{1};
% % Set computation
% disp('Compute stochastic reach sets from the optimal value functions');
% timer_DP_set = tic;
% stoch_reach_set_dp = getDynProgLevelSets2D(cell_of_xvecs, prob_x, ...
%     prob_thresh, safety_tube);
% elapsed_time_DP_set = toc(timer_DP_set);
% elapsed_time_DP_total = elapsed_time_DP_recursion + elapsed_time_DP_set;    
% len_grid = length(cell_of_xvecs{1});
% figure(1);
% surf(cell_of_xvecs{2}, cell_of_xvecs{1}, reshape(prob_x, len_grid, []))
    

% %% Chance-open
% set_of_dir_vecs = spreadPointsOnUnitSphere(2, 8);
% options_cco = SReachSetOptions('term', 'chance-open', 'set_of_dir_vecs', ...
%     set_of_dir_vecs, 'compute_style', 'cheby', 'verbose', 1);
% stoch_reach_set_cco = SReachSet('term', 'chance-open', sys, ...
%         prob_thresh, safety_tube, options_cco);
    
%% Lagrangian-under
% template_polytope = Polyhedron('V', spreadPointsOnUnitSphere(2, 16)');
template_polytope = safe_set;
options_lag_under = SReachSetOptions('term', 'lag-under', 'verbose', 1, ...
    'bound_set_method', 'polytope', 'template_polytope', template_polytope, ...
    'compute_style', 'vfmethod', 'vf_enum_method', 'lrs');
%     'bound_set_method', 'ellipsoid', ...                           % Smaller
%     'compute_style', 'support', 'system', sys, 'n_vertices', 50);  % Slower
try
    stoch_reach_set_lag_under = SReachSet('term', 'lag-under', sys, ...
            prob_thresh, safety_tube, options_lag_under);
catch
    fprintf('\n\nEmpty Lagrangian set\n\n');
    stoch_reach_set_lag_under = Polyhedron();
end

% %% Lagrangian-over    
% options_lag_over = SReachSetOptions('term', 'lag-over', ...
%     'bound_set_method', 'polytope', 'template_polytope', safe_set, ...
%     'compute_style', 'vfmethod', 'vf_enum_method', 'lrs');
% %     'compute_style', 'support', 'system', sys, 'n_vertices', 100); % Slower
% %     'bound_set_method', 'ellipsoid', ...                           % Smaller
% stoch_reach_set_lag_over = SReachSet('term', 'lag-over', sys, ...
%         prob_thresh, safety_tube, options_lag_over);

%% Plotting
figure(2);    
% clf;    
plot(safe_set, 'color', 'w','alpha',0);
hold on;
% plot(stoch_reach_set_lag_over, 'color', 'm', 'alpha', 0.3);
% plot(stoch_reach_set_dp(1), 'color', 'c' , 'alpha', 0.5);
% plot(stoch_reach_set_cco, 'color', 'r');
plot(stoch_reach_set_lag_under, 'color', 'g');
box on;
grid on;
axis square;

%% Problem transfer for FAUST^2
A = sys.state_mat;
B = sys.input_mat;
Sigma = gaussian_dist.cov;