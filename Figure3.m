%
% Name        : Figure3.m
% Authors     : Joseph D. Gleason and Abraham P. Vinod
% Date        : 2018-10-12
%
% Description : Generate Figure 3 from submitted work; verification of satellite
%               rendezvous-docking problem using Clohessy-Wiltshire-Hill
%               dynamics
% 

close all;
clearvars;

%% System definition
umax = 0.1;
mean_disturbance = zeros(4,1);
covariance_disturbance = diag([1e-4, 1e-4, 5e-8, 5e-8]);
% Define the CWH (planar) dynamics of the deputy spacecraft relative to the
% chief spacecraft as a LtiSystem object
sys = getCwhLtiSystem(4, Polyhedron('lb', -umax*ones(2,1),...
                                    'ub',  umax*ones(2,1)),...
       RandomVector('Gaussian', mean_disturbance,covariance_disturbance));

%% Target tube construction --- reach-avoid specification
time_horizon=5;          % Stay within a line of sight cone for 4 time steps and 
                         % reach the target at t=5% Safe Set --- LoS cone
slice_at_vx_vy = 0*ones(2,1);
% Safe set definition --- LoS cone |x|<=y and y\in[0,ymax] and |vx|<=vxmax and 
% |vy|<=vymax
ymax = 2;
vxmax = 0.5;
vymax = 0.5;
A_safe_set = [1, 1, 0, 0;           
             -1, 1, 0, 0; 
              0, -1, 0, 0;
              0, 0, 1,0;
              0, 0,-1,0;
              0, 0, 0,1;
              0, 0, 0,-1];
b_safe_set = [0;
              0;
              ymax;
              vxmax;
              vxmax;
              vymax;
              vymax];
safe_set = Polyhedron(A_safe_set, b_safe_set);
safe_set_init = safe_set.intersect(Polyhedron('He', ...
    [zeros(2) eye(2) slice_at_vx_vy]));
% Target set --- Box [-0.1,0.1]x[-0.1,0]x[-0.01,0.01]x[-0.01,0.01]
target_set = Polyhedron('lb', [-0.1; -0.1; -0.01; -0.01],...
                        'ub', [0.1; 0; 0.01; 0.01]);
% target_tube = Tube('reach-avoid',safe_set, target_set, time_horizon);                    
target_tube = Tube(safe_set_init, safe_set, safe_set, safe_set, safe_set, ...
    target_set);

%% Preparation for set computation
prob_thresh = 0.8;

vecs_per_orth = [1]%[4, 8, 32, 64];
lag_comptimes = zeros(size(vecs_per_orth));
elapsed_time_cc_open = zeros(size(vecs_per_orth));
elapsed_time_genzps = zeros(size(vecs_per_orth));
lag_polys = [];
n_dim = sys.state_dim + sys.input_dim;

for lv = 1:length(vecs_per_orth)
    % Option for Lagrangian underapproximation
    timer_lagunder_options = tic;
    lagunder_options = SReachSetOptions('term', 'lag-under',...
        'bound_set_method', 'ellipsoid', 'compute_style','support',...
        'system', sys, 'vf_enum_method', 'lrs', 'verbose', 1,...
        'n_vertices', 2^n_dim * vecs_per_orth(lv) + 2*n_dim);
    elapsed_time_lagunder_options = toc(timer_lagunder_options);
    
    % Option for chance-const and genzps
    equi_dir_vecs_over_state = spreadPointsOnUnitSphere(sys.state_dim, ...
                        lagunder_options.n_vertices, lagunder_options.verbose);
    rand_index = randperm(size(equi_dir_vecs_over_state,2));
    equi_dir_vecs_over_state = equi_dir_vecs_over_state(:, ...
        rand_index(1:lagunder_options.n_vertices));
    
    ccopen_options = SReachSetOptions('term', 'chance-open', 'verbose', 1, ...
        'compute_style', 'mve', 'set_of_dir_vecs', equi_dir_vecs_over_state);

    genzps_options = SReachSetOptions('term', 'genzps-open', 'verbose', 1, ...
        'desired_accuracy', 5e-2, ...
        'set_of_dir_vecs', ccopen_options.set_of_dir_vecs);

    %% Lagrangian underapproximation
    fprintf('Lagrangian approximation with %d vectors per orthant\n', ...
        vecs_per_orth(lv));
    fprintf('------------------------------------------------------------\n');
    timer_lagunder = tic;
    [polytope_lagunder, extra_info_under] = SReachSet('term', 'lag-under', ...
        sys, prob_thresh, target_tube, lagunder_options);
    lag_comptimes(lv) = toc(timer_lagunder);
    lag_polys = [lag_polys; polytope_lagunder];

    %% CC (Linear program approach)    
    timer_cc_open = tic;
    [polytope_cc_open, extra_info] = SReachSet('term','chance-open', sys,...
        prob_thresh, target_tube, ccopen_options);  
    elapsed_time_cc_open(lv) = toc(timer_cc_open);
    
    %% Fourier transform (Genz's algorithm and MATLAB's patternsearch)
    timer_genzps = tic;
    [polytope_ft, ~] = SReachSet('term','genzps-open', sys, prob_thresh,...
        target_tube, genzps_options);  
    elapsed_time_genzps(lv) = toc(timer_genzps);
end

%% Plotting and Monte-Carlo simulation-based validation
n_direction_vectors_sv = 60;
hf = figure();
box on;
hold on;
safe_set_2D = safe_set_init.slice([3,4], slice_at_vx_vy);
% safe_set_2D = support_vector_based_slice(safe_set, n_direction_vectors_sv, slice_at_vx_vy);
plot(safe_set_2D, 'color', [0.95, 0.95, 0]);
ha = gca;
ha.Children(1).DisplayName = 'Safe Set';
target_set_2D = target_set.slice([3,4], slice_at_vx_vy);
% target_set_2D = support_vector_based_slice(target_set, n_direction_vectors_sv, slice_at_vx_vy);
plot(target_set_2D, 'color', [0, 0, 0]);
ha.Children(1).DisplayName = 'Target Set';
% polytope_cc_open_2D = polytope_cc_open.slice([3,4], slice_at_vx_vy)
polytope_cc_open_2D = support_vector_based_slice( ...
    Polyhedron('V', polytope_cc_open.V), n_direction_vectors_sv, slice_at_vx_vy);
plot(polytope_cc_open_2D, 'color',[1, 0.6, 0],'alpha', 1);
ha.Children(1).DisplayName = sprintf('Chance Constraint: %d Directions', ...
         2^n_dim * vecs_per_orth(lv) + 2*n_dim);
% polytope_ft_2D = polytope_ft.slice([3,4], slice_at_vx_vy)
polytope_ft_2D = support_vector_based_slice( ...
    Polyhedron('V', polytope_ft.V), n_direction_vectors_sv, slice_at_vx_vy);
plot(polytope_ft_2D, 'color',[0, 0.6, 1],'alpha',1);
ha.Children(1).DisplayName = 'Fourier Transforms';
cl = [0, 0.5, 0.5];
for lv = length(lag_polys):-1:1
    poly = lag_polys(lv);
    cl = cl + [0, 0.5/4, -0.5/4];
%     poly_2D = poly.slice([3,4], slice_at_vx_vy);
    poly_2D = support_vector_based_slice(Polyhedron('V', poly.V), n_direction_vectors_sv, ...
        slice_at_vx_vy);
    plot(poly_2D, 'color', cl, 'alpha',1);
    ha.Children(1).DisplayName = sprintf('Lagrangian Approximation: %d Directions', ...
         2^n_dim * vecs_per_orth(lv) + 2*n_dim);
end
hl = legend();
hl.Location = 'South';
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
box on;
grid on;

% Formatting for paper
ha.FontSize = 9;
ha.TickLabelInterpreter = 'latex';
hl.FontSize = 9;
hl.Interpreter = 'latex';

hf.Units = 'inches';
hf.Position(3) = 3.2;

ha.Position = [0.1497, 0.4677, 0.7650, 0.4892];
hl.Position = [0.0398, 0.0202, 0.9213, 0.3298];
hl.Box = 'off';

print(gcf, '-r200', '-dpng', sprintf('~/Dropbox/cwh-ex-%s.png', datestr(now, 'yyyy-mm-dd-HHMMSS')));

%% Compute time
fprintf('Elapsed time: \n');
fprintf('    Fourier Transforms (genzps-open) %1.3f\n', elapsed_time_genzps);
fprintf('    Chance Constraints (chance-open) with %d Directions %1.3f\n', ...
    size(equi_dir_vecs_over_state,2), elapsed_time_cc_open);
for lv = 1:length(vecs_per_orth)
    fprintf('    Lagrangian Approximation         with %d Directions %1.3f\n', ... 
        lagunder_options.n_vertices, lag_comptimes(lv)); 
end

%% Helper functions
% Plotting function
function [legend_cell] = plotMonteCarlo(method_str, mcarlo_result,...
    concat_state_realization, n_mcarlo_sims, n_sims_to_plot, state_dim,...
    initial_state, legend_cell)
% Plots a selection of Monte-Carlo simulations on top of the plot

    green_legend_updated = 0;
    red_legend_updated = 0;
    traj_indices = floor(n_mcarlo_sims*rand(1,n_sims_to_plot));
    for realization_index = traj_indices
        % Check if the trajectory satisfies the reach-avoid objective
        if mcarlo_result(realization_index)
            % Assign green triangle as the marker
            plotOptions = {'Color', 'g', 'Marker', '^', ...
                'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g' 'MarkerSize', 5};
            markerString = 'g^';
        else
            % Assign red asterisk as the marker
            plotOptions = {'Color', 'r', 'Marker', 'x', ...
                'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r' 'MarkerSize', 5};
            markerString = 'rx';
        end

        % Create [x(t_1) x(t_2)... x(t_N)]
        reshaped_X_vector = reshape(...
            concat_state_realization(:,realization_index), state_dim,[]);

        % This realization is to be plotted
        h = plot([initial_state(1), reshaped_X_vector(1,:)], ...
                 [initial_state(2), reshaped_X_vector(2,:)], ...
                 plotOptions{:});
    end
end
