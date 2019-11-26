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

% Target set --- Box [-0.1,0.1]x[-0.1,0]x[-0.01,0.01]x[-0.01,0.01]
target_set = Polyhedron('lb', [-0.1; -0.1; -0.01; -0.01],...
                        'ub', [0.1; 0; 0.01; 0.01]);
target_tube = Tube('reach-avoid',safe_set, target_set, time_horizon);                    
slice_at_vx_vy = zeros(2,1);

%% Preparation for set computation
prob_thresh = 0.8;

% n_dir_vecs = 16;
% theta_vec = linspace(0, 2*pi, n_dir_vecs);
% set_of_dir_vecs_ft = [cos(theta_vec);
%                       sin(theta_vec);
%                       zeros(2,n_dir_vecs)];
% n_dir_vecs = 40;
% theta_vec = linspace(0, 2*pi, n_dir_vecs);
% set_of_dir_vecs_cc_open = [cos(theta_vec);
%                            sin(theta_vec);
%                            zeros(2,n_dir_vecs)];
% 
% %% CC (Linear program approach)
% options = SReachSetOptions('term', 'chance-open', 'verbose', 1, ...
%     'compute_style', 'max_safe_init','set_of_dir_vecs',set_of_dir_vecs_cc_open);
% timer_cc_open = tic;
% [polytope_cc_open, extra_info] = SReachSet('term','chance-open', sys,...
%     prob_thresh, target_tube, options);  
% elapsed_time_cc_open = toc(timer_cc_open);
% 
% %% Fourier transform (Genz's algorithm and MATLAB's patternsearch)
% options = SReachSetOptions('term', 'genzps-open', 'verbose', 1, ...
%     'desired_accuracy', 5e-2, 'set_of_dir_vecs', set_of_dir_vecs_ft);
% timer_ft = tic;
% [polytope_ft, ~] = SReachSet('term','genzps-open', sys, prob_thresh,...
%     target_tube, options);  
% elapsed_time_ft = toc(timer_ft);

% %% Lagrangian underapproximation
vecs_per_orth = [6]%[4, 8, 32, 64];
lag_comptimes = zeros(size(vecs_per_orth));
elapsed_time_cc_open = zeros(size(vecs_per_orth));
elapsed_time_genzps = zeros(size(vecs_per_orth));
lag_polys = [];
n_dim = sys.state_dim + sys.input_dim;
for lv = 1:length(vecs_per_orth)
    %% Lagrangian underapproximation
    fprintf('Lagrangian approximation with %d vectors per orthant\n', ...
        vecs_per_orth(lv));
    fprintf('------------------------------------------------------------\n');
    timer_lagunder_options = tic;
    lagunder_options = SReachSetOptions('term', 'lag-under',...
        'bound_set_method', 'ellipsoid', 'compute_style','support',...
        'system', sys, 'vf_enum_method', 'lrs', 'verbose', 1,...
        'n_vertices', 2^n_dim * vecs_per_orth(lv) + 2*n_dim);
    elapsed_time_lagunder_options = toc(timer_lagunder_options);
%     disp(elapsed_time_lagunder_options)
    
    timer_lagunder = tic;
    [polytope_lagunder, extra_info_under] = SReachSet('term', 'lag-under', ...
        sys, prob_thresh, target_tube, lagunder_options);
    lag_comptimes(lv) = toc(timer_lagunder);
    lag_polys = [lag_polys; polytope_lagunder];

    %% CC (Linear program approach)
%     equi_dir_vecs_over_state=lagunder_options.equi_dir_vecs(1:sys.state_dim,:);
    equi_dir_vecs_over_state=spreadPointsOnUnitSphere(sys.state_dim,...
                        lagunder_options.n_vertices, lagunder_options.verbose);
    ccopen_options = SReachSetOptions('term', 'chance-open', 'verbose', 1, ...
        'compute_style', 'max_safe_init', ...
        'set_of_dir_vecs', equi_dir_vecs_over_state);
    timer_cc_open = tic;
    [polytope_cc_open, extra_info] = SReachSet('term','chance-open', sys,...
        prob_thresh, target_tube, ccopen_options);  
    elapsed_time_cc_open(lv) = toc(timer_cc_open);
    
%     %% Fourier transform (Genz's algorithm and MATLAB's patternsearch)
%     genzps_options = SReachSetOptions('term', 'genzps-open', 'verbose', 1, ...
%         'desired_accuracy', 5e-2, ...
%         'set_of_dir_vecs', equi_dir_vecs_over_state);
%     timer_genzps = tic;
%     [polytope_ft, ~] = SReachSet('term','genzps-open', sys, prob_thresh,...
%         target_tube, genzps_options);  
%     elapsed_time_genzps(lv) = toc(timer_genzps);
end

%% Plotting and Monte-Carlo simulation-based validation
hf = figure();
box on;
hold on;
plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', [0.95, 0.95, 0]);
ha = gca;
ha.Children(1).DisplayName = 'Safe Set';
plot(target_set.slice([3,4], slice_at_vx_vy), 'color', [0, 0, 0]);
ha.Children(1).DisplayName = 'Target Set';
plot(polytope_cc_open.slice([3,4], slice_at_vx_vy), 'color',[1, 0.6, 0],'alpha', 1);
ha.Children(1).DisplayName = 'Chance Constraint';
% plot(polytope_ft.slice([3,4], slice_at_vx_vy), 'color',[0, 0.6, 1],'alpha',1);
% ha.Children(1).DisplayName = 'Fourier Transforms';
cl = [0, 0.5, 0.5];
for lv = length(lag_polys):-1:1
    poly = lag_polys(lv);
    cl = cl + [0, 0.5/4, -0.5/4];
    plot(poly.slice([3,4], slice_at_vx_vy), 'color', cl, 'alpha',1);
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
fprintf('    Chance Constraints (chance-open) %1.3f\n', elapsed_time_cc_open);
for lv = 1:length(vecs_per_orth)
    fprintf('    Lagrangian Approximation with %d Directions %1.3f\n', ... 
        2^n_dim * vecs_per_orth(lv) + 2*n_dim, lag_comptimes(lv)); 
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
