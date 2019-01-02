% Prescript running: Initializing srtinit, if it already hasn't been initialized
close all;clearvars;srtinit;srtinit --version;

time_horizon = 50;
init_heading = pi/10;
sampling_time = 0.1;
box_halflength = 4;
omega = pi/time_horizon/sampling_time;
turning_rate = omega*ones(time_horizon,1);
umax = 10;
input_space = Polyhedron('lb',0,'ub',umax);
n_mcarlo_sims = 1e3;                        % Monte-Carlo simulation particles
n_mcarlo_sims_affine = 1e3;                 % For affine controllers

%% LTV system definition
% Disturbance matrix and random vector definition
dist_matrix = eye(2);
eta_dist_nongauss = RandomVector('UserDefined', @(N) ...
    0.1*(betarnd(2, 2, 2, N) - 0.5));

[sys_nongauss, heading_vec] = getDubinsCarLtv('add-dist', turning_rate, ...
    init_heading, sampling_time, input_space, dist_matrix, eta_dist_nongauss);
% Compute the mean trajectory of the concatenated disturbance vector
muW_nongauss = sys_nongauss.dist.concat(time_horizon).mean();


% Probability threshold of interest
prob_thresh = 0.8;

v_nominal = umax * 2/3;                 % Nominal trajectory's heading velocity
                                        % TODO: Gaussian was 3/2
% Construct the nominal trajectory
[~,H,~] = sys_nongauss.getConcatMats(time_horizon);
center_box_X = [zeros(2,1);H * (v_nominal * ones(time_horizon,1))];
center_box = reshape(center_box_X,2,[]);
% Box sizes
box_halflength_at_0 = 4;                % Box half-length at t=0
time_const = 1/2*time_horizon;          % Time constant characterize the
                                        % exponentially decaying box half-length

% Target tube definition as well as plotting
target_tube_cell = cell(time_horizon + 1,1); % Vector to store target sets
figure(100);clf;hold on
a= gca;
for itt = 0:time_horizon
    % Define the target set at time itt
    target_tube_cell{itt+1} = Polyhedron(...
        'lb',center_box(:, itt+1) -box_halflength_at_0*exp(- itt/time_const),...
        'ub', center_box(:, itt+1) + box_halflength_at_0*exp(- itt/time_const));
    if itt==0
        % Remember the first the tube
        h_target_tube = plot(target_tube_cell{1},'alpha',0.5,'color','y');
    else
        plot(target_tube_cell{itt+1},'alpha',0.08,'LineStyle',':','color','y');
        a.Children(1).Annotation.LegendInformation.IconDisplayStyle='off';
    end            
end
axis equal        
h_nominal_traj = scatter(center_box(1,:), center_box(2,:), 50,'ks','filled');        
h_vec = [h_target_tube, h_nominal_traj];
legend_cell = {'Target tube', 'Nominal trajectory'};
legend(h_vec, legend_cell, 'Location','EastOutside', 'interpreter','latex');
xlabel('x');
ylabel('y');
axis equal
box on;
grid on;
drawnow;

% Target tube definition
target_tube = Tube(target_tube_cell{:});

%% Lagrangian approach
%% Underapproximation of the stochastic reach set
% 1. We use ellipsoidal set to describe the Gaussian disturbance
% 2. We use support vector-based underapproximation of the one step
%    backward reach set (We use 2^n_dim * 6 + 2*n_dim well-separated
%    vectors; n_dim = sys_nongauss.state_dim + sys_nongauss.input_dim)
% 3. We use LRS library (via GeoCalcLib) for vertex-facet enumeration
n_dim = sys_nongauss.state_dim + sys_nongauss.input_dim;
timer_lagunder_options = tic;
theta_polytope_vec = linspace(0,2*pi,10)';
lagunder_options = SReachSetOptions('term', 'lag-under', ...
    'bound_set_method', 'polytope', 'template_polytope', ...
    Polyhedron('V',[cos(theta_polytope_vec),sin(theta_polytope_vec)]), ...
    'compute_style', 'vfmethod', 'vf_enum_method', 'lrs', 'verbose', 2);
    %'compute_style','support', 'n_particles', 1e3, ...
    %'system', sys_nongauss, 'n_underapprox_vertices', 2^n_dim * 6 + 2*n_dim,...
    

elapsed_time_lagunder_options = toc(timer_lagunder_options);

timer_lagunder = tic;
[polytope_lagunder, extra_info_under] = SReachSet('term', 'lag-under', ...
    sys_nongauss, prob_thresh, target_tube, lagunder_options);
elapsed_time_lagunder = toc(timer_lagunder);

%% Overapproximation of the stochastic reach set
% 1. We use ellipsoidal set to describe the Gaussian disturbance
% 2. We use support vector-based underapproximation of the overapproximation of 
%    the stochastic reach set (We use 2^n_dim * 6 + 2*n_dim well-separated
%    vectors; n_dim = sys_nongauss.state_dim)
% 3. We use LRS library (via GeoCalcLib) for vertex-facet enumeration
n_dim = sys_nongauss.state_dim;
timer_lagover_options = tic;
lagover_options = SReachSetOptions('term', 'lag-over', ...
    'bound_set_method', 'polytope', 'template_polytope', ...
    Polyhedron('V',[cos(theta_polytope_vec),sin(theta_polytope_vec)]), ...
    'compute_style', 'vfmethod', 'vf_enum_method', 'lrs', 'verbose', 2);
elapsed_time_lagover_options = toc(timer_lagover_options);
    
timer_lagover = tic;
polytope_lagover = SReachSet('term', 'lag-over', sys_nongauss,  prob_thresh, ...
    target_tube, lagover_options);
elapsed_time_lagover = toc(timer_lagover);

%% Plot the set
figure(101);
clf;
hold on;
plot(target_tube(1));
% plot(underapproximate_stochastic_reach_avoid_polytope_ccc,'color','m');
plot(polytope_lagover,'color','g', 'alpha', 0.3);
plot(polytope_lagunder,'color','b', 'alpha', 0.6);
axis equal
axis (1.2*[-box_halflength box_halflength -box_halflength box_halflength]);
box on;
legend('Target set at t=0','Stochastic reach set (Outer-approx)', ...
    'Stochastic reach set (Inner-approx)', 'Location','SouthEast');
set(gca,'FontSize',20);
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');

%% Compute a far-away safe initial
cvx_begin quiet
    variable initial_state(sys_nongauss.state_dim, 1)
    minimize ([1 1]*initial_state)
    subject to
        polytope_lagunder.A*initial_state <= polytope_lagunder.b;
        target_tube(1).A*initial_state <= target_tube(1).b;
cvx_end
switch cvx_status
    case 'Solved'
        fprintf('Testing initial state: ');
        disp(initial_state');
        
        % Create a controller based on the underapproximation
        srlcontrol = SReachLagController(sys_nongauss, ... 
            extra_info_under.bounded_dist_set, ...
            extra_info_under.stoch_reach_tube);
        % Generate Monte-Carlo simulations using the srlcontrol and
        % generateMonteCarloSims
        timer_mcarlo = tic;
        [X,U,W] = generateMonteCarloSims(n_mcarlo_sims, sys_nongauss, ...
            initial_state, time_horizon, srlcontrol, [], ...
            lagunder_options.verbose);
        elapsed_time_mcarlo = toc(timer_mcarlo);
        figure(100);
        hold on;
        % Plot the convex hull of the spread of the points
        polytopesFromMonteCarloSims(X, 4, [1,2], {'color','k','alpha',0});
        a = gca;            
        for tindx = 1:time_horizon-1
            a.Children(tindx).Annotation.LegendInformation.IconDisplayStyle='off';
        end
        a.Children(1).Annotation.LegendInformation.IconDisplayStyle='on';
        a.Children(1).DisplayName = 'Trajectory spread at various time steps';           
        % Plot the initial state
        scatter(initial_state(1), initial_state(2), 200, 'ko', 'filled', ...
            'DisplayName','Initial state');
    otherwise        
end


%% Reporting solution times
fprintf('Elapsed time (Online): (lag-under) %1.3f | (lag-over) %1.3f s\n', ...
    elapsed_time_lagunder, elapsed_time_lagover);
fprintf('Elapsed time (Offline): (lag-under) %1.3f | (lag-over) %1.3f s\n', ...
    elapsed_time_lagunder_options, elapsed_time_lagover_options);
fprintf('Monte-Carlo simulation: Success prob: %1.3f | Time: %1.3f s\n', ...
    size(X,2)/n_mcarlo_sims, elapsed_time_mcarlo);

