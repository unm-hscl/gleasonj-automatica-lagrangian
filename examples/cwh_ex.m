close all;
clearvars;
cd('../SReachTools-private');
srtinit;
srtinit --version;
cd('../examples');

% Monte-Carlo simulation parameters
n_mcarlo_sims = 1e3;

%% System definition
umax = 0.1;
mean_disturbance = zeros(4,1);
covariance_disturbance = diag([1e-4, 1e-4, 5e-8, 5e-8]);
% Define the CWH (planar) dynamics of the deputy spacecraft relative to the
% chief spacecraft as a LtiSystem object
sys = getCwhLtiSystem(4, Polyhedron('lb', -umax*ones(2,1), ...
                                    'ub',  umax*ones(2,1)), ...
       RandomVector('Gaussian', mean_disturbance,covariance_disturbance));


%% Target tube construction --- reach-avoid specification
time_horizon = 5;          % Stay within a line of sight cone for 4 time steps and 
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
target_set = Polyhedron('lb', [-0.1; -0.1; -0.01; -0.01], ...
                        'ub', [0.1; 0; 0.01; 0.01]);
target_tube = Tube('reach-avoid',safe_set, target_set, time_horizon);                    
slice_at_vx_vy = zeros(2,1);
%%

%% Preparation for set computation
prob_thresh = 0.8;

n_dir_vecs = 10;
theta_vec = linspace(0, 2*pi, n_dir_vecs);
set_of_dir_vecs_ft = [cos(theta_vec);
                      sin(theta_vec);
                      zeros(2,n_dir_vecs)];
n_dir_vecs = 40;
theta_vec = linspace(0, 2*pi, n_dir_vecs);
set_of_dir_vecs_cc_open = [cos(theta_vec);
                           sin(theta_vec);
                           zeros(2,n_dir_vecs)];
init_safe_set_affine = Polyhedron('He',[zeros(2,2) eye(2,2) slice_at_vx_vy]);

%% Lagrangian approach
%% Underapproximation of the stochastic reach set
% 1. We use ellipsoidal set to describe the Gaussian disturbance
% 2. We use support vector-based underapproximation of the one step
%    backward reach set (We use 2^n_dim * 6 + 2*n_dim well-separated
%    vectors; n_dim = sys.state_dim + sys.input_dim)
% 3. We use LRS library (via GeoCalcLib) for vertex-facet enumeration
n_dim = sys.state_dim + sys.input_dim;
timer_lagunder_options = tic;
lagunder_options = SReachSetOptions('term', 'lag-under',...
    'bound_set_method', 'ellipsoid', 'compute_style','support',...
    'system', sys, 'vf_enum_method', 'lrs', 'verbose', 2,...
    'n_vertices', 2^n_dim * 6 + 2*n_dim);
elapsed_time_lagunder_options = toc(timer_lagunder_options);

timer_lagunder = tic;
[polytope_lagunder, extra_info_under] = SReachSet('term', 'lag-under', ...
    sys, prob_thresh, target_tube, lagunder_options);
elapsed_time_lagunder = toc(timer_lagunder);

%% Overapproximation of the stochastic reach set
% 1. We use ellipsoidal set to describe the Gaussian disturbance
% 2. We use support vector-based underapproximation of the overapproximation of 
%    the stochastic reach set (We use 2^n_dim * 6 + 2*n_dim well-separated
%    vectors; n_dim = sys.state_dim)
% 3. We use LRS library (via GeoCalcLib) for vertex-facet enumeration
n_dim = sys.state_dim;
timer_lagover_options = tic;
lagover_options = SReachSetOptions('term', 'lag-over', ...
        'bound_set_method', 'ellipsoid', 'compute_style','support', ...
        'system', sys, 'vf_enum_method', 'lrs', 'verbose', 1, ...
        'n_vertices', 2^n_dim * 6 + 2*n_dim);
elapsed_time_lagover_options = toc(timer_lagover_options);
    
timer_lagover = tic;
polytope_lagover = SReachSet('term', 'lag-over', sys,  prob_thresh, ...
    target_tube, lagover_options);
elapsed_time_lagover = toc(timer_lagover);

%% Plot the sets
figure(101);
clf
box on;
hold on;
plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y');
plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'g');
legend_cell = {'Safe set','Target set'};
plot(Polyhedron('V',polytope_lagunder.V(:,1:2)), 'color','r','alpha',0.5);    
legend_cell{end+1} = 'Underapprox. polytope (lag-under)';
plot(Polyhedron('V',polytope_lagover.V(:,1:2)), 'color','m','alpha',0.5);    
legend_cell{end+1} = 'Overapprox. polytope (lag-over)';

%% Compute a far-away safe initial
cvx_begin quiet
    variable initial_state(sys.state_dim, 1)
    minimize ([-1 1 0 0]*initial_state)
    subject to
        polytope_lagunder.A*initial_state <= polytope_lagunder.b;
        target_tube(1).A*initial_state <= target_tube(1).b;
cvx_end
switch cvx_status
    case 'Solved'
        fprintf('Testing initial state: ');
        disp(initial_state');
        
        % Create a controller based on the underapproximation
        srlcontrol = SReachLagController(sys, ... 
            extra_info_under.bounded_dist_set, ...
            extra_info_under.stoch_reach_tube);
        % Generate Monte-Carlo simulations using the srlcontrol and
        % generateMonteCarloSims
        timer_mcarlo = tic;
        [X,U,W] = generateMonteCarloSims(n_mcarlo_sims, sys, initial_state, ...
            time_horizon, srlcontrol);
        elapsed_time_mcarlo = toc(timer_mcarlo);
        % Plot the convex hull of the spread of the points
        polytopesFromMonteCarloSims(X, 4, [1,2], {'color','k','alpha',0});
        a = gca;            
        for tindx = 1:time_horizon - 1
            a.Children(tindx).Annotation.LegendInformation.IconDisplayStyle='off';
        end
        legend_cell{end+1} = 'Trajectory spread at various time steps';           
        % Plot the initial state
        scatter(initial_state(1), initial_state(2), 200, 'ko', 'filled', ...
            'DisplayName','Initial state');
    otherwise
        
end
legend(legend_cell, 'Location','South');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');

%% Reporting solution times
fprintf('Elapsed time (Online): (lag-under) %1.3f | (lag-over) %1.3f s\n', ...
    elapsed_time_lagunder, elapsed_time_lagover);
fprintf('Elapsed time (Offline): (lag-under) %1.3f | (lag-over) %1.3f s\n', ...
    elapsed_time_lagunder_options, elapsed_time_lagover_options);
fprintf('Monte-Carlo simulation: Success prob: %1.3f | Time: %1.3f s\n', ...
    size(X,2)/n_mcarlo_sims, elapsed_time_mcarlo);

