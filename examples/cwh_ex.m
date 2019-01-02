


close all;
clearvars;
cd('../SReachTools-private');
srtinit;
srtinit --version;
cd('../examples');


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
n_dim = sys.state_dim + sys.input_dim;
lagunder_options = SReachSetOptions('term', 'lag-under',...
    'bound_set_method', 'ellipsoid', 'compute_style','support',...
    'system', sys, 'verbose', 2,...
    'n_underapprox_vertices', 2^n_dim * 6 + 2*n_dim);

timer_lagunder = tic;
[polytope_lagunder, extra_info_under] = SReachSet('term', 'lag-under', ...
    sys, prob_thresh, target_tube, lagunder_options);
elapsed_time_lagunder = toc(timer_lagunder);

n_dim = sys.state_dim;
lagover_options = SReachSetOptions('term', 'lag-over', ...
    'bound_set_method', 'ellipsoid', 'compute_style','support', ...
    'system', sys, 'verbose', 1, ...
    'n_underapprox_vertices', 2^n_dim * 6 + 2*n_dim);

timer_lagover = tic;
polytope_lagover = SReachSet('term', 'lag-over', sys,  prob_thresh, ...
    target_tube, lagover_options);
elapsed_time_lagover = toc(timer_lagover);
