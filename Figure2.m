% 
% Filename    : Figure3.m
% Author      : Joseph Gleason
% Description : Script to generate Figure 3 of Automatica paper "Lagrangian
%               Approximations for Stochastic Reachability of a Target Tube."
%
%	        Dubins Car Fiure
%

% Prescript running: Initializing srtinit, if it already hasn't been initialized
clearvars;
close all;

% Time horizon
time_horizon = 21;
init_heading = pi/10;
sampling_time = 0.1;
box_halflength = 4;
omega = pi/time_horizon/sampling_time;
turning_rate = omega*ones(time_horizon,1);
umax = 10;
input_space = Polyhedron('lb',0,'ub',umax);
n_mcarlo_sims = 1e3;                        % Monte-Carlo simulation particles
n_mcarlo_sims_affine = 1e3;                 % For affine controllers

% LTV system definition
% Disturbance matrix and random vector definition
dist_matrix = eye(2);
eta_dist_nongauss = RandomVector('UserDefined', @(N) ...
    0.1*(betarnd(2, 2, 2, N) - 0.5));

[sys, heading_vec] = getDubinsCarLtv('add-dist', turning_rate, ...
    init_heading, sampling_time, input_space, dist_matrix, eta_dist_nongauss);

% Compute the mean trajectory of the concatenated disturbance vector
muW_nongauss = sys.dist.concat(time_horizon).mean();
mu = reshape(muW_nongauss, 2, []);



% Probability threshold of interest
prob_thresh = 0.8;

v_nominal = umax * 2/3;                 % Nominal trajectory's heading velocity
                                        % TODO: Gaussian was 3/2
% Construct the nominal trajectory
[~,H,~] = sys.getConcatMats(time_horizon);
center_box_X = [zeros(2,1);H * (v_nominal * ones(time_horizon,1))];
center_box = reshape(center_box_X,2,[]);
% Box sizes
box_halflength_at_0 = 4;                % Box half-length at t=0
time_const = 1/2*time_horizon;          % Time constant characterize the
                                        % exponentially decaying box half-length

% Target tube definition as well as plotting
target_tube_cell = cell(time_horizon + 1,1);
for itt = 0:time_horizon
    % Define the target set at time itt
    target_tube_cell{itt+1} = Polyhedron(...
        'lb',center_box(:, itt+1) -box_halflength_at_0*exp(- itt/time_const),...
        'ub', center_box(:, itt+1) + box_halflength_at_0*exp(- itt/time_const));
end

target_tube = Tube(target_tube_cell{:});

% Underapproximation of the stochastic reach set
n_dim = sys.state_dim + sys.input_dim;
timer_lagunder_options = tic;
theta_polytope_vec = linspace(0,2*pi,10)';
template_poly = Polyhedron('V',[cos(theta_polytope_vec),sin(theta_polytope_vec)]);
template_poly.computeHRep();
lagunder_options = SReachSetOptions('term', 'lag-under', ...
    'bound_set_method', 'polytope', 'template_polytope', ...
    template_poly, ...
    'compute_style', 'vfmethod', 'vf_enum_method', 'lrs', 'verbose', 2);
    

elapsed_time_lagunder_options = toc(timer_lagunder_options);

timer_lagunder = tic;
[polytope_lagunder, extra_info_under] = SReachSet('term', 'lag-under', ...
    sys, prob_thresh, target_tube, lagunder_options);
elapsed_time_lagunder = toc(timer_lagunder);

% Overapproximation of the stochastic reach set
n_dim = sys.state_dim;
timer_lagover_options = tic;
lagover_options = SReachSetOptions('term', 'lag-over', ...
    'bound_set_method', 'polytope', 'template_polytope', ...
    template_poly, ...
    'compute_style', 'vfmethod', 'vf_enum_method', 'lrs', 'verbose', 2);
elapsed_time_lagover_options = toc(timer_lagover_options);
    
timer_lagover = tic;
polytope_lagover = SReachSet('term', 'lag-over', sys,  prob_thresh, ...
    target_tube, lagover_options);
elapsed_time_lagover = toc(timer_lagover);

% Compute a far-away safe initial
cvx_begin quiet
    variable initial_state(sys.state_dim, 1)
    minimize ([1 1]*initial_state)
    subject to
        polytope_lagunder.A*initial_state <= polytope_lagunder.b;
        target_tube(1).A*initial_state <= target_tube(1).b;
cvx_end

srp_timer = tic;
point_opts = SReachPointOptions('term', 'particle-open', ...
    'n_particles', 100, 'verbose', 1, 'max_particles', 200);

[p, u, K] = SReachPoint('term', 'particle-open', sys, initial_state, ...
    target_tube, point_opts);
srp_time = toc(srp_timer);

% Create a controller based on the underapproximation
srp_cont_timer = tic;
srlcontrol = SReachLagController(sys, ... 
    extra_info_under.bounded_dist_set, ...
    extra_info_under.stoch_reach_tube);
srl_time = toc(srp_cont_timer);

mX_particle   = zeros(2, time_horizon+1);
mX_lagrangian = zeros(2, time_horizon+1);
mX_particle(:, 1) = initial_state;
mX_lagrangian(:, 1) = initial_state;
for lv = 1:length(u)
    mX_particle(:, lv+1) = sys.state_mat(lv-1) * mX_particle(:, lv) + ...
        sys.input_mat(lv-1) * u(lv) + sys.dist_mat(lv-1) * mu(:, lv);
    mX_lagrangian(:, lv+1) = sys.state_mat(lv-1) * mX_lagrangian(:, lv) + ...
        sys.input_mat(lv-1) * srlcontrol.getInput(mX_lagrangian(:, lv), lv-1) + ...
        sys.dist_mat(lv-1) * mu(:, lv);
end

%% Reporting solution times
fprintf('Elapsed time (Lagrangian Underapproximation)      : %1.3f s\n', ...
    elapsed_time_lagunder);
fprintf('Elapsed time (Lagrangian Overapproximation)       : %1.3f s\n', ...
    elapsed_time_lagover);
% fprintf('Elapsed time (Offline): (lag-under) %1.3f | (lag-over) %1.3f s\n', ...
%     elapsed_time_lagunder_options, elapsed_time_lagover_options);
fprintf('Elapsed time (SReachPoint Particle Approximation) : %1.3f s\n', ...
    srp_time);

% Plotting
hf = figure();
hold on;
for lt = length(target_tube):-1:1
    plot(target_tube(lt), 'Color', 'y');
end

for lv = 2:length(hf.Children.Children)
    ptch = hf.Children.Children(lv);
    ptch.FaceAlpha = 0.2;
    ptch.EdgeAlpha = 0.2;
end

plot(polytope_lagover, 'color', 'r');
plot(polytope_lagunder, 'color', 'g');

xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');

plot(mX_particle(1, :), mX_particle(2, :), '-o', ...
    'Color', 'm', ...
    'DisplayName', 'Particle Approximation Mean Trajectory');
plot(mX_lagrangian(1, :), mX_lagrangian(2, :), '-^', ...
    'Color', 'k', ...
    'DisplayName', 'Lagrangian Mean Trajectory');

hf.Units = 'inches';
hf.Position = [hf.Position(1:2), 3.1979, 4.7917];

box on;

ha = gca;
axis equal;
ha.XLim = [-4.2000, 5.5000];
ha.YLim = [-4.2000, 10];
ha.Position = [0.1628, 0.2496, 0.7158, 0.7287];
ha.TickLabelInterpreter = 'latex';

ha.Children(1).MarkerFaceColor = ha.Children(1).Color;
ha.Children(2).MarkerFaceColor = ha.Children(2).Color;

ha.Children(3).DisplayName = 'Lagrangian Underapproximation';
ha.Children(4).DisplayName = 'Lagrangian Overapproximation';
ha.Children(5).DisplayName = 'Target Tube';

hl = legend(ha.Children(1:5));
hl.Interpreter = 'latex';
hl.Box = 'off';
hl.Position = [0.0977, 0.0531, 0.7983, 0.1247];


print(hf, '-r200', '-dpng', sprintf('var/figs/dubincs-car-%s.png', ...
    datestr(now, 'yyyy-mm-dd-HHMMSS')));

% fprintf('Monte-Carlo simulation: Success prob: %1.3f | Time: %1.3f s\n', ...
%     size(X,2)/n_mcarlo_sims, mc_time);

% savefig(hf, sprintf('~/Dropbox/dubins-automatica-%s.fig', datestr(now, 'yyyy-mm-dd-HHMMSS')));
% save(sprintf('~/Dropbox/dubins-automatica-%s.mat', datestr(now, 'yyyy-mm-dd-HHMMSS')));
