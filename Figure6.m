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
%     Define the target set at time itt
    target_tube_cell{itt+1} = Polyhedron(...
        'lb',center_box(:, itt+1) -box_halflength_at_0*exp(- itt/time_const),...
        'ub', center_box(:, itt+1) + box_halflength_at_0*exp(- itt/time_const));
end

target_tube = Tube(target_tube_cell{:});

% Underapproximation of the stochastic reach set
n_dim = sys.state_dim + sys.input_dim;
timer_lagunder_options = tic;
theta_polytope_vec = linspace(0,2*pi,50)';
% template_poly = Polyhedron('V',[cos(theta_polytope_vec),sin(theta_polytope_vec)]);
template_poly = Polyhedron('V', [-1, 0; 0, -1; 1,0; 0, 1]);
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

% Underapproximation of the stochastic reach set with different template
n_dim = sys.state_dim + sys.input_dim;
timer_lagunder_options = tic;
theta_polytope_vec = linspace(0,2*pi,10)';
template_poly = Polyhedron('lb', [-1, -1], 'ub', [1, 1]);
template_poly.computeHRep();
lagunder_options2 = SReachSetOptions('term', 'lag-under', ...
    'bound_set_method', 'polytope', 'template_polytope', ...
    template_poly, ...
    'compute_style', 'vfmethod', 'vf_enum_method', 'lrs', 'verbose', 2);
    

elapsed_time_lagunder_options = toc(timer_lagunder_options);

timer_lagunder = tic;
[polytope_lagunder2, extra_info_under2] = SReachSet('term', 'lag-under', ...
    sys, prob_thresh, target_tube, lagunder_options2);
elapsed_time_lagunder2 = toc(timer_lagunder);

figure(1)
plot(target_tube(1), 'Color', 'y');
hold on;
plot(polytope_lagunder2, 'Color', 'c');
plot(polytope_lagunder, 'Color', 'g');
hold off;

figure(2)
plot(extra_info_under.bounded_dist_set, 'Color', 'g');
hold on;
plot(extra_info_under2.bounded_dist_set, 'Color', 'c');
hold off;

hf = figure(1);
hf.Units = 'inches';
hf.Position(3) = 3.2;
hf.Position(4) = 3.2;

box on;
xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');

ha = gca;
axis equal;
ha.XLim = [-4.2000, 4.2];
ha.YLim = [-4.2000, 4.2];
% ha.Position = [0.1628, 0.2496, 0.7158, 0.7287];
ha.TickLabelInterpreter = 'latex';

hf = figure(2);
hf.Units = 'inches';
hf.Position(3) = 3.2;
hf.Position(4) = 2.8;

ha = gca;

box on;
xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');
ha.TickLabelInterpreter = 'latex';

axis equal;
