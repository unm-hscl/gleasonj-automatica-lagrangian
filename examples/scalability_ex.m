% Prescript running: Initializing srtinit, if it already hasn't been initialized
close all;
clearvars;
srtinit;
srtinit --version;

% example parameters
T = 0.25;

time_horizon = 5;

% probability threshold desired
beta = 0.8;

fprintf('Computing inner and outer approximation for chain of integrators...\n');

MAX_DIM = 12;
under_runtimes = zeros(1, MAX_DIM-1);
over_runtimes = zeros(1, MAX_DIM-1);
for lv = 2:MAX_DIM
    fprintf('    Dimension       : %d\n', lv);

    tm = tic;

    % safe set definition
    safe_set = Polyhedron('lb', -1 * ones(1, lv), 'ub', ones(1, lv));
    % target tube definition
    target_tube = Tube('viability', safe_set, time_horizon);

    % define the system
    sys = getChainOfIntegLtiSystem(lv, ...
        T, ...
        Polyhedron('lb', -0.1, 'ub', 0.1), ...
        RandomVector('Gaussian', zeros(lv,1), 0.001*eye(lv)));

    %inner approximation
    ni_dim = sys.state_dim + sys.input_dim;
    no_dim = sys.state_dim;
    lagunder_options = SReachSetOptions('term', 'lag-under',...
        'bound_set_method', 'ellipsoid', 'compute_style','support',...
        'system', sys, 'vf_enum_method', 'lrs', 'verbose', 2,...
        'n_underapprox_vertices', 2^ni_dim * 6 + 2*ni_dim);

    lagover_options = SReachSetOptions('term', 'lag-over', ...
        'bound_set_method', 'ellipsoid', 'compute_style','support', ...
        'system', sys, 'vf_enum_method', 'lrs', 'verbose', 2, ...
        'n_underapprox_vertices', 2^no_dim * 6 + 2*no_dim);


    setup_time = toc(tm);

    if lv < 6
        tm = tic;

        [polytope_lagunder, extra_info_under] = SReachSet('term', 'lag-under', ...
            sys, beta, target_tube, lagunder_options);

        comp_time_under = toc(tm);
    else
        comp_time_under = NaN;
    end

    tm = tic;

    polytope_lagover = SReachSet('term', 'lag-over', sys,  beta, ...
        target_tube, lagover_options);

    comp_time_over = toc(tm);

    fprintf('    Setup Time      : %.5f\n', setup_time);
    fprintf('        Number of Apporximation Vertices:\n');
    fprintf('            Under: %d\n', 2^ni_dim * 6 + 2*ni_dim);
    fprintf('            Over : %d\n', 2^no_dim * 6 + 2*no_dim);
    fprintf('\n')
    fprintf('    Computation Times:\n');
    fprintf('        Under: %.5f\n', comp_time_under);
    fprintf('        Over : %.5f\n', comp_time_over);
    fprintf('        Total: %.5f\n', comp_time_under + comp_time_over)
    fprintf('--------------------------------------------------------------\n');
    fprintf('\n');

    under_runtimes(lv-1) = comp_time_under;
    over_runtimes(lv-1) = comp_time_over;
end