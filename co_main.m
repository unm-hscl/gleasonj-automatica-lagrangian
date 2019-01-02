function co_main()
% CodeOcean main executable file for running Lagrangian method examples for 
% stochastic reachbility of a target tube from submitted work:
% 
% J. D. Gleason, A. P. Vinod, M. M. K. Oishi, "Lagrangian Approximations for 
% Stochastic Reachability of a Target Tube," submitted to _Automatica_, 2019.
% 


fprintf(['CodeOcean Module for Lagrangian methods for the stochastic \n', ...
         'reachbility of a target tube. This code will run several \n', ...
         'examples that generate the figures in the submitted paper:\n\n']);
fprintf(['    J. D. Gleason, A. P. Vinod, M. M. K. Oishi, "Lagrangian \n', ...
         '    Approximations for Stochastic Reachability of a Target \n', ...
         '    Tube," submitted to _Automatica_, 2019.\n\n']);

% choosing which examples to run
%   - double integrator
%   - scalability 
%   - CWH (Spacecraft Rendezvous-docking)
run_dblint = 1;

% For the code ocean module we do not automatically perform the scalability
% example because it is very time-consuming. If you would like to perform this
% example please set the run flag to '1'; we suggest that this is only done
% if you download the full module for offline use or if you have a rather large
% amount of time available for running code with CodeOcean
run_scalability = 0;

run_cwh = 1;

if run_dblint
    run('examples/dblint_ex.m');
end

if run_scalability
    run('examples/scalability_ex.m');
end

if run_cwh
    run('examples/cwh_ex.m');
end

end