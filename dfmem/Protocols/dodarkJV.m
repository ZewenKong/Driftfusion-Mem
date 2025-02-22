function darkJVsol = dodarkJV(sol_ini, JV_scan_rate, JV_scan_pnts, mobseti, Vstart, Vend)
% sol_ini   	= an initial solution
% JV_scan_rate
% JVscan_pnts   = no. of points in time mesh

% mobseti       = determines whether ion mobility is on or off 
%                 (0, ion mobility off; 1, ion mobility on)

% Vstart        = scan start voltage
% Vend          = scan end voltage

%% Function Start
tic % Start stopwatch timer
disp('CURRENT-VOLTAGE SCAN IN DARK');

par = sol_ini.par; % read parameters structure into structure par
par.int1 = 0;
par.mobseti = mobseti;

%% Settings
par.tmesh_type = 1;
par.t0 = 0;
par.tmax = abs(Vend - Vstart)/JV_scan_rate;
par.tpoints = JV_scan_pnts;
par.V_fun_type = 'sweep';
par.V_fun_arg(1) = Vstart;
par.V_fun_arg(2) = Vend;
par.V_fun_arg(3) = par.tmax;

disp('>');
disp('> FORWARD SCAN');
darkJVsol.dk.f = df(sol_ini, par);
disp('> FORWARD COMPLETE');
disp('>');
disp('> REVERSE SCAN');
par.V_fun_arg(1) = Vend;
par.V_fun_arg(2) = Vstart;
par.V_fun_arg(3) = par.tmax;
darkJVsol.dk.r = df(darkJVsol.dk.f, par);
disp('> REVERSE COMPLETE');
disp('>');
toc
end