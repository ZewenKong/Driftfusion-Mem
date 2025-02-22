function darkCVsol = dodarkCV(sol_ini, V0, Vmax, Vmin, CV_scan_rate, cycles, tpnts)

% SOL_INI = solution containing intitial conditions
% V0 = Starting voltage (V)
% VMAX = Maximum voltage point (V)
% VMIN = Minimum voltage point (V)
% SCAN_RATE = Scan rate (Vs-1)
% CYCLES = No. of scan cycles
% TPOINTS = No. of points in output time array

tic
disp('> C-V SIMULATION STARTING');
par = sol_ini.par;
sol = sol_ini;

if V0 ~= par.Vapp % if V0 not equal to par.Vapp
    sol = genVappStructs(sol, V0, 0);
end

deltaV = abs(Vmax - V0) + abs(Vmin - Vmax) + abs(V0 - Vmin);
tmax = cycles * (deltaV/CV_scan_rate);

disp('>');
disp('> PERFORMING C-V SIMULATION'); % cyclic voltammogram
darkCVsol = VappFunction(sol, 'tri', [V0, Vmax, Vmin, cycles, tmax/cycles], tmax, tpnts, 0);
disp('> C-V SIMULATION COMPLETE');
disp('>');
toc
end