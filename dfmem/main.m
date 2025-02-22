% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- SIMULATION BEGIN --------------

% prepare solutions and indirectly test equilibrate and genIntStructs func
initialise_df

% input = 'Input_files/pedotpss_mapi_pcbm.csv';
input = 'Input_files/1_layer_single_carrier.csv';

% build a parameters object
par = pc(input);
par.prob_distro_function = 'Boltz';
soleq = equilibrate(par);

%% Current-Voltage
scan_rt = 5e-2; % use for both J-V and C-V
pnts = 500; % use for both J-V and C-V

% Simulation in dark (J-V)
mobseti = true; Vstart = 0; Vend = 1.5;
JVsol = dodarkJV(soleq.ion, scan_rt, pnts, mobseti, Vstart, Vend);
% Plot
% dfmem.JVmem(JVsol);

% Simulation in dark (C-V)
V0 = 0; Vmax = 1.5; Vmin = -1.5; cycles = 1;
% CVsol = dodarkCV(soleq.ion, V0, Vmax, Vmin, scan_rt, cycles, pnts);
% Plot
% dfmem.CVmem(CVsol, 0);

%% Energy lvl v.s. position

dfmem.ELJVfr(JVsol); % J-V (at a special time point)
% figure(200);
% h1 = dfmem.ELJVff(JVsol);
% hold on
% h2 = dfmem.ELJVrr(JVsol);
% hold off

% dfmem.ELCV(CVsol, tp); % C-V

%% Charge carriers density v.s. position

% dfmem.npxJVfr(JVsol); % J-V; electrons, holes
% dfmem.acxJVfr(JVsol); % J-V; cations, anions

% dfplot.npx(CVsol, tp); % C-V; electrons, holes
% dfplot.acx(CVsol, tp); % C-V; cations, anions

%% Impedance Spectroscopy

%% Chronoamperometry


