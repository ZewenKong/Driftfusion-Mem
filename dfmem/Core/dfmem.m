classdef dfmem
    %% LICENSE
    % Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
    % Imperial College London
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU Affero General Public License as published
    % by the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    %
    methods(Static)
        %% (J-V simulation) Plot
        function JVmem(sol)       
            % Variable define
            J_f = dfana.calcJ(sol.dk.f);
            J_fl = J_f.tot(:, end); % the last coloumn (slice) of total forward current
            V_f = dfana.calcVapp(sol.dk.f); % forward sweep voltage

            J_r = dfana.calcJ(sol.dk.r);
            J_rl = J_r.tot(:, end);
            J_rlf = flipud(J_rl); % flip the last coloumn of total reverse current
            V_r = dfana.calcVapp(sol.dk.r);

            J_diff = J_fl - J_rlf; % curve of r/f difference
            
            [J_maxDiff, idx] = max(abs(J_diff)); % find the max abs difference value and the index to find the corr. voltage value
            V_corr = V_f(idx);
            J_fcorr = J_fl(idx);
            J_rcorr = J_rlf(idx);
            
            % Graph plot
            figure(001);
            hold on
            plot([V_corr, V_corr], [J_fcorr, J_rcorr], 'k-'); % max diff plot
            bar_length = 0.015; % Adjust the width of the top and bottom horizontal bars
            plot([V_corr - bar_length, V_corr + bar_length], [J_fcorr, J_fcorr], 'k-');
            plot([V_corr - bar_length, V_corr + bar_length], [J_rcorr, J_rcorr], 'k-');
            plot(V_f, J_f.tot(:,end), '--', V_r, J_r.tot(:,end)); % f/r J-V plot
            xlabel('Applied voltage [V]')
            ylabel('Current density [A cm^{-2}]');
            hold off
        end
        
        %% (J-V simulation) Max hysteresis
        function idx = JVcalcIdx(sol)
            % Find the index where has the largest hystersis
            % J_f, J_r; matrix [no.of.pnts x 220]
            % J_fl, J_rl, J_rlf; matrix [no.of.pnts x 1]
            % V_r; matrix [1 x no.of.pnts]
            % V_r_corr; matrix [1 x 1]; a value
            % - - - - - - - - - -
            J_f = dfana.calcJ(sol.dk.f);
            J_fl = J_f.tot(:, end); % choose the last coloumn of the J
            
            J_r = dfana.calcJ(sol.dk.r);
            J_rl = J_r.tot(:, end); % choose the last coloumn of the J
            J_rlf = flipud(J_rl); % flip the coloumn/array
            
            V_f = dfana.calcVapp(sol.dk.f); % calculate the voltage, both f/r can be used
            V_r = dfana.calcVapp(sol.dk.r);

            J_diff = J_fl - J_rlf; % the difference curve between f/r
            [J_maxDiff, J_maxDiff_idx] = max(abs(J_diff)); % find the max abs difference value and the index
            V_corr = V_f(J_maxDiff_idx); % the corresponding V point (time point)
            idx = J_maxDiff_idx; % make the variable name shorter
        end
        
        %% (J-V simulation) Energy lvl diagram at special time point
        function ELJV(sol, dirc)
            idx = dfmem.JVcalcIdx(sol);
            if dirc == "f"
                sol = sol.dk.f;
                tp = sol.t(idx);
                dfmem.ELxJVff(sol, tp);
            elseif dirc == "r"
                sol = sol.dk.r;
                tp = sol.t(length(sol.t) - idx + 1);
                dfmem.ELxJVrr(sol, tp);
            end
            % dfplot.ELx(sol, tp);
        end

        % - - - - - - - - - - - - - - - - - - - -
        function ELxJVf(varargin)
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.calcEnergies(sol);
    
            % figure(22);
            dfplot.x2d(sol, x, {Efn, Efp, Ecb, Evb}, {'E_{fn_f}', 'E_{fp_f}', 'E_{CB_f}', 'E_{VB_f}'},...
                {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0)
        end

        function ELxJVr(varargin)
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.calcEnergies(sol);
    
            % figure(22);
            dfplot.x2d(sol, x, {Efn, Efp, Ecb, Evb}, {'E_{fn_r}', 'E_{fp_r}', 'E_{CB_r}', 'E_{VB_r}'},...
                {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0)
        end
        % - - - - - - - - - - - - - - - - - - - -

        % - - - - - - - - - - - - - - - - - - - -
        function ELJVfr(sol)
            figure(002);
            dfmem.ELJV(sol, "f");
            figure(003);
            dfmem.ELJV(sol, "r");
        end

        function ELJVf(sol)
            dfmem.ELJV(sol, "f");
        end

        function ELJVr(sol)
            dfmem.ELJV(sol, "r");
        end
        % - - - - - - - - - - - - - - - - - - - -

        %% (J-V simulation) e-/h+ density at special time point
        function npxJV(sol, dirc)
            idx = dfmem.JVcalcIdx(sol);
            if dirc == "f"
                sol = sol.dk.f;
                tp = sol.t(idx);
            elseif dirc == "r"
                sol = sol.dk.r;
                tp = sol.t(length(sol.t) - idx + 1);
            end
            dfplot.npx(sol, tp);
        end

        function npxJVfr(sol)
            figure(101);
            dfmem.npxJV(sol, "f");
            figure(102);
            dfmem.npxJV(sol, "r");
        end
        
        %% (J-V simulation) c/a density at special time point
        function acxJV(sol, dirc)
            idx = dfmem.JVcalcIdx(sol);
            if dirc == "f"
                sol = sol.dk.f;
                tp = sol.t(idx);
            elseif dirc == "r"
                sol = sol.dk.r;
                tp = sol.t(length(sol.t) - idx + 1);
            end
            dfplot.acx(sol, tp);
        end

        function acxJVfr(sol)
            figure(103);
            dfmem.acxJV(sol, "f");
            figure(104);
            dfmem.acxJV(sol, "r");
        end

        %% (C-V simulation) Plot
        function CVmem(sol, xpos)
            xmesh = sol.x;
            ppos =getpointpos(xpos, xmesh);

            J = dfana.calcJ(sol);
            Vapp = dfana.calcVapp(sol);
            t = sol.t;

            hold on
            figure(004);
            % Main plot
            plot(Vapp, abs(J.tot(:, ppos)));
            xlabel('Applied Voltage, Vapp [V]');
            ylabel('Current Density, J [A cm^{-2}]');
            % Sub plot
            ax1 = axes('Position', [0.75, 0.2, 0.15, 0.15]);
            set(ax1, 'Box', 'on');
            plot(t, Vapp);
            xlabel('t');
            ylabel('Vapp');
            
            hold off
        end
        
        %% (C-V simulation) Time array
        function tp = CVtarr(v, vmax, scan_rt)
            if v >= 0
                tp = v/scan_rt;
            elseif v < 0
                tp = ((2*vmax)+abs(v))/scan_rt;
            end
        end
        
        %% (C-V simulation) Energy lvl diagram
        function ELCV(sol, tp)
            figure(100);
            dfplot.ELx(sol, tp);
        end
        

        
        % NOT IN USE
        % function EL_JVf(sol)
        % 
        %     % - - - - - - - - - -
        %     % idx, the time idx where has the largest hystersis
        %     % - - - - - - - - - -
        %     idx = dfmem.JVcalcJDiff(sol);
        % 
        %     % - - - - - - - - - -
        %     % u; matrix [no.of.pnts x 221 x 5]
        %     % t; matrix [1 x no.of.pnts]
        %     % x; matrix [1 x 221]
        %     % par; parameters
        %     % dev; device related parameters
        %     % tarr; time array (tplot)
        %     % pointtype = t
        %     % xrange; 2 element array with [xmin, xmax], in nm
        %     % - - - - - - - - - -
        %     sol = sol.dk.f;
        %     u = sol.u; % voltage, and charge carriers conc.
        %     t = sol.t(1:size(u, 1)); % 从 t 列表里提取前 [size(u, 1)/(矩阵 u 的行数)] 个元素
        %     x = sol.x;
        %     par = sol.par;
        %     dev = sol.par.dev; % device related parameters
        %     tarr = sol.t(size(u,1));
        %     pointtype = 't';
        %     xrange = [x(1), x(end)]*1e7;
        % 
        %     % - - - - - - - - - -
        %     % V_sol; matrix [no.of.pnts x 221]
        %     % n_sol; matrix [no.of.pnts x 221]
        %     % p_sol; matrix [no.of.pnts x 221]
        %     % V_sol_idx; matrix [1 x 221]
        %     % n_sol_idx; matrix [1 x 221]
        %     % p_sol_idx; matrix [1 x 221]
        %     % - - - - - - - - - -
        %     V_sol = sol.u(:,:,1); % u 的第一层 => voltage
        %     n_sol = sol.u(:,:,2); % u 的第二层 => e- conc.
        %     p_sol = sol.u(:,:,3); % u 的第三层 => h+ conc.
        %     c_sol = sol.u(:,:,4); % cation
        %     a_sol = sol.u(:,:,5); % anion
        % 
        %     V_idx = V_sol(idx, :); % choose the idx th coloumn
        %     n_idx = n_sol(idx, :);
        %     p_idx = p_sol(idx, :);
        %     c_idx = c_sol(idx, :);
        %     a_idx = a_sol(idx, :);
        % 
        %     % Calculate
        %     % - - - - - - - - - -
        %     % EA; matrix [1 x 221]
        %     % IP; matrix [1 x 221]
        %     % Ecb; matrix [1 x 221]
        %     % Evb; matrix [1 x 221]
        %     % Efn; matrix [1 x 221]
        %     % Efp; matrix [1 x 221]
        %     % - - - - - - - - - -
        %     EA = sol.par.dev.EA;
        %     IP = sol.par.dev.IP;
        % 
        %     Ecb = EA - V_idx; % minus the voltage at special time point
        %     Evb = IP - V_idx;
        % 
        %     Efn = zeros(size(n_idx,1), size(n_idx,2)); % 创建一个 n 的全零矩阵
        %     Efp = zeros(size(p_idx,1), size(p_idx,2)); % 创建一个 p 的全零矩阵 !!!
        % 
        %     % Boltzmann
        %     Efn = real(Ecb + (par.kB*par.T/par.q)*log((n_idx)./dev.Nc)); % Electron quasi-Fermi level
        %     Efp = real(Evb - (par.kB*par.T/par.q)*log((p_idx)./dev.Nv)); % Hole quasi-Fermi level
        % 
        %     figure(1);
        %     dfplot.x2d(sol, x, {Efn, Efp, Ecb, Evb}, {'E_{fn}', 'E_{fp}', 'E_{CB}', 'E_{VB}'},...
        %         {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0)
        % end
        % 
        % function EL_JVr(sol)
        %     [idx] = dfmem.JVcalcJDiff(sol);
        % 
        %     sol = sol.dk.r;
        % 
        %     u = sol.u;
        %     t = sol.t(1:size(u, 1));
        %     x = sol.x;
        %     par = sol.par;
        %     dev = sol.par.dev;
        %     tarr = sol.t(size(u,1));
        %     pointtype = 't';
        %     xrange = [x(1), x(end)]*1e7;
        % 
        %     V_sol = sol.u(:,:,1);
        %     n_sol = sol.u(:,:,2);
        %     p_sol = sol.u(:,:,3);
        %     c_sol = sol.u(:,:,4);
        %     a_sol = sol.u(:,:,5);
        % 
        %     V_idx = V_sol(length(sol.t) - idx, :);
        %     n_idx = n_sol(length(sol.t) - idx, :);
        %     p_idx = p_sol(length(sol.t) - idx, :);
        %     c_idx = c_sol(length(sol.t) - idx, :);
        %     a_idx = a_sol(length(sol.t) - idx, :);
        % 
        %     EA = sol.par.dev.EA;
        %     IP = sol.par.dev.IP;
        % 
        %     Ecb = EA - V_idx;
        %     Evb = IP - V_idx;
        % 
        %     Efn = zeros(size(n_idx,1), size(n_idx,2));
        %     Efp = zeros(size(p_idx,1), size(p_idx,2));
        % 
        %     % Boltzmann
        %     Efn = real(Ecb + (par.kB*par.T/par.q)*log((n_idx)./dev.Nc));
        %     Efp = real(Evb - (par.kB*par.T/par.q)*log((p_idx)./dev.Nv));
        % 
        %     figure(2);
        %     dfplot.x2d(sol, x, {Efn, Efp, Ecb, Evb}, {'E_{fn}', 'E_{fp}', 'E_{CB}', 'E_{VB}'},...
        %         {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0)
        % end
     end
end

