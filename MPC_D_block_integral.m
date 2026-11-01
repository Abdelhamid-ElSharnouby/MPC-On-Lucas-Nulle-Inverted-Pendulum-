function a = MPC_D(xhat, setpos_mm, enable, mpcCore, mpcState0, mpcOnline, zRefD, nMPC)
%#codegen
% MPC (candidate D) with integral action on the cart position (offset-free).
% Paste this whole function into the MATLAB Function block "MPC (candidate D)".
% Inputs/parameters are unchanged -> no change in the Ports and Data Manager.
%
% Integral action = reference shifting: the MPC tracks r + xi, where xi integrates
% the position error (r - x). Removes any constant offset (angle-sensor offset
% after Reset, drive friction / y1 offset). MPC data (mpcCore...) is unchanged.
Ti    = single(2);      % integral time [s]  (1..5 s tested in sim; 2 s chosen)
Tmpc  = single(0.01);   % MPC sample time [s] (= Ts in Calculations.m)
eWin  = single(0.02);   % integrate only while |error| < 20 mm (no windup during steps/pushes)
xiMax = single(0.03);   % limit of the correction: +-30 mm

persistent sd cnt aHold xi
if isempty(sd), sd = mpcState0; cnt = 0; aHold = single(0); xi = single(0); end
if enable                                   % enable = "State" (balance flag) = SSC Reset inport
    if cnt == 0
        r = single(setpos_mm)/single(1000);
        e = r - single(xhat(1));
        if abs(e) < eWin
            xi = xi + Tmpc/Ti*e;
        end
        xi = min(max(xi, -xiMax), xiMax);
        re = r + xi;                        % shifted reference
        od = mpcOnline;
        od.signals.ym  = single([xhat(1); xhat(3)]);
        od.signals.ref = [re single(0) single(0) re*single(zRefD)];
        sd.Plant = single(xhat(:));         % custom state estimation (LN observer)
        [aHold, sd] = mpcmoveCodeGeneration(mpcCore, sd, od);
    end
    cnt = cnt + 1; if cnt >= nMPC, cnt = 0; end
else
    cnt = 0; aHold = single(0); sd = mpcState0; xi = single(0);
end
a = double(aHold);
