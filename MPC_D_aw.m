function [a, qpIter, qpIterMax] = MPC_D_aw(xhat, setpos_mm, enable, mpcCore, mpcState0, mpcOnline, zRefD, nMPC, Tisr, awOn)
%#codegen
% MATLAB Function block "MPC (candidate D)" with v_set anti-windup.
% Paste this into the block in A62_SwingUp_Control_MPC.slx (and A62_MPC_Sim.slx).
%
% Changes against the current block:
%   * new non-tunable parameters Tisr (= 1/fISR) and awOn (1 = anti-windup on)
%   * vInt copies the course "Discrete-Time Integrator" (a -> v_set, Forward
%     Euler, reset to 0 when the balance flag rises). The output a is clipped
%     so that v_set stays within +-0.5 m/s, the drive's speed limit. Without
%     it, v_set winds up to -1.3 m/s in a hard push and braking comes too late.
%   * only valid while useMPC = 1 (with useMPC = 0 the integrator gets the LQR
%     output and vInt is just a shadow value - harmless).
%
% Block parameters (Scope = Parameter, non-tunable):
%   mpcCore, mpcState0, mpcOnline, zRefD, nMPC, Tisr, awOn
persistent sd cnt aHold vInt itHold itMax
if isempty(sd)
    sd = mpcState0; cnt = 0; aHold = single(0); vInt = single(0);
    itHold = single(0); itMax = single(0);
end
vmax = single(0.5);
T    = single(Tisr);
if enable
    if cnt == 0
        r  = single(setpos_mm)/single(1000);
        od = mpcOnline;
        od.signals.ym  = single([xhat(1); xhat(3)]);
        od.signals.ref = [r single(0) single(0) r*single(zRefD)];   % outputs [x phi xdot z1..z4]
        sd.Plant = single(xhat(:));                                  % custom estimation (LN observer)
        [aHold, sd, info] = mpcmoveCodeGeneration(mpcCore, sd, od);
        itHold = single(info.Iterations);       % >0 ok, 0 = max iter, -1 infeasible, -2 numerical
        if itHold > itMax, itMax = itHold; end
    end
    cnt = cnt + 1; if cnt >= nMPC, cnt = 0; end
    aOut = aHold;
    if awOn
        aOut = min(max(aOut, (-vmax - vInt)/T), (vmax - vInt)/T);   % keep v_set within +-0.5 m/s
    end
    vInt = vInt + T*aOut;                                            % same as the course integrator
else
    cnt = 0; aHold = single(0); sd = mpcState0; vInt = single(0); itHold = single(0);
    aOut = single(0);
end
a = double(aOut);  qpIter = double(itHold);  qpIterMax = double(itMax);
