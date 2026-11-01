function [mpcD, LcD, info] = buildMPC(cfg, A, B)
%BUILDMPC  MPC Toolbox object for the LN pendulum, from one config struct.
%   [mpcD, LcD, info] = buildMPC(cfg, A, B)
%
%   Same design as Calculations.m Step 5 (terminal cost, custom estimator),
%   generalised so experiments.m and Calculations.m build exactly the same
%   controller. Works on R2022a (lab) and R2025b (home).
%
%   cfg fields (defaults in brackets):
%     Ts     MPC sample time [s]                         (0.01)
%     Np     prediction horizon (steps)                  (20)
%     moves  Nc (scalar) or blocking vector, sum <= Np   (5)
%     q      [q_x q_v q_phi]  (q_phidot = 0), R = r      ([5000 100 1000])
%     r      input weight                                (1)
%     xLim   soft cart limit [m], Inf = none             (Inf)
%     xFrom  first horizon step with the x limit         (1)
%     vLim   soft speed limit [m/s], Inf = none          (Inf)
%     ecr    ECR of the soft limits                      (1)
%     rho    Weights.ECR (slack penalty)                 (1e5)
%     du     Weights.MVRate (Delta-a weight)             (0)
%
%   Outputs of the plant model: [x phi xdot z1..z4], z = LcD*x (terminal cost).
%   Reference for the MPC block: [r 0 0 r*LcD(:,1)'],  i.e. zRefD = LcD(:,1)'.

d = struct('Ts',0.01,'Np',20,'moves',5,'q',[5000 100 1000],'r',1, ...
           'xLim',Inf,'xFrom',1,'vLim',Inf,'ecr',1,'rho',1e5,'du',0);
fn = fieldnames(d);
for k = 1:numel(fn), if ~isfield(cfg,fn{k}) || isempty(cfg.(fn{k})), cfg.(fn{k}) = d.(fn{k}); end, end

Ts = cfg.Ts;  Np = cfg.Np;
Cm = [1 0 0 0; 0 0 1 0; 0 1 0 0];                 % x, phi (measured), xdot (unmeasured)
pd = c2d(ss(A, B, Cm, 0), Ts, 'zoh');
Qd = diag([cfg.q(1) cfg.q(2) cfg.q(3) 0])/90;   Rd = cfg.r/90;
PD = idare(pd.A, pd.B, Qd, Rd);                   % terminal cost = infinite-horizon cost-to-go
LcD = chol(PD);

pD = ss(pd.A, pd.B, [Cm; LcD], 0, Ts);
pD = setmpcsignals(pD, 'MV',1, 'MO',[1 2], 'UO',3:7);

old = mpcverbosity('off');
mpcD = mpc(pD, Ts, Np, cfg.moves);
wu   = 20*sqrt(Rd);
setterminal(mpcD, struct('Weight',[0 0 0 1 1 1 1]), struct('Weight',wu));   % FIRST (R2022a quirk)

% --- limits AFTER setterminal ---
mpcD.MV.Min = -10;  mpcD.MV.Max = 10;  mpcD.MV.ScaleFactor = 20;
mpcD.MV.MinECR = 0; mpcD.MV.MaxECR = 0;                                      % a: hard
mpcD.OV(1).ScaleFactor = 0.3;  mpcD.OV(2).ScaleFactor = 0.7;  mpcD.OV(3).ScaleFactor = 1;
for k = 1:7, mpcD.OV(k).Min = -Inf; mpcD.OV(k).Max = Inf; end
if isfinite(cfg.xLim)
    lim = cfg.xLim*ones(1,Np);  lim(1:cfg.xFrom-1) = Inf;                   % time-varying: off on early steps
    mpcD.OV(1).Min = -lim;  mpcD.OV(1).Max = lim;
    mpcD.OV(1).MinECR = cfg.ecr;  mpcD.OV(1).MaxECR = cfg.ecr;
end
if isfinite(cfg.vLim)
    mpcD.OV(3).Min = -cfg.vLim;  mpcD.OV(3).Max = cfg.vLim;
    mpcD.OV(3).MinECR = cfg.ecr; mpcD.OV(3).MaxECR = cfg.ecr;
end

% --- weights: exactly Np rows, last row = terminal row ---
mpcD.Weights.MV     = wu;
mpcD.Weights.MVRate = cfg.du;
mpcD.Weights.ECR    = cfg.rho;
Wov = repmat([0.3*sqrt(Qd(1,1)) 0.7*sqrt(Qd(3,3)) sqrt(Qd(2,2)) 0 0 0 0], Np, 1);
Wov(end,:) = [0 0 0 1 1 1 1];
mpcD.Weights.OV = Wov;
assert(size(mpcD.Weights.OV,1) == Np && isequal(mpcD.Weights.OV(end,:), [0 0 0 1 1 1 1]), ...
    'buildMPC: terminal weight row missing');
assert(mpcD.MV.Min(1) == -10 && mpcD.MV.MinECR(1) == 0, 'buildMPC: MV limit not hard +-10');

% --- custom state estimation: MPC state = [x xdot phi phidot] from the LN observer ---
setoutdist(mpcD, 'model', tf(zeros(7,1)));
setEstimator(mpcD, 'custom');
mpcverbosity(old);

info.zRefD = LcD(:,1)';
info.Qd = Qd;  info.Rd = Rd;  info.PD = PD;  info.Ad = pd.A;  info.Bd = pd.B;
info.cfg = cfg;
end
