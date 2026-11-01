%% experiments.m - systematic MPC tuning against the course LQR (simulation only)
%
%  Goal: find the MPC parameter combination that beats the course LQR on the
%  CO3620-2G, in simulation, and hand the best 2-3 candidates to the lab.
%  Needs: Control System Toolbox, MPC Toolbox, buildMPC.m (same folder).
%  No Simulink, no LN blocks -> runs on the home PC (R2025b).
%
%  ----------------------------------------------------------------------
%  APPROACH (and why it is not a plain grid over everything)
%  ----------------------------------------------------------------------
%  A full factorial grid (Ts x Np x blocking x limits x 3 weights x ...)
%  is thousands of runs, and most of them answer questions nobody asked.
%  The search is done in stages, each with a clear question:
%
%   0  verify    Is the fast QP used here identical to the toolbox MPC?
%   1  baseline  Does the simulator reproduce what the hardware showed?
%                (plot11: MPC with x/xdot limits drops catches;
%                 plot13: MPC without limits holds; LQR holds easy catches.)
%                If not, every ranking below is worthless.
%   2  structure Ts / Np / blocking / soft limits / Delta-a weight /
%                anti-windup, at fixed weights (candidate D). These choices
%                change the QP and the memory, so they are settled first.
%   3  weights   Log grid over q_x, q_v, q_phi on the best 2 structures.
%   4  refine    fminsearch in log space from the grid optimum.
%   5  robust    Best candidates + LQR: pendulum length +-20 %, 0.2 deg
%                angle offset, 10 ms sensor delay.
%   6  report    Choose 2-3 finalists, plots, paste-ready configuration.
%
%  The scenarios copy the hardware protocol: catches at given xi = phidot +
%  omega*phi (with the cart standing, or already moving at +-0.25 m/s), a
%  Set_Pos 0->100->0 step, and pushes of 1.2/1.5/1.8 rad/s. The plant includes everything the
%  earlier simulations missed: the drive's 0.5 m/s / 10 m/s^2 limits, the
%  v_set integrator that is RESET TO 0 at the catch, the LN observer with
%  its own drive model, quantisation and the 2 ms Control-Routine.
%
%  Score J (lower = better, LQR = 1.0 by construction) - relative to the
%  simulated course LQR, per scenario:
%    0.50 * catches  (peak x, t(|phi|<1 deg), t(|x|<5 mm); a drop costs 4,
%                     holding a catch that the LQR drops counts 0.5)
%    0.25 * step     (5 mm settling time, overshoot + undershoot)
%    0.25 * pushes   (peak x, t(|x|<5 mm); same drop rules as catches)
%    + 0.2 * max(0, aRMS/aRMS_LQR - 1)        (do not buy speed with noise)
%    + 1 per held scenario with peak |x| > 140 mm (10 mm from the drop limit)
%
%  Usage:  experiments                        % all stages (about 5-10 min)
%          RUN = struct('report',1); experiments   % only selected stages
%  Results are saved in experiments_results.mat after every stage.
%
%  FINDINGS THAT SHAPED THIS SCRIPT (2026-11-01)
%   * The course v_set integrator (a -> v_set -> y1) has NO limits, while the
%     drive caps the speed at 0.5 m/s. When the MPC keeps asking for a, v_set
%     winds up (-1.3 m/s in the push case) and braking comes far too late.
%     -> option aw: the MPC block tracks v_set itself and clips a so that
%        v_set stays within +-0.5 (anti-windup inside the MPC block only).
%   * Near an active limit, the MPC alternates +10/-10 every sample (the
%     pattern of plot10/plot11): the prediction model has no drive lag.
%     A Delta-a weight (Weights.MVRate, option du) damps it.
%   * A soft x limit only helps with a horizon of at least 0.2 s at Ts = 5 ms
%     or 0.3 s at Ts = 10 ms (blocking keeps the QP small).
%   * Single-precision QP is NOT a problem: u_single - u_double < 0.003 m/s^2.
%  ----------------------------------------------------------------------

% (no clear/clearvars here: RUN may be set before calling the script)
if ~exist('RUN','var') || ~isstruct(RUN)
    RUN = struct('verify',1,'baseline',1,'structure',1,'weights',1,'refine',1,'robust',1,'report',1);
end
for s = ["verify" "baseline" "structure" "weights" "refine" "robust" "report"]
    if ~isfield(RUN, s), RUN.(s) = 0; end
end
here    = fileparts(mfilename('fullpath'));  if isempty(here), here = pwd; end
resFile = fullfile(here, 'experiments_results.mat');
if isfile(resFile), S = load(resFile); else, S = struct(); end
mpcverbosity('off');

%% ---------------- Plant, observer and course LQR (as Calculations.m) ----
P.m = 0.3;  P.l = 0.165;  P.g = 9.81;  P.I = P.m*P.l^2/3;  P.fr = 0.3;
P.A = [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 P.m*P.g*P.l/(P.I+P.m*P.l^2) 0];
P.B = [0; 1; 0; P.m*P.l/(P.I+P.m*P.l^2)];
P.C = [1 0 0 0; 0 0 1 0];
TL  = 0.005;
P.L = place(P.A', P.C', [-1/TL -1/TL/1.01 -1/TL/1.02 -1/TL/1.03])';
P.K = lqr(P.A, P.B, diag([1000 0 1000 0]), 1);
P.Nbar = -1/([1 0 0 0]*((P.A - P.B*P.K)\P.B));      % = rscale(...) = -31.62
P.omega = sqrt(P.A(4,3));                          % 6.68 rad/s
P.qPos = 0.054/16000;  P.qAng = 2*pi/65535;
P.TpMax = 2.5e-4;                                  % plant integration step <= 0.25 ms

nominal = struct('lScale',1,'angOff',0,'delay',0,'vsetInit',0,'vsetClamp',0);
SC = makeScenarios(P);
fprintf('%d scenarios: %s\n', numel(SC), strjoin(string({SC.name}), ', '));

%% ---------------- Course LQR reference (500 Hz, as built in the lab) ----
LQR = struct('type','LQR','fISR',500,'name',"course LQR");
S.lqrM = runAll(LQR, SC, P, nominal);
fprintf('\nCourse LQR (sim):\n');  disp(metricTable(S.lqrM, SC));

%% ========================================================================
%% Stage 0: verify the fast QP against the toolbox object
if RUN.verify
    fprintf('\n=== Stage 0: QP vs MPC Toolbox ===\n');
    vc = { baseCfg(), ...
           setf(baseCfg(), 'moves',[1 1 1 2 3 4 8], 'xLim',0.14, 'vLim',0.5, 'du',10), ...
           setf(baseCfg(), 'Ts',0.005,'fISR',1000,'Np',40,'moves',[1 1 2 2 4 6 8 16],'xLim',0.14,'xFrom',10) };
    rng(1);
    for i = 1:numel(vc)
        [d, nAct] = verifyQP(vc{i}, P);
        fprintf('  %-34s max|u_QP - u_toolbox| = %.2e m/s^2  (%d of 300 states with active limits)\n', ...
            cfgName(vc{i}), d, nAct);
    end
end

%% ========================================================================
%% Stage 1: baseline - does the simulator reproduce the hardware?
if RUN.baseline
    fprintf('\n=== Stage 1: baseline ===\n');
    b1 = setf(baseCfg(), 'q',[1000 0 1000], 'xLim',0.15, 'vLim',0.5);   % plot10/11 build
    b2 = setf(baseCfg(), 'q',[1000 0 1000]);                            % plot13 build (noOVlimits)
    b3 = baseCfg();                                                     % candidate D, no limits
    S.base = evalCfgs({b1,b2,b3}, SC, P, nominal, S.lqrM);
    showRank(S.base, 3);
    fprintf(['  Hardware: the x0.15/v0.5 build chattered +-10 and dropped easy catches (plot10/11).\n' ...
             '  The simulator shows the chattering at active limits (push case), but not the\n' ...
             '  drop of the plot11 catch itself - keep that in mind when reading the ranking.\n']);
    save(resFile, '-struct', 'S');
end

%% ========================================================================
%% Stage 2: structure (weights fixed at candidate D)
if RUN.structure
    fprintf('\n=== Stage 2: structure (q = [5000 100 1000]) ===\n');
    hz = { struct('Ts',0.010,'fISR', 500,'Np',20,'blk',[1 1 1 2 3 4 8]), ...       % 0.2 s, 500 Hz
           struct('Ts',0.010,'fISR', 500,'Np',30,'blk',[1 1 1 2 3 4 8 10]), ...    % 0.3 s, 500 Hz
           struct('Ts',0.005,'fISR',1000,'Np',30,'blk',[1 1 2 2 4 6 14]), ...      % 0.15 s, 1 kHz
           struct('Ts',0.005,'fISR',1000,'Np',40,'blk',[1 1 2 2 4 6 8 16]) };      % 0.2 s, 1 kHz (memory tight)
    lims = { struct('xLim',Inf, 'vLim',Inf), struct('xLim',0.13,'vLim',Inf), ...
             struct('xLim',Inf, 'vLim',0.5), struct('xLim',0.13,'vLim',0.5) };
    cf = {};
    for h = 1:numel(hz)
        for li = 1:numel(lims)
            for du = [0 10]
                for aw = 0:1
                    c = baseCfg();
                    c.Ts = hz{h}.Ts;  c.fISR = hz{h}.fISR;  c.Np = hz{h}.Np;  c.moves = hz{h}.blk;
                    c.xLim = lims{li}.xLim;  c.vLim = lims{li}.vLim;  c.du = du;  c.aw = aw;
                    cf{end+1} = c; %#ok<SAGROW>
                end
            end
        end
    end
    cf{end+1} = setf(baseCfg(), 'aw',1);    % current build structure (Nc = 5) + anti-windup
    S.struct = evalCfgs(cf, SC, P, nominal, S.lqrM, true);
    showRank(S.struct, 15);
    save(resFile, '-struct', 'S');
end

%% ========================================================================
%% Stage 3: weights on the best two structures
if RUN.weights
    fprintf('\n=== Stage 3: weights ===\n');
    top = pickStructures(S.struct, 2);
    [qx, qv, qp] = ndgrid([1000 2000 5000 10000 20000], [0 100 300 1000], [300 1000 3000]);
    cf = {};
    for t = 1:numel(top)
        for k = 1:numel(qx)
            cf{end+1} = setf(top{t}, 'q', [qx(k) qv(k) qp(k)]); %#ok<SAGROW>
        end
    end
    % resumable: checkpoint every 10 configs (a long run can be interrupted)
    key = strjoin(cellfun(@cfgName, top, 'UniformOutput', false), ' | ');
    if ~isfield(S,'weightsKey') || ~strcmp(S.weightsKey, key) || ~isfield(S,'weights')
        S.weights = evalCfgs({}, SC, P, nominal, S.lqrM);  S.weightsKey = key;
    end
    for k = numel(S.weights)+1:numel(cf)
        S.weights(k) = evalCfgs(cf(k), SC, P, nominal, S.lqrM);
        if mod(k,10) == 0, fprintf('  %d/%d\n', k, numel(cf)); save(resFile, '-struct', 'S'); end
    end
    showRank(S.weights, 15);
    save(resFile, '-struct', 'S');
end

%% ========================================================================
%% Stage 4: local refinement of the weights (fminsearch in log space)
if RUN.refine
    fprintf('\n=== Stage 4: refine ===\n');
    [~, i0] = min([S.weights.J]);  c0 = S.weights(i0).cfg;
    obj = @(z) evalCfgs({setf(c0,'q',[10^z(1) max(0,10^z(2)-1) 10^z(3)])}, SC, P, nominal, S.lqrM).J;
    z0  = [log10(c0.q(1)) log10(c0.q(2)+1) log10(c0.q(3))];
    zo  = fminsearch(obj, z0, optimset('MaxFunEvals',40,'TolX',0.02,'Display','iter'));
    cR  = setf(c0, 'q', round([10^zo(1) max(0,10^zo(2)-1) 10^zo(3)], -1));
    S.refined = evalCfgs({cR}, SC, P, nominal, S.lqrM);
    fprintf('  grid best  %s  J = %.3f\n  refined    %s  J = %.3f\n', ...
        cfgName(c0), S.weights(i0).J, cfgName(cR), S.refined.J);
    save(resFile, '-struct', 'S');
end

%% ========================================================================
%% Stage 5: robustness of the best candidates (and of the LQR)
if RUN.robust
    fprintf('\n=== Stage 5: robustness ===\n');
    pool = S.weights;  if isfield(S,'refined'), pool = [pool S.refined]; end
    pool = [pool S.struct];
    pool = uniqueBest(pool, 6);
    cases = { setf(nominal,'lScale',0.8), setf(nominal,'lScale',1.2), ...
              setf(nominal,'angOff',0.2*pi/180), setf(nominal,'delay',0.010) };
    caseNames = ["l -20%","l +20%","offset 0.2deg","delay 10ms"];
    R = table();
    for c = 1:numel(cases)
        lq = runAll(LQR, SC, P, cases{c});
        R = [R; table(caseNames(c), "course LQR", sum(~[lq.held]), 1.0, ...
             'VariableNames',{'Case','Controller','Drops','J'})]; %#ok<AGROW>
        for i = 1:numel(pool)
            e = evalCfgs({pool(i).cfg}, SC, P, cases{c}, lq);
            R = [R; table(caseNames(c), string(cfgName(pool(i).cfg)), e.drops, e.J, ...
                 'VariableNames',{'Case','Controller','Drops','J'})]; %#ok<AGROW>
            pool(i).robJ(c) = e.J;  pool(i).robDrops(c) = e.drops - sum(~[lq.held]);
        end
    end
    disp(R);
    S.robust = pool;  S.robustTable = R;
    save(resFile, '-struct', 'S');
end

%% ========================================================================
%% Stage 6: report - finalists, course-model fixes, plots, paste-ready code
if RUN.report
    fprintf('\n=== Stage 6: report ===\n');
    pool = S.robust;
    ok   = arrayfun(@(p) all(p.robDrops <= 0), pool);       % never drops more than the LQR
    jAll = arrayfun(@(p) p.J + mean(p.robJ), pool);         % nominal + mean robust score
    jAll(~ok) = jAll(~ok) + 10;
    [~, ord] = sort(jAll);
    fin = pool(ord(1));                                      % F1: best overall
    is500 = arrayfun(@(p) p.cfg.fISR == 500, pool(ord));     % F2: best at the current 500 Hz build
    if fin.cfg.fISR ~= 500 && any(is500), fin(end+1) = pool(ord(find(is500,1))); end
    for k = ord(:)'                                          % F3: best with other limit setting
        if numel(fin) >= 3, break, end
        if ~any(arrayfun(@(f) isequal(limTag(f.cfg), limTag(pool(k).cfg)), fin)) && ok(k)
            fin(end+1) = pool(k); %#ok<AGROW>
        end
    end
    b0 = setf(baseCfg(),'q',[1000 0 1000],'xLim',0.15,'vLim',0.5);
    bBase = memBytes(b0, P);  bBig = memBytes(setf(b0,'Np',50), P);
    for i = 1:numel(fin), fin(i).memWords = memWordsEst(memBytes(fin(i).cfg, P), bBase, bBig); end
    S.finalists = fin;

    fprintf('\nFINALISTS (score J: LQR = 1, lower is better)\n');
    T = table();
    for i = 1:numel(fin)
        f = fin(i);  m = f.mets;  ca = [SC.group] == "catch";  lq = S.lqrM;
        T = [T; table("F"+i, string(cfgName(f.cfg)), f.J, mean(f.robJ), f.drops, ...
            mean([m(ca & [lq.held]).xPk]), mean([m(ca & [lq.held]).tPhi1]), mean([m(ca & [lq.held]).tX5]), ...
            m([SC.group]=="step").ts5, m(strcmp({SC.name},'push 1.5rad/s')).xPk, f.aRMS, f.memWords, f.maxIt, ...
            'VariableNames',{'ID','Config','J','Jrobust','Drops','catch_xPk_mm','catch_tPhi1_s', ...
            'catch_tX5_s','step_ts5_s','push_xPk_mm','aRMS','memWords','maxQPit'})]; %#ok<AGROW>
    end
    ca = [SC.group] == "catch";  lq = S.lqrM;
    T = [table("LQR","course LQR",1,1,sum(~[lq.held]), mean([lq(ca&[lq.held]).xPk]), ...
        mean([lq(ca&[lq.held]).tPhi1]), mean([lq(ca&[lq.held]).tX5]), lq([SC.group]=="step").ts5, ...
        lq(strcmp({SC.name},'push 1.5rad/s')).xPk, mean([lq([lq.held]).aRMS]), 0, 0, 'VariableNames',T.Properties.VariableNames); T];
    disp(T);  S.finalTable = T;

    % --- course-model fixes, shared by both controllers (separate question) ---
    fx = { nominal, setf(nominal,'vsetInit',1), setf(nominal,'vsetClamp',1), setf(nominal,'vsetInit',1,'vsetClamp',1) };
    fxn = ["as built","v_set IC = xdot_hat","v_set clamp +-0.5","both"];
    FX = table();
    for k = 1:numel(fx)
        lq = runAll(LQR, SC, P, fx{k});
        e  = evalCfgs({fin(1).cfg}, SC, P, fx{k}, S.lqrM);
        FX = [FX; table(fxn(k), sum(~[lq.held]), scoreVs(lq, S.lqrM, SC), e.drops, e.J, ...
            'VariableNames',{'CourseModel','LQR_drops','LQR_J','F1_drops','F1_J'})]; %#ok<AGROW>
    end
    fprintf('\nCourse-model fixes (J relative to the as-built LQR):\n');  disp(FX);  S.fixTable = FX;

    % --- plots: hardest held catch, step, push ---
    pick = [find(strcmp({SC.name},'catch xi-1.2 v0'),1), find([SC.group]=="step",1), find(strcmp({SC.name},'push 1.5rad/s'),1)];
    figure('Name','experiments: finalists vs course LQR','Color','w');
    for j = 1:3
        sc = SC(pick(j));
        o  = simRun(LQR, sc, P, nominal);
        subplot(3,3,j);   plot(o.t, o.x*1e3,'k--','LineWidth',1.2); hold on; title(sc.name); ylabel('x [mm]'); grid on
        subplot(3,3,3+j); plot(o.t, o.phi*180/pi,'k--','LineWidth',1.2); hold on; ylabel('\phi [deg]'); grid on
        subplot(3,3,6+j); plot(o.t, o.v,'k--','LineWidth',1.2); hold on; ylabel('xdot [m/s]'); xlabel('t [s]'); grid on
        for i = 1:numel(fin)
            o = simRun(makeCtl(fin(i).cfg, P), sc, P, nominal);
            subplot(3,3,j);   plot(o.t, o.x*1e3);
            subplot(3,3,3+j); plot(o.t, o.phi*180/pi);
            subplot(3,3,6+j); plot(o.t, o.v);
        end
    end
    subplot(3,3,1); legend(["LQR" "F"+(1:numel(fin))], 'Location','best');

    % --- paste-ready configuration + toolbox check (single precision, memory) ---
    fprintf('\n---- paste into Calculations.m (Step 5), one finalist at a time ----\n');
    for i = 1:numel(fin)
        c = fin(i).cfg;
        [mpcF, ~, bi] = buildMPC(c, P.A, P.B);
        [core, st0] = getCodeGenerationData(mpcF, 'DataType','single'); %#ok<ASGLU>
        fprintf(['%% F%d: %s\n' ...
                 'cfgD = struct(''Ts'',%g, ''Np'',%d, ''moves'',%s, ''q'',%s, ''r'',1, ...\n' ...
                 '              ''xLim'',%g, ''xFrom'',%d, ''vLim'',%g, ''ecr'',%g, ''du'',%g);\n' ...
                 'fISR = %d;   %% anti-windup in the MPC block: %d (1 = use MPC_D_aw block code)\n' ...
                 '[mpcD, LcD] = buildMPC(cfgD, A, B);  Ts = cfgD.Ts;  NpD = cfgD.Np;  qD = cfgD.q;\n' ...
                 '%% est. MPC memory %.1fk words (current build 7.0k, overflow above ~12.4k)\n\n'], ...
            i, cfgName(c), c.Ts, c.Np, mat2str(c.moves), mat2str(c.q), c.xLim, c.xFrom, c.vLim, c.ecr, c.du, ...
            c.fISR, c.aw, fin(i).memWords/1000);
        assert(isequal(class(st0.Plant),'single'));
        S.finalists(i).zRefD = bi.zRefD;
    end
    save(resFile, '-struct', 'S');
end

%% ########################################################################
%% Local functions
%% ########################################################################

function c = baseCfg()
% Candidate D, current hardware structure, no x/xdot limits (= plot13 build with D weights)
c = struct('Ts',0.01,'fISR',500,'Np',20,'moves',5,'q',[5000 100 1000],'r',1, ...
           'xLim',Inf,'xFrom',1,'vLim',Inf,'ecr',1,'rho',1e5,'obsClamp',0,'du',0,'aw',0);
end

function s = setf(s, varargin)
for k = 1:2:numel(varargin), s.(varargin{k}) = varargin{k+1}; end
end

function t = limTag(c)
t = sprintf('x%g/%d v%g du%g aw%d', c.xLim, c.xFrom, c.vLim, c.du, c.aw);
end

function n = cfgName(c)
if isscalar(c.moves), mv = sprintf('Nc%d', c.moves); else, mv = 'blk'; end
n = sprintf('q%g-%g-%g Ts%g Np%d %s', c.q(1), c.q(2), c.q(3), c.Ts*1e3, c.Np, mv);
if isfinite(c.xLim), n = [n sprintf(' x%g', c.xLim)]; if c.xFrom > 1, n = [n sprintf('@%d', c.xFrom)]; end, end
if isfinite(c.vLim), n = [n sprintf(' v%g', c.vLim)]; end
if c.du > 0, n = [n sprintf(' du%g', c.du)]; end
if c.rho ~= 1e5, n = [n sprintf(' rho%g', c.rho)]; end
if c.obsClamp, n = [n ' obsClamp']; end
if c.aw, n = [n ' AW']; end
end

function bv = blockVec(moves, Np)
if isscalar(moves), bv = [ones(1,moves-1) Np-moves+1];
else, bv = moves(:)'; bv(end) = bv(end) + Np - sum(bv); end
assert(all(bv >= 1) && sum(bv) == Np, 'blocking vector must sum to Np');
end

%% ---------------- condensed QP (identical cost/constraints to buildMPC) --
function qp = buildQP(c, P)
Ts = c.Ts;  Np = c.Np;
sd = c2d(ss(P.A, P.B, eye(4), 0), Ts, 'zoh');  Ad = sd.A;  Bd = sd.B;
Q  = diag([c.q(1) c.q(2) c.q(3) 0])/90;  R = c.r/90;
Pt = idare(Ad, Bd, Q, R);
bv = blockVec(c.moves, Np);  nb = numel(bv);
M  = zeros(Np, nb);  k0 = 0;
for j = 1:nb, M(k0+(1:bv(j)), j) = 1;  k0 = k0 + bv(j); end
Sx = zeros(4*Np, 4);  Su = zeros(4*Np, Np);  Ai = eye(4);
for i = 1:Np
    Ai = Ad*Ai;  Sx(4*i-3:4*i, :) = Ai;
    for j = 1:i, Su(4*i-3:4*i, j) = Ad^(i-j)*Bd; end
end
Qb = kron(eye(Np), Q);  Qb(end-3:end, end-3:end) = Pt;     % last step: terminal cost only
Sub = Su*M;  E = repmat([1;0;0;0], Np, 1);
H  = 2*(Sub'*Qb*Sub + R*(M'*M));
F  = 2*Sub'*Qb*Sx;  Fr = -2*Sub'*Qb*E;
Rdu = (c.du/20)^2;                                        % Weights.MVRate = du (MV scale 20)
DM  = (eye(Np) - diag(ones(Np-1,1),-1))*M;                % Delta-u = DM*z - e1*u_prev
H   = H + 2*Rdu*(DM'*DM);  Fu = -2*Rdu*DM(1,:)';
Ain = [eye(nb); -eye(nb)];  b0 = 10*ones(2*nb,1);  bx = zeros(2*nb, 4);
soft = isfinite(c.xLim) || isfinite(c.vLim);
if soft
    H = blkdiag(H, 2*c.rho);  F = [F; zeros(1,4)];  Fr = [Fr; 0];  Fu = [Fu; 0];
    Ain = [Ain zeros(2*nb,1); zeros(1,nb) -1];  b0 = [b0; 0];  bx = [bx; zeros(1,4)];
    if isfinite(c.xLim)
        ix = 4*(max(1,c.xFrom):Np) - 3;
        G = Sub(ix,:);  Gx = Sx(ix,:);  n = numel(ix);  e = -0.3*c.ecr*ones(n,1);
        Ain = [Ain; G e; -G e];  b0 = [b0; c.xLim*ones(2*n,1)];  bx = [bx; -Gx; Gx];
    end
    if isfinite(c.vLim)
        iv = 4*(1:Np) - 2;
        G = Sub(iv,:);  Gx = Sx(iv,:);  n = numel(iv);  e = -1.0*c.ecr*ones(n,1);
        Ain = [Ain; G e; -G e];  b0 = [b0; c.vLim*ones(2*n,1)];  bx = [bx; -Gx; Gx];
    end
end
qp.H = (H+H')/2;  qp.F = F;  qp.Fr = Fr;  qp.Fu = Fu;  qp.A = Ain;  qp.b0 = b0;  qp.bx = bx;
qp.nz = size(H,1);  qp.opt = mpcActiveSetOptions;  qp.opt.IntegrityChecks = false;
end

function [u, iA, it] = qpStep(qp, x, r, iA, up)
f = qp.F*x + qp.Fr*r + qp.Fu*up;                     % up = last applied MPC move
b = qp.b0 + qp.bx*x;
[z, ef, iA] = mpcActiveSetSolver(qp.H, f, qp.A, b, zeros(0,qp.nz), zeros(0,1), iA, qp.opt);
u = z(1);  it = max(ef, 0);
if ef <= 0, iA = false(size(b)); end                 % 0 = max iter, <0 infeasible: reset warm start
end

function [dmax, nAct] = verifyQP(c, P)
[o, ~, bi] = buildMPC(c, P.A, P.B);  qp = buildQP(c, P);
xc = mpcstate(o);  dmax = 0;  nAct = 0;
for k = 1:300
    x = [0.12*randn; 0.4*randn; 0.15*randn; 1.5*randn];  r = 0.1*(rand > 0.5);
    up = 4*randn;  xc.Plant = x;  xc.LastMove = up;
    ut = mpcmove(o, xc, x([1 3])', [r 0 0 r*bi.zRefD]);
    [uq, iA] = qpStep(qp, x, r, false(size(qp.b0)), up);
    dmax = max(dmax, abs(ut - uq));  nAct = nAct + any(iA(2*numel(blockVec(c.moves,c.Np))+2:end));   % skip MV rows and the slack row
end
end

%% ---------------- scenarios ---------------------------------------------
function SC = makeScenarios(P)
SC = struct('name',{},'group',{},'x0',{},'T',{},'rt',{},'rv',{},'xi',{});
d2r = pi/180;
% catches: phi0 = 10 deg coming through upright (phidot from xi), cart at +20 mm
for xi = [-0.6 -1.2 -1.7]
    for v0 = [-0.25 0 0.25]
        phi0 = 10*d2r;  pd0 = xi - P.omega*phi0;
        SC(end+1) = struct('name',sprintf('catch xi%.1f v%s', xi, vtag(v0)),'group',"catch", ...
            'x0',[0.02; v0; phi0; pd0],'T',3,'rt',0,'rv',0,'xi',xi); %#ok<AGROW>
    end
end
% catches falling away from upright (phi0 = -5 deg), cart standing
for xi = [-0.6 -1.2 -1.7]
    phi0 = -5*d2r;  pd0 = xi - P.omega*phi0;
    SC(end+1) = struct('name',sprintf('catch xi%.1f fall', xi),'group',"catch", ...
        'x0',[0.02; 0; phi0; pd0],'T',3,'rt',0,'rv',0,'xi',xi); %#ok<AGROW>
end
SC(end+1) = struct('name','step 0-100-0','group',"step",'x0',[0;0;0;0],'T',5,'rt',[0 0.1 2.6],'rv',[0 0.1 0],'xi',0);
for pv = [1.2 1.5 1.8]
    SC(end+1) = struct('name',sprintf('push %.1frad/s',pv),'group',"push",'x0',[0;0;0;pv],'T',3,'rt',0,'rv',0,'xi',pv); %#ok<AGROW>
end
end

function s = vtag(v)
if v == 0, s = '0'; elseif v > 0, s = sprintf('+%.2g', v); else, s = sprintf('%.2g', v); end
end

%% ---------------- closed-loop simulation (replica of A62_MPC_Sim) ------
function ctl = makeCtl(c, P)
ctl = struct('type','MPC','fISR',c.fISR,'nMPC',round(c.Ts*c.fISR),'qp',buildQP(c,P), ...
             'obsClamp',c.obsClamp,'aw',c.aw,'name',cfgName(c));
end

function out = simRun(ctl, sc, P, opt)
Tisr = 1/ctl.fISR;  nSub = ceil(Tisr/P.TpMax - 1e-9);  Tp = Tisr/nSub;
N = round(sc.T/Tisr);  l = P.l*opt.lScale;  kph = 3/(4*l);  fr4 = P.fr/4;
xc = sc.x0(1); vc = sc.x0(2); ph = sc.x0(3); phd = sc.x0(4); iP = 0;      % plant + drive PI
xh = sc.x0(:); vm = sc.x0(2); iPm = 0;                                   % LN observer (converged)
if opt.vsetInit, vset = xh(2); else, vset = 0; end                       % integrator reset at catch
isMPC = strcmp(ctl.type,'MPC');
if isMPC, iA = false(size(ctl.qp.b0)); end
cnt = 0; aHold = 0; itMax = 0;
nDel = round(opt.delay/Tisr);  buf = repmat([xc; ph+opt.angOff], 1, nDel+1);
dec = max(1, round(0.01/Tisr));  nL = floor(N/dec);
out.t = zeros(nL,1); out.x = out.t; out.phi = out.t; out.v = out.t; out.a = out.t; out.r = out.t;
iL = 0;  dropped = false;  aSq = 0;  nA = 0;
for n = 0:N-1
    t = n*Tisr;  r = sc.rv(find(t >= sc.rt - 1e-12, 1, 'last'));
    ym = [round(xc/P.qPos)*P.qPos; round((ph+opt.angOff)/P.qAng)*P.qAng];
    buf = [buf(:,2:end) ym];  y = buf(:,1);
    if isMPC
        if cnt == 0
            xs = xh;  if ctl.obsClamp, xs(2) = min(max(xs(2),-0.5),0.5); end
            [aHold, iA, it] = qpStep(ctl.qp, xs, r, iA, aHold);  itMax = max(itMax, it);
        end
        cnt = cnt + 1;  if cnt >= ctl.nMPC, cnt = 0; end
        a = aHold;
        if ctl.aw, a = min(max(a, (-0.5 - vset)/Tisr), (0.5 - vset)/Tisr); end   % v_set anti-windup
    else
        a = P.Nbar*r - P.K*xh;
    end
    y1 = 20*vset;  vcmd = min(max(0.05*y1,-0.5),0.5);
    % LN observer (Forward Euler, drive model V->m/s, sat 0.5, PI 100/30 +-10 back-calc)
    e = vcmd - vm;  u = 100*e + iPm;  am = min(max(u,-10),10);
    iPm = iPm + Tisr*(30*e + (am - u));
    xh  = xh + Tisr*(P.A*xh + P.B*am + P.L*(y - P.C*xh));  vm = vm + Tisr*am;
    % v_set integrator (Forward Euler)
    vset = vset + Tisr*a;  if opt.vsetClamp, vset = min(max(vset,-0.5),0.5); end
    % plant: drive speed loop + nonlinear pendulum
    for s = 1:nSub
        e = vcmd - vc;  u = 100*e + iP;  ac = min(max(u,-10),10);  iP = iP + Tp*(30*e + (ac - u));
        phdd = kph*(P.g*sin(ph) + ac*cos(ph)) - fr4*phd;
        vc = vc + Tp*ac;  xc = xc + Tp*vc;  phd = phd + Tp*phdd;  ph = ph + Tp*phd;
    end
    aSq = aSq + a^2;  nA = nA + 1;
    if mod(n, dec) == 0 && iL < nL
        iL = iL + 1;  out.t(iL) = t;  out.x(iL) = xc;  out.phi(iL) = ph;  out.v(iL) = vc;  out.a(iL) = a;  out.r(iL) = r;
    end
    if abs(ph) > 20*pi/180 || abs(xc) > 0.150, dropped = true; break, end   % course hold condition
end
fn = {'t','x','phi','v','a','r'};
for k = 1:numel(fn), out.(fn{k}) = out.(fn{k})(1:iL); end
out.dropped = dropped;  out.aRMS = sqrt(aSq/max(nA,1));  out.itMax = itMax;
end

%% ---------------- metrics and score -------------------------------------
function m = scenMetrics(o, sc)
t = o.t;  T = sc.T;
m = struct('held',~o.dropped,'xPk',max(abs(o.x))*1e3,'phiPk',max(abs(o.phi))*180/pi, ...
           'tPhi1',T,'tX5',T,'ts5',T,'ovun',100,'aRMS',o.aRMS,'it',o.itMax);
if o.dropped, return, end
m.tPhi1 = settleT(t, abs(o.phi), pi/180, 0, T);
m.tX5   = settleT(t, abs(o.x - o.r), 0.005, 0, T);
if sc.group == "step"
    ts = 0;  ou = 0;
    for k = 2:numel(sc.rt)
        t0 = sc.rt(k);  t1 = T;  if k < numel(sc.rt), t1 = sc.rt(k+1); end
        w = t >= t0 & t < t1;  sg = sign(sc.rv(k) - sc.rv(k-1));
        ts = ts + settleT(t(w), abs(o.x(w) - sc.rv(k)), 0.005, t0, t1);
        ou = ou + max(0, max(sg*(o.x(w) - sc.rv(k))))*1e3 + max(0, max(-sg*(o.x(w) - sc.rv(k-1))))*1e3;
    end
    m.ts5 = ts/(numel(sc.rt)-1);  m.ovun = ou/(numel(sc.rt)-1);
end
end

function ts = settleT(t, e, tol, t0, t1)
k = find(e > tol, 1, 'last');
if isempty(k), ts = 0; elseif k == numel(e), ts = t1 - t0; else, ts = t(k+1) - t0; end
end

function M = runAll(ctl, SC, P, opt)
for i = 1:numel(SC)
    M(i) = scenMetrics(simRun(ctl, SC(i), P, opt), SC(i)); %#ok<AGROW>
end
end

function [J, parts] = scoreVs(m, L, SC)
ca = find([SC.group] == "catch");  pc = [];
for i = ca
    if L(i).held
        if m(i).held, pc(end+1) = mean([m(i).xPk/L(i).xPk, (m(i).tPhi1+0.05)/(L(i).tPhi1+0.05), (m(i).tX5+0.05)/(L(i).tX5+0.05)]); %#ok<AGROW>
        else, pc(end+1) = 4; end %#ok<AGROW>
    else
        pc(end+1) = 1 - 0.5*m(i).held; %#ok<AGROW>
    end
end
i = find([SC.group] == "step");
if m(i).held, ps = mean([(m(i).ts5+0.05)/(L(i).ts5+0.05), (m(i).ovun+1)/(L(i).ovun+1)]); else, ps = 4; end
pq = [];
for i = find([SC.group] == "push")
    if L(i).held
        if m(i).held, pq(end+1) = mean([m(i).xPk/L(i).xPk, (m(i).tX5+0.05)/(L(i).tX5+0.05)]); %#ok<AGROW>
        else, pq(end+1) = 4; end %#ok<AGROW>
    else
        pq(end+1) = 1 - 0.5*m(i).held; %#ok<AGROW>
    end
end
pp = mean(pq);
aR = mean([m([m.held]).aRMS]) / mean([L([L.held]).aRMS]);
nearLimit = sum([m.held] & [m.xPk] > 140);
J = 0.5*mean(pc) + 0.25*ps + 0.25*pp + 0.2*max(0, aR-1) + nearLimit;
parts = [mean(pc) ps pp aR nearLimit];
end

function E = evalCfgs(cfgs, SC, P, opt, L, withMem)
if nargin < 6, withMem = false; end
E = struct('cfg',{},'name',{},'J',{},'parts',{},'drops',{},'aRMS',{},'maxIt',{},'memWords',{},'mets',{},'robJ',{},'robDrops',{});
if withMem
    base = setf(baseCfg(), 'q',[1000 0 1000], 'xLim',0.15, 'vLim',0.5);   % current lab build = 7.0k words
    bBase = memBytes(base, P);
    bBig  = memBytes(setf(base,'Np',50), P);       % Np = 50 build = 14.4k words (overflowed)
end
for i = 1:numel(cfgs)
    c = cfgs{i};
    m = runAll(makeCtl(c, P), SC, P, opt);
    [J, parts] = scoreVs(m, L, SC);
    mw = NaN;  if withMem, mw = memWordsEst(memBytes(c, P), bBase, bBig); end
    E(i) = struct('cfg',c,'name',cfgName(c),'J',J,'parts',parts,'drops',sum(~[m.held]), ...
        'aRMS',mean([m([m.held]).aRMS]),'maxIt',max([m.it]),'memWords',mw,'mets',{m}, ...
        'robJ',[],'robDrops',[]);
    if numel(cfgs) > 5 && mod(i, 10) == 0, fprintf('  %d/%d\n', i, numel(cfgs)); end
end
end

function w = memWordsEst(b, bBase, bBig)
% Linear fit through the two builds whose linker sizes are known:
% current lab build 7.0k words, Np = 50 build 14.4k words (RAML4L7 overflow above ~12.4k).
w = 7000 + (b - bBase)*(14400 - 7000)/(bBig - bBase);
end

function b = memBytes(c, P)
o = buildMPC(c, P.A, P.B);
[core, st0, on] = getCodeGenerationData(o, 'DataType','single'); %#ok<ASGLU>
w = whos('core','st0','on');  b = sum([w.bytes]);
end

function showRank(E, n)
[~, o] = sort([E.J]);  o = o(1:min(n,numel(o)));
T = table(string({E(o).name})', [E(o).J]', vertcat(E(o).parts), [E(o).drops]', [E(o).aRMS]', [E(o).maxIt]', [E(o).memWords]', ...
    'VariableNames',{'Config','J','parts_catch_step_push_aR_near','Drops','aRMS','maxQPit','memWords'});
disp(T);
end

function top = pickStructures(E, n)
[~, o] = sort([E.J]);  top = {};  seen = strings(0);
for k = o
    if E(k).memWords > 10500, continue, end          % keep a margin to the ~12.4k-word RAML4L7 limit
    key = string(E(k).cfg.Ts) + "/" + E(k).cfg.Np;    % one structure per horizon setup
    if ~any(seen == key), top{end+1} = E(k).cfg; seen(end+1) = key; end %#ok<AGROW>
    if numel(top) >= n, break, end
end
end

function U = uniqueBest(E, n)
[~, o] = sort([E.J]);  U = E([]);  seen = strings(0);
for k = o
    if ~any(seen == E(k).name), U(end+1) = E(k); seen(end+1) = E(k).name; end %#ok<AGROW>
    if numel(U) >= n, break, end
end
end

function T = metricTable(M, SC)
T = table(string({SC.name})', [M.held]', [M.xPk]', [M.phiPk]', [M.tPhi1]', [M.tX5]', [M.ts5]', [M.ovun]', [M.aRMS]', ...
    'VariableNames',{'Scenario','held','xPk_mm','phiPk_deg','tPhi1_s','tX5_s','ts5_s','ovun_mm','aRMS'});
end
