%% CalcExp.m - Calculations.m with the MPC tuned by experiments.m (2026-11-01)
%  Differences to Calculations.m:
%    * finalist selects the candidate: "F1" (default), "F2", "F3"
%    * Step 2: Ts / fISR follow the finalist (F1/F2: 5 ms, 1 kHz; F3: 10 ms, 500 Hz)
%    * Step 5: mpcD is built by buildMPC(cfgD, A, B) (needs buildMPC.m in this folder)
%    * Step 6: awOn = 1 for the MPC block with anti-windup (MPC_D_aw.m)
%  Run:  clear; CalcExp;     then set the LucasNuelleHardware frequency to fISR, Ctrl+B.
%  Details: claude/Tuning_Handover.md
%% Inverted Pendulum State Space Model
m = 0.3;
l = 0.165;
I = m*l^2/3;
g = 9.81;

A=[ 0 1 0 0;
    0 0 0 0
    0 0 0 1
    0 0 m*g*l/(I+m*l^2) 0];
B=[0;1;0;m*l/(I+m*l^2)];
C=[ 1 0 0 0
    0 0 1 0];
%% Inverted Pendulum State Space Observer
TL = 0.005;

EW =[-1/TL -1/TL/1.01 -1/TL/1.02 -1/TL/1.03];

L = (place(A',C',EW))';
%% Inverted Pendulum State Space Controller
Q = diag([1000 0 1000 0]);

R = 1;

K = lqr(A,B,Q,R);

Cn=[1 0 0 0];

sys_ss = ss(A,B,Cn,0);

Nbar = rscale(sys_ss,K);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Container Control
% m = 0.3;
% l = 0.165;
% I = m*l^2/3;
% g = 9.81;
% A=[ 0 1 0 0;
%     0 0 0 0
%     0 0 0 1
%     0 0 -m*g*l/(I+m*l^2) 0];
% B=[0;1;0;-m*l/(I+m*l^2)];
% C=[ 1 0 0 0
%     0 0 1 0];
% % Observer
% TL = 0.005;
% 
% EW =[-1/TL -1/TL/1.01 -1/TL/1.02 -1/TL/1.03];
% 
% L = (place(A',C',EW))';
% % Controller
% Q = diag([1000 0 1000 0]);
% 
% R = 1;
% 
% K = lqr(A,B,Q,R);
% 
% Cn=[1 0 0 0];
% 
% sys_ss = ss(A,B,Cn,0);
% 
% Nbar = rscale(sys_ss,K);
finalist = "F1";       % finalist from experiments.m: "F1" | "F2" | "F3"
switch finalist
    case {"F1","F2"}
        Ts   = 0.005;  % MPC sample time [s]
        fISR = 1000;   % Control-Routine rate -> set the LucasNuelleHardware block to 1000 Hz
    case "F3"
        Ts   = 0.01;
        fISR = 500;    % current lab rate
    otherwise
        error('finalist must be "F1", "F2" or "F3"');
end
nMPC = Ts*fISR;        % ISR calls per MPC step (= 5 for all finalists)

Cm = [1 0 0 0;         % x
      0 0 1 0;         % phi
      0 1 0 0];        % xdot (unmeasured, for the speed constraint)
plant_c = ss(A, B, Cm, 0);
plant_c.InputName  = {'a'};
plant_c.OutputName = {'x','phi','xdot'};
plant_c.StateName  = {'x','xdot','phi','phidot'};

plant_d = c2d(plant_c, Ts, 'zoh');
plant_d = setmpcsignals(plant_d, 'MV',1, 'MO',[1 2], 'UO',3);

Np = 50;  Nc = 5;                      % prediction / control horizon (NOT m!)
mpcobj = mpc(plant_d, Ts, Np, Nc);
% Manipulated variable: cart acceleration a [m/s^2]
mpcobj.MV.Min = -10;   mpcobj.MV.Max = 10;
mpcobj.MV.ScaleFactor = 20;

% Output constraints (soft by default)
mpcobj.OV(1).Min = -0.15;  mpcobj.OV(1).Max = 0.15;   % cart travel [m]
mpcobj.OV(3).Min = -0.5;   mpcobj.OV(3).Max = 0.5;    % cart speed  [m/s]

% Scale factors
mpcobj.OV(1).ScaleFactor = 0.3;
mpcobj.OV(2).ScaleFactor = 0.7;
mpcobj.OV(3).ScaleFactor = 1;

% First-guess weights
mpcobj.Weights.MV     = 0;
mpcobj.Weights.MVRate = 0.1;
mpcobj.Weights.OV     = [1 1 0];

r  = [0.1 0 0];                 % references for [x phi xdot]
[y, t, u] = sim(mpcobj, 300, r);   % 300 steps = 3 s

% figure
% subplot(4,1,1); plot(t, y(:,1)*1000); ylabel('x [mm]');      grid on
% subplot(4,1,2); plot(t, y(:,2)*180/pi); ylabel('\phi [deg]'); grid on
% subplot(4,1,3); plot(t, y(:,3)); ylabel('xdot [m/s]');        grid on
% subplot(4,1,4); stairs(t, u); ylabel('a [m/s^2]'); xlabel('t [s]'); grid on
%% 4a: LQR-equivalent stage weights (scaled units)
wx   = 1;
wphi = wx*0.7/0.3;                   % 2.333
wu   = 20*sqrt((wx/0.3)^2/1000);     % 2.108  -> same Q/R ratio as Calculations.m
Qd = diag([(wx/0.3)^2, 0, (wphi/0.7)^2, 0]);
Rd = (wu/20)^2;

%% 4b: terminal weight = infinite-horizon cost-to-go
[P, Kd] = idare(plant_d.A, plant_d.B, Qd, Rd);
Lc = chol(P);                        % Lc'*Lc = P

%% MPC with 4 extra outputs z = Lc*x (used only at the last step)
plantT = ss(plant_d.A, plant_d.B, [Cm; Lc], 0, Ts);
plantT = setmpcsignals(plantT, 'MV',1, 'MO',[1 2], 'UO',3:7);

mpcT = mpc(plantT, Ts, Np, Nc);
Yt.Weight = [0 0 0 1 1 1 1];         % last step: only z'z = x'Px
Ut.Weight = wu;
setterminal(mpcT, Yt, Ut);           % FIRST setterminal ...

mpcT.MV.Min = -10;  mpcT.MV.Max = 10;  mpcT.MV.ScaleFactor = 20;   % ... THEN limits
mpcT.OV(1).Min = -0.15; mpcT.OV(1).Max = 0.15; mpcT.OV(1).ScaleFactor = 0.3;
mpcT.OV(2).ScaleFactor = 0.7;
mpcT.OV(3).Min = -0.5;  mpcT.OV(3).Max = 0.5;  mpcT.OV(3).ScaleFactor = 1;
mpcT.MV.MinECR = 0;  mpcT.MV.MaxECR = 0;      % hard MV limit (setterminal makes it soft)

mpcT.Weights.MV     = wu;
mpcT.Weights.MVRate = 0;
Wov = repmat([wx wphi 0 0 0 0 0], Np, 1);
Wov(end,:) = [0 0 0 1 1 1 1];               % last step: z'z = x'Px
mpcT.Weights.OV = Wov;
%% Simulate the 100 mm step
xr = [0.1; 0; 0; 0];
rT = [0.1 0 0 (Lc*xr)'];             % z must track Lc*x_ref, or there is an offset
[yT, tT, uT] = sim(mpcT, 300, rT);

%% Course LQR for comparison (continuous, no limits)
clLQR = ss(A - B*K, B*Nbar*0.1, eye(4), 0);
[xL, tL] = step(clLQR, 3);
aL = Nbar*0.1 - xL*K';

% figure
% subplot(3,1,1); plot(t,y(:,1)*1e3, tT,yT(:,1)*1e3, tL,xL(:,1)*1e3); ylabel('x [mm]'); grid on
% legend('MPC Step 3','MPC LQR-equiv + terminal','LQR (course)','Location','southeast')
% subplot(3,1,2); plot(t,y(:,2)*180/pi, tT,yT(:,2)*180/pi, tL,xL(:,3)*180/pi); ylabel('\phi [deg]'); grid on
% subplot(3,1,3); stairs(t,u); hold on; stairs(tT,uT); plot(tL,aL); ylabel('a [m/s^2]'); xlabel('t [s]'); grid on
opts = {'SettlingTimeThreshold',0.02, 'RiseTimeLimits',[0.1 0.9]};
s3 = stepinfo(y(:,1),  t,  0.1, opts{:});    % MPC Step 3
sT = stepinfo(yT(:,1), tT, 0.1, opts{:});    % MPC LQR-equiv + terminal
sL = stepinfo(xL(:,1), tL, 0.1, opts{:});    % LQR (course)
S  = [s3; sT; sL];

% extra figures that stepinfo does not give
phiMax = [max(abs(y(:,2)))  max(abs(yT(:,2)))  max(abs(xL(:,3)))]'*180/pi;
vMax   = [max(abs(y(:,3)))  max(abs(yT(:,3)))  max(abs(xL(:,2)))]';
aMax   = [max(abs(u))       max(abs(uT))       max(abs(aL))]';
aRMS   = [sqrt(trapz(t,u.^2)/t(end)) ...
          sqrt(trapz(tT,uT.^2)/tT(end)) ...
          sqrt(trapz(tL,aL.^2)/tL(end))]';

Controller = ["MPC Step 3"; "MPC LQR-eq + terminal"; "LQR (course)"];
Results = table(Controller, [S.RiseTime]', [S.SettlingTime]', ...
    [S.Overshoot]', [S.Undershoot]', phiMax, vMax, aMax, aRMS, ...
    'VariableNames', {'Controller','Rise_s','Settle_2pct_s', ...
    'Overshoot_mm','Undershoot_mm','phiMax_deg','vMax_mps','aMax_mps2','aRMS_mps2'});
mpcverbosity('off');                  % silence the default-value messages
opts  = {'SettlingTimeThreshold',0.02, 'RiseTimeLimits',[0.1 0.9]};
cand  = [ 1000   0 1000;              % [q_x  q_v  q_phi],  R = 1
          3000   0 1000;
          5000 100 1000;
         10000   0 1000];
names = ["LQR-eq"; "A qx3000"; "D qx5000 qv100"; "B qx10000"];
res = table();  Tc = {};  Yc = {};  Uc = {};  mpcs = {};

for i = 1:size(cand,1)
    Qd = diag([cand(i,1) cand(i,2) cand(i,3) 0])/90;   Rd = 1/90;
    P  = idare(plant_d.A, plant_d.B, Qd, Rd);
    Lc = chol(P);

    pT  = ss(plant_d.A, plant_d.B, [Cm; Lc], 0, Ts);
    pT  = setmpcsignals(pT, 'MV',1, 'MO',[1 2], 'UO',3:7);
    obj = mpc(pT, Ts, 50, 5);
    wu  = 20*sqrt(Rd);
    setterminal(obj, struct('Weight',[0 0 0 1 1 1 1]), struct('Weight',wu));   % FIRST

    obj.MV.Min = -10;  obj.MV.Max = 10;  obj.MV.ScaleFactor = 20;            % THEN limits
    obj.OV(1).Min = -0.15; obj.OV(1).Max = 0.15; obj.OV(1).ScaleFactor = 0.3;
    obj.OV(2).ScaleFactor = 0.7;
    obj.OV(3).Min = -0.5;  obj.OV(3).Max = 0.5;  obj.OV(3).ScaleFactor = 1;
    obj.MV.MinECR = 0;  obj.MV.MaxECR = 0;        % hard MV limit (setterminal makes it soft)

    obj.Weights.MV     = wu;
    obj.Weights.MVRate = 0;
    Wov = repmat([0.3*sqrt(Qd(1,1))  0.7*sqrt(Qd(3,3))  sqrt(Qd(2,2))  0 0 0 0], 50, 1);
    Wov(end,:) = [0 0 0 1 1 1 1];                 % step 50 = terminal step: z'z = x'Px
    obj.Weights.OV = Wov;
    rr = [0.1 0 0 (Lc*[0.1;0;0;0])'];
    [yy, tt, uu] = sim(obj, 300, rr);
    s = stepinfo(yy(:,1), tt, 0.1, opts{:});
    res = [res; table(names(i), s.RiseTime, s.SettlingTime, s.Overshoot, s.Undershoot, ...
           max(abs(yy(:,2)))*180/pi, max(abs(yy(:,3))), max(abs(uu)), ...
           sqrt(trapz(tt,uu.^2)/tt(end)))];
    Tc{i} = tt;  Yc{i} = yy;  Uc{i} = uu;  mpcs{i} = obj;
end
res.Properties.VariableNames = {'Controller','Rise_s','Settle_2pct_s','Overshoot_mm', ...
    'Undershoot_mm','phiMax_deg','vMax_mps','aMax_mps2','aRMS_mps2'};

% figure
% for i = 1:4
%     subplot(3,1,1); plot(Tc{i}, Yc{i}(:,1)*1e3); hold on
%     subplot(3,1,2); plot(Tc{i}, Yc{i}(:,2)*180/pi); hold on
%     subplot(3,1,3); stairs(Tc{i}, Uc{i}); hold on
% end
% subplot(3,1,1); plot(tL, xL(:,1)*1e3, 'k--'); ylabel('x [mm]'); grid on
% legend([names; "LQR (course)"], 'Location','southeast')
% subplot(3,1,2); plot(tL, xL(:,3)*180/pi, 'k--'); ylabel('\phi [deg]'); grid on
% subplot(3,1,3); plot(tL, aL, 'k--'); ylabel('a [m/s^2]'); xlabel('t [s]'); grid on
% mpcverbosity('on');

%% Check: limits and weight matrices of the candidate MPCs
for i = 1:numel(mpcs)
    W = mpcs{i}.Weights.OV;
    fprintf('%-16s MV [%g %g]  xdot [%g %g]  W size %dx%d  row1 %s  row50 %s\n', names(i), ...
        mpcs{i}.MV.Min(1), mpcs{i}.MV.Max(1), mpcs{i}.OV(3).Min(1), mpcs{i}.OV(3).Max(1), ...
        size(W,1), size(W,2), mat2str(W(1,:),3), mat2str(W(end,:),3));
end

%% ========================================================================
%% Step 5: final MPC (candidate D) with custom state estimation
%  The MPC receives x_hat = [x; xdot; phi; phidot] from the LN observer
%  ("Inverted Pendulum State Space Observer" in A62), so the MPC Toolbox's
%  own Kalman filter and its output-disturbance state are removed.
%  Toolbox quirk (R2022a): setterminal changes the limits and ECRs, so it
%  is called FIRST and all limits/ECRs are set AFTER it.
mpcverbosity('off');
% Finalists of experiments.m: blocking, soft x limit +-0.13 m, no xdot limit,
% a = +-10 hard, terminal cost, custom estimation (all inside buildMPC.m).
switch finalist
    case "F1", qCand = [5000 100 1000];  mvCand = [1 1 2 2 4 6 14];      % candidate D weights
    case "F2", qCand = [4010 130 350];   mvCand = [1 1 2 2 4 6 14];      % fminsearch optimum
    case "F3", qCand = [5000 100 1000];  mvCand = [1 1 1 2 3 4 8 10];    % 500 Hz fallback
end
cfgD = struct('Ts',Ts, 'Np',30, 'moves',mvCand, 'q',qCand, 'r',1, ...
              'xLim',0.13, 'xFrom',1, 'vLim',Inf, 'ecr',1, 'du',0);
[mpcD, LcD] = buildMPC(cfgD, A, B);           % same builder as experiments.m
NpD = cfgD.Np;  qD = cfgD.q;
xcD = mpcstate(mpcD);
fprintf('\nStep 5: MPC state size = %d (must be 4: [x xdot phi phidot])\n', numel(xcD.Plant));

% reference vector for the 7 outputs [x phi xdot z1..z4], x_ref in metres
mpcRef = @(xref) [xref 0 0 (LcD*[xref;0;0;0])'];

%% Step 5 check: closed loop with the state fed in directly (perfect x_hat)
%  1) 100 mm step   2) push on the pendulum (phidot0 = 1.5 rad/s)   3) 140 mm step
scen5 = {'100 mm step', [0;0;0;0],   0.10;
         'push 1.5 rad/s', [0;0;0;1.5], 0;
         '140 mm step', [0;0;0;0],   0.14};
N5 = 300;  t5 = (1:N5)'*Ts;
chk5 = table();
for j = 1:size(scen5,1)
    x = scen5{j,2};  rr = mpcRef(scen5{j,3});  xc = mpcstate(mpcD);
    xs = zeros(N5,4);  us = zeros(N5,1);  it = zeros(N5,1);  ok = true;
    for k = 1:N5
        xc.Plant = x;                                   % <- LN observer output on the DSP
        [uk, info] = mpcmove(mpcD, xc, x([1 3])', rr);
        ok = ok && any(strcmp(info.QPCode, {'feasible','optimal'}));
        x = plant_d.A*x + plant_d.B*uk;
        xs(k,:) = x';  us(k) = uk;  it(k) = info.Iterations;
    end
    chk5 = [chk5; table(string(scen5{j,1}), min(xs(:,1))*1e3, max(xs(:,1))*1e3, ...
        max(abs(xs(:,3)))*180/pi, max(abs(xs(:,2))), max(abs(us)), max(it), ok)];
    X5{j} = xs;  U5{j} = us;
end
chk5.Properties.VariableNames = {'Scenario','xMin_mm','xMax_mm','phiMax_deg', ...
    'vMax_mps','aMax_mps2','maxQPiter','allOptimal'};
s5 = stepinfo(X5{1}(:,1), t5, 0.1, opts{:});
fprintf('100 mm step: settle %.3f s, overshoot %.2f mm, undershoot %.2f mm (must match res row D)\n', ...
    s5.SettlingTime, s5.Overshoot, s5.Undershoot);


% figure('Name','Step 5: final MPC (D), custom state')
% lbl = string(scen5(:,1));
% for j = 1:size(scen5,1)
%     subplot(4,1,1); plot(t5, X5{j}(:,1)*1e3); hold on
%     subplot(4,1,2); plot(t5, X5{j}(:,3)*180/pi); hold on
%     subplot(4,1,3); plot(t5, X5{j}(:,2)); hold on
%     subplot(4,1,4); stairs(t5, U5{j}); hold on
% end
% subplot(4,1,1); ylabel('x [mm]'); grid on; yline([-150 150],'r:'); legend(lbl,'Location','best')
% subplot(4,1,2); ylabel('\phi [deg]'); grid on
% subplot(4,1,3); ylabel('xdot [m/s]'); grid on; yline([-0.5 0.5],'r:')
% subplot(4,1,4); ylabel('a [m/s^2]'); xlabel('t [s]'); grid on; yline([-10 10],'r:')
% mpcverbosity('on');
%% ========================================================================
%% Step 6: workspace data for the Simulink test model A62_MPC_Sim.slx
%  (no LN blocks -> runs on any PC without the dongle)
mpcverbosity('off');
[mpcCore, mpcState0, mpcOnline] = getCodeGenerationData(mpcD, 'DataType', 'single');
mpcverbosity('on');
zRefD  = LcD(:,1)';          % reference of the 4 terminal outputs per metre of x_ref
nMPC   = round(Ts*fISR);     % MPC runs every nMPC-th Control-Routine call (5 at 500 Hz)
Tisr   = 1/fISR;             % Control-Routine sample time (125 us at 8 kHz, 1 ms at 1 kHz)
Tplant = 1/8000;             % plant simulation step in A62_MPC_Sim (fixed, independent of fISR)
% plant model (same equations as lucasnuelle/Simulation/Position and Angle Control System)
fr     = 0.3;                % pendulum friction of the LN plant block
qPos   = 0.054/16000;        % encoder resolution [m]      (3.375 um)
qAng   = 2*pi/65535;         % CAN-angle resolution [rad]
theta0 = pi + 1.5*pi/180;    % plant start angle: pi = upright, 0 = hanging
useMPC = 1;                  % 1 = MPC (candidate D), 0 = course LQR (K, Nbar) in A62_MPC_Sim
awOn   = 1;                  % MPC block (MPC_D_aw.m): 1 = v_set anti-windup on
spStepTime = 3;              % Set_Pos step time [s]  (use ~12 s for a swing-up test: catch needs |x| < 30 mm)
spStepSize = 100;            % Set_Pos step size [mm]

%% Lab check (printed every run)
fprintf('\nLAB CHECK: finalist=%s  fISR=%d  nMPC=%d  Tisr=%g  NpD=%d  W rows=%d  last row=%s  qD=%s  moves=%s  xLim=%g  data=%s  useMPC=%d  awOn=%d\n', ...
    finalist, fISR, nMPC, Tisr, mpcD.PredictionHorizon, size(mpcD.Weights.OV,1), mat2str(mpcD.Weights.OV(end,:)), ...
    mat2str(qD), mat2str(mpcD.ControlHorizon), mpcD.OV(1).Max(end), class(mpcState0.Plant), useMPC, awOn);
