# Handover: MPC tuning vs the course LQR (2026-11-01)

**Read first:** `claude/MPC_LN_Pendulum_Context.md`, `claude/Recovery_Handover.md`, `claude/A62_MPC_lab_checklist.md`.
This file covers the systematic tuning in simulation (`experiments.m`) and the hardware plan for the finalists.

New files in the repo `MPC-On-Lucas-Nulle-Inverted-Pendulum-` (written 2026-11-01, **not committed yet**):

| File | Purpose |
|---|---|
| `experiments.m` | Staged search (verify → baseline → structure → weights → refine → robustness → report). Home PC, no Simulink, no LN blocks. Results in `experiments_results.mat`, log in `experiments_log.txt`. |
| `buildMPC.m` | Builds the toolbox MPC from one `cfg` struct (Ts, Np, blocking, q, soft x/ẋ limits, ECR, ΔMV weight). Same design as `Calculations.m` Step 5, R2022a-compatible. `Calculations.m` should call it so the lab builds exactly what was tuned. |
| `MPC_D_aw.m` | New code for the MATLAB Function block "MPC (candidate D)": v_set anti-windup + qpIter outputs. |

---

## 1. Method (and why staged, not one big grid)

- **Simulator** = replica of `A62_MPC_Sim`: nonlinear pendulum, drive PI with **±0.5 m/s and ±10 m/s²**, the course **v_set integrator (reset to 0 at the catch, no limits)**, the LN observer with its own drive model, quantisation, 500 Hz / 1 kHz Control-Routine, MPC held between solves.
- **Fast QP** (`mpcActiveSetSolver`, condensed) instead of `mpcmoveCodeGeneration` (7 ms/call interpreted). **Stage 0 checks it against the toolbox:** max |Δu| ≤ 1e-10 m/s², with active soft limits, blocking and ΔMV weight.
- **Scenarios (16):** 12 catches (ξ = −0.6/−1.2/−1.7, cart at 0 or ±0.25 m/s, plus 3 "falling away"), Set_Pos 0→100→0, pushes 1.2/1.5/1.8 rad/s.
- **Score J** relative to the simulated course LQR (LQR = 1): 0.5·catches + 0.25·step + 0.25·pushes + smoothness penalty + 1 per scenario with peak |x| > 140 mm.
- Stages: structure first (changes QP and memory), then weights on the best structures, then robustness (l ±20 %, 0.2° offset, 10 ms delay).
- **Why staged:** a full grid is thousands of runs and mixes questions. The first structural sweep showed that **all** structures dropped the same scenarios, which pointed to a model problem (v_set windup, below) that no weight grid could have fixed.

## 2. Findings

1. **v_set integrator windup (main new cause).** The course `Discrete-Time Integrator` (a → v_set) has no limits; the drive caps at 0.5 m/s. In a hard push, the MPC kept commanding a ≈ −3 and **v_set wound up to −1.3 m/s** (LQR: −1.07). When braking was needed, v_set first had to unwind → cart into the 150 mm limit. **Fix inside the MPC block:** track v_set (the block knows its own output and the reset) and clip `a` so that v_set stays within ±0.5 (`MPC_D_aw.m`).
2. **±10 chattering at active limits** reproduced in simulation (as in plot10/plot11): the prediction model has no drive lag (a → v_set → drive PI ≈ 10 ms). With anti-windup it mostly disappears; a ΔMV weight (du = 10) damps it further but **costs the 1.5 rad/s push** → not used.
3. **Soft x limit needs a long enough horizon:** ≥ 0.15 s at Ts = 5 ms or 0.3 s at 10 ms (with blocking). At Ts = 10 ms / Np = 20 it does nothing.
4. **ẋ ≤ 0.5 soft limit is unnecessary once anti-windup is in** (identical results). Observer speed clamp ±0.5: no effect in simulation.
5. **Single precision is not a problem:** u_single − u_double ≤ 0.003 m/s² over 200 catch-like states (hypothesis rejected).
6. **Blocking `[1 1 1 2 3 4 8]` alone** (instead of Nc = 5) changes little; it matters together with longer horizons.
7. **Not reproduced:** the drop of the specific plot11 catch (14°, −120°/s). The limited MPC holds it in simulation. The chattering mechanism is reproduced, the drop is not.
8. **5 catches + the 1.8 rad/s push are dropped by every controller** (cart already at speed, no reserve to lean the pendulum back within 150 mm). They only count as "holds what the LQR can't".
9. **LQR at 1 kHz = LQR at 500 Hz** (scores identical), so changing fISR doesn't favour either controller.

## 3. Finalists (simulation, nominal plant)

Structure for all: soft x limit ±0.13 m (ECR 1), no ẋ limit, MV ±10 hard, terminal cost, **anti-windup on**.

| | Weights q = [x v φ] | Ts / fISR | Np, blocking | Est. memory |
|---|---|---|---|---|
| **F1** | [5000 100 1000] | 5 ms / **1 kHz** | 30, `[1 1 2 2 4 6 14]` | ≈ 9k words |
| **F2** | [10000 300 1000] | 5 ms / **1 kHz** | 30, `[1 1 2 2 4 6 14]` | ≈ 9k words |
| **F3** | [5000 100 1000] | 10 ms / 500 Hz | 30, `[1 1 1 2 3 4 8 10]` | ≈ 8.9k words |

Memory: linear fit through the two known builds (7.0k words now; Np = 50 → 14.4k, overflowed above ≈ 12.4k). F1/F2 with Np = 40 would be ≈ 11.1k (tight); Np = 30 performs the same.

| Controller | catch peak x | catch \|φ\|<1° | catch \|x\|<5 mm | step ts(5 mm) | step over+under | push 1.2 peak / settle | push 1.5 | aRMS |
|---|---|---|---|---|---|---|---|---|
| Course LQR | 67 mm | 0.51 s | 1.24 s | 1.12 s | 11.9 mm | 90 mm / 1.36 s | 149 mm, 1.25 s | 1.11 |
| **F1** | 67 | 0.59 | **0.96** | **0.96** | 14.9 | 94 / **1.11** | **133 mm**, 1.72 s | 1.26 |
| **F2** | 67 | 0.57 | **0.93** | **0.96** | 15.9 | 94 / **1.10** | **136 mm**, 1.55 s | 1.37 |
| **F3** | 67 | 0.59 | **0.96** | **0.96** | 15.1 | 94 / **1.11** | **138 mm**, 1.77 s | 1.37 |
| Lab MPC now (D, Nc5, no AW) | 71 | 0.56 | 0.91 | 0.98 | 14.6 | 101 / 1.03 | **drops** | 1.58 |

**Where the MPC wins:** x returns 20–25 % faster after catches, small pushes and steps; under a hard push it stops 11–16 mm further from the drop limit (LQR: 1 mm margin).
**Where it loses:** φ settles ~0.07 s later, step over/undershoot +3–4 mm, 15–25 % more control activity, and after a hard push it parks at the soft limit for 0.3–0.6 s (slower return).
**Robustness** (l ±20 %, 0.2° offset, 10 ms delay): the finalists never drop more than the LQR, and hold more in three of the four cases.

## 4. Hardware plan (lab PC, R2022a)

1. Copy `buildMPC.m` next to `Calculations.m`. In Step 5 replace the hand-built `mpcD` with:
   ```matlab
   cfgD = struct('Ts',0.005,'Np',30,'moves',[1 1 2 2 4 6 14],'q',[5000 100 1000],'r',1, ...
                 'xLim',0.13,'xFrom',1,'vLim',Inf,'ecr',1,'du',0);           % F1
   fISR = 1000;                                                              % F1/F2: 1 kHz (F3: 500, Ts 0.01, Np 30, moves [1 1 1 2 3 4 8 10])
   [mpcD, LcD] = buildMPC(cfgD, A, B);  Ts = cfgD.Ts;  NpD = cfgD.Np;  qD = cfgD.q;
   [mpcCore, mpcState0, mpcOnline] = getCodeGenerationData(mpcD,'DataType','single');
   zRefD = LcD(:,1)';  nMPC = round(Ts*fISR);  Tisr = 1/fISR;  awOn = 1;
   ```
   and set the LucasNuelleHardware block's Control-Routine frequency to `fISR`. The LAB CHECK line must show `W rows=30`, the last row `[0 0 0 1 1 1 1]`, `nMPC=5`.
2. MPC block: paste `MPC_D_aw.m`; add the parameters `Tisr` and `awOn` (Scope = Parameter, non-tunable). Outputs stay `a, qpIter, qpIterMax`.
3. Build: watch the **linker** (`.ebss`/`.econst` in RAML4L7) and the **Processor Load** display at 1 kHz.
4. Fallback if memory/time fails: F3 (500 Hz, nothing else changes).
5. **Comparison protocol (same day, after homing):** Set_Pos 0→100→0, 3 taps, 3 catches. Record Balance, Pos, Phi, y1, a_MPC, a_LQR; note ξ = φ̇ + 6.68·φ for each catch. Order: LQR (useMPC = 0) → F1 → F2 → F3 → LQR again.
6. **What to compare:** catch peak x, t(|φ| < 1°), t(|x| < 5 mm), step settling and dip, push peak x. Expected from simulation: MPC faster on x, more margin on hard pushes, slightly more φ activity.

## 5. Open points

- Run `experiments.m` fully on the home PC and compare its report with the table above (the Python replica gives identical LQR numbers).
- If F1/F2 chatter at the x limit on hardware: try `xLim 0.14`, or `du` 3–5, or an x limit on later steps only (`xFrom`).
- Anti-windup would also help the course LQR (integrator limits ±0.5). It's a course-model change; keep the LQR as built for the comparison and mention it in the thesis.
- Commit the new files.
