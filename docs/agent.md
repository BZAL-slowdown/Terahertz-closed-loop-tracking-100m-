# 220 GHz THz Reproduction Checkpoint

Updated: 2026-11-02

## [System Status]

- Current workspace source: `太赫兹/main_THZ_220GHz_deterministic_mechanical_tracking.py`.
- This continuation corrected the capture-time semantics, made the joint-experiment registry an optional dependency, and added eight focused tests in `太赫兹/test_capture_window.py`.
- Verification: 8/8 unit tests passed; both Python files passed `py_compile`; the representative simulation completed at 10 s and at the configured 100 s.
- Full representative output is isolated in `太赫兹/outputs_capture_semantics_full/representative/`. It reports first 3 dB entry `0.00 s`, sustained-window start `16.62 s`, sustained-window confirmation `21.43 s`, `convergence_success=1`, `tracking_pass=1`, and accepted-update fraction `1.0`. The 10 s smoke run is under `outputs_capture_semantics_smoke/`; it correctly has no confirmed sustained lock.
- The separate Stage 1-2 100 m baseline is `outputs_stage1_2_100m_optimized_recheck/stage1_2_summary.json`: final error norm `0.054135 deg`, final-20-step RMS `0.043232 deg`, IEKF accepted `200/200`, minimum covariance eigenvalue `4.7381e-6`, Jacobian check maximum error `2.2918e-9`, peak SNR `67.0454 dB`, and P.676 gas attenuation `0.24868 dB`. This is a distinct script/result and must not be conflated with the 100 s representative result above.
- The base workspace `D:\太赫兹跟瞄` is not itself a Git repository. The separate organized repository exists at `D:\太赫兹跟瞄_副本_20261028\Terahertz-closed-loop-tracking`, clean at commit `f435edf` (`Add organized MATLAB stage 3-4 THz reproduction`). No files were changed there.
- The other repository has a dedicated `code/stage3_4_thz.py`, MATLAB Stage 3-4 implementation, static/dynamic CSV+JSON outputs, figures and Word reports. Its `agent.md` records historical results; this continuation only inspected those artifacts and did not rerun Stage 3-4.
- Stage 3-4 comparison notes: its passive branch freezes states/covariance blocks with a `1e-12` assertion. Its radius law actually uses `2.2 * sqrt(mean(recent lambda_max))` in THz units; `K_DELTA_AUTHOR=1/(7000*NS)` is retained as a reference constant rather than used directly by that law. `FEEDBACK_GAIN` is set to `1.0`, while nearby comments discuss risks and mention `0.6`; the implementation/comment mismatch should be reconciled before relying on that explanation. These are findings from source inspection, not a fresh benchmark comparison.
- The other repository's main tracking script still labels the first one-sample 3 dB entry as `capture_time_to_3db_lock_s`; it does separately compute a continuous convergence window. The current workspace now distinguishes these meanings.

## [Mathematical Anchors]

- Six-state order: `x=[q, q_dot, e_H, e_E, omega_H, omega_E]^T`.
- Measurement model: `z_i=q-q_dot*tau_i+g(e-omega*tau_i-s_i)+v_i`; H/V antenna cuts are based on measured patterns with PCHIP interpolation. Link-budget outputs use free-space loss plus ITU-R P.676-12 oxygen and water-vapour attenuation.
- A 3 dB lock sample is `pointing_loss_db <= 3`. The configured sustained dwell is `convergence_consecutive_cycles * mean_terminal_update_interval`; sample count is `ceil(dwell/dt)+1`, so the first-to-last sample timestamps span the required time. The zero-duration initialization snapshot is excluded.
- Capture fields: `capture_time_first_3db_entry_s` is diagnostic; `capture_time_sustained_lock_start_s` is the start of the first qualifying uninterrupted window; `capture_time_sustained_lock_confirmed_s` and legacy `capture_time_to_3db_lock_s` are the last sample timestamp that confirms it. `convergence_time_s` remains a backward-compatible alias for the window start.
- `_first_true_window_start()` finds the first uninterrupted boolean window. `_dwell_window_sample_count()` converts physical dwell time to a correctly elapsed sample span and validates its inputs.
- The representative simulation's 21.43 s sustained confirmation is a local simulation metric, not a hardware lock certification. The receiver noise figure, detector noise assumptions, and simulated measurement models are not receiver-calibrated.

## [Remaining Todo & Traps]

- Reconcile the inspected Stage 3-4 script's `FEEDBACK_GAIN=1.0` with its nearby explanatory comments, then rerun its static/dynamic cases and regenerate its figures/report before citing those outputs as current.
- Confirm the standalone Stage 3-4 MATLAB/Python results in the organized repository against the actual CSV/JSON values; do not substitute its README/agent summary for artifact inspection.
- Calibrate receiver noise/readout statistics versus input power, integration time, bandwidth, AGC and detector mode before making physical claims from the Gaussian or finite-sample/Gamma simulation assumptions.
- The `10 dB` Delta-method SNR validity condition applies only to the approximation; it is not a validity gate for the raw Gamma EnergyMAP path.
- Capture/dwell success remains a local screening definition. Extend it with innovation consistency and no-divergence dwell criteria, and sample more initial-error phases around observed failure boundaries.
- Do not claim Stage 3-4 was revalidated in this continuation. The full Stage 3-4 speed sweep and a high-frequency inertial/IMU feed-forward assessment remain open.
