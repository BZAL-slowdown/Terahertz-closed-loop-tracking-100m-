# 220 GHz Terahertz Closed Loop Tracking

This repository contains the current Stage 1-2 reproduction work for a 220 GHz active mechanical pointing loop. The work starts from measured H/V antenna cuts, builds a 220 GHz link budget, then runs a local six-state IEKF and PI control loop. It also includes the follow-up noise, range, capture-boundary, and capture-time checks used to understand where the present model is useful and where it is not.

The repository is deliberately scoped. Stage 3-4 results are not included as current validated results because they were not rerun in this workspace.

## What is here

- `src/`: the three simulation entry points and focused capture-window tests.
- `data/`: the measured antenna-cut workbook used by the simulation.
- `results/`: compact summaries and figures for the 100 m baseline, supplementary experiments, and the representative capture-time run. Full step-by-step traces are intentionally excluded because they are large and can be regenerated.
- `docs/`: the technical notes and the latest Chinese Word report.

## Main results

The 100 m Stage 1-2 baseline is a single fixed-seed, static, one-end local-loop simulation. It runs for 100 s at 0.5 s intervals. The final pointing-error norm is 0.054135 degrees, and the RMS over the final 20 samples is 0.043232 degrees. All 200 IEKF batches were accepted; the minimum covariance eigenvalue was `4.7381e-6`; the analytic PCHIP Jacobian differed from the central finite-difference check by at most `2.2918e-9`.

At 100 m, the measured-cut model gives H/V peak gains of 44.955/44.691 dBi and H/V 3 dB beamwidths of 0.7026/0.8328 degrees. The link budget gives 119.296 dB free-space loss, 0.249 dB P.676-12 gaseous attenuation, and 67.05 dB peak SNR under the current assumed receiver parameters.

These are simulation results, not antenna pointing measurements, communications-range certification, or hardware lock certification. In particular, the current receiver noise, bandwidth, detector, AGC, and quantization assumptions have not been calibrated against hardware.

## Reproduce

Tested with Python 3.11. Install the Python dependencies:

```bash
python -m pip install -r requirements.txt
```

Run the 100 m baseline from the repository root:

```bash
python src/stage1_2_thz.py --distance-m 100 --output-dir results/stage1_2_100m
```

This reads `data/LWB_CRS-220-100MM-2406-2.xlsx` and writes a JSON summary, a time-series CSV, and two figures.

Run the supplementary sensitivity experiments:

```bash
python src/stage1_2_supplementary_experiments.py
```

The supplementary driver writes its configured outputs to `results/supplementary/`. It compares fixed-noise and SNR-dependent noise assumptions, performs paired-seed noise trials, and scans local initial offsets. These experiments take longer than the deterministic baseline.

Check the capture-window semantics:

```bash
python -m unittest src/test_capture_window.py
```

Run the representative continuous-tracking case:

```bash
python src/main_THZ_220GHz_deterministic_mechanical_tracking.py --experiment representative --runs 1 --output results/capture_semantics
```

The representative run distinguishes a first sample entering the 3 dB region from a sustained-lock confirmation. The saved representative result reports first entry at 0.00 s, sustained-window start at 16.62 s, and sustained confirmation at 21.43 s. The confirmation time is a simulation metric defined by the configured dwell rule, not a hardware result.

## Reading the supplementary results

The noise sweep retains outliers. At 0.1, 0.3, and 0.6 dB assumed measurement noise, the local screening pass counts were 19/20, 18/20, and 19/20 respectively. A batch update being accepted does not by itself establish closed-loop convergence.

The range sweep is mainly a model-sensitivity exercise. A fixed 0.1 dB measurement-noise model makes range look artificially benign because the noise does not change with received power. With the SNR-driven delta-method assumption, the 5 km case has 20.88 dB peak SNR and 10/10 local screening passes, while the 10 km and 20 km cases have 2.43 dB and -28.46 dB peak SNR and 0/10 passes. The low-SNR points are outside the approximation's intended domain and must not be reported as physical maximum-range predictions.

The initial-offset stress sweep is non-monotonic. For example, 0.50 HPBW yields 3/20 passes, while 0.65 and 1.00 HPBW yield 20/20 for the particular tested offset direction and seeds. It is therefore not valid to claim a general capture radius from this sweep.

## Limits and next work

The present Stage 1-2 work is static, single-ended, and locally initialized. It does not validate bidirectional time slots, passive-state freezing, moving platforms, platform vibration, high-rate scan mismatch, or receiver-calibrated observation statistics. The next practical step is to calibrate receiver/readout noise versus input power and integration settings, then rerun the range and capture studies with that measurement model. After that, the bidirectional state machine and motion cases can be checked on the same evidence standard.

## Files not committed

Generated per-step traces, Python caches, Word lock files, and local rendering QA files are ignored. They are reproducible intermediates rather than source artifacts.
