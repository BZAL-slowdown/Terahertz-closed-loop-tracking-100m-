"""独立补充实验：距离/噪声敏感性、多随机种子统计与局部捕获范围。"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import stage1_2_thz as sim
import main_THZ_220GHz_deterministic_mechanical_tracking as thz


ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent / "results" / "supplementary"
SEEDS = tuple(range(20261101, 20261121))
DISTANCES_M = (10.0, 25.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 5000.0, 10000.0, 20000.0)
NOISE_FLOORS_DB = (0.10, 0.30, 0.60)
COMPARISON_SEEDS = tuple(range(20261121, 20261131))
CAPTURE_OFFSETS_HPBW = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.65, 0.80, 1.00)
EFFECTIVE_SAMPLES = 256
NOISE_SIGMA_CAP_DB = 3.0
DELTA_METHOD_MIN_SNR_DB = 10.0
LOCK_DWELL_STEPS = 10
LOCK_LOSS_THRESHOLD_DB = 3.0


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_case(pattern, link, distance, sigma_db, seed, offset=(0.12, -0.10),
             noise_model="fixed", effective_samples=EFFECTIVE_SAMPLES):
    rows = sim.stage2(
        pattern, link, distance_m=distance, seed=seed,
        measurement_sigma_db=sigma_db, noise_model=noise_model,
        effective_samples=effective_samples,
        noise_sigma_cap_db=NOISE_SIGMA_CAP_DB,
        delta_method_min_snr_db=DELTA_METHOD_MIN_SNR_DB,
        design_sigma_db=max(0.18, sigma_db),
        initial_offset_hpbw=offset,
    )
    errors = np.array([r["error_norm_deg"] for r in rows], dtype=float)
    pointing_loss = np.array([r["pointing_loss_db"] for r in rows], dtype=float)
    initial = math.hypot(offset[0] * pattern.hpbw(0), offset[1] * pattern.hpbw(1))
    final_rms = float(np.sqrt(np.mean(errors[-20:] ** 2)))
    accepted = int(sum(r["update_accepted"] for r in rows))
    min_eig = float(min(r["minimum_covariance_eigenvalue"] for r in rows))
    sigma_trace = np.array([r["measurement_sigma_db"] for r in rows], dtype=float)
    noise_model_valid = bool(all(r["noise_model_valid"] for r in rows))
    dwell = min(LOCK_DWELL_STEPS, len(rows))
    lock_window = pointing_loss[-dwell:]
    sustained_lock = bool(np.isfinite(lock_window).all() and np.all(lock_window <= LOCK_LOSS_THRESHOLD_DB))
    finite = bool(np.isfinite(errors).all() and np.isfinite(pointing_loss).all()
                  and np.isfinite(sigma_trace).all() and np.isfinite(min_eig))
    # 成功必须同时满足数值稳定、误差改善和末端持续位于 3 dB 主瓣区。
    # 这仍是静态局部实验判据，不等价于双端通信系统的捕获认证。
    success = (finite and noise_model_valid and accepted == len(rows)
               and final_rms < 0.8 * initial and sustained_lock)
    return rows, {
        "distance_m": distance, "noise_model": noise_model,
        "measurement_sigma_db": sigma_db, "noise_floor_db": sigma_db,
        "energy_log_sigma_db": (10.0 / math.log(10.0)
                                * math.sqrt(float(sim.polygamma(1, effective_samples)))
                                if noise_model == "gamma_energy" else math.nan),
        "minimum_sample_snr_db": float(min(r["minimum_sample_snr_db"] for r in rows)),
        "mean_scan_sigma_db": float(np.mean(sigma_trace)),
        "max_scan_sigma_db": float(np.max(sigma_trace)),
        "noise_model_valid": noise_model_valid,
        "effective_samples": effective_samples, "seed": seed,
        "initial_error_norm_deg": initial, "final_error_norm_deg": float(errors[-1]),
        "final_20_step_rms_deg": final_rms, "accepted_batches": accepted,
        "total_batches": len(rows), "acceptance_rate": accepted / len(rows),
        "minimum_covariance_eigenvalue": min_eig, "finite": finite,
        "sustained_lock": sustained_lock, "lock_dwell_steps": dwell,
        "lock_loss_threshold_db": LOCK_LOSS_THRESHOLD_DB,
        "success": bool(success), "initial_h_offset_hpbw": offset[0],
        "initial_v_offset_hpbw": offset[1],
    }


def noise_sigma_from_snr(snr_db: float, floor_db: float, effective_samples: float) -> float:
    """delta-method 功率估计噪声 + 明示噪声底；effective_samples 是敏感性假设而非硬件参数。"""
    snr_linear = 10.0 ** (snr_db / 10.0)
    thermal_estimate_sigma_db = 10.0 / math.log(10.0) / math.sqrt(effective_samples * snr_linear)
    return math.hypot(floor_db, thermal_estimate_sigma_db)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pattern = sim.WorkbookPattern(sim.WORKBOOK)
    link = thz.LinkBudget()

    # 先固定 100 m 只改变观测噪声，隔离统计稳定性；保留每次运行原始逐步轨迹。
    monte_carlo = []
    for sigma in NOISE_FLOORS_DB:
        for seed in SEEDS:
            rows, result = run_case(pattern, link, 100.0, sigma, seed)
            monte_carlo.append(result)
            sim._write_csv(OUT / f"trace_mc_sigma{sigma:.2f}_seed{seed}.csv", rows)
    write_csv(OUT / "monte_carlo_trials.csv", monte_carlo)
    monte_summary = []
    for sigma in NOISE_FLOORS_DB:
        group = [r for r in monte_carlo if r["measurement_sigma_db"] == sigma]
        values = np.array([r["final_20_step_rms_deg"] for r in group])
        monte_summary.append({
            "measurement_sigma_db": sigma, "trials": len(group),
            "successes": sum(r["success"] for r in group),
            "success_rate": float(np.mean([r["success"] for r in group])),
            "rms_mean_deg": float(np.mean(values)), "rms_std_deg": float(np.std(values, ddof=1)),
            "rms_p05_deg": float(np.quantile(values, 0.05)),
            "rms_median_deg": float(np.median(values)), "rms_p95_deg": float(np.quantile(values, 0.95)),
            "acceptance_rate_mean": float(np.mean([r["acceptance_rate"] for r in group])),
            "minimum_covariance_eigenvalue": float(min(r["minimum_covariance_eigenvalue"] for r in group)),
        })
    write_csv(OUT / "monte_carlo_summary.csv", monte_summary)

    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    distributions = [[r["final_20_step_rms_deg"] for r in monte_carlo
                      if r["measurement_sigma_db"] == sigma] for sigma in NOISE_FLOORS_DB]
    ax.boxplot(distributions, tick_labels=[f"{s:.1f}" for s in NOISE_FLOORS_DB], showfliers=True)
    ax.set(xlabel="Measurement noise sigma (dB)", ylabel="Final 20-step RMS (deg)",
           title="Paired-seed noise sensitivity (outliers retained)")
    ax.grid(axis="y", alpha=0.25); fig.tight_layout()
    fig.savefig(OUT / "monte_carlo_noise_sensitivity.png", dpi=180); plt.close(fig)

    # 距离扫描同时给出固定噪声对照与逐扫描点 SNR 变化的工程敏感性模型。
    # 后者在 stage2 内部按每个样本的接收功率重新计算 sigma，而不是只用峰值 SNR 生成常数。
    distance_rows = []
    distance_trials = []
    for distance in DISTANCES_M:
        peak_dbm = link.received_power_dbm(distance, pattern.h_peak, pattern.e_peak)
        snr_db = link.snr_db(peak_dbm)
        models = (("fixed_0.1dB", "fixed", 0.10),
                  ("snr_delta_method_floor_0.1dB", "snr_delta", 0.10))
        for model, noise_model, sigma in models:
            trial_rms = []
            trial_records = []
            for seed in SEEDS[:10]:
                rows, result = run_case(
                    pattern, link, distance, sigma, seed,
                    noise_model=noise_model, effective_samples=EFFECTIVE_SAMPLES,
                )
                sim._write_csv(OUT / f"trace_distance_{model}_d{distance:g}m_seed{seed}.csv", rows)
                result["noise_model"] = model
                result["peak_received_power_dbm"] = peak_dbm
                result["peak_snr_db"] = snr_db
                distance_trials.append(result)
                trial_records.append(result)
                trial_rms.append(result["final_20_step_rms_deg"])
            distance_rows.append({
                "distance_m": distance, "noise_model": model,
                "peak_received_power_dbm": peak_dbm, "peak_snr_db": snr_db,
                "noise_floor_db": sigma, "effective_samples": EFFECTIVE_SAMPLES,
                "delta_method_min_snr_db": DELTA_METHOD_MIN_SNR_DB,
                "noise_sigma_cap_db": NOISE_SIGMA_CAP_DB,
                "mean_scan_sigma_db": float(np.mean([r["mean_scan_sigma_db"] for r in trial_records])),
                "max_scan_sigma_db": float(np.max([r["max_scan_sigma_db"] for r in trial_records])),
                "trials": len(trial_rms),
                "success_rate": float(np.mean([r["success"] for r in trial_records])),
                "valid_model_trial_rate": float(np.mean([r["noise_model_valid"] for r in trial_records])),
                "lock_rate": float(np.mean([r["sustained_lock"] for r in trial_records])),
                "rms_mean_deg": float(np.mean(trial_rms)),
                "rms_std_deg": float(np.std(trial_rms, ddof=1)),
                "rms_p05_deg": float(np.quantile(trial_rms, 0.05)),
                "rms_p95_deg": float(np.quantile(trial_rms, 0.95)),
                "rms_max_deg": float(np.max(trial_rms)),
            })
    write_csv(OUT / "distance_sweep_summary.csv", distance_rows)
    write_csv(OUT / "distance_sweep_trials.csv", distance_trials)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for model in ("fixed_0.1dB", "snr_delta_method_floor_0.1dB"):
        group = [r for r in distance_rows if r["noise_model"] == model]
        x = np.array([r["distance_m"] for r in group])
        y = np.array([r["rms_mean_deg"] for r in group])
        sd = np.array([r["rms_std_deg"] for r in group])
        axes[0].errorbar(x, y, yerr=sd, marker="o", capsize=3, label=model)
    axes[0].set(xscale="log", xlabel="Distance (m)", ylabel="Final 20-step RMS (deg)", title="Range and noise sensitivity")
    axes[0].grid(alpha=0.25); axes[0].legend(fontsize=8)
    snr_rows = [r for r in distance_rows if r["noise_model"] == "fixed_0.1dB"]
    sigma_rows = [r for r in distance_rows if r["noise_model"] == "snr_delta_method_floor_0.1dB"]
    axes[1].plot([r["distance_m"] for r in snr_rows],
                 [r["peak_snr_db"] for r in snr_rows],
                 marker="s", color="#d17a22", label="Peak link SNR")
    axes[1].plot([r["distance_m"] for r in sigma_rows],
                 [r["mean_scan_sigma_db"] for r in sigma_rows],
                 marker="^", color="#8b4e75", label="Mean sample sigma")
    axes[1].set(xscale="log", xlabel="Distance (m)", ylabel="dB", title="Link SNR and sample-noise model")
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "distance_noise_sensitivity.png", dpi=180); plt.close(fig)

    # 同一随机种子配对比较固定 dB 高斯观测与有限样本 Gamma 能量似然。
    # 这只比较两种统计假设在当前链路参数下的行为，不代表接收机噪声已标定。
    model_comparison = []
    for seed in COMPARISON_SEEDS:
        for model, noise_model in (("fixed_gaussian_0.1dB", "fixed"),
                                   ("gamma_energy_M256", "gamma_energy")):
            rows, result = run_case(pattern, link, 100.0, 0.10, seed,
                                    noise_model=noise_model, effective_samples=EFFECTIVE_SAMPLES)
            result["comparison_model"] = model
            result["noise_model_valid"] = bool(
                result["noise_model_valid"] and result["minimum_sample_snr_db"] >= DELTA_METHOD_MIN_SNR_DB
            )
            result["success"] = bool(result["success"] and result["noise_model_valid"])
            sim._write_csv(OUT / f"trace_model_{model}_seed{seed}.csv", rows)
            model_comparison.append(result)
    write_csv(OUT / "noise_model_comparison_trials.csv", model_comparison)
    comparison_summary = []
    for model in ("fixed_gaussian_0.1dB", "gamma_energy_M256"):
        group = [r for r in model_comparison if r["comparison_model"] == model]
        values = np.array([r["final_20_step_rms_deg"] for r in group])
        comparison_summary.append({
            "model": model, "trials": len(group),
            "valid_domain_trials": sum(r["noise_model_valid"] for r in group),
            "successes": sum(r["success"] for r in group),
            "success_rate_all_trials": float(np.mean([r["success"] for r in group])),
            "rms_median_deg": float(np.median(values)),
            "rms_p95_deg": float(np.quantile(values, 0.95)),
            "rms_max_deg": float(np.max(values)),
            "acceptance_rate_mean": float(np.mean([r["acceptance_rate"] for r in group])),
        })
    write_csv(OUT / "noise_model_comparison_summary.csv", comparison_summary)
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    groups = [[r["final_20_step_rms_deg"] for r in model_comparison
               if r["comparison_model"] == model] for model in
              ("fixed_gaussian_0.1dB", "gamma_energy_M256")]
    ax.boxplot(groups, tick_labels=["Fixed Gaussian", "Gamma energy, M=256"], showfliers=True)
    ax.set(ylabel="Final 20-step RMS (deg)", title="Paired-seed observation-model comparison at 100 m")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout(); fig.savefig(OUT / "noise_model_comparison.png", dpi=180); plt.close(fig)

    # 沿 H 轴逐渐增大初始误差，固定 100 m/0.1 dB/20 个配对种子。
    capture_trials = []
    for multiplier in CAPTURE_OFFSETS_HPBW:
        for seed in SEEDS:
            rows, result = run_case(pattern, link, 100.0, 0.10, seed,
                                    offset=(multiplier, -0.8 * multiplier))
            sim._write_csv(OUT / f"trace_capture_{multiplier:.2f}hpbw_seed{seed}.csv", rows)
            result["offset_multiplier_hpbw"] = multiplier
            capture_trials.append(result)
    write_csv(OUT / "capture_boundary_trials.csv", capture_trials)
    capture_summary = []
    for multiplier in CAPTURE_OFFSETS_HPBW:
        group = [r for r in capture_trials if r["offset_multiplier_hpbw"] == multiplier]
        rms = np.array([r["final_20_step_rms_deg"] for r in group])
        capture_summary.append({
            "offset_multiplier_hpbw": multiplier, "trials": len(group),
            "successes": sum(r["success"] for r in group),
            "success_rate": float(np.mean([r["success"] for r in group])),
            "rms_mean_deg": float(np.mean(rms)), "rms_std_deg": float(np.std(rms, ddof=1)),
            "rms_p95_deg": float(np.quantile(rms, 0.95)),
            "rms_max_deg": float(np.max(rms)),
        })
    write_csv(OUT / "capture_boundary_summary.csv", capture_summary)

    fig, ax = plt.subplots(figsize=(7.6, 4.5))
    ax.plot([r["offset_multiplier_hpbw"] for r in capture_summary],
            [r["success_rate"] for r in capture_summary], marker="o", color="#137c8b")
    ax.set(xlabel="Initial offset (HPBW multiplier per axis)", ylabel="Local-loop screening pass rate",
           ylim=(-0.05, 1.05), title="Initial-offset stress test (not a capture boundary)")
    ax.grid(alpha=0.25); fig.tight_layout(); fig.savefig(OUT / "capture_boundary.png", dpi=180); plt.close(fig)

    payload = {
        "scope": "Stage 1-2 only; Stage 3-4 not run",
        "noise_model_note": "Fixed Gaussian dB noise, per-sample delta-method noise, and Gamma finite-sample energy likelihood are alternative simulation assumptions. M=256, 0.1 dB floor, 3 dB cap, and 10 dB validity threshold are not receiver calibration. Gamma energy_log_sigma_db is the exact log-energy standard deviation from trigamma(M), recorded descriptively only; the Gamma estimator uses the raw Gamma likelihood, not this value as a Gaussian IEKF measurement covariance.",
        "effective_samples": EFFECTIVE_SAMPLES,
        "noise_sigma_cap_db": NOISE_SIGMA_CAP_DB,
        "delta_method_min_snr_db": DELTA_METHOD_MIN_SNR_DB,
        "distance_m": list(DISTANCES_M), "monte_carlo_seeds": list(SEEDS),
        "capture_offsets_hpbw": list(CAPTURE_OFFSETS_HPBW),
        "success_definition": f"all finite, noise model valid throughout, all filter batches accepted, final 20-step RMS < 0.8 * initial error norm, and the final {LOCK_DWELL_STEPS} steps remain within {LOCK_LOSS_THRESHOLD_DB} dB pointing loss; screening criterion only, not communications lock certification",
        "monte_carlo_summary": monte_summary,
        "distance_sweep_summary": distance_rows,
        "capture_boundary_summary": capture_summary,
        "noise_model_comparison_seeds": list(COMPARISON_SEEDS),
        "noise_model_comparison_summary": comparison_summary,
    }
    (OUT / "supplementary_experiments_summary.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Supplementary Stage 1-2 experiments complete: {OUT}")
    for row in monte_summary:
        print(f"MC sigma={row['measurement_sigma_db']:.2f} dB: success={row['successes']}/{row['trials']}, "
              f"RMS={row['rms_mean_deg']:.6f} +/- {row['rms_std_deg']:.6f} deg")
    print("Distance models: fixed 0.1 dB and SNR delta-method + 0.1 dB floor (assumption only)")
    print("Capture sweep: " + ", ".join(f"{r['offset_multiplier_hpbw']:.2f} HPBW={r['successes']}/{r['trials']}"
                                    for r in capture_summary))


if __name__ == "__main__":
    main()
