# 802.11a/g L-STF Start-of-Packet Detector Evaluation

A self-contained MATLAB framework (zero toolbox dependencies; base MATLAB only)
that evaluates three L-STF Start-of-Packet (SoP) detectors on a realistic
**end-to-end RF + ADC** receive chain under AWGN and 3GPP TDL-D Rician
multipath, with four adversarial jammer families:

| Detector | Description |
|---|---|
| `AC_BASELINE`    | Schmidl–Cox lag-`L` autocorrelator + dual-CFAR gate (relative + absolute) |
| `XC_SINGLE_MF`   | Single 16-tap matched filter at zero CFO + 4-tap stride-`L` coherent comb + dual-CFAR gate |
| `XC_BANK` (M=3)  | 3-hypothesis bank `{-145, 0, +145}` kHz covering the IEEE ±233 kHz CFO envelope, with `max_k |v_k|²` + dual-CFAR gate |

Each detector runs with **calibrated thresholds frozen in the runner**
(joint Pd–Pfa grid-search on AWGN, transfers to TDL-D without re-tuning at
K ≥ 10 dB). The CFAR pre-block runs both **mean** and **median** floor
estimators in parallel so the report covers both.

---

## Why this matters

Receivers that re-use the 802.11a/g legacy preamble — including tactical
UAV / mesh OFDM PHYs that share the L-STF — inherit the same SoP stage.
This study quantifies how each detector behaves under four adversarial
jammer families (random-carrier CW, random-carrier 1 MHz filtered noise,
wideband Gaussian, linear chirp) at realistic UAV CFO envelopes
(±233 kHz at 5.825 GHz, 40 ppm combined TX+RX) and Rician multipath
(TDL-D, K = 13.3 dB).

The headline finding: **the autocorrelator alone is uniquely vulnerable
to in-band 1 MHz filtered noise** (Pfa ≈ 0.36 in AWGN, 0.37 in TDL-D)
while the XC bank holds Pfa ≤ 0.042 across every jammer family in both
channels.

---

## Repo layout

```
802.11a_sop_det_eval/
├── eval_80211_sop_e2e_run.m   % SINGLE entry point — run this in MATLAB
├── eval/
│   └── eval_80211_sop_kernel.m  % Per-trial CRN-optimised kernel (5 sweeps)
├── config/
│   └── sim_params.m             % All knobs in one place
├── phy_tx/                       % TX side
│   ├── tx_frame.m
│   ├── make_lstf.m              % 802.11a L-STF (10× tiled 16-sample STS)
│   └── make_lltf.m              % 802.11a L-LTF (waveform realism)
├── rffe_tx/
│   └── resample_int.m            % Integer-ratio resampler
├── channel/                      % Wireless channel + jammer models
│   ├── channel_tdl_d_gen.m      % 3GPP TR 38.901 TDL-D tap generator (Jakes-FFT)
│   ├── channel_tdl_d_apply.m    % Apply the cached realisation
│   └── make_jammer.m            % CW / 1 MHz noise / wideband / chirp jammers
├── rffe_rx/                      % RX analog/mixed-signal model
│   ├── apply_rf_adc.m           % Combined RX gain + ADC quantize + clip
│   ├── decimate_hb_2to1.m       % First-principles half-band 2:1 decimator
│   ├── design_lpf.m              % Windowed-sinc LPF (first-principles FIR)
│   └── fractional_delay.m       % Kaiser-windowed sinc fractional delay
├── phy_rx/                       % SoP detection — the methods under test
│   ├── ac_dual_stat.m           % Schmidl–Cox autocorrelator (T_abs, T_rel)
│   ├── xc_dual_stat.m           % Matched-filter cross-correlator
│   ├── rx_lstf_xcorr_bank.m     % Frequency-shifted MF bank (M hypotheses)
│   ├── estimate_floor.m         % CFAR floor (mean + median, parallel)
│   ├── detector_passes_window.m % Dual-gate detection event
│   └── compute_stf_window_offsets.m % L-STF inclusion window
├── helpers/
│   └── plot_eval_results.m       % Figure rendering from results.mat
└── results/
    ├── sop_e2e_baseline/         % AWGN production: results.mat + 7 PNGs
    └── sop_e2e_baseline_tdl/     % TDL-D production: results.mat + 7 PNGs
```

---

## How to run

In MATLAB R2018+ (no toolboxes required):

```matlab
cd 802.11a_sop_det_eval
eval_80211_sop_e2e_run    % AWGN by default; toggle cfg.channel_id = 1 for TDL-D
```

The runner self-pathes its 8 domain folders, builds the TX preamble at
`fs_tx = 80 MHz`, calls the per-trial kernel for 5 sweeps × `n_trials`
realisations, and writes `results.mat` + a 7-PNG figure set to
`results/sop_e2e_baseline/` (AWGN) or `results/sop_e2e_baseline_tdl/`
(TDL-D).

Default Monte Carlo budget is **5000 trials per sweep point** — clean
confidence intervals at every operating point. The shipped
`results/sop_e2e_baseline*/results.mat` are production-grade 5000-trial
runs; you can drop the budget by editing `cfg.n_trials_pd` / `cfg.n_trials_pfa`
at the top of the runner for fast smoke tests.

---

## Calibrated thresholds (frozen in the runner)

```matlab
cfg.alpha_ac           = 0.50;   % AC relative gate
cfg.K_abs_ac           = 4.40;   % AC absolute gate (units of √L·σ²)
cfg.beta_sq_single     = 0.020;  % XC single relative gate (β²)
cfg.K_abs_xc_single    = 3.90;   % XC single absolute gate (units of σ)
cfg.beta_sq_bank       = 0.040;  % XC bank relative gate
cfg.K_abs_xc_bank      = 4.70;   % XC bank absolute gate
```

These were tuned offline by a joint grid search over a 100 k-pool of
1-SIFS listen-window false-alarm samples (per-sample Pfa < 9.4 × 10⁻⁸,
λ_FA < 1.9 FA/sec) at the strictest Pfa = 0 target, then verified to
transfer cleanly to TDL-D at K ≥ 10 dB.

The calibration scripts themselves are **not shipped** (`calibration/`
stays local); the thresholds above are the frozen production values.
No re-calibration is required for AWGN or TDL-D LOS-dominated channels.

---

## Architecture

```
TX (80 MHz)            CHANNEL                RFFE_RX (80 → 20 MHz)        DETECTOR (20 MHz)
┌──────────┐    ┌───────────────────────┐    ┌─────────────────────┐    ┌──────────────────────┐
│ L-STF +  │───▶│ × P_in (SNR), CFO,    │───▶│ apply_rf_adc:       │───▶│ AC: Schmidl–Cox      │
│ L-LTF    │    │ fractional delay,     │    │   80 dB RX gain     │    │ XC_s: 16-tap MF +   │
│ tx_frame │    │ TDL-D [optional],     │    │   12-bit ADC clip   │    │       4-tap comb    │
└──────────┘    │ + AWGN + jammer       │    │ + 2-stage half-band │    │ XC_b: 3-hyp bank    │
                └───────────────────────┘    │   decim (11+23)     │    │                      │
                                              └─────────────────────┘    │ + CFAR floor (mean   │
                                                                          │   AND median)        │
                                                                          │ + dual-gate decision │
                                                                          └──────────────────────┘
```

* **Common Random Numbers (CRN)** restructuring: per-trial channel,
  noise, and jammer-base waveform drawn ONCE per trial and reused across
  every SNR / JNR / detector cell of the 5 sweeps — drops TDL-D
  wall-clock by ~10×.
* **Both CFAR estimators in parallel** (one-shot mean over the listen
  buffer, and median-via-histogram tracker); the trailing dim-2 axis of
  every `arr_*_sweepK` is `{mean, median}`.

---

## Sweeps

| # | Sweep | Detectors per cell |
|---|---|---|
| 0 | Pd vs SNR (random CFO, no jammer) | AC, XC_s, XC_b |
| 1 | Pd vs CFO, 4 panels (incl. CW jammer panel) | AC, XC_s, XC_b |
| 2 | Pd vs CFO at on-comb / off-comb CW jammer offsets | AC, XC_s, XC_b |
| 3 | Pd vs residual JNR for 4 jammer types (CW, 1 MHz noise, wideband, chirp) | AC, XC_s, XC_b |
| 4 | Pd vs 1 MHz-noise jammer carrier offset | XC_s, XC_b |

The 7 PNGs in each `results/sop_e2e_baseline*/` correspond to the 5
sweeps (sweep 3 fans out to two SNR panels = 2 figures, hence 7).

---

## Algorithmic notes

* **AC's traditional role** is the *coarse CFO estimator* in deployed
  802.11 receivers (the argument of `M[n]` recovers an unbiased CFO
  estimate unambiguous over ±`fs_rx`/(2L) ≈ ±625 kHz). This study asks
  the orthogonal question of whether AC is *also* adequate as the SoP
  *trigger* under realistic jamming — it is not. AC and the matched-filter
  detectors are not mutually exclusive in a real receiver.

* **XC bank vs single**: the bank covers the full ±233 kHz CFO envelope
  by frequency-shifting the template into 3 hypotheses pre-baked in
  ROM; the single-MF lobe is ±78 kHz at 3 dB — useful as the comb-vs-MF
  baseline only, not a production detector on its own.

* **CFAR**: both mean and median estimators run in parallel; mean
  over-counts impulsive interference into σ̂² and shuts out valid
  detections, while median is robust. The shipped headline numbers use
  the median path; the mean path is available in every `arr_*_sweepK`
  for ablation.

* **No toolboxes**: half-band decimation (`decimate_hb_2to1`), windowed-sinc
  LPF (`design_lpf`), integer-ratio resampler (`resample_int`), Kaiser
  windows (`fractional_delay` internals), and TDL-D Jakes generation
  (`channel_tdl_d_gen`) are all first-principles implementations using
  only base MATLAB (`conv`, `fft`, `filter`, `besseli`, `randn`, etc.).

---

## What's **not** in this repo

These directories are kept in the author's local but **excluded from the
public release**; they are not on the runtime path and the script runs
without them.

| Omitted | Why |
|---|---|
| `calibration/` | Threshold optimisers (joint AC + XC grid search). The frozen production thresholds are already in the runner. |
| `diagnostics/` | Verification testbenches (SNR/JNR, decimator lead-in, threshold sensitivity, high-SNR failure isolation). |
| `doc/`         | LaTeX evaluation report (.tex + .pdf + figures). The full 50-page report is available on request. |
| `sandbox/`     | Archived legacy code (baseline runner predating the E2E pipeline, MEX experiments, smoke tests). Not part of the production tree. |

---

## Coding conventions

* `snake_case` only. Type-indicating suffixes mandatory: `_iq`, `_fd`,
  `_freq`, `_pwr`, `_db`, `_mag`, `_idx`, `_sym`, `_bit`, `_byte`. FIR
  coefficient vectors carry the `hir_` prefix; array constants carry `arr_`.
* Every function has a docstring header with Author / License / IEEE
  reference / Algorithm / Inputs (with sizes) / Outputs.
* Function name = filename. Domain prefixes: `tx_*`, `rx_*`,
  `channel_*`, plus single-word names inside the `phy_rx/`, `phy_tx/`,
  `rffe_rx/`, `rffe_tx/` folders.

---

## Author

Panos Alevizos &lt;bigpan27@gmail.com&gt;
AI assistance: Claude (Anthropic) — code drafting and verification.
© 2026 Panos Alevizos. Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
