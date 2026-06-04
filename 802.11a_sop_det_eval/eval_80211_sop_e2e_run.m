%% EVAL_80211_SOP_E2E_RUN  802.11a/g SoP detector evaluation, end-to-end
%%                          variant: includes combined RX-chain gain and
%%                          ADC quantization + clipping in front of the
%%                          detector. Adds the RF + ADC modelling block
%%                          (\S "RF front-end + ADC") over the pure-IQ
%%                          chain shared with the per-trial kernel.
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
%   IEEE Ref: derived from IEEE Std 802.11-2016 Sec. 17.3.3 (OFDM PHY
%   preamble) and IEEE Std 802.11-2020 Sec. 17.3.9.7 (CFO worst-case
%   tolerance: +/- 20 ppm TX + 20 ppm RX = 40 ppm combined at the
%   U-NII-3 top edge 5.825 GHz -> +/- 233 kHz).
%
%   Driver responsibilities:
%     1. Populate cfg with EVERY knob the kernel will read (sample
%        rates, sweep grids, thresholds, Monte-Carlo trial counts,
%        channel model parameters).
%     2. One-time setup: build the sim_params struct, synthesise the
%        unit-power TX preamble at fs_tx, build N_h = 1 and N_h = 3
%        sim_params variants for the two XC detectors, compute the
%        L-STF window offsets in rx_iq (and widen the trailing edge
%        by the TDL-D max delay when cfg.channel_id = 1).
%     3. Assemble the per-detector threshold matrix:
%             [alpha_ac,        K_abs_ac;
%              beta_sq_single,  K_abs_xc_single;
%              beta_sq_bank,    K_abs_xc_bank]
%     4. Call the CRN-optimised kernel eval_80211_sop_kernel(...).
%     5. Save results.mat and render every sweep to PNG. Output goes
%        to results/sop_e2e_baseline/      (cfg.channel_id = 0, AWGN), or
%             results/sop_e2e_baseline_tdl/ (cfg.channel_id = 1, TDL-D).
%
%   Detectors under test (calibrated jointly on AWGN, then frozen):
%     AC_BASELINE      alpha       = 0.50,   K_abs = 4.40
%     XC_SINGLE_MF     beta_sq     = 0.020,  K_abs = 3.90
%     XC_BANK (N_h=3)  beta_sq     = 0.040,  K_abs = 4.70
%   Sources: optimize_ac_dual_joint.m, optimize_xc_dual_joint.m
%   (joint grid search vs the AWGN Pd-Pfa front).
%
%   Channel selector:
%     cfg.channel_id = 0  -> AWGN  (no multipath; default)
%     cfg.channel_id = 1  -> TDL-D (3GPP TR 38.901 Table 7.7.2-4;
%                                   K = 13.3 dB, 13 taps, LOS at tap 0,
%                                   delay spread = cfg.tdl_delay_spread_s,
%                                   max Doppler = cfg.tdl_f_doppler_hz)
%
%   See also: eval/eval_80211_sop_kernel.m,
%             phy_rx/detector_passes_window.m, phy_rx/estimate_floor.m,
%             phy_rx/compute_stf_window_offsets.m,
%             channel/channel_tdl_d_gen.m, channel/channel_tdl_d_apply.m,
%             channel/make_jammer.m,
%             config/sim_params.m, phy_tx/tx_frame.m,
%             helpers/plot_eval_results.m.

clear; clc; close all;
% ---- Domain-folder path setup ---------------------------------------
% Production layout: each signal-chain stage lives in its own folder.
% The driver puts every domain on the MATLAB path so the kernel and the
% domain functions can call each other by plain function name (no
% package prefix). Order doesn't matter -- addpath collisions are
% impossible since every function name is unique across folders.
%
% Anchor every path to the SCRIPT'S directory via mfilename('fullpath')
% so the runner works regardless of the MATLAB current-directory pwd
% (otherwise running it from a parent or sibling repo triggers
% "Name is nonexistent or not a directory" warnings).
repo_root = fileparts(mfilename('fullpath'));
addpath(repo_root);                                     % runner
addpath(fullfile(repo_root, 'eval'));                   % runner + kernel split into eval/
addpath(fullfile(repo_root, 'config'));                 % sim_params
addpath(fullfile(repo_root, 'phy_tx'));                 % tx_frame, make_lstf, make_lltf
addpath(fullfile(repo_root, 'rffe_tx'));                % resample_int
addpath(fullfile(repo_root, 'channel'));                % channel_tdl_d_*, make_jammer
addpath(fullfile(repo_root, 'rffe_rx'));                % decimate_hb_2to1, fractional_delay, design_lpf
addpath(fullfile(repo_root, 'phy_rx'));                 % ac/xc_dual_stat, detector_passes_window,
                                                        % estimate_floor, rx_lstf_xcorr_bank,
                                                        % compute_stf_window_offsets
addpath(fullfile(repo_root, 'helpers'));                % plot_eval_results

%% ====================================================================
%% USER KNOBS (cfg)
%% Every field of cfg is read by eval_80211_sop_kernel; no kernel-internal
%% knobs exist. Group ordering below mirrors the physical pipeline:
%%   (1) RF front-end & noise floor
%%   (2) Detector smoothing & loading
%%   (3) Random impairment envelopes
%%   (4) Monte-Carlo sample budget
%%   (5) AWGN-calibrated thresholds
%%   (6) Channel selector
%%   (7) RNG control
%%   (8) Sweep grids (one block per sweep 0..4)
%% ====================================================================

% -------- (1) RF front-end & noise floor -----------------------------
cfg.fs_bb_freq  = 20e6;             % baseband / RX sample rate [Hz] (IEEE 802.11a)
cfg.fs_tx_freq  = 80e6;             % oversampled TX sample rate [Hz] (4x baseband)
cfg.fs_rx_freq  = 20e6;             % RX sample rate after halfband cascade [Hz]
cfg.nf_db       = 6;                % RX noise figure [dB] (typical SDR front-end)
cfg.pn_dbm      = -174 + 10*log10(cfg.fs_bb_freq) + cfg.nf_db;
                                    % thermal noise in 20 MHz at NF [dBm]
                                    % = kTB + NF   (=> ~ -95 dBm)
cfg.pn_pwr      = 10^((cfg.pn_dbm - 30) / 10);
                                    % thermal noise power, linear [W]
                                    % = reference for SNR/JNR conversions

% -------- (2) Detector smoothing & loading ---------------------------
cfg.L_lag_samples = 16;             % AC lag length L [samples @ fs_rx]
                                    % = one short-symbol period (STS) of L-STF
cfg.S_ac          = 16;             % AC post-correlator MA window [samples @ fs_rx]
                                    % latency-bounded to one STS (~0.8 us)
cfg.S_xc          = 16;             % XC post-correlator MA window [samples @ fs_rx]
cfg.eps_load      = 1e-15;          % regulariser added to (P + eps) inside the
                                    % rel-gate denominator to avoid div-by-zero
                                    % in the leading guard region

% -------- (3) Random impairment envelopes ----------------------------
cfg.random_cfo_max_freq = 233e3;    % +/- max CFO drawn per trial [Hz]
                                    % IEEE 802.11a worst case @ 5.825 GHz
                                    % (top of U-NII-3): 20 ppm TX + 20 ppm RX
                                    % = 40 ppm combined => +/- 233 kHz
cfg.chirp_f_start_freq  = -5e6;     % chirp jammer sweep start frequency [Hz]
cfg.chirp_f_end_freq    = +5e6;     % chirp jammer sweep end   frequency [Hz]

cfg.guard_min_samples_rx = 160;     % min pre-packet quiet interval [samples @ fs_rx]
                                    % = 8 us @ fs_rx = 20 MHz. When the
                                    % receiver is already in continuous Rx
                                    % mode, the protocol-minimum inter-frame
                                    % gap before the next L-STF can shrink
                                    % below SIFS (e.g., DCF EDCA TXOP
                                    % bursting, A-MPDU subframes, opportunistic
                                    % piggybacking). This 8-us floor is
                                    % also EXACTLY the number of samples the
                                    % CFAR one-shot estimator averages over
                                    % (cfg.cfar_n_samples = 160): the first
                                    % 160 samples of rx_iq are guaranteed
                                    % pre-packet noise+jammer in every trial.
cfg.guard_max_samples_rx = 320;     % max pre-packet quiet interval [samples @ fs_rx]
                                    % = 16 us = 1 SIFS @ 5 GHz (IEEE 802.11a
                                    % Sec. 17.3.4 short inter-frame space).
                                    % Cap above SIFS would push the L-STF
                                    % well outside the detector's CFAR-warm
                                    % region and is not representative of
                                    % a receiver that is actively listening.

% -------- (4) Monte-Carlo sample budget ------------------------------
cfg.n_trials_pd               = 100;
                                    % independent signal-frame trials per Pd cell
                                    % 5000 was used for the published curves;
                                    % 100 is a fast first-light setting for E2E
cfg.n_trials_pfa              = 100;
                                    % independent noise+jam-only trials per Pfa cell
cfg.n_pfa_observation_samples = 320;
                                    % Pfa observation window length [samples @ fs_rx]
                                    % = 16 us = 1 SIFS @ 5 GHz. Matches the
                                    % maximum pre-packet idle the receiver
                                    % expects (cfg.guard_max_samples_rx).
                                    % A false alarm anywhere in these 320
                                    % samples of noise commits the PHY chain
                                    % to a wrong timing reference.
cfg.n_pfa_cfar_pre_samples    = 160;
                                    % extra noise+jam samples reserved BEFORE the
                                    % Pfa window so the CFAR estimator window
                                    % is fully populated at the first Pfa
                                    % test sample (= 8 us = guard_min)
cfg.cfar_excl_lead_n          = 32;
                                    % leading samples dropped from the start of
                                    % rx_iq (halfband decimator transient)
cfg.cfar_excl_trail_n         = 32;
                                    % LEGACY: was the trailing slab cut from
                                    % the oracle one-shot estimator. The
                                    % streaming CFAR ignores this field; kept
                                    % only for back-compat with the cfg struct.
cfg.cfar_n_samples            = 160;

% -------- (3b) RF front-end + ADC model (E2E variant ONLY) -----------
% Models the analog RX chain (combined gain at MAX AGC setting) and
% the ADC (uniform mid-tread quantization + symmetric clipping at
% fullscale) in front of the detector. AGC is FIXED at maximum gain.
% This captures the dominant impairments in the sensitivity regime
% (per \S limitations of the report):
%   1. Quantization noise from the ADC LSB
%   2. ADC clipping under strong-signal / strong-jammer conditions
%   3. Combined RX-chain gain (LNA + mixer + IF VGA + baseband amp)
%   4. Constant LNA NF stays as cfg.nf_db = 6 (unchanged)
cfg.rx_total_gain_db   = 80;
                                    % Combined RX-chain gain at MAX AGC
                                    % setting [dB]. Typical UAV C/L-band
                                    % digital receiver:
                                    %   LNA (1st stage)         ~15-18 dB
                                    %   1st mixer / converter   ~10 dB
                                    %   IF amplifier + VGA      ~30-40 dB
                                    %   baseband amp + ch-sel   ~10-15 dB
                                    %   ----------------------------
                                    %   Total at max AGC        ~70-90 dB
                                    % 80 dB is a representative midpoint.
cfg.adc_n_bits         = 12;
                                    % ADC resolution [bits], signed
                                    % two's-complement (so 11 effective
                                    % bits per I/Q axis, LSB range
                                    % [-2048, +2047]).
                                    % Typical UAV digital receiver:
                                    % AD9361/AD9363, LMS7002M family.
                                    % apply_rf_adc.m returns the I/Q
                                    % samples as LSBs (integers); the
                                    % downstream pipeline operates
                                    % directly on LSB-valued samples.
cfg.adc_fullscale_dbm  = 10;
                                    % ADC fullscale sine-wave power [dBm]
                                    % into 50 Ohm reference impedance.
                                    % Typical 12-bit ADC, 1 Vpp diff:
                                    % FS = +10 dBm => 1 V peak amplitude
                                    % (used inside apply_rf_adc.m only to
                                    % set the LSB step on the analogue
                                    % side; the returned samples are
                                    % already in LSB units).
                                    % number of leading rx_iq samples
                                    % averaged for the one-shot CFAR floor
                                    % estimate [samples @ fs_rx]. 160 = 8 us
                                    % @ fs_rx = 20 MHz: guaranteed pure
                                    % noise+jammer because guard_min_samples_rx
                                    % >= 160. Matches the minimum inter-frame
                                    % quiet a receiver in continuous Rx mode
                                    % can expect; a real receiver computes
                                    % the same estimate once during this
                                    % quiet period and holds it for the next
                                    % listen interval. Report Sec. 9.3.4
                                    % discusses the adaptive streaming variant.

% -------- (5) AWGN-calibrated thresholds -----------------------------
% Joint (alpha or beta_sq, K_abs) values from grid search on the AWGN
% Pd-Pfa front. Source helper scripts (now archived in sandbox/):
%   sandbox/optimize_ac_dual_joint.m,
%   sandbox/optimize_xc_dual_joint.m
% These values are FROZEN across every sweep, every channel, every
% jammer scenario -- they are NOT re-calibrated per scenario (which
% would be cheating since the rel-gate would absorb the jammer it is
% trying to characterise). Per-trial CFAR adaptation handles in-band
% noise+jammer floor via estimate_floor() inside the kernel.
cfg.alpha_ac        = 0.50;         % AC rel-gate threshold (Cauchy-Schwarz ratio)
                                    % Re-calibrated under first-firing Pd +
                                    % 1-SIFS Pfa semantics (optimize_ac_refine_joint.m)
cfg.beta_sq_single  = 0.020;        % XC_SINGLE_MF rel-gate threshold (in [0, 1])
                                    % Re-calibrated (optimize_xc_dual_joint.m)
cfg.beta_sq_bank    = 0.040;        % XC_BANK rel-gate threshold (in [0, 1])
                                    % Re-calibrated (optimize_xc_dual_joint.m)
cfg.K_abs_ac        = 4.40;         % AC abs-gate multiplier:
                                    % abs_floor = K * sqrt(L) * sigma_sq_hat
                                    % Strict Pfa = 0 / 100k SIFS windows
                                    % at (alpha=0.50, K=4.40). SNR@Pd>=0.99
                                    % = -1.00 dB in AWGN with CFO U[-233,
                                    % +233] kHz and guard U[160, 320].
cfg.K_abs_xc_single = 3.90;         % XC_SINGLE_MF abs-gate multiplier:
                                    % abs_floor = K * sqrt(sigma_sq_hat)
                                    % Strict Pfa = 0 / 100k at (beta^2=0.020,
                                    % K=3.90). SNR@Pd>=0.99 = +3.25 dB.
cfg.K_abs_xc_bank   = 4.70;         % XC_BANK abs-gate multiplier:
                                    % Strict Pfa = 0 / 100k at (beta^2=0.040,
                                    % K=4.70). SNR@Pd>=0.99 = -1.75 dB --
                                    % the most sensitive of the three.

% -------- (6) Channel selector ---------------------------------------
% =====================================================================
% CHANNEL SELECTOR    (the single knob to flip between the two campaigns)
%   0 = AWGN                       -> results/sop_e2e_baseline/
%   1 = TDL-D (3GPP TR 38.901)     -> results/sop_e2e_baseline_tdl/
% =====================================================================
cfg.channel_id         = 0;         % 0=AWGN, 1=TDL-D
cfg.tdl_delay_spread_s = 30e-9;     % TDL-D delay spread [s]
                                    % 30 ns = UAV LOS air-to-ground regime
                                    % used only when channel_id == 1
cfg.tdl_f_doppler_hz   = 500;       % TDL-D max Doppler [Hz]
                                    % ~100 km/h relative motion @ 5.25 GHz
                                    % used only when channel_id == 1

% -------- (7) RNG control --------------------------------------------
cfg.rng_seed = -1;                  % -1 = let RNG advance naturally (production)
                                    % >=0 = seed Twister at kernel entry
                                    %       (used only in dev smoke tests)

% -------- (8) Sweep grids --------------------------------------------
% Sweep 0:  Pd vs SNR, random CFO ~ U[-233, +233] kHz, no jammer.
cfg.snr_grid_sweep0_db = -5 : 1 : 15;       % SNR axis [dB], 21 points

% Sweep 1:  Pd vs CFO, 4 panels (last panel adds a CW jammer at DC).
cfo_grid_inner_khz       = -200 : 40 : 200;             % dense interior [kHz]
cfo_grid_outer_khz       = [-400 -320 -240, 240 320 400]; % coarse exterior [kHz]
cfg.cfo_grid_sweep1_khz  = sort([cfo_grid_inner_khz, cfo_grid_outer_khz]);
                                                        % 17 CFO points [kHz]
cfg.snr_sweep1_panels_db    = [-5, 0, 5, 10];           % SNR per panel [dB]
cfg.jammer_on_sweep1_panel  = [0, 0, 0, 1];             % 0/1 jammer flag per panel
cfg.jnr_sweep1_panel4_db    = 5;                        % CW JNR on panel 4 [dB]
cfg.jam_offs_sweep1_panel4_mhz = 0;                     % CW carrier offset [MHz]

% Sweep 2:  Pd vs CFO with stronger CW jammer, on-comb / off-comb test.
cfg.cfo_grid_sweep2_khz = cfg.cfo_grid_sweep1_khz;      % reuse Sweep-1 CFO axis
cfg.snr_sweep2_db       = [10, 20];                     % 2 SNR panels [dB]
cfg.jnr_sweep2_db       = 10;                           % stronger JNR [dB]
cfg.jam_offs_sweep2_mhz = [1.25, 0.625];                % on-comb / off-comb [MHz]

% Sweep 3:  Pd vs residual JNR, 4 jammer types.
cfg.res_jnr_grid_sweep3_db   = 0 : 2.5 : 30;            % residual JNR axis [dB], 13 pts
cfg.snr_sweep3_db            = [10, 20];                % 2 SNR panels [dB]
cfg.jam_types_sweep3_codes   = [1, 2, 3, 4];            % 1=CW, 2=1MHz, 3=WB, 4=chirp
cfg.jam_random_freq_grid_mhz = -20 : 1 : 20;            % random carrier draws [MHz]

% Sweep 4:  Pd vs 1-MHz-noise jammer carrier offset (spectral fingerprint).
cfg.jam_offs_grid_sweep4_mhz = -20 : 1 : 20;            % carrier offset axis [MHz], 41 pts
cfg.snr_sweep4_db            = [0, 10];                 % 2 SNR panels [dB]
cfg.jnr_sweep4_db            = 10;                      % 1-MHz noise JNR [dB]
cfg.sweep4_det_idx           = [1, 2, 3];               % detectors plotted in Sweep 4
                                                        % (1=AC, 2=XC_SINGLE, 3=XC_BANK)

%% ====================================================================
%% ONE-TIME SETUP (cheap, runs once before the kernel)
%% ====================================================================

% ---- Simulator parameter structs ------------------------------------
% p_init      : N_h = 3 cross-correlator bank at {-145, 0, +145} kHz
%               (minimax-optimal coverage of the +/- 233 kHz CFO envelope).
%               Also used for tx_frame and as p_xc_bank.
% p_xc_single : degenerate N_h = 1 bank at f = 0 (single matched filter).
% p_xc_bank   : alias to p_init for code clarity in the kernel call.
xcb_freq_min  = -145e3;
xcb_freq_max  = 145e3;
xcb_freq_step = 145e3;
p_init = sim_params(cfg.fs_bb_freq, cfg.fs_tx_freq, cfg.fs_rx_freq, ...
    xcb_freq_min, xcb_freq_max, xcb_freq_step);
p_xc_single = sim_params(cfg.fs_bb_freq, cfg.fs_tx_freq, cfg.fs_rx_freq, ...
    0, 0, 1);
p_xc_bank = p_init;

% ---- Unit-power TX baseband waveform at fs_tx -----------------------
% tx_iq_base_full : raw L-STF+L-LTF+SIGNAL fields synthesised per IEEE
%                   802.11-2016 Sec. 17.3.3 at fs_tx (oversampled 4x).
% tx_iq_base      : same waveform normalised so mean|tx|^2 = 1, so that
%                   downstream P_in scaling cleanly equals the per-sample
%                   signal power and the SNR axis is exact.
[tx_iq_base_full, info] = tx_frame(p_init);
% Normalise tx_iq_base so the ACTIVE preamble region (L-STF + L-LTF)
% has unit mean per-sample power. Averaging over the whole frame
% would include the zero pre/post guards and bias the SNR convention
% by ~1.5 dB (cf. diag_verify_snr_jnr.m).
idx_active   = info.preamble_start_tx : info.preamble_end_tx;
active_pwr   = mean(abs(tx_iq_base_full(idx_active)).^2);
tx_iq_base   = tx_iq_base_full / sqrt(active_pwr);

% ---- L-STF window offsets in rx_iq ----------------------------------
% stf_offset_start              : index of the first L-STF sample in
%                                 rx_iq, relative to the post-guard start.
% stf_offset_end_no_channel     : index of the last  L-STF sample, AWGN
%                                 case (no multipath extension).
% Cross-checked against an empirical envelope in compute_stf_window_offsets.
[stf_offset_start, stf_offset_end_no_channel] = ...
    compute_stf_window_offsets(p_init, tx_iq_base, ...
        cfg.fs_tx_freq, cfg.fs_bb_freq, cfg.fs_rx_freq);

% ---- Widen trailing edge under TDL-D (multipath spread) -------------
% Under TDL-D the first arrival stays at the AWGN index (LOS tap at
% delay 0), but trailing taps spread the L-STF by max_delay samples.
% We widen stf_offset_end by EXACTLY this much so multipath copies of
% the L-STF inside the window still count as valid detection samples.
% tdl_max_norm = 9.708 -> last tap of TDL-D (TR 38.901 Table 7.7.2-4).
if cfg.channel_id == 1
    tdl_max_norm   = 9.708;
    max_delay_tx   = round(tdl_max_norm * cfg.tdl_delay_spread_s * cfg.fs_tx_freq);
                                    % max channel delay at fs_tx [samples]
    max_delay_rx   = ceil(max_delay_tx / round(cfg.fs_tx_freq / cfg.fs_rx_freq));
                                    % same, at fs_rx (after halfband cascade)
    stf_offset_end = stf_offset_end_no_channel + max_delay_rx;
else
    stf_offset_end = stf_offset_end_no_channel;
end

% ---- Per-detector threshold matrix ----------------------------------
% Row layout (one row per detector, indexed 1..3 matching det_id):
%   col 1: rel-gate threshold (alpha for AC, beta_sq for XC family)
%   col 2: K_abs multiplier on the CFAR-derived abs floor
%   col 3: legacy hard-coded abs floor at the nominal pn_pwr -- kept for
%          back-compat with non-CFAR diagnostics; ignored by the kernel,
%          which always uses (col 1, col 2) + per-trial sigma_sq_hat.
thresh_mat = [ ...
    cfg.alpha_ac,       cfg.K_abs_ac,        cfg.K_abs_ac        * sqrt(cfg.L_lag_samples) * cfg.pn_pwr;
    cfg.beta_sq_single, cfg.K_abs_xc_single, cfg.K_abs_xc_single * sqrt(cfg.pn_pwr);
    cfg.beta_sq_bank,   cfg.K_abs_xc_bank,   cfg.K_abs_xc_bank   * sqrt(cfg.pn_pwr)];

%% ====================================================================
%% RUN (call the CRN-optimised evaluation kernel)
%% ====================================================================
fprintf('=== eval_80211_sop_e2e_run ===\n');
if cfg.channel_id == 0
    fprintf(' Channel: AWGN\n');
else
    fprintf(' Channel: TDL-D DS=%.0f ns Fd=%.0f Hz\n', ...
        cfg.tdl_delay_spread_s*1e9, cfg.tdl_f_doppler_hz);
end
fprintf(' Trials: %d Pd / %d Pfa per sweep point\n', ...
    cfg.n_trials_pd, cfg.n_trials_pfa);

t_start        = tic;
results        = eval_80211_sop_kernel(cfg, tx_iq_base, p_xc_single, p_xc_bank, ...
                                       thresh_mat, stf_offset_start, stf_offset_end);
compute_time_s = toc(t_start);                          % wall-clock kernel time [s]
fprintf('Done. Compute time: %.1f s\n', compute_time_s);

% ---- Cosmetic side-info attached to cfg for plot_eval_results --------
cfg.detectors_str        = {'AC_BASELINE', 'XC_SINGLE_MF', 'XC_BANK'};
cfg.jam_types_sweep3_str = {'cw_random', 'noise_1mhz_random', ...
                            'wideband_gauss', 'chirp'};
cfg.sweep4_det_str       = cfg.detectors_str(cfg.sweep4_det_idx);

%% ====================================================================
%% SAVE + PLOT
%% ====================================================================

% ---- Output directory (one per channel campaign) --------------------
% channel_descr : human-readable channel tag stamped into results.mat
%                 so consumers do not need to re-parse cfg fields.
if cfg.channel_id == 0
    out_dir       = fullfile(fileparts(mfilename('fullpath')), ...
                             'results', 'sop_e2e_baseline');
    channel_descr = 'AWGN, E2E (RF+ADC)';
else
    out_dir       = fullfile(fileparts(mfilename('fullpath')), ...
                             'results', 'sop_e2e_baseline_tdl');
    channel_descr = sprintf('TDL-D, DS=%.0fns, Fd=%.0fHz, E2E (RF+ADC)', ...
        cfg.tdl_delay_spread_s*1e9, cfg.tdl_f_doppler_hz);
end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% ---- Stamp the saved artefact with first-class provenance metadata --
% Every field below is written TWICE:
%   (a) as a field of the results struct (travels with the data)
%   (b) as a loose top-level variable in the .mat
%       (load('results.mat', 'n_trials_pd') works without unpacking)
% These are the fields any downstream consumer (plotting, report
% regeneration, archival audit) should read first to understand what
% the rest of the arrays mean.
results.n_trials_pd    = cfg.n_trials_pd;       % Monte-Carlo Pd trial count
results.n_trials_pfa   = cfg.n_trials_pfa;      % Monte-Carlo Pfa trial count
results.channel_id     = cfg.channel_id;        % 0=AWGN, 1=TDL-D
results.channel_descr  = channel_descr;         % human-readable channel tag
results.run_date       = datestr(now, 'yyyy-mm-dd HH:MM:SS');
                                                % timestamp of this run
results.compute_time_s = compute_time_s;        % wall-clock kernel runtime [s]

% Loose top-level mirrors (for direct workspace access from load()):
n_trials_pd    = cfg.n_trials_pd;     %#ok<NASGU>
n_trials_pfa   = cfg.n_trials_pfa;    %#ok<NASGU>
channel_id     = cfg.channel_id;      %#ok<NASGU>
run_date       = results.run_date;    %#ok<NASGU>

save(fullfile(out_dir, 'results.mat'), 'cfg', 'results', ...
    'n_trials_pd', 'n_trials_pfa', 'channel_id', 'channel_descr', ...
    'run_date', 'compute_time_s');

% ---- Render every sweep to PNG (results/sop_e2e_baseline*/*.png) ---
plot_eval_results(cfg, results, out_dir);
