function results = eval_80211_sop_kernel(cfg, tx_iq_base, p_xc_single, p_xc_bank, ...
    thresh_mat, stf_offset_start, stf_offset_end)
%EVAL_80211_SOP_KERNEL  CRN-optimised SoP detector evaluation kernel,
%                 E2E variant (3 dual-gate detectors, 5 sweeps, AWGN or
%                 TDL-D channel, RX gain + ADC quantize + clip block).
%
%   The kernel interposes a combined RX gain + ADC quantize + clip step
%   between the noise+jammer addition and the half-band decimator
%   cascade. The cfg fields rx_total_gain_db, adc_n_bits,
%   adc_fullscale_dbm parameterise that block (see apply_rf_adc.m).
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
%   IEEE Ref: derived from IEEE Std 802.11-2016 Sec. 17.3.3 (OFDM PHY
%   preamble: L-STF / L-LTF subcarrier and timing definitions).
%
%   Inner-loop pseudocode of one Pd cell (signal present):
%
%     tx_unit  = tx_iq_base + CFO + fractional delay + [TDL-D]
%     rx_iq    = decimate( P_in * tx_unit + noise [+ jammer], 11) -> 23
%     sigma_sq = estimate_floor(rx_iq, guard region)
%     det      = (T_abs > K * sigma_sq) AND (T_rel > threshold)
%     count detection iff trigger sample lies inside L-STF window
%
%   Detectors evaluated (det_id):
%     1 = AC_BASELINE   Schmidl-Cox autocorrelator, abs-power gate
%     2 = XC_SINGLE_MF  matched-filter cross-correlator at f = 0
%     3 = XC_BANK       3-hypothesis coherent comb at {-145, 0, +145} kHz
%   See phy_rx/detector_passes_window.m for the per-detector test.
%
%   Sweeps performed (each populates one arr_pd_sweepK / arr_pfa_sweepK):
%     0  Pd vs SNR, random CFO, no jammer       -> baseline sensitivity
%     1  Pd vs CFO, 4 panels (last w/ CW jam)   -> CFO robustness
%     2  Pd vs CFO, on-comb / off-comb CW jam   -> spectral selectivity
%     3  Pd vs residual JNR (4 jammer types)    -> jammer robustness
%     4  Pd vs 1-MHz-noise jammer carrier offs. -> 1-MHz spectral map
%
%   Variance reduction: *common random numbers* (CRN). Outer = trial,
%   inner = sweep parameter axis. Per trial, draw ONCE: CFO (when
%   random), pre-packet guard length, fractional delay tau, TX-channel
%   realisation, jammer-channel realisation, AWGN noise sequence, and
%   the random parts of the jammer waveform. These draws are reused
%   across every SNR / JNR / detector point of the sweep, dropping
%   channel-generation cost by an order of magnitude on TDL channels
%   while leaving the marginal Pd estimator statistically unbiased.
%
%   CFAR floor (per trial, both estimators are computed in parallel
%   on the same pre-packet guard samples, then fed to all 3 detectors):
%     sigma_sq_mean   = mean(|r|^2)
%     sigma_sq_median = median(|r|^2) / ln(2)
%   The trailing dim-2 axis of every arr_*_sweepK is mean / median.
%
%   Inputs:
%     cfg                struct, every sweep grid / threshold / RNG knob
%                        (see eval_80211_sop_e2e_run.m for the full set)
%     tx_iq_base         unit-power TX baseband at fs_tx_freq [1 x N_tx]
%     p_xc_single        sim_params variant with N_h = 1 (XC_SINGLE_MF)
%     p_xc_bank          sim_params variant with N_h = 3 (XC_BANK)
%     thresh_mat         3 x 2 = [rel_thresh, K_abs] per detector
%     stf_offset_start   first sample of L-STF in rx_iq (post-guard)
%     stf_offset_end     last sample of L-STF (widened by TDL spread)
%
%   Output:
%     results struct with arr_pd_sweep{0..4}, arr_pfa_sweep{0..4} and
%     the {stf_offset_start, stf_offset_end} used. The trailing axis
%     of every arr_* is {mean, median} CFAR estimator.
%
%   See also: eval_80211_sop_e2e_run.m  (driver script that fills cfg,
%             builds tx_iq_base, calls this kernel, then renders the
%             results to PNG / .mat),
%             phy_rx/detector_passes_window.m, phy_rx/estimate_floor.m,
%             channel/channel_tdl_d_gen.m, channel/channel_tdl_d_apply.m,
%             channel/make_jammer.m
%
if isfield(cfg, 'rng_seed') && cfg.rng_seed >= 0
    rng(cfg.rng_seed, 'twister');
end

%% ---- Cached scalars ---------------------------------------------
osr_tx_bb = round(cfg.fs_tx_freq / cfg.fs_bb_freq);
osr_tx_rx = round(cfg.fs_tx_freq / cfg.fs_rx_freq);
n_tx_base = numel(tx_iq_base);
n_dets = 3;
pn_pwr = cfg.pn_pwr;
noise_pwr_tx = pn_pwr * (cfg.fs_tx_freq / cfg.fs_bb_freq);
sqrt_noise = sqrt(noise_pwr_tx / 2);

sqrt_pin_sweep0 = sqrt(pn_pwr * 10.^(cfg.snr_grid_sweep0_db / 10));
sqrt_pin_sweep1 = sqrt(pn_pwr * 10.^(cfg.snr_sweep1_panels_db / 10));
sqrt_pin_sweep2 = sqrt(pn_pwr * 10.^(cfg.snr_sweep2_db / 10));
sqrt_pin_sweep3 = sqrt(pn_pwr * 10.^(cfg.snr_sweep3_db / 10));
sqrt_pin_sweep4 = sqrt(pn_pwr * 10.^(cfg.snr_sweep4_db / 10));
sqrt_jam_sweep1 = sqrt(pn_pwr * 10^(cfg.jnr_sweep1_panel4_db / 10));
sqrt_jam_sweep2 = sqrt(pn_pwr * 10^(cfg.jnr_sweep2_db / 10));
sqrt_jam_sweep4 = sqrt(pn_pwr * 10^(cfg.jnr_sweep4_db / 10));
sqrt_jam_sweep3 = sqrt(pn_pwr * 10.^(cfg.res_jnr_grid_sweep3_db / 10));

n_snr0 = numel(cfg.snr_grid_sweep0_db);
n_cfo1 = numel(cfg.cfo_grid_sweep1_khz);
n_panel1 = numel(cfg.snr_sweep1_panels_db);
n_cfo2 = numel(cfg.cfo_grid_sweep2_khz);
n_snr2 = numel(cfg.snr_sweep2_db);
n_pos2 = numel(cfg.jam_offs_sweep2_mhz);
n_res3 = numel(cfg.res_jnr_grid_sweep3_db);
n_snr3 = numel(cfg.snr_sweep3_db);
n_jt3 = numel(cfg.jam_types_sweep3_codes);
n_off4 = numel(cfg.jam_offs_grid_sweep4_mhz);
n_snr4 = numel(cfg.snr_sweep4_db);
n_sweep4_det = numel(cfg.sweep4_det_idx);
n_random_grid = numel(cfg.jam_random_freq_grid_mhz);

arr_pd_sweep0 = zeros(n_dets, n_snr0, 2);
arr_pfa_sweep0 = zeros(n_dets, 1, 2);
arr_pd_sweep1 = zeros(n_dets, n_panel1, n_cfo1, 2);
arr_pfa_sweep1 = zeros(n_dets, n_panel1, 2);
arr_pd_sweep2 = zeros(n_dets, n_snr2, n_pos2, n_cfo2, 2);
arr_pfa_sweep2 = zeros(n_dets, n_snr2, n_pos2, 2);
arr_pd_sweep3 = zeros(n_dets, n_snr3, n_jt3, n_res3, 2);
arr_pfa_sweep3 = zeros(n_dets, n_snr3, n_jt3, n_res3, 2);
arr_pd_sweep4 = zeros(n_sweep4_det, n_snr4, n_off4, 2);
arr_pfa_sweep4 = zeros(n_sweep4_det, n_snr4, n_off4, 2);

%% ==================================================================
%% SWEEP 0: P_d vs SNR (random CFO, no jammer)
%% ==================================================================
fprintf(' sweep 0 (P_d vs SNR, random CFO, no jammer)...\n');
nmn = zeros(n_dets, n_snr0); nmd = zeros(n_dets, n_snr0);
for k_trial = 1:cfg.n_trials_pd
    cfo_k = (2*rand()- 1) * cfg.random_cfo_max_freq;
    guard_rx = randi([cfg.guard_min_samples_rx, cfg.guard_max_samples_rx]);
    guard_tx = guard_rx * osr_tx_rx;
    tau_k = rand()* osr_tx_bb;
    n_total = guard_tx + n_tx_base;

    tx_unit = build_tx_through_channel(tx_iq_base, guard_tx, cfo_k, tau_k, ...
        cfg.fs_tx_freq, cfg.channel_id, cfg.tdl_delay_spread_s, cfg.tdl_f_doppler_hz);
    noise_seq = sqrt_noise * (randn(1, n_total) + 1j*randn(1, n_total));

    for k_snr = 1:n_snr0
        rx_full = sqrt_pin_sweep0(k_snr) * tx_unit + noise_seq;
        rx_full = apply_rf_adc(rx_full, cfg.rx_total_gain_db, ...
                               cfg.adc_n_bits, cfg.adc_fullscale_dbm);
        rx_iq = decimate_hb_2to1(rx_full, 11);
        rx_iq = decimate_hb_2to1(rx_iq, 23);
        [sm, sM] = estimate_floor(rx_iq, cfg.cfar_n_samples);
        for k_det = 1:n_dets
            det = detector_passes_window(k_det, rx_iq, thresh_mat, ...
                p_xc_single, p_xc_bank, ...
                cfg.L_lag_samples, cfg.S_ac, cfg.S_xc, cfg.eps_load, ...
                guard_rx + stf_offset_start, guard_rx + stf_offset_end, sm, sM);
            nmn(k_det, k_snr) = nmn(k_det, k_snr) + det(1);
            nmd(k_det, k_snr) = nmd(k_det, k_snr) + det(2);
        end
    end
end
arr_pd_sweep0 = cat(3, nmn, nmd) / cfg.n_trials_pd;

arr_pfa_sweep0 = reshape(run_pfa_crn(cfg, tx_iq_base, p_xc_single, p_xc_bank, ...
    thresh_mat, cfg.random_cfo_max_freq, 0, 0, 0, false), n_dets, 1, 2);

%% ==================================================================
%% SWEEP 1: P_d vs CFO (4 panels)
%% ==================================================================
fprintf(' sweep 1 (P_d vs CFO, 4 panels)...\n');
for k_panel = 1:n_panel1
    jammer_on = cfg.jammer_on_sweep1_panel(k_panel) ~= 0;
    if jammer_on
        jam_pwr_k = sqrt_jam_sweep1^2;
        jam_fhz_k = cfg.jam_offs_sweep1_panel4_mhz * 1e6;
        jam_id = 1;
    else
        jam_pwr_k = 0; jam_fhz_k = 0; jam_id = 0;
    end
    sqrt_pin_p = sqrt_pin_sweep1(k_panel);

    nmn = zeros(n_dets, n_cfo1); nmd = zeros(n_dets, n_cfo1);
    for k_trial = 1:cfg.n_trials_pd
        guard_rx = randi([cfg.guard_min_samples_rx, cfg.guard_max_samples_rx]);
        guard_tx = guard_rx * osr_tx_rx;
        tau_k = rand()* osr_tx_bb;
        n_total = guard_tx + n_tx_base;

        if cfg.channel_id == 1
            [fading_tx, delays_tx] = channel_tdl_d_gen(n_total, ...
                cfg.fs_tx_freq, cfg.tdl_delay_spread_s, cfg.tdl_f_doppler_hz);
        else
            fading_tx = []; delays_tx = [];
        end
        noise_seq = sqrt_noise * (randn(1, n_total) + 1j*randn(1, n_total));

        if jammer_on
            jam_base = make_jammer(jam_id, jam_pwr_k, jam_fhz_k, n_total, ...
                cfg.fs_tx_freq, cfg.fs_bb_freq, ...
                cfg.chirp_f_start_freq, cfg.chirp_f_end_freq);
            if cfg.channel_id == 1
                [fading_jam, delays_jam] = channel_tdl_d_gen(n_total, ...
                    cfg.fs_tx_freq, cfg.tdl_delay_spread_s, cfg.tdl_f_doppler_hz);
                jam_iq = channel_tdl_d_apply(jam_base, fading_jam, delays_jam);
            else
                jam_iq = jam_base;
            end
        else
            jam_iq = complex(zeros(1, n_total));
        end

        for k_cfo = 1:n_cfo1
            cfo_k = cfg.cfo_grid_sweep1_khz(k_cfo) * 1e3;
            tx_unit = build_tx_with_precomp_channel(tx_iq_base, guard_tx, ...
                cfo_k, tau_k, cfg.fs_tx_freq, cfg.channel_id, fading_tx, delays_tx);
            rx_full = sqrt_pin_p * tx_unit + noise_seq + jam_iq;
            rx_full = apply_rf_adc(rx_full, cfg.rx_total_gain_db, ...
                                   cfg.adc_n_bits, cfg.adc_fullscale_dbm);
            rx_iq = decimate_hb_2to1(rx_full, 11);
            rx_iq = decimate_hb_2to1(rx_iq, 23);
            [sm, sM] = estimate_floor(rx_iq, cfg.cfar_n_samples);
            for k_det = 1:n_dets
                det = detector_passes_window(k_det, rx_iq, thresh_mat, ...
                    p_xc_single, p_xc_bank, ...
                    cfg.L_lag_samples, cfg.S_ac, cfg.S_xc, cfg.eps_load, ...
                    guard_rx + stf_offset_start, guard_rx + stf_offset_end, sm, sM);
                nmn(k_det, k_cfo) = nmn(k_det, k_cfo) + det(1);
                nmd(k_det, k_cfo) = nmd(k_det, k_cfo) + det(2);
            end
        end
    end
    arr_pd_sweep1(:, k_panel, :, 1) = nmn / cfg.n_trials_pd;
    arr_pd_sweep1(:, k_panel, :, 2) = nmd / cfg.n_trials_pd;

    pfa_panel = run_pfa_crn(cfg, tx_iq_base, p_xc_single, p_xc_bank, ...
        thresh_mat, 0, jam_pwr_k, jam_fhz_k, jam_id, false);
    arr_pfa_sweep1(:, k_panel, :) = reshape(pfa_panel, n_dets, 1, 2);
end

%% ==================================================================
%% SWEEP 2: P_d vs CFO with CW jammer (on-comb / off-comb)
%% ==================================================================
fprintf(' sweep 2 (P_d vs CFO, CW jammer comb/off-comb)...\n');
for k_snr = 1:n_snr2
    sqrt_pin = sqrt_pin_sweep2(k_snr);
    jam_pwr_k = sqrt_jam_sweep2^2;
    for k_pos = 1:n_pos2
        jam_fhz_k = cfg.jam_offs_sweep2_mhz(k_pos) * 1e6;

        nmn = zeros(n_dets, n_cfo2); nmd = zeros(n_dets, n_cfo2);
        for k_trial = 1:cfg.n_trials_pd
            guard_rx = randi([cfg.guard_min_samples_rx, cfg.guard_max_samples_rx]);
            guard_tx = guard_rx * osr_tx_rx;
            tau_k = rand()* osr_tx_bb;
            n_total = guard_tx + n_tx_base;

            if cfg.channel_id == 1
                [fading_tx, delays_tx] = channel_tdl_d_gen(n_total, ...
                    cfg.fs_tx_freq, cfg.tdl_delay_spread_s, cfg.tdl_f_doppler_hz);
                [fading_jam, delays_jam] = channel_tdl_d_gen(n_total, ...
                    cfg.fs_tx_freq, cfg.tdl_delay_spread_s, cfg.tdl_f_doppler_hz);
            else
                fading_tx = []; delays_tx = [];
                fading_jam = []; delays_jam = [];
            end
            noise_seq = sqrt_noise * (randn(1, n_total) + 1j*randn(1, n_total));

            jam_base = make_jammer(1, jam_pwr_k, jam_fhz_k, n_total, ...
                cfg.fs_tx_freq, cfg.fs_bb_freq, ...
                cfg.chirp_f_start_freq, cfg.chirp_f_end_freq);
            if cfg.channel_id == 1
                jam_iq = channel_tdl_d_apply(jam_base, fading_jam, delays_jam);
            else
                jam_iq = jam_base;
            end

            for k_cfo = 1:n_cfo2
                cfo_k = cfg.cfo_grid_sweep2_khz(k_cfo) * 1e3;
                tx_unit = build_tx_with_precomp_channel(tx_iq_base, guard_tx, ...
                    cfo_k, tau_k, cfg.fs_tx_freq, cfg.channel_id, fading_tx, delays_tx);
                rx_full = sqrt_pin * tx_unit + noise_seq + jam_iq;
                rx_full = apply_rf_adc(rx_full, cfg.rx_total_gain_db, ...
                                       cfg.adc_n_bits, cfg.adc_fullscale_dbm);
                rx_iq = decimate_hb_2to1(rx_full, 11);
                rx_iq = decimate_hb_2to1(rx_iq, 23);
                [sm, sM] = estimate_floor(rx_iq, cfg.cfar_n_samples);
                for k_det = 1:n_dets
                    det = detector_passes_window(k_det, rx_iq, thresh_mat, ...
                        p_xc_single, p_xc_bank, ...
                        cfg.L_lag_samples, cfg.S_ac, cfg.S_xc, cfg.eps_load, ...
                        guard_rx + stf_offset_start, guard_rx + stf_offset_end, sm, sM);
                    nmn(k_det, k_cfo) = nmn(k_det, k_cfo) + det(1);
                    nmd(k_det, k_cfo) = nmd(k_det, k_cfo) + det(2);
                end
            end
        end
        arr_pd_sweep2(:, k_snr, k_pos, :, 1) = nmn / cfg.n_trials_pd;
        arr_pd_sweep2(:, k_snr, k_pos, :, 2) = nmd / cfg.n_trials_pd;

        pfa = run_pfa_crn(cfg, tx_iq_base, p_xc_single, p_xc_bank, thresh_mat, ...
            0, jam_pwr_k, jam_fhz_k, 1, false);
        arr_pfa_sweep2(:, k_snr, k_pos, :) = reshape(pfa, n_dets, 1, 1, 2);
    end
end

%% ==================================================================
%% SWEEP 3: P_d vs residual JNR (4 jam types, random CFO)
%% ==================================================================
fprintf(' sweep 3 (P_d vs residual JNR)...\n');
for k_snr = 1:n_snr3
    sqrt_pin = sqrt_pin_sweep3(k_snr);
    for k_jt = 1:n_jt3
        switch cfg.jam_types_sweep3_codes(k_jt)
        case 1, jam_id_s3 = 1; jam_random = true; jam_fhz_fixed = 0;
        case 2, jam_id_s3 = 2; jam_random = true; jam_fhz_fixed = 0;
        case 3, jam_id_s3 = 3; jam_random = false; jam_fhz_fixed = 0;
        case 4, jam_id_s3 = 4; jam_random = false; jam_fhz_fixed = 0;
        otherwise, jam_id_s3 = 0; jam_random = false; jam_fhz_fixed = 0;
        end

        nmn = zeros(n_dets, n_res3); nmd = zeros(n_dets, n_res3);
        for k_trial = 1:cfg.n_trials_pd
            cfo_k = (2*rand()- 1) * cfg.random_cfo_max_freq;
            guard_rx = randi([cfg.guard_min_samples_rx, cfg.guard_max_samples_rx]);
            guard_tx = guard_rx * osr_tx_rx;
            tau_k = rand()* osr_tx_bb;
            n_total = guard_tx + n_tx_base;

            tx_unit = build_tx_through_channel(tx_iq_base, guard_tx, cfo_k, tau_k, ...
                cfg.fs_tx_freq, cfg.channel_id, cfg.tdl_delay_spread_s, cfg.tdl_f_doppler_hz);
            noise_seq = sqrt_noise * (randn(1, n_total) + 1j*randn(1, n_total));

            if jam_random
                jam_fhz_k = cfg.jam_random_freq_grid_mhz(randi(n_random_grid)) * 1e6;
            else
                jam_fhz_k = jam_fhz_fixed;
            end

            % Unit-power jammer (scale by sqrt(jam_pwr) per JNR point)
            jam_unit = make_jammer(jam_id_s3, 1.0, jam_fhz_k, n_total, ...
                cfg.fs_tx_freq, cfg.fs_bb_freq, ...
                cfg.chirp_f_start_freq, cfg.chirp_f_end_freq);
            if cfg.channel_id == 1
                [fading_jam, delays_jam] = channel_tdl_d_gen(n_total, ...
                    cfg.fs_tx_freq, cfg.tdl_delay_spread_s, cfg.tdl_f_doppler_hz);
                jam_unit_post = channel_tdl_d_apply(jam_unit, fading_jam, delays_jam);
            else
                jam_unit_post = jam_unit;
            end

            for k_res = 1:n_res3
                jam_iq = sqrt_jam_sweep3(k_res) * jam_unit_post;
                rx_full = sqrt_pin * tx_unit + noise_seq + jam_iq;
                rx_full = apply_rf_adc(rx_full, cfg.rx_total_gain_db, ...
                                       cfg.adc_n_bits, cfg.adc_fullscale_dbm);
                rx_iq = decimate_hb_2to1(rx_full, 11);
                rx_iq = decimate_hb_2to1(rx_iq, 23);
                [sm, sM] = estimate_floor(rx_iq, cfg.cfar_n_samples);
                for k_det = 1:n_dets
                    det = detector_passes_window(k_det, rx_iq, thresh_mat, ...
                        p_xc_single, p_xc_bank, ...
                        cfg.L_lag_samples, cfg.S_ac, cfg.S_xc, cfg.eps_load, ...
                        guard_rx + stf_offset_start, guard_rx + stf_offset_end, sm, sM);
                    nmn(k_det, k_res) = nmn(k_det, k_res) + det(1);
                    nmd(k_det, k_res) = nmd(k_det, k_res) + det(2);
                end
            end
        end
        arr_pd_sweep3(:, k_snr, k_jt, :, 1) = nmn / cfg.n_trials_pd;
        arr_pd_sweep3(:, k_snr, k_jt, :, 2) = nmd / cfg.n_trials_pd;

        for k_res = 1:n_res3
            pfa = run_pfa_crn(cfg, tx_iq_base, p_xc_single, p_xc_bank, thresh_mat, ...
                cfg.random_cfo_max_freq, sqrt_jam_sweep3(k_res)^2, ...
                jam_fhz_fixed, jam_id_s3, jam_random);
            arr_pfa_sweep3(:, k_snr, k_jt, k_res, :) = reshape(pfa, n_dets, 1, 1, 1, 2);
        end
    end
end

%% ==================================================================
%% SWEEP 4: P_d vs 1-MHz-noise carrier offset (random CFO)
%% ==================================================================
fprintf(' sweep 4 (P_d vs 1-MHz-noise carrier offset)...\n');
for k_snr = 1:n_snr4
    sqrt_pin = sqrt_pin_sweep4(k_snr);
    jam_pwr_k = sqrt_jam_sweep4^2;

    nmn = zeros(n_sweep4_det, n_off4); nmd = zeros(n_sweep4_det, n_off4);
    for k_trial = 1:cfg.n_trials_pd
        cfo_k = (2*rand()- 1) * cfg.random_cfo_max_freq;
        guard_rx = randi([cfg.guard_min_samples_rx, cfg.guard_max_samples_rx]);
        guard_tx = guard_rx * osr_tx_rx;
        tau_k = rand()* osr_tx_bb;
        n_total = guard_tx + n_tx_base;

        tx_unit = build_tx_through_channel(tx_iq_base, guard_tx, cfo_k, tau_k, ...
            cfg.fs_tx_freq, cfg.channel_id, cfg.tdl_delay_spread_s, cfg.tdl_f_doppler_hz);
        noise_seq = sqrt_noise * (randn(1, n_total) + 1j*randn(1, n_total));

        if cfg.channel_id == 1
            [fading_jam, delays_jam] = channel_tdl_d_gen(n_total, ...
                cfg.fs_tx_freq, cfg.tdl_delay_spread_s, cfg.tdl_f_doppler_hz);
        else
            fading_jam = []; delays_jam = [];
        end

        for k_off = 1:n_off4
            jam_fhz_k = cfg.jam_offs_grid_sweep4_mhz(k_off) * 1e6;
            jam_base = make_jammer(2, jam_pwr_k, jam_fhz_k, n_total, ...
                cfg.fs_tx_freq, cfg.fs_bb_freq, ...
                cfg.chirp_f_start_freq, cfg.chirp_f_end_freq);
            if cfg.channel_id == 1
                jam_iq = channel_tdl_d_apply(jam_base, fading_jam, delays_jam);
            else
                jam_iq = jam_base;
            end
            rx_full = sqrt_pin * tx_unit + noise_seq + jam_iq;
            rx_full = apply_rf_adc(rx_full, cfg.rx_total_gain_db, ...
                                   cfg.adc_n_bits, cfg.adc_fullscale_dbm);
            rx_iq = decimate_hb_2to1(rx_full, 11);
            rx_iq = decimate_hb_2to1(rx_iq, 23);
            [sm, sM] = estimate_floor(rx_iq, cfg.cfar_n_samples);
            for k_local = 1:n_sweep4_det
                k_det = cfg.sweep4_det_idx(k_local);
                det = detector_passes_window(k_det, rx_iq, thresh_mat, ...
                    p_xc_single, p_xc_bank, ...
                    cfg.L_lag_samples, cfg.S_ac, cfg.S_xc, cfg.eps_load, ...
                    guard_rx + stf_offset_start, guard_rx + stf_offset_end, sm, sM);
                nmn(k_local, k_off) = nmn(k_local, k_off) + det(1);
                nmd(k_local, k_off) = nmd(k_local, k_off) + det(2);
            end
        end
    end
    arr_pd_sweep4(:, k_snr, :, 1) = nmn / cfg.n_trials_pd;
    arr_pd_sweep4(:, k_snr, :, 2) = nmd / cfg.n_trials_pd;

    for k_off = 1:n_off4
        jam_fhz_k = cfg.jam_offs_grid_sweep4_mhz(k_off) * 1e6;
        pfa_all = run_pfa_crn(cfg, tx_iq_base, p_xc_single, p_xc_bank, ...
            thresh_mat, cfg.random_cfo_max_freq, ...
            jam_pwr_k, jam_fhz_k, 2, false);
        arr_pfa_sweep4(:, k_snr, k_off, :) = reshape( ...
            pfa_all(cfg.sweep4_det_idx, :), n_sweep4_det, 1, 1, 2);
    end
end

%% ---- Pack --------------------------------------------------------
results = struct( ...
    'arr_pd_sweep0', arr_pd_sweep0, ...
    'arr_pfa_sweep0', arr_pfa_sweep0, ...
    'arr_pd_sweep1', arr_pd_sweep1, ...
    'arr_pfa_sweep1', arr_pfa_sweep1, ...
    'arr_pd_sweep2', arr_pd_sweep2, ...
    'arr_pfa_sweep2', arr_pfa_sweep2, ...
    'arr_pd_sweep3', arr_pd_sweep3, ...
    'arr_pfa_sweep3', arr_pfa_sweep3, ...
    'arr_pd_sweep4', arr_pd_sweep4, ...
    'arr_pfa_sweep4', arr_pfa_sweep4, ...
    'stf_offset_start', stf_offset_start, ...
    'stf_offset_end', stf_offset_end);
end

%% ======================================================================
%% Local helpers (kernel-private)
%% ======================================================================

function tx_unit = build_tx_through_channel(tx_iq_base, guard_tx, cfo_freq, ...
    tau, fs_tx_freq, channel_id, tdl_ds_s, tdl_fd_hz)
    n_total = guard_tx + numel(tx_iq_base);
    tx_unit = complex([zeros(1, guard_tx), tx_iq_base]);
    if cfo_freq ~= 0
        tx_unit = tx_unit .* exp(1j*2*pi*cfo_freq/fs_tx_freq * (0:n_total-1));
    end
    tx_unit = fractional_delay(tx_unit, tau);
    if channel_id == 1
        [fading_tx, delays_tx] = channel_tdl_d_gen(n_total, fs_tx_freq, ...
            tdl_ds_s, tdl_fd_hz);
        tx_unit = channel_tdl_d_apply(tx_unit, fading_tx, delays_tx);
    end
end

function tx_unit = build_tx_with_precomp_channel(tx_iq_base, guard_tx, ...
    cfo_freq, tau, fs_tx_freq, channel_id, fading_tx, delays_tx)
    n_total = guard_tx + numel(tx_iq_base);
    tx_unit = complex([zeros(1, guard_tx), tx_iq_base]);
    if cfo_freq ~= 0
        tx_unit = tx_unit .* exp(1j*2*pi*cfo_freq/fs_tx_freq * (0:n_total-1));
    end
    tx_unit = fractional_delay(tx_unit, tau);
    if channel_id == 1
        tx_unit = channel_tdl_d_apply(tx_unit, fading_tx, delays_tx);
    end
end

function pfa = run_pfa_crn(cfg, tx_iq_base, p_xc_single, p_xc_bank, thresh_mat, ...
    cfo_max_random, jam_pwr, jam_fhz_fixed, jam_id, jam_random)
%RUN_PFA_CRN  N_pfa-window false-alarm rate (signal = 0, channel = no-op).
% Per trial, draw ONCE: noise, jammer carrier offset, jammer waveform,
% jammer-channel realisation. Run all 3 detectors on the same rx_iq.

    osr_tx_rx = round(cfg.fs_tx_freq / cfg.fs_rx_freq);
    % Pfa observation window: 1 SIFS = cfg.n_pfa_observation_samples
    % (= 320 samples @ fs_rx = 20 MHz, = 16 us, IEEE 802.11 Short
    % Inter-Frame Space listen interval). Starts at the first digital
    % sample (sample 1).
    n_pfa_start = 1;
    n_pfa_end   = cfg.n_pfa_observation_samples;
    n_tx_iq_rx = ceil(numel(tx_iq_base) * cfg.fs_rx_freq / cfg.fs_tx_freq);
    guard_pfa = max(0, n_pfa_end - n_tx_iq_rx) + 32;
    guard_tx_pfa = guard_pfa * osr_tx_rx;
    n_total = guard_tx_pfa + numel(tx_iq_base);
    noise_pwr_tx = cfg.pn_pwr * (cfg.fs_tx_freq / cfg.fs_bb_freq);
    sqrt_noise = sqrt(noise_pwr_tx/2);
    n_random_grid = numel(cfg.jam_random_freq_grid_mhz);
    n_dets = 3;

    n_fa = zeros(n_dets, 2);
    for k = 1:cfg.n_trials_pfa
        rx_full = sqrt_noise * (randn(1, n_total) + 1j*randn(1, n_total));
        if jam_pwr > 0 && jam_id ~= 0
            if jam_random && n_random_grid > 0
                jam_fhz_k = cfg.jam_random_freq_grid_mhz(randi(n_random_grid)) * 1e6;
            else
                jam_fhz_k = jam_fhz_fixed;
            end
            jam_iq = make_jammer(jam_id, jam_pwr, jam_fhz_k, n_total, ...
                cfg.fs_tx_freq, cfg.fs_bb_freq, ...
                cfg.chirp_f_start_freq, cfg.chirp_f_end_freq);
            if cfg.channel_id == 1
                [fading_jam, delays_jam] = channel_tdl_d_gen(n_total, ...
                    cfg.fs_tx_freq, cfg.tdl_delay_spread_s, cfg.tdl_f_doppler_hz);
                jam_iq = channel_tdl_d_apply(jam_iq, fading_jam, delays_jam);
            end
            rx_full = rx_full + jam_iq;
        end
        if cfo_max_random > 0
            cfo_freq_k = (2*rand()- 1) * cfo_max_random;
            rx_full = rx_full .* exp(1j*2*pi*cfo_freq_k/cfg.fs_tx_freq * (0:n_total-1));
        end
        rx_full = apply_rf_adc(rx_full, cfg.rx_total_gain_db, ...
                               cfg.adc_n_bits, cfg.adc_fullscale_dbm);
        rx_iq = decimate_hb_2to1(rx_full, 11);
        rx_iq = decimate_hb_2to1(rx_iq, 23);
        [sm, sM] = estimate_floor(rx_iq, cfg.cfar_n_samples);
        % Pfa: any dual-gate firing in the 1-SIFS listen window
        % [n_pfa_start, n_pfa_end] = [1, 320] counts as a false alarm.
        % The buffer contains NO signal (zero TX); the window starts
        % at sample 1 and runs for 1 SIFS = 320 samples.
        for k_det = 1:n_dets
            det = detector_passes_window(k_det, rx_iq, thresh_mat, ...
                p_xc_single, p_xc_bank, ...
                cfg.L_lag_samples, cfg.S_ac, cfg.S_xc, cfg.eps_load, ...
                n_pfa_start, n_pfa_end, sm, sM, 'pfa');
            n_fa(k_det, :) = n_fa(k_det, :) + det;
        end
    end
    pfa = n_fa / cfg.n_trials_pfa;
end
