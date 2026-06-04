function [metric_pwr, freq_argmax, per_freq_pwr, xc_mag_max_per_sample] = ...
    rx_lstf_xcorr_bank(rx_iq, p)
%RX_LSTF_XCORR_BANK  L-STF matched-filter bank with COHERENT comb integration
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
% across STS periods. Optimised for jammer robustness.
%
% IEEE Ref : Operates on the L-STF defined in IEEE Std 802.11-2016, Sec. 17.3.3.
% No direct sample-based equivalent in wifi_e2e_sim -- the reference
% repo uses a 96-tap L-LTF cross-correlator (rx_lltf_crosscorr_sb.m);
% this is its L-STF analogue with an explicit CFO hypothesis grid
% AND a coherent comb across STS periods (the new optimisation).
%
% Algorithm:
% For each frequency hypothesis f_k in
% f_k in [p.xcb_freq_min : p.xcb_freq_step : p.xcb_freq_max]
% compute:
% rxd[n] = rx_iq[n] * exp(-j*2*pi*f_k*n/fs_rx) (CFO derot)
% xc[n] = sum_{m=0..L_t-1} rxd[n-m] * conj(t[L_t-1-m]) (matched filter)
% v[n] = sum_{q=0..P-1} xc[n - q*L_period] (P-tap COHERENT comb)
% metric[n] = |v[n]|^2 / (P * sum_{k=0..P*L_period-1} |rxd[n-k]|^2)
% where L_period = p.n_stf_period * p.os_rx (one STS short-symbol period
% at fs_rx), L_t = template length, and P = p.xcb_coh_periods (default 4).
% The bank metric is max_k metric_k[n] (greedy over the frequency bank).
%
% ---- WHY COHERENT COMB ? (the optimisation) -------------------------
% The L-STF is ~10 repetitions of a 16-sample STS short symbol. After CFO
% derotation, matched-filter peaks land at multiples of L_period samples
% apart and IN PHASE (because each period contains the same template).
% Summing P consecutive peaks COHERENTLY (i.e. sum complex, then take
% |.|^2) builds the signal amplitude by P (signal power by P^2), whereas
% the legacy "sum |xc|^2 over P*L_period samples" formulation builds it
% by only sqrt(P) in amplitude (incoherent sum). Net SNR gain at the
% detector input ≈ factor P, i.e. ~6 dB for P=4.
%
% ---- WHY JAMMER ROBUST ? (the other side of the same coin) ----------
% The P-tap comb {1, z^-L_period, z^-2L_period, ..., z^-(P-1)L_period} is
% itself a digital filter with frequency response
% |H_comb(f)| = |sin(pi*f*P*L_period/fs) / sin(pi*f*L_period/fs)|.
% Passband peaks land at f = m * fs/L_period (m integer); deep nulls at
% the intermediate frequencies. At fs_bb=20 MHz, os_rx=1, L_period=16,
% the passband peaks sit at multiples of 1.25 MHz -- which is EXACTLY
% the L-STF subcarrier comb (the STS is built from subcarriers
% ±4,±8,±12,±16,±20,±24 of the 64-pt FFT, i.e. ±1.25, ±2.5, ..., ±7.5 MHz).
% So this comb passes everything the L-STF has spectral support on, and
% rejects everything in between -- a narrowband interferer at any other
% frequency lands in a comb null and is suppressed by ~P (~12 dB for P=4).
% This is a pure-digital echo of what the reference RTL achieves with
% analogue notch filters + energy gating.
%
% ---- WHY THE METRIC STILL LIES IN [0,1] -----------------------------
% By Cauchy-Schwarz applied to the effective length-(P*L_period) filter
% that is the template tiled P times:
% |v[n]|^2 <= ||h_eff||^2 * sum_{k=0..P*L_period-1} |rxd[n-k]|^2
% ||h_eff||^2 = P * ||template||^2 = P (template is unit-energy).
% Dividing by P * sum|r|^2 thus keeps metric in [0,1], like the legacy
% formulation, and preserves threshold semantics across migrations.
%
% ---- COST -----------------------------------------------------------
% The dominant cost is unchanged: O(n_iq * n_freq * L_t) for the matched
% filter bank. The comb adds (P-1) shift-and-adds per hypothesis (free).
% The 1-D denominator sliding sum is O(n_iq).
%
% Inputs:
% rx_iq [1 x n_iq] complex baseband at fs_rx
% p [struct] sim_params (.xcb_*, .n_stf_period, .os_rx, .fs_rx)
%
% Outputs:
% metric_pwr [1 x n_iq] max-over-frequency normalised
% coherent-comb metric (in [0,1])
% freq_argmax [1 x n_iq] winning frequency hypothesis [Hz]
% per_freq_pwr [n_freq x n_iq] full bank metric (diagnostic)
% xc_mag_max_per_sample [1 x n_iq] per-sample peak |v| (combed magnitude)
% for the abs-floor gate; the natural
% level above noise scales as P for
% STS-aligned signals.
%
% Role : RX DSP / preamble detection (jammer-resilient alternative)
% Phase : 1 (Floating-Point, frame-based)
%
% See also: make_lstf, phy_rx/xc_dual_stat.m

% --- template at fs_rx (upsampled 16-sample STS short symbol) ----------
    [~, short_sym_iq] = make_lstf(p); % 16 samples at fs_bb (20 MHz)
    short_sym_iq = short_sym_iq / sqrt(sum(abs(short_sym_iq).^2));

    template_iq = resample_int(short_sym_iq, p.os_rx, 1);
    template_iq = template_iq(1 : p.xcb_template_len_bb * p.os_rx);
    template_iq = template_iq / sqrt(sum(abs(template_iq).^2));

    % Optional template quantization (per-I/Q rail). When p.xcb_template_max_level
    % is finite, each I and Q component is rounded to the symmetric integer
    % level set {-N, ..., -1, 0, +1, ..., +N} where N = xcb_template_max_level
    % (so total levels = 2N+1). Re-normalises to unit energy so the metric
    % stays bounded.
    if isfinite(p.xcb_template_max_level) && p.xcb_template_max_level >= 1
        n_lvl = p.xcb_template_max_level;
        max_abs = max( max(abs(real(template_iq))), max(abs(imag(template_iq))) );
        if max_abs > 0
            scale = n_lvl / max_abs;
            ti_re = max( min( round(scale * real(template_iq)), n_lvl), -n_lvl);
            ti_im = max( min( round(scale * imag(template_iq)), n_lvl), -n_lvl);
            template_iq = complex(ti_re, ti_im);
            template_iq = template_iq / sqrt(sum(abs(template_iq).^2));
        end
    end

    n_template = numel(template_iq);

    % --- frequency grid ----------------------------------------------------
    freq_grid = p.xcb_freq_min : p.xcb_freq_step : p.xcb_freq_max;
    n_freq = numel(freq_grid);
    n_iq = numel(rx_iq);
    fs = p.fs_rx;

    % --- coherent comb sizing ---------------------------------------------
    % Comb stride = one STS short-symbol period (in samples at fs_rx).
    % Number of taps P = p.xcb_coh_periods (default 4 -> ~4 short symbols
    % of coherent integration). The effective integration length is
    % P * L_period_idx samples (= the support of the tiled template filter).
    l_period_idx = p.n_stf_period * p.os_rx;
    n_coh = max(p.xcb_coh_periods, 1);
    n_den_win = n_coh * l_period_idx;

    % --- precompute per-hypothesis matched-filter then comb ---------------
    per_freq_v_pwr = zeros(n_freq, n_iq);
    n_idx = (0:n_iq-1);

    mag2_iq = abs(rx_iq).^2;
    hir_template_rev = conj(fliplr(template_iq)); % matched filter (time-reversed conj)

    for k = 1:n_freq
        f_k = freq_grid(k);
        rot_iq = exp(-1j*2*pi*f_k/fs * n_idx);
        rxd_iq = rx_iq .* rot_iq;
        xc_full = conv(rxd_iq, hir_template_rev);

        % Align matched-filter output so xc_iq(n) corresponds to the response
        % to the most recent template_len samples of rxd_iq ending at index n.
        xc_iq = xc_full(n_template : n_template + n_iq - 1);

        % --- COHERENT COMB SUM ---------------------------------------------
        % v[n] = sum_{q=0..P-1} xc[n - q*L_period].
        % Add P complex-valued matched-filter outputs spaced exactly one STS
        % period apart. For an STS-aligned signal each tap contributes the
        % same complex peak -> amplitude grows by P. For an off-comb narrowband
        % interferer the P phasors land on a length-P DFT bin away from a
        % passband peak -> destructive interference (jammer rejected by ~P).
        v_iq = xc_iq;
        for q_idx = 1 : (n_coh - 1)
            shift_idx = q_idx * l_period_idx;
            if shift_idx >= n_iq
                break
            end
            v_iq(shift_idx + 1 : end) = v_iq(shift_idx + 1 : end) ...
                + xc_iq(1 : n_iq - shift_idx);
        end

        per_freq_v_pwr(k, :) = abs(v_iq).^2;
    end

    % --- Denominator: input energy over the combed support ----------------
    % The effective filter (template tiled P times at stride L_period) has
    % support P*L_period samples, so the right CS-tight denominator is the
    % input energy summed over the SAME window.
    den_int = sliding_sum_1d(mag2_iq, n_den_win);

    eps_floor = 1e-20;
    % Cauchy-Schwarz: |v[n]|^2 <= P * sum|r|^2; normalise by P so metric in [0,1].
    per_freq_pwr = per_freq_v_pwr ./ (n_coh * den_int + eps_floor);

    % Boundary: the first (n_den_win-1) samples have a partial window in BOTH
    % numerator (comb not fully filled) and denominator (sliding sum not yet
    % saturated). Zero them so the startup transient cannot masquerade as a
    % peak (small den_int can otherwise blow the ratio up).
    n_zero_head = min(n_den_win - 1, n_iq);
    per_freq_pwr(:, 1:n_zero_head) = 0;
    per_freq_v_pwr(:, 1:n_zero_head) = 0;

    [metric_pwr, k_argmax] = max(per_freq_pwr, [], 1);
    freq_argmax = freq_grid(k_argmax);

    % Per-sample peak |v| across hypotheses (combed magnitude, un-normalised).
    % This is the natural counterpart of the autocorr's |M_hat| floor: in
    % noise-only frames its mean^2 is ~ P*L_template*sigma^4, in signal frames
    % it climbs to ~ P*L_template*signal_pwr*sigma^2. Downstream code uses it
    % as the absolute-level gate (calibrated empirically in Phase A).
    xc_mag_max_per_sample = sqrt(max(per_freq_v_pwr, [], 1));

end

function y = sliding_sum_1d(x, n_win)
% Causal sliding sum of length n_win over a row vector.
% Same-length output. For k < n_win the window is shorter (k samples).
    n_cols = numel(x);
    cs = [0, cumsum(x)];
    idx_start_0 = max(0, (0:n_cols-1) - n_win + 1);
    y = cs(2:end) - cs(idx_start_0 + 1);
end
