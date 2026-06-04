function [stf_offset_start, stf_offset_end, diag] = ...
    compute_stf_window_offsets(p_init, tx_iq_base, fs_tx_freq, fs_bb_freq, fs_rx_freq)
%COMPUTE_STF_WINDOW_OFFSETS  Locate the L-STF in rx_iq buffer (analytic + empirical).
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
%
% For a TX waveform built as
% tx_full = [zeros(1, guard_samples_tx), tx_iq_base]
% then optionally fractionally delayed by tau in [0, osr_tx) samples,
% passed through CFO + AWGN/jammer add, and decimated by 2:1 twice
% (80 -> 40 -> 20 MHz half-band cascade), this function returns offsets
% such that the L-STF lives in rx_iq at fs_rx indices
%
% n_stf_start_rx = guard_samples_rx + stf_offset_start
% n_stf_end_rx = guard_samples_rx + stf_offset_end
%
% The window is widened by 1 sample on each end to absorb the [0, 1)
% fs_rx sample jitter from the random fractional delay tau.
%
% The analytic prediction (derived from the halfband cascade trim that
% advances the output by (n_taps-1)/2 per stage) is cross-checked against
% an envelope measurement on a noise-free, jammer-free, CFO-free, tau=0
% reference pass. If the empirical leading-edge index differs from the
% prediction by more than 5 samples, a warning is issued.
%
% Inputs:
% p_init sim_params struct (used to grab L-STF info from tx_frame)
% tx_iq_base [1 x N] complex baseband at fs_tx, unit-power normalised
% fs_tx_freq TX rate (Hz)
% fs_bb_freq baseband rate (Hz)
% fs_rx_freq RX rate (Hz)
%
% Outputs:
% stf_offset_start L-STF leading offset (fs_rx samples) vs guard_samples_rx
% stf_offset_end L-STF trailing offset (fs_rx samples) vs guard_samples_rx
% diag struct with diagnostic fields:
% .p_stf_start_rx_pred analytic L-STF start (fs_rx, may be fractional)
% .p_stf_end_rx_pred analytic L-STF end (fs_rx, may be fractional)
% .n_start_meas empirical envelope leading-edge sample
% .n_envelope_above empirical envelope width (L-STF + L-LTF)
% .guard_test_rx the guard used for the empirical run
%
% The returned offsets are integer-valued (floor on start, ceil on end)
% and include +/- 1 sample of slack for the [0, 1) fs_rx fractional-delay
% jitter on the leading and trailing edge.

% --- L-STF location in tx_iq_base at fs_tx (from tx_frame info) ---
    [~, info] = tx_frame(p_init);
    stf_start_tx_in_base = info.stf_start_tx;
    stf_end_tx_in_base = info.stf_end_tx;

    % --- Reference pass: noise-free, jammer-free, no CFO, tau = 0 -----
    guard_test_rx = 500;
    osr_tx = round(fs_tx_freq / fs_bb_freq);
    osr_tx_to_rx = round(fs_tx_freq / fs_rx_freq);
    guard_test_tx = guard_test_rx * osr_tx_to_rx;
    tx_full = [zeros(1, guard_test_tx), tx_iq_base];
    rx_iq = decimate_hb_2to1(tx_full, 11);
    rx_iq = decimate_hb_2to1(rx_iq, 23);

    % --- Analytic L-STF position at fs_rx -----------------------------
    % The two-stage halfband cascade is the composition of the per-stage
    % map n_in -> (n_in + 1 - (n_taps-1)/2) / 2 (in 1-indexed samples
    % at the input rate). For n_taps_1 = 11, n_taps_2 = 23, decimating
    % 80 -> 40 -> 20 MHz, the cascade maps an event at fs_tx position
    % n_tx to fs_rx position (n_tx - 24) / 4.
    p_stf_start_tx_full = guard_test_tx + stf_start_tx_in_base;
    p_stf_end_tx_full = guard_test_tx + stf_end_tx_in_base;
    p_stf_start_rx_pred = (p_stf_start_tx_full - 24) / 4;
    p_stf_end_rx_pred = (p_stf_end_tx_full - 24) / 4;

    % --- Empirical envelope cross-check ------------------------------
    % The envelope spans the WHOLE preamble (L-STF + L-LTF), so we only
    % validate the leading edge (which is the L-STF start).
    env = filter(ones(1, 8) / 8, 1, abs(rx_iq).^2);
    thresh_env = 0.05 * max(env);
    idx_above = find(env > thresh_env);
    if isempty(idx_above)
        warning('compute_stf_window_offsets:NoEnvelope', ...
            'Could not find any rx_iq envelope above 5%% of peak.');
        n_start_meas = NaN;
        n_envelope_above = 0;
    else
        n_start_meas = idx_above(1);
        n_envelope_above = idx_above(end) - idx_above(1) + 1;
        if abs(n_start_meas - p_stf_start_rx_pred) > 5
            warning('compute_stf_window_offsets:LeadingEdgeMismatch', ...
                ['Envelope leading edge at fs_rx index %d differs ' ...
                'from analytic %.2f by more than 5 samples'], ...
                n_start_meas, p_stf_start_rx_pred);
        end
    end

    % --- Final offsets with +/-1 sample of slack (covers tau in [0,1) fs_rx)
    %     The leading-edge offset is anchored to the cascade-FIR PEAK
    %     position (analytic GD-removed mapping), which has been
    %     independently verified to equal the L-STF SoP position obtained
    %     by matched-filtering the noiseless rx_iq with the 160-sample
    %     canonical L-STF template: SoP_emp = argmax|xc| - 160 + 1.
    %     See diag_sop_from_xc.m for the verification.
    %     The -1 absorbs the integer-part jitter of the fractional-delay
    %     tau in [0, OSR_tx) -> [0, 1) fs_rx.
    stf_offset_start = floor(p_stf_start_rx_pred) - guard_test_rx - 1;
    stf_offset_end = ceil(p_stf_end_rx_pred) - guard_test_rx + 1;

    diag = struct( ...
        'p_stf_start_rx_pred', p_stf_start_rx_pred, ...
        'p_stf_end_rx_pred', p_stf_end_rx_pred, ...
        'n_start_meas', n_start_meas, ...
        'n_envelope_above', n_envelope_above, ...
        'guard_test_rx', guard_test_rx);
end
