function p = sim_params(fs_bb, fs_tx, fs_rx, ...
    xcb_freq_min, xcb_freq_max, xcb_freq_step)
%SIM_PARAMS  Build the simulator parameter struct.
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
% IEEE Ref: derived from IEEE Std 802.11-2016 Sec. 17.3.3 (OFDM PHY
% preamble: L-STF / L-LTF subcarrier and timing definitions).
%
% Returns a flat, immutable struct that fixes every constant the rest
% of the simulator reads: sample rates, FFT / preamble sizing, the
% pre-packet zero-pad length, and the cross-correlator-bank knobs.
% All fields are set in a single struct(...) literal so the field set
% is fixed at first assignment and callers cannot accidentally add new
% names that the rest of the pipeline does not know about.
%
% Inputs:
% fs_bb [Hz] baseband sample rate (= 20e6 for 802.11a/g).
% fs_tx [Hz] upsampled TX rate (oversampling ratio
% fs_tx / fs_bb must be a positive integer;
% 4x oversampling at fs_tx = 80e6 is standard).
% fs_rx [Hz] RX sample rate after the decimator cascade
% (= fs_bb in this study).
% xcb_freq_min [Hz] XC-bank lowest frequency hypothesis.
% xcb_freq_max [Hz] XC-bank highest frequency hypothesis.
% xcb_freq_step [Hz] XC-bank hypothesis spacing.
% (Use (-145e3, +145e3, 145e3) for the
% 3-hypothesis 802.11a/g default; pass
% (0, 0, 1) to collapse the bank to a single
% matched filter at DC.)
%
% Output:
% p struct with the fields:
% fs_bb / fs_tx / fs_rx sample rates [Hz]
% os_tx / os_rx integer oversampling ratios
% n_fft = 64 OFDM FFT size (IEEE Eq. 17-6)
% n_tone = 52 active subcarriers per L-LTF
% n_stf = 160 L-STF length at fs_bb [samples]
% n_ltf = 160 L-LTF length at fs_bb [samples]
% n_stf_period = 16 short-symbol period at fs_bb
% guard_samples_bb = 64 pre/post-packet zero-pad
% xcb_template_len_bb = 16 MF template length at fs_bb
% xcb_template_max_level = 3 3-bit symmetric template quant
% xcb_coh_periods = 4 coherent-comb depth P
% xcb_freq_{min,max,step} hypothesis-grid bounds [Hz]
%
% Role : Simulator setup / configuration block.
%
% See also: tx_frame, rx_lstf_xcorr_bank.

% Pre-compute derived values OUTSIDE the struct(...) call so no
% field of p is read before construction.
    os_tx = round(fs_tx / fs_bb);
    os_rx = round(fs_rx / fs_bb);

    p = struct( ...
        'fs_bb', fs_bb, ...
        'fs_tx', fs_tx, ...
        'fs_rx', fs_rx, ...
        'os_tx', os_tx, ...
        'os_rx', os_rx, ...
        'n_fft', 64, ...
        'n_tone', 52, ...
        'n_stf', 160, ...
        'n_ltf', 160, ...
        'n_stf_period', 16, ...
        'guard_samples_bb', 64, ...
        'xcb_template_len_bb', 16, ...
        'xcb_template_max_level', 3, ...
        'xcb_coh_periods', 4, ...
        'xcb_freq_min', xcb_freq_min, ...
        'xcb_freq_max', xcb_freq_max, ...
        'xcb_freq_step', xcb_freq_step);
end
