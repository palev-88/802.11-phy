function [tx_iq, info] = tx_frame(p)
%TX_FRAME  Build a frame-based TX waveform at fs_tx.
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
%
% IEEE Ref : Frame layout follows IEEE Std 802.11-2016, Sec. 17.3.3
% (preamble) -- L-STF + L-LTF, no L-SIG / DATA in this study.
%
% Algorithm:
% 1. Generate L-STF and L-LTF at fs_bb (20 MHz).
% 2. Normalise both so that mean(|lstf_iq|^2) = 1 (unit-power reference
% for the SNR / JNR definitions in channel_top).
% 3. Concatenate [zero guard | L-STF | L-LTF | zero guard].
% 4. Upsample to fs_tx (default 80 MHz, 4x) using a polyphase LPF.
%
% Inputs:
% p [struct] sim_params (.fs_tx, .fs_bb, .os_tx, .n_stf,
% .guard_samples_bb, .n_ltf)
%
% Outputs:
% tx_iq [1 x n_iq] complex baseband at fs_tx
% info [struct] frame layout info:
% .scale_dc unit-power normalisation factor
% .stf_start_tx first L-STF sample index (1-based) in tx_iq
% .stf_end_tx last L-STF sample index in tx_iq
% .preamble_start_tx alias of stf_start_tx
% .preamble_end_tx last preamble (LTF) sample index in tx_iq
% .fs fs_tx
% .frame_len numel(tx_iq)
%
% Role : TX waveform / framing
% Phase : 1 (Floating-Point, frame-based)
%
% See also: make_lstf, make_lltf, resample_int

% --- preamble at baseband ----------------------------------------------
    lstf_iq = make_lstf(p);
    lltf_iq = make_lltf(p);

    % Normalise to mean(|lstf|^2) == 1 (SNR/JNR reference).
    stf_pwr = mean(abs(lstf_iq).^2);
    scale_dc = 1 / sqrt(stf_pwr);
    lstf_iq = lstf_iq * scale_dc;
    lltf_iq = lltf_iq * scale_dc;

    guard_iq = zeros(1, p.guard_samples_bb);
    frame_iq_bb = [guard_iq, lstf_iq, lltf_iq, guard_iq];

    stf_start_bb = p.guard_samples_bb + 1;
    stf_end_bb = stf_start_bb + p.n_stf - 1;

    % --- upsample baseband -> fs_tx ----------------------------------------
    tx_iq = resample_int(frame_iq_bb, p.os_tx, 1);

    % Indices in the upsampled stream. Build the info struct in a single
    % struct(...) call so the struct's field set is fixed at first assignment
    % (field-by-field construction with intermediate reads is rejected).
    stf_start_tx_val = (stf_start_bb - 1) * p.os_tx + 1;
    stf_end_tx_val = stf_end_bb * p.os_tx;
    preamble_end_tx_val = stf_end_bb * p.os_tx + p.n_ltf * p.os_tx;
    info = struct( ...
        'scale_dc', scale_dc, ...
        'stf_start_tx', stf_start_tx_val, ...
        'stf_end_tx', stf_end_tx_val, ...
        'preamble_start_tx', stf_start_tx_val, ...
        'preamble_end_tx', preamble_end_tx_val, ...
        'fs', p.fs_tx, ...
        'frame_len', numel(tx_iq));

end
