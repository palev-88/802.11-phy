function [lstf_iq, short_sym_iq] = make_lstf(p)
%MAKE_LSTF  Generate the 802.11a/g/n L-STF time-domain waveform at fs_bb.
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
%
% IEEE Ref : IEEE Std 802.11-2016, Eq. 17-10 (L-STF subcarriers and Eq. 17-6
% for the IFFT normalisation).
%
% Algorithm:
% 1. Build the frequency-domain L-STF on subcarriers -26..+26: 12 active
% subcarriers spaced every 4 in frequency, normalised by sqrt(13/6).
% 2. Take the 64-point IFFT, scaled by n_fft / sqrt(n_tone) per Eq. 17-6.
% 3. Take one short-symbol period (n_stf_period = 16 samples) and tile
% 10x to form the 160-sample L-STF.
%
% Inputs:
% p [struct] sim_params (.n_fft, .n_tone, .n_stf_period)
%
% Outputs:
% lstf_iq [1 x 160] full L-STF time-domain at fs_bb (20 MHz)
% short_sym_iq [1 x 16] one short-symbol period at fs_bb (20 MHz)
%
% Role : TX waveform / preamble
% Phase : 1 (Floating-Point, frame-based)
%
% See also: make_lltf, tx_frame

    n_fft = p.n_fft;
    n_tone = p.n_tone;

    % Frequency-domain L-STF on subcarriers -26..+26 (mapped into 64-pt FFT bins).
    arr_stf_fd_vals = sqrt(13/6) * ...
        [0, 0, 1+1j, 0, 0, 0, -1-1j, 0, 0, 0, 1+1j, 0, 0, 0, -1-1j, 0, ...
        0, 0, -1-1j, 0, 0, 0, 1+1j, 0, 0, 0, 0, 0, 0, 0, -1-1j, 0, ...
        0, 0, -1-1j, 0, 0, 0, 1+1j, 0, 0, 0, 1+1j, 0, 0, 0, 1+1j, 0, ...
        0, 0, 1+1j, 0, 0];

    lstf_fd = complex(zeros(1, n_fft)); % must be complex (assigned to below)
    sc_idx = -26:26;
    sc_bin_idx = mod(sc_idx + n_fft, n_fft);
    for k = 1:length(sc_idx)
        lstf_fd(sc_bin_idx(k) + 1) = arr_stf_fd_vals(k);
    end

    % IFFT with the IEEE 802.11-2016 Eq. 17-6 normalisation.
    td_one_period_iq = ifft(lstf_fd, n_fft) * n_fft / sqrt(n_tone);

    % One short-symbol period (16 samples), tiled 10x => 160 samples.
    short_sym_iq = td_one_period_iq(1:p.n_stf_period);
    lstf_iq = repmat(short_sym_iq, 1, 10);

end
