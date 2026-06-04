function lltf_iq = make_lltf(p)
%MAKE_LLTF  Generate the 802.11a/g/n L-LTF time-domain waveform at fs_bb.
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
%
% IEEE Ref : IEEE Std 802.11-2016, Eq. 17-11.
%
% Algorithm:
% 1. Build the frequency-domain L-LTF on subcarriers -26..+26: 52 active
% BPSK +/-1 (DC = 0).
% 2. IFFT scaled by n_fft / sqrt(n_tone) per Eq. 17-6.
% 3. Pre-pend a 32-sample double-guard interval (GI2 = second half of T1)
% and append two 64-sample long-symbol periods (T1, T2) -> 160 samples.
%
% Inputs:
% p [struct] sim_params (.n_fft, .n_tone)
%
% Outputs:
% lltf_iq [1 x 160] full L-LTF time-domain at fs_bb (20 MHz)
%
% Notes:
% The L-LTF is not used by the SoP detectors evaluated in this repo,
% but is included so that a realistic frame can be transmitted for
% visual inspection and future extension.
%
% Role : TX waveform / preamble
% Phase : 1 (Floating-Point, frame-based)
%
% See also: make_lstf, tx_frame

    n_fft = p.n_fft;
    n_tone = p.n_tone;

    arr_lltf_fd_vals = [ 1, 1,-1,-1, 1, 1,-1, 1,-1, 1, 1, 1, 1, 1, 1,-1,-1, 1, 1, ...
        -1, 1,-1, 1, 1, 1, 1, 0, 1,-1,-1, 1, 1,-1, 1,-1, 1,-1,-1,-1, ...
        -1,-1, 1, 1,-1,-1, 1,-1, 1,-1, 1, 1, 1, 1];

    lltf_fd = zeros(1, n_fft);
    sc_idx = -26:26;
    sc_bin_idx = mod(sc_idx + n_fft, n_fft);
    for k = 1:length(sc_idx)
        lltf_fd(sc_bin_idx(k) + 1) = arr_lltf_fd_vals(k);
    end

    td_one_period_iq = ifft(lltf_fd, n_fft) * n_fft / sqrt(n_tone);
    gi2_iq = td_one_period_iq(n_fft/2+1:end);
    lltf_iq = [gi2_iq, td_one_period_iq, td_one_period_iq]; % 160 samples

end
