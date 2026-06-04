function jam_iq = make_jammer(jam_type_id, jam_pwr, jam_freq_hz, ...
    n_samples, fs_tx_freq, fs_bb_freq, ...
    chirp_f_start_freq, chirp_f_end_freq)
%MAKE_JAMMER  jammer waveform synthesis @ fs_tx.
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
% jam_type_id: 1 = cw, 2 = noise_1mhz (filtered Gaussian, 1 MHz BW),
% 3 = wideband_gauss (white over fs_tx),
% 4 = chirp (linear, [f_start..f_end]).
% jam_pwr is referenced to fs_bb (matches the SNR/JNR baseband convention).

    n_idx_v = 0 : n_samples - 1;
    switch jam_type_id
    case 1 % cw
        jam_iq = sqrt(jam_pwr) * ...
            exp(1j*2*pi*jam_freq_hz/fs_tx_freq * n_idx_v);
    case 2 % noise_1mhz
        ble_bw_hz = 1e6;
        n_lpf_taps = 51;
        h_lpf = design_lpf(n_lpf_taps, ble_bw_hz/2 / fs_tx_freq);
        h_lpf = h_lpf / sum(h_lpf);
        sigma_sq = jam_pwr / sum(abs(h_lpf).^2);
        bb_noise = sqrt(sigma_sq/2) * ...
            (randn(1, n_samples + n_lpf_taps) + ...
            1j*randn(1, n_samples + n_lpf_taps));
        bb_filt = filter(h_lpf, 1, bb_noise);
        bb_filt = bb_filt(n_lpf_taps + (1:n_samples));
        jam_iq = bb_filt .* exp(1j*2*pi*jam_freq_hz/fs_tx_freq * n_idx_v);
    case 3 % wideband_gauss
        jam_pwr_tx = jam_pwr * (fs_tx_freq / fs_bb_freq);
        jam_iq = sqrt(jam_pwr_tx/2) * ...
            (randn(1, n_samples) + 1j*randn(1, n_samples));
    case 4 % chirp
        f0 = chirp_f_start_freq;
        f1 = chirp_f_end_freq;
        t_v = n_idx_v / fs_tx_freq;
        t_sweep = n_samples / fs_tx_freq;
        phi_v = 2*pi * (f0 * t_v + 0.5 * (f1 - f0)/t_sweep * t_v.^2);
        jam_iq = sqrt(jam_pwr) * exp(1j * (phi_v + 2*pi*rand()));
    otherwise
        jam_iq = complex(zeros(1, n_samples));
    end
end
