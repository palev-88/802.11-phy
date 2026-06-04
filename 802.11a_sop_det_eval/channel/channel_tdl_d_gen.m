function [fading, delays_samples, max_delay] = channel_tdl_d_gen( ...
    N, fs_hz, delay_spread_s, f_doppler_hz)
%CHANNEL_TDL_D_GEN  Generate TDL-D fading taps + delays (no application).
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
% Returns:
% fading [n_taps x N] complex per-tap time-varying gains
% delays_samples [1xn_t] integer per-tap delays at fs_hz
% max_delay scalar = max(delays_samples)
%
% Splits cleanly from channel_tdl_d_apply so the generation cost
% (random draws + FFT-based Jakes shaping for each tap) is paid
% ONCE per trial and the same realisation is applied across all
% SNR / JNR / detector combinations in the sweep (Common Random
% Numbers variance reduction + huge compute saving).

    delays_norm = [0.0000, 0.0000, 0.0350, 0.6120, 1.3630, ...
        1.4050, 1.8040, 2.5960, 1.7750, 4.0420, ...
        7.9370, 9.4240, 9.7080];
    powers_dB = [-0.2, -13.5, -18.8, -21.0, -22.8, ...
        -17.9, -20.1, -21.9, -22.9, -27.8, ...
        -23.6, -24.8, -30.0];
    K_dB = 13.3;
    n_taps = numel(delays_norm);

    delays_s = delays_norm * delay_spread_s;
    delays_samples = round(delays_s * fs_hz);
    max_delay = max(delays_samples);

    powers_lin = 10.^(powers_dB / 10);
    powers_lin = powers_lin / sum(powers_lin);
    amplitudes = sqrt(powers_lin);

    K_lin = 10^(K_dB / 10);
    fading = complex(zeros(n_taps, N));
    t = (0 : N-1) / fs_hz;

    if f_doppler_hz == 0
        for i = 1:n_taps
            if i == 1
                phi_los = 2*pi*rand();
                h_complex = (randn()+ 1j*randn()) / sqrt(2);
                h_los = sqrt(K_lin/(K_lin+1)) * exp(1j*phi_los);
                h_nlos = sqrt(1/(K_lin+1)) * h_complex;
                fading(i, :) = amplitudes(i) * (h_los + h_nlos);
            else
                phi = 2*pi*rand();
                fading(i, :) = amplitudes(i) * exp(1j*phi);
            end
        end
    else
        df = fs_hz / N;
        k_idx = 0 : N-1;
        f_axis = k_idx * df;
        f_axis(f_axis > fs_hz/2) = f_axis(f_axis > fs_hz/2) - fs_hz;
        mask_in = abs(f_axis) < f_doppler_hz;
        H_jakes = zeros(1, N);
        H_jakes(mask_in) = 1 ./ sqrt(1 - (f_axis(mask_in) / f_doppler_hz).^2);
        H_jakes = H_jakes / sqrt(sum(H_jakes.^2) / N);

        phi_los_outer = 2*pi*rand();
        h_los_carrier = exp(1j*(2*pi*f_doppler_hz*t + phi_los_outer));

        for i = 1:n_taps
            W = (randn(1, N) + 1j*randn(1, N)) / sqrt(2);
            h_complex = ifft(fft(W) .* H_jakes);

            if i == 1
                h_los_w = sqrt(K_lin/(K_lin+1)) * h_los_carrier;
                h_nlos_w = sqrt(1/(K_lin+1)) * h_complex;
                fading(i,:) = amplitudes(i) * (h_los_w + h_nlos_w);
            else
                fading(i,:) = amplitudes(i) * h_complex;
            end
        end
    end
end
