function y_iq = fractional_delay(x_iq, tau_samples, n_taps)
%FRACTIONAL_DELAY  Apply a fractional sample-delay via windowed-sinc FIR.
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
%
% IEEE Ref : Not normative -- classical fractional-delay interpolator,
% see Laakso et al., "Splitting the unit delay: tools for
% fractional delay filter design," IEEE Signal Processing
% Magazine, vol. 13, no. 1, pp. 30-60, Jan. 1996.
%
% Algorithm:
% 1. Design a length-n_taps Kaiser-windowed sinc with delay center at
% index ( (n_taps-1)/2 + tau_samples ). This centres the main lobe
% inside the filter span when tau_samples is small (e.g., < 5).
% 2. Filter x_iq with this FIR. Output sample n approximates
% x_iq[n - ((n_taps-1)/2 + tau_samples)].
% 3. Shift the output left by the constant integer delay (n_taps-1)/2
% so the kept output represents EXACTLY x_iq delayed by tau_samples.
%
% A tau_samples = 0 reduces to a near-identity filter (slight low-pass
% from the Kaiser window). For OFDM signals occupying < 25% of the
% Nyquist band, the passband ripple is < 0.05 dB and irrelevant.
%
% Inputs:
% x_iq [1 x n_iq] complex baseband samples (any rate).
% tau_samples scalar fractional delay in SAMPLES (>= 0). For the
% SPO use case, tau in [0, osr) at the TX rate.
% n_taps scalar FIR length; must be odd. Default 31 (gives
% < 0.05 dB amplitude error across the 802.11
% 20 MHz BB at fs_tx = 80 MHz).
%
% Outputs:
% y_iq [1 x n_iq] delayed signal, same length as input. Leading
% samples that would require x_iq from before
% index 1 are zero (filter startup transient
% region).
%
% Verification (see verify_fractional_delay.m):
% For x_iq[n] = exp(j*2*pi*f0*n/fs), the output should obey
% y_iq[n] = x_iq[n - tau_samples] = exp(j*phi) * x_iq[n] with
% phi = -2*pi*f0*tau_samples/fs. The script measures phi empirically
% and compares with this analytical value across a frequency grid.
%
% Role : RF/clock impairment model (sampling phase offset, SPO).
% Phase : 1 (Floating-Point, frame-based).
%
% See also: decimate_hb_2to1, design_lpf

    if nargin < 3 || isempty(n_taps)
        n_taps = 31;
    end
    if mod(n_taps, 2) == 0
        n_taps = n_taps + 1; % force odd length
    end
    if tau_samples < 0
        error('fractional_delay:NegativeTau', 'tau_samples must be >= 0');
    end

    x_iq = x_iq(:).';
    n_iq_idx = numel(x_iq);

    % ---- Design the windowed-sinc filter centred at (center + tau_samples).
    n_center_idx = (n_taps - 1) / 2;
    k_idx = (0:n_taps - 1);
    tau_total = n_center_idx + tau_samples;

    % Ideal fractional-delay impulse response = sinc shifted by tau_total.
    % sinc(x) = sin(pi*x) / (pi*x), with sinc(0) = 1.
    hir_sinc = sinc_local(k_idx - tau_total);

    % Kaiser window: beta = 8 gives stop-band ripple < -80 dB.
    beta_kaiser = 8;
    hir_kaiser = kaiser_local(n_taps, beta_kaiser);

    hir_total = hir_sinc .* hir_kaiser;

    % Normalise to unit DC gain so a static-DC input passes through unchanged.
    hir_total = hir_total / sum(hir_total);

    % ---- Apply the filter.
    y_full_iq = filter(hir_total, 1, x_iq);

    % ---- Compensate the constant integer delay (n_center_idx) introduced by
    % the filter so the kept output reflects ONLY tau_samples of delay.
    y_iq = [y_full_iq(n_center_idx + 1 : end), zeros(1, n_center_idx)];
end

function y = sinc_local(x)
% sinc: sin(pi*x)/(pi*x), with sinc(0) = 1 by convention.
    y = ones(size(x));
    mask_idx = (x ~= 0);
    y(mask_idx) = sin(pi * x(mask_idx)) ./ (pi * x(mask_idx));
end

function arr_w = kaiser_local(n_taps, beta_kaiser)
% Kaiser window:
% w[k] = besseli(0, beta * sqrt(1 - ((k - alpha) / alpha)^2)) / besseli(0, beta),
% alpha = (n_taps - 1) / 2, k = 0, ..., n_taps - 1.
    alpha_idx = (n_taps - 1) / 2;
    k_idx = (0 : n_taps - 1);
    arg_in = beta_kaiser * sqrt( max(0, 1 - ((k_idx - alpha_idx) / alpha_idx).^2) );
    arr_w = besseli(0, arg_in) / besseli(0, beta_kaiser);
end
