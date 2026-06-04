function y_iq = decimate_hb_2to1(x_iq, n_taps)
%DECIMATE_HB_2TO1  2:1 decimator built around a half-band FIR.
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
%
% IEEE Ref : Not normative -- standard half-band decimator design.
%
% Algorithm:
% 1. Design a half-band FIR by windowing the ideal sinc at fc = fs/4.
% For a sinc with cutoff exactly 1/4 cycles/sample, the impulse
% response has mathematical zeros at every even index away from
% the centre (the sinc nulls land on integer multiples of 2),
% which is the half-band property: roughly half the taps are zero
% so a real implementation does ~N/2 multiplies per output sample.
% 2. PRE-PAD x_iq with (n_taps-1) zeros so the filter's startup
% transient happens INSIDE the zero-padded region. Without this,
% a large stopband input (e.g. an out-of-band jammer) produces a
% visible spurious peak at the first kept output sample because
% the filter takes ~n_taps samples to build full stopband rejection.
% 3. Skip the transient portion of the convolution and keep n_iq_in
% steady-state samples. The output is NOT group-delay-compensated:
% each output sample is delayed by (n_taps-1) input samples from
% the "centered" alignment. The caller is responsible for adding
% a constant offset (either computed analytically or measured
% empirically from a clean test signal) to align downstream
% sample-index references (e.g. STF window position).
% 4. Decimate by 2.
%
% Inputs:
% x_iq [1 x n_iq_in] input samples (real or complex)
% n_taps [1 x 1] optional FIR length (default 23). Forced to
% satisfy mod(n_taps, 4) == 3 so the half-band
% zero-coefficient symmetry lines up exactly.
%
% Outputs:
% y_iq [1 x ceil(n_iq_in/2)] decimated samples (uncompensated:
% sample j corresponds to the filter's response
% AFTER it has seen input sample (2j - 1 + n_taps - 1)
% at the pre-decimation rate).
%
% Notes:
% - This is the building block of a cascaded multi-stage decimator:
% to drop by 4 (e.g. 80 -> 20 MHz), call decimate_hb_2to1 twice; to
% drop by 8, three times. Each stage operates at progressively lower
% rates, which is what makes the cascade compute-efficient vs. a
% single high-rate FIR.
% - The first-stage stopband at the operating rate is [fs/4, fs/2];
% this is exactly what the cascade requires (everything in the
% upper half of the input Nyquist band is rejected before decim).
% - Cumulative delay through a K-stage cascade in seconds:
% t_delay = (n_taps - 1) * sum_{k=1..K} (1/fs_stage_k)
% At a final rate fs_out this is (t_delay * fs_out) samples.
% Use this to compute the constant offset for the alignment.
%
% Role : Common DSP utility (anti-alias decimation, half-band)
% Phase : 1 (Floating-Point, frame-based)
%
% See also: design_lpf, resample_int

    if nargin < 2 || isempty(n_taps)
        n_taps = 23;
    end
    % Force mod(n_taps,4)==3 so that with centre at (n_taps-1)/2 (even),
    % the taps at +/-2, +/-4, ... land on the sinc zeros exactly. This is the
    % half-band property: ~half the coefficients are mathematically zero.
    n_taps = max(n_taps, 11);
    if mod(n_taps, 4) ~= 3
        n_taps = n_taps + (3 - mod(n_taps, 4));
        if n_taps < 11
            n_taps = n_taps + 4;
        end
    end

    % Half-band FIR: windowed sinc at fc = 1/4 cycles/sample (the natural
    % half-band cutoff). design_lpf already DC-normalises so sum(hir_hb) = 1.
    hir_hb = design_lpf(n_taps, 0.25);

    x_iq = x_iq(:).';
    n_iq_in = numel(x_iq);
    n_pad_idx = n_taps - 1; % minimal pad to absorb conv startup transient

    % --- Option 2: pre-pad with (n_taps-1) zeros. The convolution startup
    % transient now happens entirely within the zero-padded region, NOT
    % in the kept output. This is what prevents the spurious leading
    % peak that a large stopband input would otherwise produce.
    arr_x_padded_iq = [zeros(1, n_pad_idx), x_iq];
    arr_xf_iq = conv(arr_x_padded_iq, hir_hb); % length n_iq_in + 2*(n_taps-1)

    % --- Trim: keep n_iq_in samples starting at conv index 2*n_taps - 1
    % (the first index at which ALL n_taps filter taps overlap with the
    % real input -- no zero pad in the window, hence no transient).
    % The output is NOT group-delay-compensated. Each output sample j
    % corresponds to the filter's response after seeing input sample
    % (j + n_taps - 1) of the original x_iq -- i.e. a built-in delay of
    % (n_taps - 1) input samples vs the centered alignment.
    arr_xf_steady_iq = arr_xf_iq(2*n_taps - 1 : 2*n_taps - 1 + n_iq_in - 1);

    % Decimate by 2
    y_iq = arr_xf_steady_iq(1:2:end);
end
