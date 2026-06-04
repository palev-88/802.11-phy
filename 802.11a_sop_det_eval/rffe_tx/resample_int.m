function y_iq = resample_int(x_iq, l_up, m_dn, n_taps)
%RESAMPLE_INT  Integer-ratio resampler (from first principles).
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
%
% Algorithm:
% 1. Zero-stuff x_iq by l_up (insert l_up-1 zeros between each sample).
% 2. LPF with a windowed-sinc filter cut at min(1/(2*l_up), 1/(2*m_dn))
% on the intermediate fs = l_up * fs_in. A 5% guard inside the
% Nyquist edge prevents passband ripple at the corner.
% 3. Compensate the linear-phase group delay (n_taps-1)/2 samples.
% 4. Downsample by m_dn.
%
% Inputs:
% x_iq [1 x n_in] input samples (real or complex)
% l_up [1 x 1] upsampling factor (>=1)
% m_dn [1 x 1] downsampling factor (>=1)
% n_taps [1 x 1] optional FIR length (default = 8*max(l_up,m_dn)+1, odd)
%
% Outputs:
% y_iq [1 x n_out] resampled samples
%
% Role : Common DSP utility (rate conversion)
% Phase : 1 (Floating-Point, frame-based)
%
% See also: design_lpf

    if nargin < 4 || isempty(n_taps)
        n_taps = 8 * max(l_up, m_dn) + 1;
        if mod(n_taps, 2) == 0
            n_taps = n_taps + 1;
        end
    end

    x_iq = x_iq(:).';

    % --- 1. zero-stuff (upsample by l_up) ---------------------------------
    if l_up > 1
        n_in = numel(x_iq);
        xs_iq = zeros(1, n_in*l_up);
        xs_iq(1:l_up:end) = x_iq;
    else
        xs_iq = x_iq;
    end

    % --- 2. anti-image / anti-alias LPF -----------------------------------
    fc_norm = 0.95 * min(1/(2*l_up), 1/(2*m_dn));
    hir_lpf = design_lpf(n_taps, fc_norm);
    hir_lpf = hir_lpf * l_up; % gain compensation

    xf_iq = conv(xs_iq, hir_lpf);
    n_gd = (n_taps - 1) / 2;
    xf_iq = xf_iq(n_gd + 1 : n_gd + numel(xs_iq));

    % --- 3. downsample by m_dn --------------------------------------------
    if m_dn > 1
        y_iq = xf_iq(1:m_dn:end);
    else
        y_iq = xf_iq;
    end

end
