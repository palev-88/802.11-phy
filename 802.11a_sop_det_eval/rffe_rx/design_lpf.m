function hir_lpf = design_lpf(n_taps, fc_norm)
%DESIGN_LPF  Linear-phase windowed-sinc low-pass filter (from first principles).
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
%
% Algorithm:
% 1. Build the centred sample index n = (0..n_taps-1) - (n_taps-1)/2.
% 2. Ideal sinc impulse response h_ideal = 2*fc_norm * sinc(2*fc_norm * n).
% 3. Apply a Hamming window.
% 4. DC-normalise: hir_lpf = h_windowed / sum(h_windowed).
%
% Inputs:
% n_taps [1 x 1] filter length (odd recommended for integer group delay)
% fc_norm [1 x 1] cutoff in cycles/sample (in (0, 0.5))
%
% Outputs:
% hir_lpf [1 x n_taps] FIR impulse response, sum(hir_lpf) = 1
%
% Role : Common DSP utility (anti-alias / interpolation LPF)
% Phase : 1 (Floating-Point, frame-based)
%
% See also: resample_int

    if nargin < 2
        error('design_lpf:nargs', 'Need n_taps and fc_norm');
    end
    if fc_norm <= 0 || fc_norm >= 0.5
        error('design_lpf:fc', 'fc_norm must lie in (0, 0.5)');
    end

    n_centered = (0:n_taps-1) - (n_taps-1)/2;
    hir_ideal = 2*fc_norm * sinc_local(2*fc_norm * n_centered);

    m_idx = 0:n_taps-1;
    hir_win = 0.54 - 0.46*cos(2*pi*m_idx/(n_taps-1));

    hir_lpf = hir_ideal .* hir_win;
    hir_lpf = hir_lpf / sum(hir_lpf); % DC-normalise

end

function y = sinc_local(x)
% sinc(x) = sin(pi*x)/(pi*x), with sinc(0) = 1.
    y = ones(size(x));
    mask_idx = (x ~= 0);
    y(mask_idx) = sin(pi * x(mask_idx)) ./ (pi * x(mask_idx));
end
