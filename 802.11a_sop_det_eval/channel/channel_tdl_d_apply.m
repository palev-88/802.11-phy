function rx_iq = channel_tdl_d_apply(tx_iq, fading, delays_samples)
%CHANNEL_TDL_D_APPLY  Apply precomputed TDL fading taps to a TX waveform.
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
% rx_iq[n] = sum_i fading(i, n) * tx_iq(n - delays_samples(i))
%
% Pure linear operation: scaling tx_iq by sqrt(P_in) commutes with this,
% so for CRN sweeps you can apply the channel to a unit-power waveform
% once and rescale by sqrt(P_in) per SNR.

    tx_iq = tx_iq(:).';
    N = numel(tx_iq);
    n_taps = size(fading, 1);
    rx_iq = complex(zeros(1, N));
    for i = 1:n_taps
        d = delays_samples(i);
        if d == 0
            rx_iq = rx_iq + fading(i, :) .* tx_iq;
        elseif d < N
            rx_iq(d+1:end) = rx_iq(d+1:end) + ...
                fading(i, d+1:end) .* tx_iq(1:N-d);
        end
    end
end
