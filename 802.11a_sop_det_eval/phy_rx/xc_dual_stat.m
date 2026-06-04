function [T_abs_xc, T_rel_xc] = xc_dual_stat(r_iq, p_bank, S)
%XC_DUAL_STAT  L-STF cross-correlator dual statistics (single MF or bank).
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
% Delegates the per-hypothesis correlation to rx_lstf_xcorr_bank,
% then applies an S-tap moving average on the metric and on |xc|^2.
%
% Intentional input delay (N_INPUT_DELAY = 16 samples):
%   The half-band decimator cascade has a SYMMETRIC impulse response of
%   ~14 fs_rx samples centered on its (GD-removed) peak, so the L-STF
%   rising edge appears in rx_iq up to ~7 samples BEFORE the rigorous
%   SoP. At high SNR the matched filter + comb + smoother respond to
%   that lead-in energy and fire BEFORE the SoP, breaking first-firing
%   monotonicity. A 16-sample delay applied at the detector input pushes
%   the entire XC response window past the cascade smear, so the first
%   firing is guaranteed to land at or after SoP + 16. Downstream timing
%   refinement (e.g. L-LTF correlation) absorbs the 16-sample offset.

    [metric_pwr, ~, ~, xc_mag_max] = rx_lstf_xcorr_bank(r_iq, p_bank);
    if S > 1
        sm_kernel = ones(1, S) / S;
        T_rel_xc = filter(sm_kernel, 1, metric_pwr);
        T_abs_xc = sqrt(filter(sm_kernel, 1, xc_mag_max.^2));
    else
        T_rel_xc = metric_pwr;
        T_abs_xc = xc_mag_max;
    end

    % Apply the intentional input delay: output at sample n now reflects
    % what the un-delayed pipeline would have produced at sample n-N_DELAY.
    % Equivalent to inserting a 16-sample shift register at rx_iq -> MF.
    N_INPUT_DELAY = 16;
    n_buf = numel(T_abs_xc);
    if n_buf > N_INPUT_DELAY
        T_abs_xc = [zeros(1, N_INPUT_DELAY), T_abs_xc(1 : n_buf - N_INPUT_DELAY)];
        T_rel_xc = [zeros(1, N_INPUT_DELAY), T_rel_xc(1 : n_buf - N_INPUT_DELAY)];
    else
        T_abs_xc = zeros(1, n_buf);
        T_rel_xc = zeros(1, n_buf);
    end
end
