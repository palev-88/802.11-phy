function [sigma_sq_mean, sigma_sq_median] = estimate_floor(rx_iq, n_samples)
%ESTIMATE_FLOOR  One-shot CFAR floor over the first n_samples of rx_iq.
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
%   The noise+jammer floor is estimated ONCE per buffer by averaging the
%   first n_samples samples of rx_iq (the post-decimation receive stream):
%
%       sigma_sq_mean   = (1/n) * sum_{k=1..n}    |r[k]|^2
%       sigma_sq_median = median(|r[1..n]|^2) / ln(2)
%
%   These first samples are GUARANTEED to be pre-packet noise+jammer in
%   our simulator because the runner draws guard_min_samples_rx >= 320 =
%   1 SIFS samples of inter-frame quiet before every L-STF. A real
%   receiver sees the same situation: standards-compliant 802.11
%   transmitters cannot place a packet closer than SIFS (16 us) after
%   the previous transmission, so the receiver always has >= 320
%   samples of pure noise+jammer at the start of every listen interval.
%
%   Why ONE-SHOT (not sample-by-sample streaming):
%     This is the minimum-complexity practical implementation. A real
%     receiver could either:
%       (a) compute (sigma_sq_mean, sigma_sq_median) ONCE during the
%           inter-frame quiet period, hold the result, and apply it as
%           the abs-gate floor throughout the next L-STF detection
%           attempt (= what this function models); OR
%       (b) run an adaptive sample-by-sample trailing-window CFAR
%           continuously, so the floor tracks any non-stationarity in
%           the jammer (HW cost in report Sec. 9.3.4, ~2.2 MAC/sample).
%     For the steady-state jammers modelled in this study (stationary
%     CW / 1-MHz noise / wideband / chirp throughout each trial),
%     forms (a) and (b) give statistically equivalent floor estimates;
%     we adopt the simpler (a) for compute speed and direct mapping
%     to the published curves. The report discusses (b) as the
%     non-stationary-jammer extension.
%
%   Why TWO estimators (mean AND median):
%     - mean    : MLE under pure complex-Gaussian noise (lowest variance
%                 of any unbiased estimator under that hypothesis).
%     - median  : robust under burst / pulsed jammers and ADC spikes
%                 that violate Gaussianity; the 1/ln(2) factor
%                 bias-corrects the median to match mean(|r|^2) in
%                 expectation when |r|^2 is exponential.
%     The kernel feeds BOTH to the detector in parallel; the trailing
%     dim-2 axis of every results array is {mean, median} so the
%     report can compare them apples-to-apples.
%
%   Inputs:
%     rx_iq      [1 x N] complex baseband at fs_rx
%     n_samples  scalar count of leading samples to average
%                (typical: 320 = 1 SIFS @ fs_rx = 20 MHz)
%
%   Outputs:
%     sigma_sq_mean    scalar mean   |r|^2 over rx_iq(1..n_samples)
%     sigma_sq_median  scalar median |r|^2 over rx_iq(1..n_samples), /ln(2)
%
%   See also: detector_passes_window.m  (consumer of these scalars),
%             eval/eval_80211_sop_kernel.m (per-trial driver).

    n_buf = numel(rx_iq);
    n_hi = min(n_buf, n_samples);
    if n_hi < 1
        sigma_sq_mean   = 0;
        sigma_sq_median = 0;
        return
    end
    x_sq = abs(rx_iq(1:n_hi)).^2;
    sigma_sq_mean   = mean(x_sq);
    sigma_sq_median = median(x_sq) / log(2);
end
