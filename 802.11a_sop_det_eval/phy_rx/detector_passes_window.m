function det_pair = detector_passes_window(det_id, rx_iq, thresh_mat, ...
    p_xc_single, p_xc_bank, ...
    L_lag_samples, S_ac, S_xc, eps_load, ...
    n_win_start, n_win_end, sigma_sq_mean, sigma_sq_median, mode)
%DETECTOR_PASSES_WINDOW  Dual-gate detection event with selectable semantics.
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
%   Two semantics are supported via the mode argument:
%
%   mode = 'pd'   (default; receiver-realistic Pd):
%     Scan the detector output from sample 1 onwards. A detection is
%     declared only when the FIRST sample at which the dual gate fires
%     lands inside the L-STF window [n_win_start, n_win_end]. If a
%     firing occurs anywhere in [1, n_win_start - 1] (pre-packet false
%     trigger), the trial is scored as a MISS, even if the detector
%     would later fire again inside the L-STF window. This models a
%     real receiver that locks to the first trigger and commits the
%     downstream PHY chain to a wrong timing reference if that trigger
%     is pre-packet.
%
%   mode = 'pfa'  (false-alarm-only frames):
%     Return 1 iff the dual gate fires at any sample inside the listen
%     window [n_win_start, n_win_end] when no signal is present. The
%     pre-window samples [1, n_win_start - 1] are CFAR-warm-up
%     territory and are not tested by the receiver during the listen
%     interval, so any activity there is ignored. This is the standard
%     window-event Pfa metric.
%
%   Per-detector test (both gates must pass at the same sample):
%     AC   :  T_abs(n) > K_abs * sqrt(L) * sigma_sq   AND  T_rel(n) > alpha
%     XC*  :  T_abs(n) > K_abs * sqrt(sigma_sq)       AND  T_rel(n) > beta_sq
%
%   Inputs:
%     det_id            1 = AC_BASELINE, 2 = XC_SINGLE_MF, 3 = XC_BANK
%     rx_iq             [1 x N] complex baseband at fs_rx
%     thresh_mat        3 x >=2, row = detector, col 1 = rel_thresh,
%                       col 2 = K_abs
%     p_xc_single       sim_params variant (N_h = 1) for XC_SINGLE_MF
%     p_xc_bank         sim_params variant (N_h = 3) for XC_BANK
%     L_lag_samples     AC lag length L [samples @ fs_rx]
%     S_ac, S_xc        post-correlator MA window lengths
%     eps_load          rel-gate denominator regulariser
%     n_win_start       first sample of the inclusion window
%                       (L-STF start for 'pd', listen-interval start for 'pfa')
%     n_win_end         last sample of the inclusion window
%     sigma_sq_mean     scalar mean   CFAR floor for this buffer
%     sigma_sq_median   scalar median CFAR floor for this buffer
%     mode              'pd' (default) or 'pfa'
%
%   Output:
%     det_pair          [det_mean, det_median] as 0/1 doubles
%
%   See also: estimate_floor.m, ac_dual_stat.m, xc_dual_stat.m,
%             eval/eval_80211_sop_kernel.m (per-trial driver).

    if nargin < 14
        mode = 'pd';
    end

    n_buf = numel(rx_iq);
    n_search_end = min(n_buf, n_win_end);
    if n_search_end < 1
        det_pair = [0, 0];
        return
    end
    n_win_lo = max(1, n_win_start);
    rel_thresh = thresh_mat(det_id, 1);
    K_abs = thresh_mat(det_id, 2);

    if det_id == 1 % AC_BASELINE
        [T_abs, T_rel] = ac_dual_stat(rx_iq, L_lag_samples, S_ac, eps_load);
        abs_eff_mean   = K_abs * sqrt(L_lag_samples) * sigma_sq_mean;
        abs_eff_median = K_abs * sqrt(L_lag_samples) * sigma_sq_median;
    elseif det_id == 2 % XC_SINGLE_MF
        [T_abs, T_rel] = xc_dual_stat(rx_iq, p_xc_single, S_xc);
        abs_eff_mean   = K_abs * sqrt(sigma_sq_mean);
        abs_eff_median = K_abs * sqrt(sigma_sq_median);
    else % XC_BANK
        [T_abs, T_rel] = xc_dual_stat(rx_iq, p_xc_bank, S_xc);
        abs_eff_mean   = K_abs * sqrt(sigma_sq_mean);
        abs_eff_median = K_abs * sqrt(sigma_sq_median);
    end

    switch mode
    case 'pd'
        % Search [1, n_search_end] for the first sample where both
        % gates fire; declare detection iff that sample lies inside
        % [n_win_lo, n_search_end].
        idx = 1 : n_search_end;
        rel_pass = T_rel(idx) > rel_thresh;
        fires_mean   = (T_abs(idx) > abs_eff_mean)   & rel_pass;
        fires_median = (T_abs(idx) > abs_eff_median) & rel_pass;

        first_mean   = find(fires_mean,   1, 'first');
        first_median = find(fires_median, 1, 'first');

        det_mean   = ~isempty(first_mean)   && first_mean   >= n_win_lo;
        det_median = ~isempty(first_median) && first_median >= n_win_lo;

    case 'pfa'
        % Standard window-event Pfa: any firing in [n_win_lo, n_win_end]
        % counts. The pre-window samples are CFAR-warm-up.
        idx = n_win_lo : n_search_end;
        rel_pass = T_rel(idx) > rel_thresh;
        det_mean   = any((T_abs(idx) > abs_eff_mean)   & rel_pass);
        det_median = any((T_abs(idx) > abs_eff_median) & rel_pass);

    otherwise
        error('detector_passes_window:bad_mode', ...
            'mode must be ''pd'' or ''pfa'' (got ''%s'')', mode);
    end

    det_pair = double([det_mean, det_median]);
end
