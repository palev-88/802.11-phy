function [T_abs, T_rel] = ac_dual_stat(r_iq, L, S, eps_load)
%AC_DUAL_STAT  Schmidl & Cox autocorrelator dual statistics.
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
% T_abs[n] = smoothed(|M[n]|)                    -- absolute level
% T_rel[n] = smoothed(|M[n]|) / smoothed(P[n] + eps) -- normalized
% where M[n] = sum_{k=0..L-1} conj(r[n-L-k]) * r[n-k]
% P[n] = sum_{k=0..L-1} |r[n-k]|^2
% S is the moving-average window length applied on top.
%
% HW-friendly variant: the absolute test uses the same MA{|M|}
% statistic as the relative test, so only one smoother delay line
% is required (instead of one for |M| and another for |M|^2). Saves
% S * w bits of state and removes the sqrt operator. The threshold
% scale K_abs is re-calibrated against MA{|M|} (typical noise level
% ~ sqrt(pi L / 4) * sigma^2, vs sqrt(L) * sigma^2 for the older
% sqrt(MA{|M|^2}) form).

    r_iq = r_iq(:).';
    N = numel(r_iq);
    if N < 2*L
        T_abs = zeros(1, N);
        T_rel = zeros(1, N);
        return
    end
    m_iq = complex(zeros(1, N)); % must be complex
    m_iq(L+1:end) = conj(r_iq(1:end-L)) .* r_iq(L+1:end);
    M = filter(ones(1, L), 1, m_iq);
    P = filter(ones(1, L), 1, abs(r_iq).^2);
    sm_kernel = ones(1, S) / S;
    if S > 1
        Mmag_smooth = filter(sm_kernel, 1, abs(M));
        P_smooth = filter(sm_kernel, 1, P);
    else
        Mmag_smooth = abs(M);
        P_smooth = P;
    end
    T_abs = Mmag_smooth;
    T_rel = Mmag_smooth ./ (P_smooth + eps_load);
end
