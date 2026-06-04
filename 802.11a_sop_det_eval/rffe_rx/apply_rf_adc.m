function rx_codes_out = apply_rf_adc(rx_iq, rx_total_gain_db, adc_n_bits, adc_fullscale_dbm)
%APPLY_RF_ADC  Combined RX-chain gain + ADC quantization + clipping.
%
%   Author:  Panos Alevizos <bigpan27@gmail.com>
%   AI assistance: Claude (Anthropic) -- code drafting and verification.
%   Copyright (c) 2026 Panos Alevizos. Licensed under CC BY 4.0.
%     https://creativecommons.org/licenses/by/4.0/
%
%   Applies the analog RX chain (combined gain at max AGC setting) and
%   the ADC (uniform mid-tread quantization + symmetric clipping at
%   fullscale) in a single block. The OUTPUT IS THE RAW ADC SAMPLES on
%   each I/Q axis expressed in LSBs (signed integer code units) --
%   exactly what a real ADC presents to the digital back-end -- not a
%   reconstructed voltage waveform.
%
%   Pipeline:
%     1. Linear gain by 10^(rx_total_gain_db/20) on amplitude.
%     2. Compute ADC fullscale amplitude in volts from the sine-power
%        fullscale spec at a 50-Ohm reference impedance.
%     3. Map the post-gain analogue amplitude to LSB units via
%        scale = code_max / FS_amp_V (1 LSB <-> FS_amp_V / code_max).
%     4. Quantize I and Q independently to signed integer LSB values
%        via mid-tread rounding, then saturate-clip to the signed range.
%
%   The downstream blocks (decimation FIRs, CFAR floor estimator,
%   dual-gate detector) operate directly on the LSB-valued samples;
%   their absolute thresholds use the CFAR-estimated noise scale from
%   the same sample stream, so the relative comparisons are unaffected
%   by the choice of representation.
%
%   AGC is FIXED at maximum gain. There is no signal-power-dependent
%   backoff. This captures the worst-case low-signal sensitivity
%   scenario where strong jammers can drive the ADC into clip while
%   weak signals are amplified to fill the ADC dynamic range.
%
%   Inputs:
%     rx_iq             [1 x N] complex baseband at fs_rx, antenna-
%                       referenced (the value mean(|rx_iq|^2) is the
%                       antenna-input power in W).
%     rx_total_gain_db  total RX chain gain at max AGC [dB]
%     adc_n_bits        ADC resolution (signed two's complement, so
%                       per-axis LSB range is
%                         [-2^(n_bits-1), 2^(n_bits-1)-1])
%     adc_fullscale_dbm ADC fullscale sine-wave power [dBm] @ 50 Ohm
%
%   Outputs:
%     rx_codes_out      [1 x N] complex samples in LSB units; real
%                       part holds the I-axis LSB value and imaginary
%                       part holds the Q-axis LSB value. Stored as
%                       double for downstream FIR arithmetic but every
%                       value is an integer in
%                       [-2^(n_bits-1), 2^(n_bits-1)-1].

    % Reference impedance (50 Ohm) for the dBm-to-amplitude conversion.
    R_ref_ohm = 50;

    % 1. Linear RX gain (amplitude factor).
    G_lin_amp = 10^(rx_total_gain_db / 20);
    rx_at_adc = rx_iq * G_lin_amp;

    % 2. ADC fullscale AMPLITUDE in volts:
    %      P_FS [W]   = 10^((adc_fullscale_dbm - 30)/10)
    %      A_FS_peak  = sqrt(2 * R_ref * P_FS)   (sine-wave equivalence)
    FS_pwr_W   = 10^((adc_fullscale_dbm - 30) / 10);
    FS_amp_V   = sqrt(2 * R_ref_ohm * FS_pwr_W);   % peak amplitude

    % 3. ADC scale: code_max LSBs per FS_amp_V, i.e. 1 LSB = FS_amp_V / code_max.
    code_max = 2^(adc_n_bits - 1) - 1;             % +full-scale LSB value
    code_min = -2^(adc_n_bits - 1);                % -full-scale LSB value
    scale    = code_max / FS_amp_V;

    % 4. Quantize I and Q axes independently (mid-tread rounding), then
    %    clip to the signed LSB range. Output stays in double precision
    %    so downstream FIR/MA arithmetic is unaffected, but the values
    %    themselves are integers in [code_min, code_max] (LSBs).
    code_i = max(code_min, min(code_max, round(real(rx_at_adc) * scale)));
    code_q = max(code_min, min(code_max, round(imag(rx_at_adc) * scale)));

    rx_codes_out = code_i + 1j * code_q;
end
