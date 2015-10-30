function [TOP, BOT] = getTopBot(sums, squares, dt, nyq)
%[TOP, BOT] = getTopBot(sums, squares, dt, nyq)
%Function which takes the signal, aligned square wave, and time
%differential and generates the top and bottom arrays.
%
%INPUTS:
%sums       = Array of pmt readings.
%squares    = Square wave.
%dt         = time difference between readings. 1 / f
%nyq        = Desired nyquist frequency (usually 1/2 of chop)
%OUTPUTS:
%TOP        = Points corresponding to the downsampled laser ON
%BOT        = Poitns corresponding to the downsampled lase OFF
    su = sums .* (squares + 1) /2;
    sb = sums .* (-squares + 1) / 2;
    [f, gu] = spec(su, dt);
    [f, gb] = spec(sb, dt);
    %change > < if returning to resample and change cut usage
    cut = find(abs(f) > nyq);
    gu(cut) = 0;
    gb(cut) = 0;
    [t, su_c] = ispec(gu, f);
    [t, sb_c] = ispec(gb, f);
    p = str2double(rat(nyq * 2))
    q = str2double(rat(1 / dt))
    [P, Q] = rat(p / q)
    sprintf('p = %f, q = %f', p, q)
    TOP = resample(su_c, P, Q);
    BOT = resample(sb_c, P, Q);
end
