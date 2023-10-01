function [hh, mm, ss] = sec2hhmmss(t)
    hh = floor(t / 3600);
    t = t - hh * 3600;
    mm = floor(t / 60);
    ss = round(t - mm * 60);
end
