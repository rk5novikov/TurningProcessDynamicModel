function rw = cut(zw, rw, zt, rt)
    r_rez = interp1(zt, rt, zw);
    hzw = zw(2) - zw(1);
    qmin = max((floor(zt(1) / hzw)) - 4, 1);
    qmax = min((ceil(zt(length(zt)) / hzw) + 4), length(zw));
    for q=qmin:qmax
        if (r_rez(q) < rw(q))
            rw(q) = r_rez(q);
        end
    end