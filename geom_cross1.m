function [zp, rp] = geom_cross1(zt1, rt1, zt2, rt2, zwc1, rwc1, zwc2, rwc2)
    kt = (rt2 - rt1) / (zt2 - zt1);
    kw = (rwc2 - rwc1) / (zwc2 - zwc1);
    zp = (rwc1 - rt1 - kw * zwc1 + kt * zt1) / (kt - kw);
    rp = rt1 + kt * (zp - zt1);