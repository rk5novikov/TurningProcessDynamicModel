function [zp, rp] = geom_cross2(zt1, rt1, zt2, rt2, zwc1, rwc1, zwc2, rwc2)
    if(rt1~=rt2)
        kn = tan((pi / 2) + atan((rt2 - rt1) / (zt2 - zt1)));
        kw = (rwc2 - rwc1) / (zwc2 - zwc1);
        zp = ( rwc1 - ((rt1 + rt2) / 2) - kw * zwc1 + kn * ((zt1 + zt2) / 2) ) / (kn - kw);
        rp = ((rt1 + rt2) / 2) + kn * (zp - ((zt1 + zt2) / 2));
    else
        kw = (rwc2 - rwc1) / (zwc2 - zwc1);
        zp = (zt1 + zt2) / 2;
        rp = rwc1 + kw * (zp - zwc1);
    end