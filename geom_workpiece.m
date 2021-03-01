function [tw, zw, rw] = geom_workpiece(L, R, ket, kel)
    tw = zeros(ket,kel+1);
    zw = zeros(ket,kel+1);
    rw = zeros(ket,kel+1);
    for p = 1:ket
        for q = 2:kel
            tw(p,q) = (p-1)*(2*pi*R/ket);
            zw(p,q) = (q-1)*(L/kel);
            rw(p,q) = R;
        end
    end
    % ends
    for p = 1:ket
            tw(p,1) = (p-1)*(2*pi*R/ket);
            zw(p,1) = 0;
            rw(p,1) = 0;
            tw(p,kel+1) = (p-1)*(2*pi*R/ket);
            zw(p,kel+1) = kel*(L/kel);
            rw(p,kel+1) = 0;
    end