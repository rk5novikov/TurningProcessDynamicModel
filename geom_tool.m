function [zt, rt, ztm, rtm] = geom_tool(R, e, hr, le1, le2, al1, al2, re, ke)
    deg = pi / 180;

    % elevation angle of the first and second sections
    alfa1 = al1*deg;
    alfa2 = al2*deg;

    % tool edge opening angle
    beta = 180 * deg - alfa1 - alfa2;

    % Here we determine auxiliary points for utting edge geometry
    r0 = R - hr + re;
    r2 = r0 - re * sin(alfa1);
    r1 = r2 + le1 * cos(alfa1);
    r3 = r0 - re * sin(alfa2);
    r4 = r3 + le2 * cos(alfa2);
    z4 = e;
    z3 = z4 - le2 * sin(alfa2);
    z0 = z3 - re * cos(alfa2);
    z2 = z0 - re * cos(alfa1);
    z1 = z2 - le1 * sin(alfa1);

    % Discretization of the cutting edge with elements approximately equal in length
    % approximate element length
    a = (le1 + le2 + re * beta) / ke;

    % number of elements on the first line
    ke1 = ceil(le1 / a);

    % number of elements on the second line
    ke2 = ceil(le2 / a);

    % number of elements on the circle
    kere = ceil(re * beta / a);

    % element length on the first line
    a1 = le1 / ke1;

    % element length on the second line
    a2 = le2 / ke2;

    % the angle of the arc corresponding to one element on the circle
    gamma = beta / kere;

    % Creating arrays storing the coordinates of the points of the cutting edge
    zt = zeros(1, ke1 + ke2 + kere + 1);
    rt = zeros(1, ke1 + ke2 + kere + 1);
    for i = 1:(ke1 + kere + ke2 + 1)
        if ((i >= 1) && (i < ke1 + 1))
            zt(i) = z1 + a1 * (i - 1) * sin(alfa1);
            rt(i) = r1 - a1 * (i - 1) * cos(alfa1);
        else if ((i >= ke1+1) && (i <= ke1 + kere + 1))
                zt(i) = z0 - re * cos(alfa1 + (i - ke1 - 1) * gamma);
                rt(i) = r0 - re * sin(alfa1 + (i - ke1 - 1) * gamma);
            else if ((i > ke1 + kere + 1) && (i <= ke1 + ke2 + kere + 1))
                    zt(i) = z3 + a2 * (i - ke1 - kere - 1) * sin(alfa2);
                    rt(i) = r3 + a2 * (i - ke1 - kere - 1) * cos(alfa2);
                end
            end
        end
    end
    ztm = zeros(1, 2);
    rtm = zeros(1, 2);
    ztm(1) = z0;
    rtm(1) = r0;
    ztm(2) = z0 + ((le1 + le2) / 2) * sin(alfa2 - alfa1);
    rtm(2) = r0 + ((le1 + le2) / 2) * cos(alfa2 - alfa1);
    end
