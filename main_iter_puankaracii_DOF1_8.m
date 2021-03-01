function YY = main_iter_puankaracii_DOF1_8(PPP)

    % Dimensionless input parameters of the problem
    % Natural frequencies of the tool excluding damping in fractions of the speed
    p_1 = PPP;
    p_2 = PPP;

    % Damping coefficient divided to the damping of the first natural frequency
    Dzitta1 = 0.05;
    Dzitta2 = 0.05;

    % Dimensionless tool stiffness
    c1 = 1.0;
    c2 = 1.0;

    % Natural frequencies of the cutting tool
    p0_21 = (2 * pi * p_1)^2;
    p0_22 = (2 * pi * p_2)^2;

    % Dimensionless damping coefficients
    n1 = 2 * pi * p_1 * Dzitta1;
    n2 = 2 * pi * p_2 * Dzitta2;

    % Natural frequencies of the cutting tool including damping effects
    p1 = sqrt(p0_21 - n1^2);
    p2 = sqrt(p0_22 - n2^2);

    % Dimensionless cutting forces coefficients.
    Kc=0.05;

    % Dimensionless geometric parameters i.e. real values divided into feed value
    % Workpiece length and radius
    L = 180;
    R = 180;

    % Radial cutting depth
    hr = 8.0;

    % Dimensionless axial feed rate
    f_ob = 1;

    % Axial coordinate of the tool point closest to the workpiece at the initial time
    e = 40;

    % Cutting tool geometry
    % Lengths of the first and second straight sections of the cutting tool
    le1 = 80;
    le2 = 80;

    % elevation angles of the first and secod section
    al1 = 45;
    al2 = 45;

    % Rounding radius between two sections
    re = 4.8;

    % Parameters of discretization of geometrical models
    % Number of elements for discretization of cutting edge
    ke = 150;

    % Parameters of discretization of the workpiece geometry
    ket = 180;
    kel = 900;

    % Number of workpiece revolutions in whole studied time period  
    oboroty = 70;
    if( ( (PPP>=2.6) && (PPP<=2.8) ) || ( (PPP>=3.5) && (PPP<=3.7) ) || ( (PPP>=4.4) && (PPP<=4.7) ) )
        oboroty = floor(0.43 * oboroty);
    end

    % The number of points to be traversed in the circumferential direction
    kmax = floor(ket * oboroty);

    % Dimensionless time step duration
    tstep = oboroty / kmax;

    % Feed value per grid step in the circumferential direction
    f_h = f_ob / ket;

    % Here we generate cutting edge geometry coordinates
    [zt, rt, ztm, rtm] = geom_tool(R, e, hr, le1, le2, al1, al2, re, ke);

    % Generating of the workpiece geometry coordinates 
    [tw, zw, rw] = geom_workpiece(L, R, ket, kel);

    % Initial values of displacements, tool speeds and forces
    y = 0.0;
    vy = 0.0;
    fy0 = 0.0;

    % Maximum permissible relative error while calculation of the cutting forces
    err = 1e-2;

    % Arrays for storing displacements, forces, radial depth of cut, number of iterations
    YY = zeros(1, kmax);

    % Calculating the constants needed
    calc_const(c1, p0_21, n1, p1, c2, p0_22, n2, p2, tstep);

    % In this cycle we solve oscillation equations and modify workpiece geometry on each step
    for pp=1:kmax
        if (mod(pp, ket) ~= 0)
            p = mod(pp, ket);
            pnext = p + 1;
        else
            p = ket;
            pnext = 1;
        end

        % We take as initial conditions at the current step the final values 
        % from the previous step
        y0 = y;
        vy0 = vy;

        % Check value to control the convergence of cutting forces
        conv = 0;

        % Prepare the input data for the next step
        fyk = fy0;

        % Iterative refinement of the cutting forces
        % Counter of the number of iterations when refining forces
        while (conv == 0)
            % Calculate a solution for displacements and speeds
            %[x,vx] = calc_difx(x0,vx0,fx0,fxk,tstep);
            [y, vy] = calc_dify(y0, vy0, fy0, fyk, tstep);

            % Move the tool iteratively using feed value and obtained dynamic displacements
            zt1 = zt + f_h;
            rt1 = rt + y - y0;
            ztm1 = ztm + f_h;
            rtm1 = rtm + y - y0;

            % Calculate the depth of cut at the end of the current step
            [s, h_z, h_r] = geom_thickness(zw(pnext,:), rw(pnext,:), zt1, rt1, ztm1, rtm1);

            % Calculate forces 
            [fxk_new, fyk_new] = calc_force(Kc, s, h_z, h_r);

            % Check the convergence condition
            if(fyk_new~=0)
                deltaFy = (abs(fyk_new - fyk)) / fyk_new;
                % Change the check value if convergence has been achieved
                if ( deltaFy < err )
                    conv = 1;
                end
            else
                conv = 1;
            end
            % Memorize the calculated forces
            fyk = fyk_new;
        end
        
        % Move the tool using final displacements values
        zt = zt1;
        rt = rt1;
        ztm = ztm1;
        rtm = rtm1;

        % Remove material from the workpiece geometrical model
        rw(p,:) = cut(zw(pnext,:), rw(pnext,:), zt, rt);

        % Prepare the input data for the next step
        fy0 = fyk;
        % Store dynamic displacements and forces to arrays
        YY(pp) = y;
    end
    disp('Last calculated frequency:');
    disp(PPP);
end