clear all
close all

% Dimensionless input parameters of the problem
% Natural frequencies of the tool excluding damping in fractions of the speed
p_1 = 1.8;
p_2 = 1.8;

% Damping coefficient divided to the damping of the first natural frequency
Dzitta1 = 0.05;
Dzitta2 = 0.05;

% Dimensionless tool stiffness
c1 = 1;
c2 = 1;

% Natural frequencies of the cutting tool
p0_21 = (2 * pi * p_1)^2;
p0_22 = (2 * pi * p_2)^2;

% Dimensionless damping coefficients
n1 = 2*pi*p_1*Dzitta1;
n2 = 2*pi*p_2*Dzitta2;

% Natural frequencies of the cutting tool including damping effects
p1 = sqrt(p0_21 - n1^2);
p2 = sqrt(p0_22 - n2^2);

% Dimensionless cutting forces coefficients.
Kc=1.0;

% Dimensionless geometric parameters i.e. real values divided into feed value
% Workpiece length and radius
L = 40;
R = 3;

% Radial cutting depth
hr = 0.2;

% Dimensionless axial feed rate
f_ob = 1;

% Axial coordinate of the tool point closest to the workpiece at the initial time
e = 0;

% Cutting tool geometry
% Lengths of the first and second straight sections of the cutting tool
le1 = 10;
le2 = 10;

% elevation angles of the first and secod section
al1 = 45;
al2 = 45;

% Rounding radius between two sections
re = 2;

% Parameters of discretization of geometrical models
% Number of elements for discretization of cutting edge
ke = 500;

% Parameters of discretization of the workpiece geometry
ket = 180;
kel = 400;

% Number of workpiece revolutions in whole studied time period  
oboroty = 30;

% Number of time steps 
kmax = floor(ket * oboroty);

% Dimensionless time step duration
tstep = oboroty / kmax;

% Feed value per grid step in the circumferential direction
f_h = f_ob / ket;

% Here we generate cutting edge geometry coordinates
[zt, rt] = geom_tool(R, e, hr, le1, le2, al1, al2, re, ke);

% Generating of the workpiece geometry coordinates 
[tw, zw, rw] = geom_workpiece(L, R, ket, kel);

% Initial values of displacements, tool speeds and forces
x = 0.0;
vx = 0.0;
y = 0.0;
vy = 0.0;
fx0 = 0.0;
fy0 = 0.0;

% Arrays for storing displacments, forces, radial depth of cut
XX = zeros(1, kmax);
YY = zeros(1, kmax);
FFX = zeros(1, kmax);
FFY = zeros(1, kmax);
hr_time = zeros(1, kmax);

% Array for timesteps storing
t = 0.0:tstep:(kmax*tstep);

% Calculating the constants needed
calc_const(c1, p0_21, n1, p1, c2, p0_22, n2, p2, tstep);

%–ешение уравнени€ колебаний и срезание материала
for pp=1:kmax
    if (mod(pp,ket) ~= 0)
        p = mod(pp,ket);
    else
        p = ket;
    end
   
    % We take as initial conditions at the current step the final values 
    % from the previous step
    x0 = x;
    y0 = y;
    vx0 = vx;
    vy0 = vy;
    
    % Calculate a solution for displacements and speeds
    [x, vx] = calc_difx(x0, vx0, fx0, fx0, tstep);
    [y, vy] = calc_dify(y0, vy0, fy0, fy0, tstep);
    
    % Move the tool iteratively using feed value and obtained dynamic displacements
    zt1 = zt + f_h + x - x0;
    rt1 = rt + y - y0;
    
    % Calculate the cutting depthat the end of current step
    if (p == ket)
        pnext = 1;
    else
        pnext = p + 1;
    end
    [s,h_z,h_r] = geom_thickness(zw(pnext,:),rw(pnext,:),zt1,rt1);
    
    % Calculate forces
    [fxk, fyk] = calc_force(Kc, s, h_z, h_r);
    
    % Calculate a solution for displacements and speeds
    [x, vx] = calc_difx(x0, vx0, fx0, fxk, tstep);
    [y, vy] = calc_dify(y0, vy0, fy0, fyk, tstep);
    
    % Remove material from the workpiece geometrical model
    rw(p,:) = cut(zw(p,:), rw(p,:), zt, rt);
    
    % Move the tool using final displacements values
    zt = zt + f_h + x - x0;
    rt = rt + y - y0;
    
    % Prepare the input data for the next step
    fx0 = fxk;
    fy0 = fyk;
    
    % Store dynamic displacements and forces into arrays
    XX(pp) = x;
    FFX(pp) = fxk;
    YY(pp) = y;
    FFY(pp) = fyk;
    hr_time(pp) = R - min(rt);
    disp('Calculation progress (%): ');
    disp(100 * pp / kmax);
end


t = 0.0:tstep:((kmax-1)*tstep);
figure(1)
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 18, 'DefaultTextFontName', 'Times New Roman'); 
plot(t, XX, 'b', 'LineWidth', 2);
grid on;
xlabel('Number of workpiece revolution');
ylabel('Nondim. x-deflection of cutter');

figure(2)
plot(t, YY, 'b', 'LineWidth', 2);
grid on;
xlabel('Number of workpiece revolution');
ylabel('Nondim. y-deflection of cutter');

figure(3)
plot(t, FFX, 'b', 'LineWidth', 2);
grid on;
xlabel('Number of workpiece revolution');
ylabel('Nondim. cutting force in x-direction');

figure(4)
plot(t, FFY, 'b', 'LineWidth', 2);
grid on;
xlabel('Number of workpiece revolution');
ylabel('Nondim. cutting force in y-direction');

figure(5)
plot(t,hr_time, 'b', 'LineWidth', 2);
grid on;
xlabel('Number of workpiece revolution');
ylabel('Nondim. radial depth of cutting');

% Nice 3d plot of the workpiece surface unfolding after cutting
rw(:, 1) = (R - hr) * ones(ket, 1);
rw(:, kel + 1) = (R - hr) * ones(ket, 1);
figure(6)
surf(zw,tw,rw);
daspect('manual');
title('Involute of workpiece surface after turning');
xlabel('z');
ylabel('t');
zlabel('r');

% Graph of cutting edge and last machined section of the workpiece
figure(7);
hold on;
plot(zt,rt,'r','LineWidth',2);
plot(zw(p,:), rw(p,:), 'b', 'LineWidth', 2);
daspect('manual');
title('Cutting edge and workpiece section');
xlabel('Nondim. z coordinate');
ylabel('Nondim. r coordinate');
grid on;
hold off;

% Preparing data for processed workpiece visualization 
xw = zeros(ket + 1, kel + 1);
yw = zeros(ket + 1, kel + 1);
for pp=1:kmax
    if (mod(pp, ket) ~= 0)
        p = mod(pp, ket);
    else
        p = ket;
    end
    for q=1:(kel + 1)
        xw(p, q) = rw(p, q) * cos(2 * (p - 1) * pi / ket);
        yw(p, q) = rw(p, q) * sin(2 * (p - 1) * pi / ket);
    end
end

% Redefine arrays of surface coordinates to get a closed
%surface
xw(ket + 1,:) = xw(1,:);
yw(ket + 1,:) = yw(1,:);
zw(ket + 1,:) = zw(1,:);

% Workpiece cylindrical surface plot
figure(8)
surf(zw, xw, yw);
title('Workpiece after cutting');
daspect('manual');
xlabel('z');
ylabel('x');
zlabel('y');

% Graph of the selected cross-section of the workpiece
figure(9)
plot(xw(:,30), yw(:,30), 'r', 'LineWidth', 2);
title('Cross section of workpiece after cutting');
daspect('manual');
xlabel('x');
ylabel('y');
grid on;
