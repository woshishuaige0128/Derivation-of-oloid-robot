%% OLOID-ROBOT ANALYSIS
% =========================================================================
% Description:
% This code implements the complete dynamic analysis of oloid robots with
% rolling locomotion, including geometric modeling, kinematic analysis, and
% derivation of equations of motion using Lagrangian formulation.
%
% Reference: Supporting Information Section - Oloid robots with rolling locomotion
% =========================================================================

%% Initialization
clc;clear;close all;
syms R th(t) phi(t) X Y Z alpha m g dis omega mp

%% Part 1: Basic Description - Fundamental geometric description (Appendix Section 1)

% Define contact point coordinates - corresponding to Appendix Equation [2]
% PA = [R sin θ, -R/2 - R cos θ, 0]  
% PB = [0, R/2 + R cos φ, R sin φ]
OA = [R*sin(th(t)) -R*cos(th(t))-R/2 0];  % Position of contact point PA in moving coordinate system
OB = [0 R*cos(phi(t))+R/2 R*sin(phi(t))]; % Position of contact point PB in moving coordinate system
AB = OB-OA;  % Vector connecting the two contact points PAPB

% Calculate tangent vectors - for geometric constraint derivation
OA_t = diff(OA,th(t));  % Partial derivative of OA with respect to θ, yielding tangent vector
OB_t = diff(OB,phi(t)); % Partial derivative of OB with respect to φ, yielding tangent vector

% Mixed product calculation - corresponding to determinant in Appendix Equation [3]
% When both tangent vectors lie within the contact plane, the mixed product equals zero
mix_product = det([AB;OA_t;OB_t]);
th_phi = simplify(mix_product);
collect_phi = collect(th_phi,cos(phi(t)));

% Geometric constraint relationship - corresponding to Appendix Equation [4]
% cos φ + cos θ + cos φ cos θ = 0
% Solving for φ in terms of θ from geometric constraints
c_phi = R^3*cos(th(t))/(-R^3*(cos(th(t)) + 1));
s_phi = -sqrt(1-c_phi^2);

% Substitute φ expression into AB vector to eliminate φ parameter
AB_nophi = subs(AB,[cos(phi(t)),sin(phi(t))],[c_phi,s_phi]);
AB_modulus_squared = simplify(AB_nophi(1)^2 + AB_nophi(2)^2 + AB_nophi(3)^2);
AB_modulus = simplify(sqrt(AB_modulus_squared));
pretty(AB_nophi(2))

%% Coordinate transformation matrix construction - corresponding to Appendix Section 3 coordinate system transformation

% xyz transformation matrix - fundamental transformation for establishing Frenet-Serret framework
xyz_trans = [cos(th(t)) -sin(th(t))*sqrt(1+2*cos(th(t)))/sqrt(2*(1+cos(th(t)))) -sin(th(t))/sqrt(2*(1+cos(th(t))));
             sin(th(t)) cos(th(t))*sqrt(1+2*cos(th(t)))/sqrt(2*(1+cos(th(t)))) cos(th(t))/sqrt(2*(1+cos(th(t))));
             0 -1/sqrt(2*(1+cos(th(t)))) sqrt(1+2*cos(th(t)))/sqrt(2*(1+cos(th(t))))];

% Additional transformation matrix
trans = [sqrt(3)/9 * (1+2*cos(th(t)))^(3/2)/(cos(th(t)/2)) -(sqrt(3)*sin(th(t)/2))/(9*cos(th(t)/2))*(5+4*cos(th(t))) 0;
         (sqrt(3)*sin(th(t)/2))/(9*cos(th(t)/2))*(5+4*cos(th(t))) sqrt(3)/9 * (1+2*cos(th(t)))^(3/2)/(cos(th(t)/2)) 0;
         0 0 1];

xyz_trans_inv = inv(xyz_trans);
xyz_trans_inv_simplified = simplify(xyz_trans_inv);
xyz = simplify(trans*xyz_trans_inv_simplified);  % Final transformation matrix
xyz_constant = xyz*[-sin(th(t))*R;R/2+R*cos(th(t));0];

%% Centroid trajectory calculation

% Trajectory equations for contact point PA - corresponding to Appendix Equation [9]
xa = (2*R*sqrt(3)/9) * (asin((2/sqrt(3))*sin(th(t)/2)) + asin((sin(th(t)/2)/sqrt(3))/cos(th(t)/2)) + 2*sin(th(t)/2)*sqrt(1+2*cos(th(t))));
ya = (8*R*sqrt(3)/9) * (sin(th(t)/2).^2) - (2*R*sqrt(3)/9) * log(cos(th(t)/2));
za = 0;

% Centroid trajectory in global coordinate system
xo = simplify(xa+xyz_constant(1));
yo = simplify(ya+xyz_constant(2));
zo = simplify(za+xyz_constant(3));

% Rotation matrix - corresponding to transformation matrix P in Appendix equations
Rotation_matrix = xyz;
simplify(Rotation_matrix.'*Rotation_matrix)  % Verify orthogonality

%% Driving point dynamics - corresponding to Appendix Equation [27] and Section 5

% Position of driving point in global coordinate system - corresponding to Appendix Equation [29]
% [Xd Yd Zd]' = [XO YO ZO]' + P[ld cos ψ; 0; ld sin ψ]
xyz_d = Rotation_matrix*[dis*cos(psi(t))-sin(th(t))*R;R/2+R*cos(th(t));dis*sin(psi(t))];
xd = simplify(xa+xyz_d(1));
yd = simplify(ya+xyz_d(2));
zd = simplify(za+xyz_d(3));

% Driving point velocity calculation
vxd = diff(xd,t);
vyd = diff(yd,t);  
vzd = diff(zd,t);

% Kinetic energy of driving point - corresponding to Appendix Equation [27]
% Ekd = (1/2)md(Ẋd² + Ẏd² + Żd²)
kd = simplify(1/2 * mp * (vxd^2 + vyd^2 + vzd^2));
k_d = collect(kd,(diff(th(t), t)));

% Potential energy of driving point - corresponding to Appendix Equation [27]
% Pd = mdgZd
potential_d = mp * g * zd;
p_d = simplify(potential_d);

%% Body angular velocity calculation - corresponding to Appendix Equation [22]

% Calculate time derivative of rotation matrix
d_Rotation_matrix = diff(Rotation_matrix,t);
tran_Rotation_matrix = Rotation_matrix.';
angular_matrix = tran_Rotation_matrix*d_Rotation_matrix;
angular = simplify(angular_matrix);
angular = simplify(subs(angular,cos(th(t)/2),sqrt((cos(th(t))+1)/2)));

% Extract angular velocity components - corresponding to Appendix Equation [22]
angular_x = simplify(angular(3,2));  % ωx component
angular_y = simplify(angular(1,3));  % ωy component  
angular_z = simplify(angular(2,1));  % ωz component

omega_x = collect(angular_x,(diff(th(t), t)));
omega_y = collect(angular_y,(diff(th(t), t)));
omega_z = collect(angular_z,(diff(th(t), t)));

%% Body center of mass kinematics - corresponding to Appendix Equation [22]

% Translational velocity of center of mass - corresponding to Appendix Equation [22]
trans_vx = simplify(diff(xo,t));
trans_vy = simplify(diff(yo,t));
trans_vz = simplify(diff(zo,t));

vx = collect(trans_vx,(diff(th(t), t)));
vy = collect(trans_vy,(diff(th(t), t))); 
vz = collect(trans_vz,(diff(th(t), t)));

%% Body kinetic energy calculation - corresponding to Appendix Equation [23]

% Translational kinetic energy - corresponding to Appendix Equation [23]
% Eb,kt = (1/2)mb(vX² + vY² + vZ²)
kinetic_trans = 1/2 * m * (vx^2 + vy^2 + vz^2);

% Rotational kinetic energy - corresponding to Appendix Equation [23] 
% Eb,kr = (1/2)Ib(ωx² + ωy² + ωz²)
% Using oloid moment of inertia coefficients 5/8, 1/2, 5/8
kinetic_angular = 1/2 * m * R^2 * (5/8*omega_x^2 + 1/2*omega_y^2 + 5/8*omega_z^2);

k_t = collect(kinetic_trans,(diff(th(t), t)));    % Translational kinetic energy
k_a = collect(kinetic_angular,(diff(th(t), t)));  % Rotational kinetic energy

% Body potential energy - corresponding to Appendix Equation [24]
% Pb = mbgZO
potential_energy = simplify(m*g*zo);

%% Total Lagrangian function - including both body and driving point
L = k_t+k_a-potential_energy;  % L = T - U
L_diff = collect(L,(diff(th(t), t)));

%% Euler-Lagrange equation solution - corresponding to Appendix Equation [25]

% d/dt(∂L/∂q̇) - ∂L/∂q = 0
dl_dqdot = diff(L_diff,(diff(th(t), t)));           % ∂L/∂θ̇
dt_d_dl_dqdot = diff(simplify(dl_dqdot),t);        % d/dt(∂L/∂θ̇)
dl_dq = diff(L_diff,th(t));                         % ∂L/∂θ

dt_coll_dt2 = collect(simplify(dt_d_dl_dqdot),diff(th(t), t, t));  
dt_coll_dt1 = collect(simplify(dt_coll_dt2),diff(th(t), t));        

dl_dq_cal = simplify(dl_dq);

%% Equation of motion construction - corresponding to Appendix Equation [26]

% Final equation of motion: A(q)q̈ + B(q)q̇² + C(q) = 0
motion_equation = dt_coll_dt1-dl_dq_cal;  % Euler-Lagrange equation
motion_equation_simp = simplify(motion_equation);

% Further simplification
motion_equation_simp2 = simplify(subs((motion_equation_simp),cos(th(t)/2),sqrt((cos(th(t))+1)/2)));

motion_d1 = collect(motion_equation_simp2,diff(th(t), t));      % Collect by θ̇
motion_d2 = collect(motion_d1,diff(th(t), t, t));              % Collect by θ̈

%% Extract coefficient matrices - corresponding to A(q), B(q), C(q) in Appendix Equation [26]

% Corresponding to A(q)θ̈ + B(q)θ̇² + C(q) = 0
[coeff_d1,name] = coeffs(motion_d2,[diff(th(t), t, t), diff(th(t), t)]);
A = simplify(coeff_d1(1,1));    % inertial term
B = simplify(coeff_d1(1,2));    % Coriolis and centrifugal force term  
C = simplify(coeff_d1(1,3));    % gravitational term