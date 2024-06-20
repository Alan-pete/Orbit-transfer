%% Part A1: Continuous Small Thrust (Continuous Burn - Trapezoidal)
close all; clear all; clc;

% Given values
a_leo = 8530;                   % semi-major axis of LEO orbit (km)
a_meo = 13200;                  % semi-major axis of MEO orbit (km)
F_e = 10;                       % electric thrust (N)
Isp_e = 2000;                   % electric specific impulse (sec)
F_k = 2000;                     % kick-motor thrust (N)
Isp_k = 270;                    % kick-motor specific impulse (sec)
mu = 398600;                    % gravitational parameter (km^3/sec^2)
g0 = 9.8067;                    % acceleration of gravity (km/s^2)
m = 1069.8;                     % Initial mass (kg)
v_final = sqrt(mu/a_meo);       % orbital velocity of final orbit (km/s)

% Initial position and velocity in orbit
r0 = a_leo;
nu_0 = 0;
w0 = 1/r0 * sqrt(mu / r0);
Vr_0 = 0;                               % radial velocity
Vnu_0 = sqrt(mu/r0);                    % orbital velocity
x0 = [Vr_0; Vnu_0; r0; nu_0; m];

dt = 1;       % sec
t = 1;          % sec

% Conditions for while loop
r_apogee = a_leo;
n = 2;
x_trap(:,1) = x0;

% Continuous burn period
while r_apogee < a_meo

    x_trap(:,n) = trapezoid(x_trap(:,n-1),mu,F_e,g0,Isp_e,dt);

    a = mu / ((2*mu/x_trap(3,n)) - (x_trap(1,n)^2 + x_trap(2,n)^2));
    e = x_trap(3,n)/mu * sqrt((x_trap(2,n)^2 - mu/x_trap(3,n))^2 + (x_trap(1,n)*x_trap(2,n))^2);

    r_apogee = a * (1 + e);
    n = n + 1;
    t = t + dt;
end

% Plotting values
n_burn1 = n;

% Coasting Period
while x_trap(3,n-1) < a_meo

    % New orbital parameters
    x_trap(:,n) = trapezoid(x_trap(:,n-1),mu,0,g0,Isp_e,dt);

    n = n + 1;
    t = t + dt;
end
n_final1 = n-1;

% Impulse at apogee of new orbit
init = x_trap(2,n-1)-Vnu_0;
delta_v = v_final - x_trap(2,n-1);
m_impulse = exp(delta_v*1000 / (g0 * Isp_k)) - 1;

% Final mass of spacecraft (should be 1000)
m_final = x_trap(5,n-1) - m_impulse;

% Initial and Final Orbits
theta = linspace(0,2*pi,360);
rho1 = zeros(length(theta),1);
rho2 = zeros(length(theta),1);
for i = 1:length(theta)
    rho1(i) = a_leo;
    rho2(i) = a_meo;
end

figure
polarplot(theta,rho1,'LineWidth',1,'color','g','LineStyle','--'), hold on
polarplot(x_trap(4,1:n_burn1), x_trap(3,1:n_burn1),'LineWidth',1,'color','b')
polarplot(x_trap(4,n_burn1:n_final1),x_trap(3,n_burn1:n_final1),'LineWidth',1,'color','r')
polarplot(theta,rho2,'LineWidth',1,'color','g','LineStyle','--')
% title('Continuous Small Thrust - Trapezoidal')
legend('Initial Orbit','Thrusting Period','Coasting Period','Final Orbit')

%% Part A2: Continuous Small Thrust (Continuous Burn - RK4)

n = 2;
t1 = 1;
x_rk(:,1) = x0;
r_apogee = a_leo;
% Continuous burn period
while r_apogee < a_meo

    x_rk(:,n) = rk4(x_rk(:,n-1),mu,F_e,g0,Isp_e,dt);

    a = mu / ((2*mu/x_rk(3,n)) - (x_rk(1,n)^2 + x_rk(2,n)^2));
    e = x_rk(3,n)/mu * sqrt((x_rk(2,n)^2 - mu/x_rk(3,n))^2 + (x_rk(1,n)*x_rk(2,n))^2);
    
    r_apogee = a * (1 + e);
    n = n + 1;
    t1 = t1 + dt;
end
n_burn2 = n;

% Coasting Period
while x_rk(3,n-1) < a_meo

    % New orbital parameters
    x_rk(:,n) = rk4(x_rk(:,n-1),mu,0,g0,Isp_e,dt);

    n = n + 1;
    t1 = t1 + dt;
end
n_final2 = n-1;

% Impulse at apogee of new orbit
init1 = x_rk(2,n-1)-Vnu_0;
delta_v1 = v_final - x_rk(2,n-1);
m_impulse1 = exp(delta_v1*1000 / (g0 * Isp_k)) - 1;

% Final mass of spacecraft (should be 1000)
m_final1 = x_rk(5,n-1) - m_impulse1;

figure
polarplot(theta,rho1,'LineWidth',1,'color','g','LineStyle','--'), hold on
polarplot(x_rk(4,1:n_burn2), x_rk(3,1:n_burn2),'LineWidth',1,'color','b')
polarplot(x_rk(4,n_burn2:n_final2),x_rk(3,n_burn2:n_final2),'LineWidth',1,'color','r')
polarplot(theta,rho2,'LineWidth',1,'color','g','LineStyle','--')
% title('Continuous Small Thrust - RK4')
legend('Initial Orbit','Thrusting Period','Coasting Period','Final Orbit')

%% Part A3 & B2: Hohmann Transfer (Impulsive Burn)

% Semi-major axis of transfer orbit
a_transfer = (a_leo + a_meo) / 2;

% Orbital velocities of initial, transfer, and final orbits
v_init = sqrt(mu/a_leo);
v_fin = sqrt(mu/a_meo);

v_perigee = sqrt(2*mu/a_leo - mu/a_transfer);
v_apogee = sqrt(2*mu/a_meo - mu/a_transfer);

% Delta v required for each impulsive maneuver
dv_perigee = v_perigee - v_init;
dv_apogee = v_fin - v_apogee;

dv_total = dv_perigee + dv_apogee;

% Mass required for Hohmann transfer
m_hohmann1 = abs(m * (1 - exp(dv_perigee*1000 / (g0 * Isp_k))));
m_hohmann2 = abs((m-m_hohmann1) * (1 - exp(dv_apogee*1000 / (g0 * Isp_k))));
m_hohmann = m_hohmann1 + m_hohmann2;

%% Part B1: Continuous Large Thrust (Continuous Burn - Trapezoidal)

n = 2;
t2 = 1;
dt = 0.1;
m_large = 1303.2;  % Initial large thrust mass
x0 = [Vr_0; Vnu_0; r0; nu_0; m_large];
x_trap2(:,1) = x0;
r_apogee = a_leo;
% Continuous burn period
while r_apogee < a_meo

    x_trap2(:,n) = trapezoid(x_trap2(:,n-1),mu,F_k,g0,Isp_k,dt);

    a = mu / ((2*mu/x_trap2(3,n)) - (x_trap2(1,n)^2 + x_trap2(2,n)^2));
    e = x_trap2(3,n)/mu * sqrt((x_trap2(2,n)^2 - mu/x_trap2(3,n))^2 + (x_trap2(1,n)*x_trap2(2,n))^2);
    
    r_apogee = a * (1 + e);
    n = n + 1;
    t2 = t2 + dt;
end
n_burn3 = n;

% Coasting Period
while x_trap2(3,n-1) < a_meo

    % New orbital parameters
    x_trap2(:,n) = rk4(x_trap2(:,n-1),mu,0,g0,Isp_k,dt);

    n = n + 1;
    t2 = t2 + dt;
end
n_final3 = n-1;

% Impulse at apogee of new orbit
init2 = x_trap2(2,n-1)-Vnu_0;
delta_v2 = v_final - x_trap2(2,n-1);
m_impulse2 = exp(delta_v2*1000 / (g0 * Isp_k)) - 1;

% Final mass of spacecraft (should be 1000)
m_final2 = x_trap2(5,n-1) - m_impulse2;

figure
polarplot(theta,rho1,'LineWidth',1,'color','g','LineStyle','--'), hold on
polarplot(x_trap2(4,1:n_burn3), x_trap2(3,1:n_burn3),'LineWidth',1,'color','b')
polarplot(x_trap2(4,n_burn3:n_final3),x_trap2(3,n_burn3:n_final3),'LineWidth',1,'color','r')
polarplot(theta,rho2,'LineWidth',1,'color','g','LineStyle','--')
% title('Continuous Large Thrust - Trapezoidal')
legend('Initial Orbit','Thrusting Period','Coasting Period','Final Orbit')

%% Part B2: Continuous Large Thrust (Continuous Burn - RK4)

n = 2;
t3 = 1;
dt = 1;
m_large = 1303.2;  % Initial large thrust mass
x0 = [Vr_0; Vnu_0; r0; nu_0; m_large];
x_rk2(:,1) = x0;
r_apogee = a_leo;

% Continuous burn period
while r_apogee < a_meo

    x_rk2(:,n) = rk4(x_rk2(:,n-1),mu,F_k,g0,Isp_k,dt);

    a = mu / ((2*mu/x_rk2(3,n)) - (x_rk2(1,n)^2 + x_rk2(2,n)^2));
    e = x_rk2(3,n)/mu * sqrt((x_rk2(2,n)^2 - mu/x_rk2(3,n))^2 + (x_rk2(1,n)*x_rk2(2,n))^2);
    
    r_apogee = a * (1 + e);
    n = n + 1;
    t3 = t3 + dt;
end
n_burn4 = n;

% Coasting Period
while x_rk2(3,n-1) < a_meo

    % New orbital parameters
    x_rk2(:,n) = rk4(x_rk2(:,n-1),mu,0,g0,Isp_k,dt);

    n = n + 1;
    t3 = t3 + dt;
end
n_final4 = n-1;

% Impulse at apogee of new orbit
init3 = x_rk2(2,n-1)-Vnu_0;
delta_v3 = v_final - x_rk2(2,n-1);
m_impulse3 = exp(delta_v3*1000 / (g0 * Isp_k)) - 1;

% Final mass of spacecraft (should be 1000)
m_final3 = x_rk2(5,n-1) - m_impulse3;

figure
polarplot(theta,rho1,'LineWidth',1,'color','g','LineStyle','--'), hold on
polarplot(x_rk2(4,1:n_burn4), x_rk2(3,1:n_burn4),'LineWidth',1,'color','b')
polarplot(x_rk2(4,n_burn4:n_final4),x_rk2(3,n_burn4:n_final4),'LineWidth',1,'color','r')
polarplot(theta,rho2,'LineWidth',1,'color','g','LineStyle','--')
% title('Continuous Large Thrust - RK4')
legend('Initial Orbit','Thrusting Period','Coasting Period','Final Orbit')

%% Part C: Continuous Large Thrust with Non-Impulsive Burns

n = 2;
t4 = 1;
m_large = 1699.75;  % Initial large thrust mass
x0 = [Vr_0; Vnu_0; r0; nu_0; m_large];
x_rk3(:,1) = x0;
r_apogee = a_leo;

% Continuous burn period
while r_apogee < a_meo

    x_rk3(:,n) = rk4(x_rk3(:,n-1),mu,F_k,g0,Isp_k,dt);

    a = mu / ((2*mu/x_rk3(3,n)) - (x_rk3(1,n)^2 + x_rk3(2,n)^2));
    e = x_rk3(3,n)/mu * sqrt((x_rk3(2,n)^2 - mu/x_rk3(3,n))^2 + (x_rk3(1,n)*x_rk3(2,n))^2);
    
    r_apogee = a * (1 + e);
    n = n + 1;
    t4 = t4 + dt;
end
n_burn5 = n;

t_split = t4/2;

% Coasting Period
while t4 < 5400

    % New orbital parameters
    x_rk3(:,n) = rk4(x_rk3(:,n-1),mu,0,g0,Isp_k,dt);

    n = n + 1;
    t4 = t4 + dt;
end
vmid5 = x_rk3(2,n-1);
n_mid5 = n;

% Circularizing Thrust
while x_rk2(3,n-1) < a_meo

    % New orbital parameters
    x_rk3(:,n) = rk4(x_rk3(:,n-1),mu,F_k,g0,Isp_k,dt);

    n = n + 1;
    t4 = t4 + dt;
end

n_final5 = n-1;

% Impulse at apogee of new orbit
delta_v4 = v_final - vmid5;
m_impulse4 = exp(delta_v4*1000 / (g0 * Isp_k)) - 1;

% Final mass of spacecraft (should be 1000)
m_final4 = x_rk3(5,n-1) - m_impulse4;

figure
polarplot(theta,rho1,'LineWidth',1,'color','g','LineStyle','--'), hold on
polarplot(x_rk3(4,1:n_burn5), x_rk3(3,1:n_burn5),'LineWidth',1,'color','b')
polarplot(x_rk3(4,n_burn5:n_mid5),x_rk3(3,n_burn5:n_mid5),'LineWidth',1,'color','r')
polarplot(x_rk3(4,n_mid5:n_final5), x_rk3(3,n_mid5:n_final5),'LineWidth',1,'color','b')
polarplot(theta,rho2,'LineWidth',1,'color','g','LineStyle','--')
% title('Continuous Large Thrust - Non-Impulsive')
legend('Initial Orbit','First Thrust','Coasting Period','Second Thrust','Final Orbit')










