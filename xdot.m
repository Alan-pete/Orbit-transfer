function X = xdot(V_r, V_nu,r, mu, F_thrust, g, Isp, m)

V_rdot = (V_nu^2 / r) - (mu / r^2);
Vnu_dot = -(V_r * V_nu)/r + F_thrust/m/1000;
rdot = V_r;
nu_dot = V_nu / r;
mdot = -F_thrust / (g * Isp);

X = [V_rdot; Vnu_dot; rdot; nu_dot; mdot];