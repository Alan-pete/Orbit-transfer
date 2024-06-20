function xnew = rk4(x, mu, F, g0, Isp, dt)

% x(1) = Vr
% x(2) = Vnu
% x(3) = r
% x(4) = nu
% x(5) = m

% K1
[x_dot] = xdot(x(1), x(2), x(3), mu, F, g0, Isp, x(5));
k_1 = x_dot;
x_new = x + dt/2*k_1;

% K2
[x_dot] = xdot(x_new(1), x_new(2), x_new(3), mu, F, g0, Isp, x_new(5));
k_2 = x_dot;
x_new = x_new + dt/2*k_2;

% K3
[x_dot] = xdot(x_new(1), x_new(2), x_new(3), mu, F, g0, Isp, x_new(5));
k_3 = x_dot;
x_new = x_new + dt*k_3;

% K4
[x_dot] = xdot(x_new(1), x_new(2), x_new(3), mu, F, g0, Isp, x_new(5));
k_4 = x_dot;

% New value
xnew = x + (dt/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);

end