function xnew = trapezoid(x, mu, F, g0, Isp, dt)

% x(1) = Vr
% x(2) = Vnu
% x(3) = r
% x(4) = nu
% x(5) = m

% Initial function value
fk = xdot(x(1), x(2), x(3), mu, F, g0, Isp, x(5));

% Prediction for x(k+1)
xguess = x + dt*fk;

% New function value
fk1 = xdot(xguess(1), xguess(2), xguess(3), mu, F, g0, Isp, xguess(5));

% New orbital parameters
xnew = x + (dt/2)*(fk + fk1);

end