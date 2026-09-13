function dxdt = CoupledPendulum_NL(t, x, cp_param, u)
% State varibles: 
% x(1) = \theta_1, x(2) = \theta_2, 
% x(3) = \dot{\theta}_1, x(4) = \dot{\theta}_2,

    g  = cp_param(1);
    m  = cp_param(2);
    L  = cp_param(3);
    k  = cp_param(4);
    d  = cp_param(5);
    xi = cp_param(6);

    tau = u(t);

    dxdt = zeros(4, 1);
    xbar = sqrt((d - L*(sin(x(1)) - sin(x(2))))^2 + ...
           (L*(cos(x(1)) - cos(x(2))))^2);

    xs = xbar - xi;

    dxbar_dtheta1 = (d*L*cos(x(1)) - sin(x(1) - x(2))*L^2)/xbar;
    dxbar_dtheta2 = (d*L*cos(x(2)) - sin(x(1) - x(2))*L^2)/xbar;

    dxdt(1) = x(3);
    dxdt(2) = x(4);
    dxdt(3) = (-m*g*L*sin(x(1)) + k*xs*dxbar_dtheta1 + tau(1))/(m*L^2);
    dxdt(4) = (-m*g*L*sin(x(2)) - k*xs*dxbar_dtheta2 + tau(2))/(m*L^2);
end