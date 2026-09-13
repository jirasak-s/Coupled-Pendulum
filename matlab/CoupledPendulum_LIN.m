function dxdt = CoupledPendulum_LIN(t, x, A, B, N)
    dxdt = A*x + B*N(t);
end