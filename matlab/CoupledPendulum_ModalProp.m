function [M, K, Lambda, wn, Phi, U, gamma] = ...
                                 CoupledPendulum_ModalProp(cp_param, F)
% This function determines the mass and stiffness matrices, 
    g  = cp_param(1);
    m  = cp_param(2);
    L  = cp_param(3);
    k  = cp_param(4);

    M = [m*L^2 0; 0 m*L^2];
    K = zeros(2, 2);
    K(1, 1) = m*g*L+k*L^2;
    K(2, 2) = K(1, 1);
    K(1, 2) = -k*L^2;
    K(2, 1) = K(1, 2);

    [Phi, Lambda] = eig(K, M);
    Phi = -Phi;
    wn = zeros(1, 2);
    U = zeros(2);
    for i = 1:2
        wn(i) = sqrt(Lambda(i, i));
        U(:,i) = Phi(:,i)/Phi(i, 1);
    end
    % Phi = Modal matrix: Phi = [phi1 phi2]
    % lambda = Eigenvalues

    gamma = Phi'*F; % Modal torque amplitude
end