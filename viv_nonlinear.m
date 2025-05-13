
function dydt = viv_nonlinear(~, y, xi, gamma, mu, delta, M, epsilon)
    y1 = y(1); dy1 = y(2); q = y(3); dq = y(4);
    damping = 2 * xi * delta + gamma / mu;
    ddy = -damping * dy1 - delta^2 * y1 + M * q;
    A_0 = 0.330;
    A_1 = 0;
    A_2 = -0.060;
    A_3 = 0.002;
    B_0 = -0.005;
    B_1 = 0.003;
    B_2 = 0;
    B_3 = 0;
    f = A_0*ddy + A_1*abs(y1)*ddy + A_2*abs(y1)^2*ddy + A_3*abs(y1)^3*ddy + B_0*dy1 + B_1*abs(y1)*dy1 + B_2*abs(y1)^2*dy1 + B_3*abs(y1)^3*dy1;
    ddq = -epsilon * (q^2 - 1) * dq - q + f;
    dydt = [dy1; ddy; dq; ddq];
end