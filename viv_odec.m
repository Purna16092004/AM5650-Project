function dydt = viv_odec(~, y, delta, d_coeff, epsilon, A, M)
    y1 = y(1); dy1 = y(2);
    q = y(3); dq = y(4);

    d2y = -d_coeff * dy1 - delta^2 * y1 + M * q;
    f = A * d2y;
    d2q = -epsilon * (q^2 - 1) * dq - q + f;

    dydt = [dy1; d2y; dq; d2q];
end