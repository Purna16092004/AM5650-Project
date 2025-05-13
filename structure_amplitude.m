function y0 = structure_amplitude(w,delta,q0)
        
        gamma = 0.8;
        xi = 3.1e-3;
        M = 2e-4;
        mu = 0.05/M;
        x = (delta^2 - w^2)^2;
        y = (2*xi*delta +gamma/mu)^2*w^2;
        y0 = q0*M*(x+y)^(-0.5);

end