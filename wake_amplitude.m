function q0 = wake_amplitude(w,delta,k)
    
    epsilon = 0.3;
    A = 12;
    gamma = 0.8;
    xi = 3.1e-3;
    M = 2e-4;
    mu = 0.05/M;
    if k == 1
        C = (-1)*(2*xi*delta + gamma/mu);
    elseif k ==2
        C = delta^2 - w^2;
    else 
        C = (2*xi*delta + gamma/mu)*w^2;
    end


    q0 = 2*(1+(A*M/epsilon)*(C/((delta^2-w^2)^2+(2*xi*delta+gamma/mu)^2*w^2)))^(0.5);


end