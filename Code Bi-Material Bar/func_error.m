function [error1,ERROR_TS] = func_error(u_ERR_full,u_ERR_red)

% Finding the sum of the RMS error in all time steps

no_N = size(u_ERR_full,1); % Number of DOF
no_T = size(u_ERR_full,2); % Number of time steps
ERROR_TS = zeros(no_T,1); % Saving error in each time step

for int1 = 1:no_T

    u_full = u_ERR_full(:,int1);
    u_red = u_ERR_red(:,int1);
    
    u_diff = u_full-u_red;
    
    sume = 0;
    for int2 = 1:no_N
        sume = sume + u_diff(int2)^2;
    end 
    
    err = (sqrt(sume))/(sqrt(no_N));
    ERROR_TS(int1) = err;
end

error1 = sum(ERROR_TS);


