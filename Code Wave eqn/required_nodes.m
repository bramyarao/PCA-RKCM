function [P] = required_nodes(x,NS,ss)

% Required nodes for constructing the moment matrix is P      
integer2 =1;
P = [];
for interger1 = 1:length(NS)
    x_m = NS(interger1);
    if (abs(x-x_m)<=ss)
        P(integer2) = x_m;
        integer2 = integer2 + 1;
    else
        continue;
    end
end