function [EE,EE_x] = E_interpolate(x,NS,basis,ss)

[P] = required_nodes(x,NS,ss);
[EE] = SF_1D(x,NS,basis,P,ss);
[EE_x] = DSF_x(x,NS,basis,P,ss);
    


 
