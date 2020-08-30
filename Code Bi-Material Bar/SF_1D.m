function [si] = SF_1D(x,NS,basis,P,ss)

%this function gives Shape function at any point (x) centered at node (xI)
% For linear basis we need 2x2 M matrix and for quadratic basis we need 3x3
% M matrix in 1D

if (basis==1)
    M=zeros(2,2);
elseif(basis==2)
    M=zeros(3,3);
elseif(basis==3)
    M=zeros(4,4);
elseif(basis==4)
    M=zeros(5,5);
elseif(basis==5)
    M=zeros(6,6);
end

%Evaluation of the Moment matrix, here we take the node positions
%as xxI and yyI for the summation process to be carried out easily

for integer_1=1:length(P)
    xxI = P(integer_1);        

    if basis==1
       HB=[1;x-xxI];
    elseif basis==2
        HB = [1;x-xxI;(x-xxI)^2];     
    elseif basis==3
        HB = [1;x-xxI;(x-xxI)^2;(x-xxI)^3]; 
    elseif basis==4
        HB = [1;x-xxI;(x-xxI)^2;(x-xxI)^3;(x-xxI)^4];
    elseif basis==5
        HB = [1;x-xxI;(x-xxI)^2;(x-xxI)^3;(x-xxI)^4;(x-xxI)^5];
    end 

    zz=(sqrt((x-xxI)^2))/ss;

    %phy = PHI is the weight function

    %USING Quintic B-SPLINE
    [phi] = phi_eval(zz);
    
    M = M + (HB*transpose(HB)*phi);
end

%After we get the Moment matrix we construct the SF
si = zeros(1,length(NS));
for int_1 = 1:length(NS)
    xI = NS(int_1);
    
    if basis==1
        Ho = [1 ;0];
        H = [1 ;x-xI];
    elseif basis==2
        Ho = [1;0;0];
        H = [1 ;x-xI;(x-xI)^2];  
    elseif basis==3
        Ho = [1 ;0;0;0];
        H = [1;x-xI;(x-xI)^2;(x-xI)^3]; 
    elseif basis==4
        Ho = [1 ;0;0;0;0];
        H = [1;x-xI;(x-xI)^2;(x-xI)^3;(x-xI)^4]; 
    elseif basis==5
        Ho = [1 ;0;0;0;0;0];
        H = [1;x-xI;(x-xI)^2;(x-xI)^3;(x-xI)^4;(x-xI)^5]; 
    end

    z=(sqrt((x-xI)^2))/ss;
    [PHI] = phi_eval(z);

    si(int_1) = transpose(Ho)*inv(M)*H*PHI;
end
%============ FUNCTIONS =================================
function [phi] = phi_eval(zz)
if zz>=0 && zz<(1/3)
    phi = (11/20)-(9/2)*zz^2 + (81/4)*zz^4 -(81/4)*zz^5;
else if zz>=1/3 && zz<2/3
        phi = (17/40) + (15/8)*zz - (63/4)*zz^2 + (135/4)*zz^3 - (243/8)*zz^4 + (81/8)*zz^5;
    else if zz>=2/3 && zz<1
            phi = (81/40) - (81/8)*zz + (81/4)*zz^2 - (81/4)*zz^3 + (81/8)*zz^4 - (81/40)*zz^5;
        else if zz>=1
                phi = 0;
            end
        end
    end
end



