function [SIxx] = DSF_xx(x,NS,basis,P,ss)

% This function gives the second derivative of Shape function at 
% any point (x) centered at nodes (xI)
% For linear basis we need 2x2 M matrix and for quadratic basis we need 3x3
% M matrix in 1D

% Evaluation of derivatives of M: M,x and M,xx
if (basis==1)
    M   =zeros(2,2);
    M_x =zeros(2,2);
    M_xx=zeros(2,2);
    
elseif(basis==2)
    M   =zeros(3,3);
    M_x =zeros(3,3);
    M_xx=zeros(3,3);
    
elseif(basis==3)
    M   =zeros(4,4);
    M_x =zeros(4,4);
    M_xx=zeros(4,4);
    
elseif(basis==4)
    M   =zeros(5,5);
    M_x =zeros(5,5);
    M_xx=zeros(5,5);
    
elseif(basis==5)
    M   =zeros(6,6);
    M_x =zeros(6,6);
    M_xx=zeros(6,6);
end

%Evaluation of the Moment matrix and its derivatives (M, M_x, M_xx), here
%we take the node positions
%as xxI and yyI for the summation process to be carried out easily
for integer_2=1:length(P)
    xxI = P(integer_2);       

    if basis==1
        h = [1 ;x-xxI];
        h_x =[0; 1];
        h_xx = [0; 0];
        
    elseif basis==2
        h = [1 ;x-xxI; (x-xxI)^2];
        h_x = [0; 1; (2*(x-xxI))];
        h_xx = [0 ; 0; 2];

    elseif basis==3
        h= [1 ;x-xxI; (x-xxI)^2;(x-xxI)^3]; 
        h_x=[0; 1; (2*(x-xxI)); (3*(x-xxI)^2)];
        h_xx = [0 ; 0; 2; 6*(x-xxI) ];
        
    elseif basis==4
        h = [1;x-xxI;(x-xxI)^2;(x-xxI)^3;(x-xxI)^4];
        h_x=[0; 1; (2*(x-xxI)); (3*(x-xxI)^2); 4*(x-xxI)^3];
        h_xx = [0 ; 0; 2; 6*(x-xxI);12*(x-xxI)^2 ];
        
    elseif basis==5
        h = [1;x-xxI;(x-xxI)^2;(x-xxI)^3;(x-xxI)^4;(x-xxI)^5];
        h_x=[0; 1; (2*(x-xxI)); (3*(x-xxI)^2); 4*(x-xxI)^3; 5*(x-xxI)^4];
        h_xx = [0 ; 0; 2; 6*(x-xxI);12*(x-xxI)^2 ;20*(x-xxI)^3];
    end    

    %For FINDING PHI
    clear zz
    zz = (sqrt((x-xxI)^2))/ss;
    [phi] = phi_eval(zz);
    [dphi] = dphi_eval(x,xxI,zz,ss); 
    [ddphi] = ddphi_eval(x,xxI,zz,ss);     
 
    M = M + (h*transpose(h)*phi);
    M_x = M_x + (h_x*transpose(h)*phi) + (h*transpose(h_x)*phi) + (h*transpose(h)*dphi);
    M_xx = M_xx + (h_xx*transpose(h)*phi) + (h_x*transpose(h_x)*phi) + (h_x*transpose(h)*dphi)  + (h_x*transpose(h_x)*phi) + (h*transpose(h_xx)*phi) + (h*transpose(h_x)*dphi)+ (h_x*transpose(h)*dphi) + (h*transpose(h_x)*dphi) + (h*transpose(h)*ddphi); 
end

%--------------------------------------------------------------------

%Evaluating inv(M),x i.e derivative of inv of M
InvM_x = -1*inv(M)*M_x*inv(M);
InvM_xx = -inv(M)*((M_xx*inv(M))+(2*M_x*InvM_x));

%--------------------------------------------------------------------

%MAIN matrices for computing Derivative of Shape function at each of the xI
SIxx = zeros(1,length(NS));
for int_1 = 1:length(NS)
    xI = NS(int_1);
    
    % Derivatives of H's
    if basis==1
        Ho = [1 ;0];
        H = [1 ;x-xI];
        H_x =[0; 1];
        H_xx = [0; 0];

    elseif basis==2
        Ho = [1 ;0;0];
        H = [1 ;x-xI; (x-xI)^2];  
        H_x = [0; 1; (2*(x-xI))];
        H_xx = [0 ; 0; 2]; 

    elseif basis==3
        Ho = [1 ;0;0;0];
        H = [1;x-xI;(x-xI)^2;(x-xI)^3]; 
        H_x=[0; 1; (2*(x-xI)); 3*(x-xI)^2];   
        H_xx = [0 ; 0 ; 2 ;6*(x-xI) ];

    elseif basis==4
        Ho = [1 ;0;0;0;0];
        H = [1;x-xI;(x-xI)^2;(x-xI)^3;(x-xI)^4];
        H_x=[0; 1; (2*(x-xI)); (3*(x-xI)^2); 4*(x-xI)^3];
        H_xx = [0 ; 0; 2; 6*(x-xI);12*(x-xI)^2 ];

    elseif basis==5
        Ho = [1 ;0;0;0;0;0];
        H = [1;x-xI;(x-xI)^2;(x-xI)^3;(x-xI)^4;(x-xI)^5];
        H_x=[0; 1; (2*(x-xI)); (3*(x-xI)^2); 4*(x-xI)^3; 5*(x-xI)^4];
        H_xx = [0 ; 0; 2; 6*(x-xI);12*(x-xI)^2 ;20*(x-xI)^3];
    end

    % For finding PHI
    clear z
    z=(sqrt((x-xI)^2))/ss;

    [PHI] = phi_eval(z);
    [DPHI] = dphi_eval(x,xI,z,ss); 
    [DDPHI] = ddphi_eval(x,xI,z,ss); 

    SIxx(int_1) = transpose(Ho)*((InvM_xx*H*PHI)+(inv(M)*H_xx*PHI)+(inv(M)*H*DDPHI)+ (2*InvM_x*H_x*PHI)+(2*InvM_x*H*DPHI)+(2*inv(M)*H_x*DPHI));
end
%================= FUNCTIONS ====================
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

function [dphi] = dphi_eval(x,xxI,zz,ss)

if zz>=0 && zz<(1/3)
    dphi = -(9*(x - xxI)*(45*((x - xxI)^2)^(3/2) - 36*ss*x^2 - 36*ss*xxI^2 + 4*ss^3 + 72*ss*x*xxI))/(4*ss^5);
else if zz>=1/3 && zz<2/3
        dphi = (3*(x - xxI)*(135*((x - xxI)^2)^(3/2) + 270*ss^2*((x - xxI)^2)^(1/2) + (5*ss^4)/((x - xxI)^2)^(1/2) - 324*ss*x^2 - 324*ss*xxI^2 - 84*ss^3 + 648*ss*x*xxI))/(8*ss^5);
    else if zz>=2/3 && zz<1
            dphi = -(81*(x - xxI)*(((x - xxI)^2)^(3/2) + 6*ss^2*((x - xxI)^2)^(1/2) + ss^4/((x - xxI)^2)^(1/2) - 4*ss*x^2 - 4*ss*xxI^2 - 4*ss^3 + 8*ss*x*xxI))/(8*ss^5);
        else if zz>=1
                dphi = 0;
            end
        end
    end
end   

function [ddphi] = ddphi_eval(x,xxI,zz,ss)

if zz>=0 && zz<(1/3)
    ddphi = -(405*((x - xxI)^2)^(3/2) - 243*ss*x^2 - 243*ss*xxI^2 + 9*ss^3 + 486*ss*x*xxI)/ss^5;
else if zz>=1/3 && zz<2/3
        ddphi = (405*(x - xxI)^6 + 405*ss^2*(x - xxI)^4 - 729*ss*((x - xxI)^2)^(5/2) - 63*ss^3*((x - xxI)^2)^(3/2))/(2*ss^5*((x - xxI)^2)^(1/2)*(x - xxI)^2);
    else if zz>=2/3 && zz<1
            ddphi = -(81*(x - xxI)^6 + 243*ss^2*(x - xxI)^4 - 243*ss*((x - xxI)^2)^(5/2) - 81*ss^3*((x - xxI)^2)^(3/2))/(2*ss^5*((x - xxI)^2)^(1/2)*(x - xxI)^2);
        else if zz>=1
                ddphi = 0;
            end
        end
    end
end

