clc
close all

%==========================================================
% SOLVING THE 1D DYNAMIC ELASTIC BAR PROBLEM WITH VARYING E
%==========================================================
yes2 = 0; % yes2 = 1 movies run
% E is interpolated using linear RK shape functions

% IC is a sine wave

% Input Parameters 
marker_size = 50;

L = 16;
rho = 1; % Density
k = (2*pi)/L; % Can be used in IC
dt = 0.2;
T = 6; % Total time for the full solution
T_red = 6; % Total time for the reduced solution
T_snap = 6; %Time upto which snapshots are taken
sp_snap = 1; % number of dt spacings between snapshots
no_reduced = 25; 

%Domain
xdim1=0;
xdim2=L;

x_interface = L/2; % Point of discontimuity of materials
E1 = 0.2; 
E2 = 1;

% No. of Source and Collocation points
Source_pts = 51;
Coll_pts = (Source_pts*2)-1;
basis = 2; 

% Nodal distance of Source and Collocation points
DX_source = (xdim2-xdim1)/(Source_pts-1);
DX_coll = (xdim2-xdim1)/(Coll_pts-1);

%---------------------------------------------------------------------
%--------------------SOURCE POINTS------------------------------------
%---------------------------------------------------------------------
% Nodal distance of source points
dx = DX_source;

% Constructing the nodes in the x direction
x_s = xdim1:dx:xdim2; % Row vector
NS = transpose(x_s);  % Column % Source points matrix
fprintf('No. of Source points = %d \n',length(x_s));

figure(1); hold on %Source points
for int1 = 1:size(NS,1)
    scatter(NS(int1),0,marker_size,'r','filled')   
    % Here the 0 is the y-value for the 2D plot
end

xlabel('x','FontSize',16)
ylabel('y','FontSize',16)
title('Source Points','FontSize',16)
axis equal
grid('on')
hold off

%---------------------------------------------------------------------
%--------------------COLLOCATION POINTS-------------------------------
%---------------------------------------------------------------------
% Nodal distance of collocation points
dx_c = DX_coll;

% Constructing the nodes in the x direction
x_c = xdim1:dx_c:xdim2;
fprintf('No. of Collocation points = %d\n',length(x_c));

N_total = transpose(x_c);   % Total Collocation points     

% Constructing the INTERIOR COLLOCATION POINTS NI_c
NI_c = N_total(2:length(x_c)-1);

%Constructing Nodes on the Essential Boundary

% In 1D there are 2 boundary points and both are EBC's for this problem
NEB = zeros(2,1);
NEB(1)=xdim1;
NEB(2)=xdim2;

fprintf('No. of EB points %d \n',length(NEB));


% FORMING THE COLLOCATION POINTS MATRIX AND PLOTTING
NC = [NI_c;NEB];
% NC Places the Collocation points in order - first the interior points
% then the points on the EB

figure(2)
hold on 
%Collocation points figure
legd  = zeros(1,2); % for Legend 
int_4 =1;
 
for int_1 = 1:length(NI_c)
    legd(int_4)=scatter(NI_c(int_1,1),0,marker_size,'b','filled');
    % 35 is the marker size
end
int_4 = int_4 + 1;
for int_1 = 1:length(NEB)
    legd(int_4) = scatter(NEB(int_1,1),0,marker_size,'g','filled') ;
    % 50 is the marker size
end

xlabel('x','FontSize',16)
ylabel('y','FontSize',16)
title('Collocation Points','FontSize',16)
axis equal
grid('on')
legend(legd,{'Interior','Essential'})
hold off


%=========================================================
%PARAMETERS

% No. of source points and total number of collocation points
no_NS = size(NS,1);
no_NC = size(NC,1);
h =(xdim2-xdim1)/(no_NS-1);
ss = (basis+1)*h; % Support size for the RK SF

NSt = transpose(NS);% to be consistent with the dimensions of P in the function

sigg = max(rho,dt^2*no_NS^2);
sq_betag = (1/no_NS)*sigg;  % Weight for the essential boundary

% For interpolating E
basis_E = 1;
ss_E = (basis_E+3)*h;

%-------------------------------------------------------------------------
% Interpolating E using the source points. Value of E at any collocation
% point can be found using this
E_NS = zeros(size(NS));
for i = 1:length(NS)
    if NS(i)<=x_interface
    E_NS(i) = E1;
    else
       E_NS(i) = E2;
    end
end

%==========================================================
% Initial conditions
%==========================================================
% Computing d0 and d1

d0 = zeros(no_NS,1); % Since u(x,0) = 0
%---------------------------------------
% Calculating Initial velocity at the collocation points
vo = zeros(no_NC,1);
for int1 = 1:length(vo)
    xc = NC(int1);
    vo(int1) = sin(k*xc);
end

%---------------------------------------
H = zeros(no_NC,no_NS); 
int_row = 1; 

for count9 = 1:length(NC) %defines row
    x = NC(count9); % Collocation point
    % Getting the nodes in the support of x for constructing the moment
    % matrix
    [P] = required_nodes(x,NSt,ss);
    [SI] = SF_1D(x,NS,basis,P,ss);
    % Evaluating SI at 'x' at all the NS points
        
    H(int_row,:) = SI(:);
    int_row = int_row + 1;
end
%---------------------------------------
d1 = H\(dt*vo+H*d0);

%==========================================================
% Constructing M, A1, A3 Matrices which remain the same and are independent
% of time
% A1 corresponds to the points satisfying the domain equation
% A3 corresponds to the points satisfying the EBC

%---------------------------------------
M = zeros(length(NI_c),no_NS);
int_row = 1; 

for count9=1:length(NI_c) %defines row
    x = NI_c(count9); % Collocation point    
    [P] = required_nodes(x,NSt,ss);
    [SI] = SF_1D(x,NS,basis,P,ss); 
        
    M(int_row,:) = SI(:);
    int_row = int_row + 1;
end
M = rho*M;

%---------------------------------------
A1 = zeros(length(NI_c), no_NS);

E_vec = zeros(length(NI_c),1);
E_vec_x = zeros(length(NI_c),1);

int_row = 1; % Row counter for A1
int_col = 1; % Column counter for A1
for int_1 = 1:length(NI_c) % Looping over the no. of collocation pts. in the domain
    % arranging rows of A
    % We need the RK shape function evaluated at the collocation point x,y 
    % and centered at the Source point

    x = NI_c(int_1);    % Collocation point    
    [P] = required_nodes(x,NSt,ss);
    % Instead of 'P' we can also send Nst which is all the nodes
    [SIx] = DSF_x(x,NS,basis,P,ss);
    [SIxx] = DSF_xx(x,NS,basis,P,ss);
    
    % Finding E(x) i.e E at the collocation point
    [EE,EE_x] = E_interpolate(x,NS,basis_E,ss_E);
    % EE and EE_x give the vector of SF
    E = EE*E_NS;
    E_x = EE_x*E_NS;
    
    E_vec(int_1) = E;
    E_vec_x(int_1) = E_x;
        
    for int_2 = 1:length(NS)  %Looping over the no. of source pts.

        % Arranging in the A matrix
        A1(int_row,int_col) = (E_x*SIx(int_2)+E*SIxx(int_2));
        int_col = int_col + 1;
    end
    int_col=1;
    int_row = int_row+1;
end

%------------------------------------------------------------------------
% Essential boundary points

A3 = zeros(length(NEB), no_NS);

int_row = 1; % Row counter for A3

for int_1 = 1:length(NEB) % Looping over the no. of collocation pts. in the domain    
    
    x = NEB(int_1);    % Collocation point     
    [P] = required_nodes(x,NSt,ss);
    [SI] = SF_1D(x,NS,basis,P,ss);
    SI = sq_betag *SI;
        
    A3(int_row,:) = SI(:);
    int_row = int_row+1;
end

%===============================================================
% ======================= Time Update===========================
%===============================================================


%---------------------
% FOR PLOTS:
%---------------------
% For finding the solution at x=L/2 we need phi vector
xmid = L/2;
[P] = required_nodes(xmid,NSt,ss);
[phi] = SF_1D(xmid,NS,basis,P,ss);

%-----------------------------------------------------
% For finding the solution at a set of points
x_plot = xdim1:0.5:xdim2;
phi_all = zeros(length(x_plot),no_NS);
dphi_all = zeros(length(x_plot),no_NS); % For strains

int_row = 1; 

for count9 = 1:length(x_plot) %defines row
    x = x_plot(count9); % Collocation point    
    [P] = required_nodes(x,NSt,ss);
    [SI] = SF_1D(x,NS,basis,P,ss); 
    [SIx] = DSF_x(x,NS,basis,P,ss); 
    
    phi_all(int_row,:) = SI(:);
    dphi_all(int_row,:) = SIx(:);
    
    int_row = int_row + 1;
end

%---------------------
% FOR ERROR
%---------------------
% For finding the RMS error at a set of points xE
xE = xdim1:0.01:xdim2;
phi_ERR = zeros(length(xE),no_NS);

int_row = 1;

for count9 = 1:length(xE) %defines row
    x = xE(count9); % Collocation point    
    [P] = required_nodes(x,NSt,ss);
    [SI] = SF_1D(x,NS,basis,P,ss);
    
    phi_ERR(int_row,:) = SI(:);
    int_row = int_row + 1;
end


%--------------------------------------------------
% Time steps
%--------------------------------------------------    
t =  0:dt:T;
len_t = length(t);

%-------------------------------------------------- 
% For finding the solution at x = L/2
%-------------------------------------------------- 
umid = zeros(len_t,1); % u at x = L/2
time = zeros(len_t,1);

% Finding the values of u in the first 2 time steps
% from the given IC's
umid(1) = phi*d0; 
umid(2) = phi*d1; 
time(1) = 0; 
time(2) = dt; 

%--------------------------------------------------
% Snapshot matrix
%-------------------------------------------------- 
len_Msnaps = length(2*dt:sp_snap*dt:T_snap);
% Getting the snapshots at different times at the Collocation points
M_snaps = len_Msnaps +2; % Number of snapshots; +2 is for the IC's
% M_snaps is same as = length(t)

S = zeros(no_NS,M_snaps); % Snapshot or Response matrix
S(:,1) = d0;
S(:,2) = d1;

%--------------------------------------------------

integer = 3; % Column number for arranging the snapshot matrix
integer_snap = 3;
timestep = 1; % Time step counter
int_snap = 1; % Time counter for snapshots
T_COUNT = zeros(len_t,1);
T_COUNT(2) = dt;
%T_COUNT(1) = 0;

%---------------------------------------
% Collecting the solution at all points at different times
u_plot_all_T = zeros(length(x_plot),M_snaps);
% Initial conditions
u_plot_all_T(:,1) = phi_all*d0;
u_plot_all_T(:,2) = phi_all*d1;

%---------------------------------------
% Collecting the strains at all points at different times
du_plot_all_T = zeros(length(x_plot),M_snaps);
% Initial conditions
du_plot_all_T(:,1) = dphi_all*d0;
du_plot_all_T(:,2) = dphi_all*d1;

% Error anlysis
% Saving the solution at xE points in every time step
u_ERR_full = zeros(length(xE),M_snaps);
% Initial conditions
u_ERR_full(:,1) = phi_ERR*d0;
u_ERR_full(:,2) = phi_ERR*d1;

%------------
% Timer 
%------------
tic;

for t_count = (2*dt):dt:T

    T_COUNT(integer) = t_count; % Getting all the time step values in a vector
    
    % Starting from 0.5 since:
    % t=0, d0
    % t=0.25, d1
    
    if timestep==1 % That is only for the first step
        dn_minus_1 = d0;
        dn = d1;
    else
        dn_minus_1 = dn;
        dn = dn_plus_1;  
    end
    
    % Calculating b1 and b3    
    b1 = 2*M*dn - M*dn_minus_1 + dt^2*A1*dn;
    b3 = zeros(length(NEB),1)*sq_betag; % u=0 on EB
    
    %-----------------------------------------------
    % Finding dn_plus_1
    A = [M;A3];
    F = [b1;b3];
    
    dn_plus_1 = A\F;    
   
    %--------------------------------------------------    
    % Finding the solution at x = L/2 
    % This is for plotting the u at x = L/2 at different times

    umid(integer) = phi*dn_plus_1;
    time(integer) = t_count;    
   
    if int_snap==timestep
        
        S(:,integer_snap) = dn_plus_1;
        int_snap = int_snap + sp_snap;
        integer_snap = integer_snap + 1;
    end      
     
    %--------------------------------------------------
    % Finding the solution at all points when T = 6
    if timestep == M_snaps-2
        u_T = phi_all * dn_plus_1;
    end    

    u_plot_all_T(:,integer) = phi_all * dn_plus_1;
    du_plot_all_T(:,integer) = dphi_all * dn_plus_1; 
        
    %--------------------------------------------------
    % For error analysis
    u_ERR_full(:,integer) = phi_ERR * dn_plus_1; 
  
    %-----------------------------------------------
    integer = integer + 1;
    timestep = timestep + 1;
end
timer_full= toc;

  
%=========================================================================
% ==================== REDUCED ORDER MODELING=============================
%=========================================================================

[U,SIG,V] = svd(S); % Here S = U*SIG*V'

% Here sizes are:
% S = nxS ; n = length of the 'd' vector = no_NS in 1D
% U = nxn 
% SIG = nxS
% V_transpose = SxS

%------------------------------------------------

Singularvalues = diag(SIG);
Eigenvalues=(Singularvalues.^2)/M_snaps;

% Printing EV to text file
fileID = fopen('EV.txt','w');
fprintf(fileID,'Eigen Values\n');
for i=1:length(Eigenvalues)
    fprintf(fileID,'%4.4e \n',Eigenvalues(i)); % Printing in E format
end
fclose(fileID);

%------------------------------------------------
% Taking the first 'r' modes

U_red = U(:,1:no_reduced); % Size of the reduced U = nxr

%-------------------------------------------------
% Percentage energy captured by the first mode
Energy = sum(Eigenvalues(1:no_reduced))/sum(Eigenvalues);

%-------------------------------------------------

fprintf('\nsize of U = %d %d \n', size(U))
fprintf('size of SIG = %d %d \n', size(SIG))
fprintf('size of V = %d %d \n', size(V))
fprintf('size of U_red = %d %d \n', size(U_red))
fprintf('Energy in r modes = %f\n', Energy)

%---------------------------------------------------------------
% Solving the problem again using the projection matrix 'U' to find the
% reduced order solution at different times at x = L/2
%--------------------------------------------------  

t_red = 0:dt:T_red;
len_t_red = length(t_red);
% len_t_red-2
umid_red = zeros(len_t_red,1); % u at x = L/2
time_reduced = zeros(len_t_red,1);

% Finding the values of u in the first 2 time steps
% from the given IC's
% umid_red = reduced solution at umid
umid_red(1) = phi*d0; 
umid_red(2) = phi*d1; 
time_reduced(1) = 0; 
time_reduced(2) = dt; 

%--------------------------------------------------
integer = 3;
timestep = 1; % Time step counter
%---------------------------------------
% Collecting the solution at all points at different times
u_plot_all_T_red = zeros(length(x_plot),M_snaps);
% Initial conditions
u_plot_all_T_red(:,1) = phi_all*d0;
u_plot_all_T_red(:,2) = phi_all*d1;
%---------------------------------------
% Collecting the solution at all points at different times
du_plot_all_T_red = zeros(length(x_plot),M_snaps);
% Initial conditions
du_plot_all_T_red(:,1) = dphi_all*d0;
du_plot_all_T_red(:,2) = dphi_all*d1;
%---------------------------------------

% Error anlysis
% Saving the solution at xE points in every time step
u_ERR_red = zeros(length(xE),M_snaps);
% Initial conditions
u_ERR_red(:,1) = phi_ERR*d0;
u_ERR_red(:,2) = phi_ERR*d1;

%--------------
% Timer 
%--------------
tic;

for t_count = (2*dt):dt:T_red

    % Starting from 0.5 since:
    % t=0, d0
    % t=0.25, d1
    
    if timestep==1 % That is only for the first step
        dn_minus_1 = d0;
        dn = d1;
    else
        dn_minus_1 = dn;
        dn = dn_plus_1_p;  
    end
    
    % Calculating b1 and b3    
    b1 = 2*M*dn - M*dn_minus_1 + dt^2*A1*dn;
    b3 = zeros(length(NEB),1)*sq_betag; % u=0 on EB
    
    %-----------------------------------------------
    % Petrov Galerkin projection, NC>NS
    %-----------------------------------------------
    % Finding dn_plus_1_reduced
    K = [M;A3];
    F = [b1;b3];
    
    % Reduced model
    K_hat = K*U_red;
    K_hat_T = transpose(K_hat);
    
    dn_plus_1_red = (K_hat_T*K_hat)\(K_hat_T*F);
    dn_plus_1_p = U_red * dn_plus_1_red; % Reduced solution, p-projected   
    
   
%     %-----------------------------------------------
%     % Galerkin projection, NC=NS
%     %-----------------------------------------------
%     % Reduced model
% 
%     dn_plus_1_red = (transpose(U_red)*K*U_red)\(transpose(U_red)*F);
%     dn_plus_1_p = U_red * dn_plus_1_red; % Reduced solution, p-projected
    
    %-----------------------------------------------
    % Finding the solution at all points when time = T 
 
    if timestep==len_t_red-2       
        u_T_red = phi_all * dn_plus_1_p;
    end
    
    u_plot_all_T_red(:,integer) = phi_all * dn_plus_1_p;    
    du_plot_all_T_red(:,integer) = dphi_all * dn_plus_1_p;    
            
    %--------------------------------------------------
    % For error analysis
    u_ERR_red(:,integer) = phi_ERR * dn_plus_1_p; 
        
    %-----------------------------------------------
    % Finding the solution at x = L/2    

    umid_red(integer) = phi*dn_plus_1_p;
    time_reduced(integer) = t_count; 
    %-----------------------------------------------

    integer = integer + 1;
    timestep = timestep + 1;
end
timer_red = toc;

% Finding the % decrease in time
per_red_time = ((timer_full-timer_red)*100)/timer_full;

%--------------------------------------------------------------------

figure(4)
plot(time,umid,'ok',time_reduced,umid_red,'+k','LineWidth',2,'MarkerSize',10)
xlabel('time t','FontSize',16,'FontName','TimesNewRoman')
ylabel('u(x) at x = L/2','FontSize',16,'FontName','TimesNewRoman')
legend('RKCM-Full',sprintf('RKCM-POD: r = %d',no_reduced))
set(gca,'FontSize',16,'FontName','TimesNewRoman')

% %--------------------------------------------------------------------
% ONLY WORKS WHEN TOTAL TIMES T = T_red
% time;
% uuu = umid-umid_red;
% 
% figure(5)
% plot(time,uuu,'-or','LineWidth',2.5)
% xlabel('time','FontSize',14)
% ylabel('u-full - u-red','FontSize',14)
% title('error: u at x = L/2','FontSize',14)

%--------------------------------------------------------------------

RRR = 1:M_snaps;
figure(6)
plot(RRR,Eigenvalues,'-ok','LineWidth',2.5)
xlabel('Mode number','FontSize',16,'FontName','TimesNewRoman')
ylabel('Eigenvalue \lambda','FontSize',16,'FontName','TimesNewRoman')
% title('Decay of Eigenvalues','FontSize',16,'FontName','TimesNewRoman')
set(gca,'FontSize',16,'FontName', 'TimesNewRoman')
%--------------------------------------------------------------------

figure(7)
plot(x_plot,u_T,'ok',x_plot,u_T_red,'+k','LineWidth',2,'MarkerSize',10)
xlabel('x','FontSize',16,'FontName','TimesNewRoman')
ylabel('u(x)','FontSize',16,'FontName','TimesNewRoman')
% title('Solution at T = 6','FontSize',20,'FontName','TimesNewRoman')
legend('RKCM-Full',sprintf('RKCM-POD: r = %d',no_reduced)) 
set(gca,'FontSize',16,'FontName', 'TimesNewRoman')
%--------------------------------------------------------------------

figure(8)
plot(NI_c,E_vec,'-ok','LineWidth',2,'MarkerSize',5)
% title('Interpolated E','FontSize',16)
xlabel('x','FontSize', 16,'FontName','TimesNewRoman')
ylabel('E','FontSize', 16,'FontName','TimesNewRoman')
% legend('RK interpolation')
set(gca,'FontSize',16,'FontName', 'TimesNewRoman')
xlim([0 16])
%--------------------------------------------------------------------

figure(9)
plot(NI_c,E_vec_x,'-ok','LineWidth',2,'MarkerSize',5)
% title('Interpolated E_x','FontSize',16)
xlabel('x','FontSize', 16,'FontName','TimesNewRoman')
ylabel('E_x','FontSize', 16,'FontName','TimesNewRoman')
% legend('RK interpolation')
set(gca,'FontSize',16,'FontName', 'TimesNewRoman')
xlim([0 16])

%--------------------------------------------------------------------
fprintf('Time full model = %4.4f \n', timer_full)
fprintf('Time reduced model = %4.4f \n', timer_red)
fprintf('Percentage Reduction in time = %4.4f \n',per_red_time)
%--------------------------------------------------------------------

if yes2 ==1 
% This plot is for comparing the reduced and the full solution at different
% times and saving figures to make a movie.
int1 = 1;
for ii = 1:size(u_plot_all_T,2)
    plot200 = figure(200);    
    plot(x_plot,u_plot_all_T(:,ii),'+g',x_plot,u_plot_all_T_red(:,ii),'ok','LineWidth',2.5)
    xlabel('x','FontSize',14)
    ylabel('u(x)','FontSize',14)
    title(sprintf('Solution at T = %4.4f',T_COUNT(int1)),'FontSize',14)
    legend('RKCM Full',sprintf('RKCM-POD: r = %d',no_reduced)) 
    xlim([xdim1 xdim2])
    ylim([-10 10])
    set(gca,'FontSize',14)
    nam = sprintf('movie snaps both/%d',int1);
    saveas(plot200,nam,'png')
    
    int1 = int1 + 1;
end

%----------------------------------------------------------------
% This plot is for comparing the reduced and the full strains at different
% times and saving figures to make a movie.
int1 = 1;
for ii = 1:size(u_plot_all_T,2)
    plot201 = figure(201);    
    plot(x_plot,du_plot_all_T(:,ii),'+g',x_plot,du_plot_all_T_red(:,ii),'ok','LineWidth',2.5)
    xlabel('x','FontSize',14)
    ylabel('du/dx','FontSize',14)
    title(sprintf('Solution at T = %4.4f',T_COUNT(int1)),'FontSize',14)
    legend('RKCM Full',sprintf('RKCM-POD: r = %d',no_reduced)) 
    xlim([xdim1 xdim2])
    ylim([-10 10])
    set(gca,'FontSize',14)
    nam = sprintf('movie snaps both strains/%d',int1);
    saveas(plot201,nam,'png')
    
    int1 = int1 + 1;
end
end

% ERROR Analysis using RMS 
% ERROR_TS = error in each time step
[error1,ERROR_TS] = func_error(u_ERR_full,u_ERR_red);

fprintf('Error in u, all time steps = %4.3E\n',error1)


% -----------------
fprintf('No. reduced = %d\n',no_reduced)

% no_reduced % of full model
per_red = (no_reduced*100)/no_NS;

fprintf('Percent of full model = %f\n',per_red)




