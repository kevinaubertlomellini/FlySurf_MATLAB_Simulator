clc, clear all, close all

% MODEL PARAMETERS
k=0.45;
k2=0.6;
c1 = 0.1;
c2 = 0.1; 
g = 9.81;
n_points=9;
n_points2=9;
l0=1/(n_points-1);
m_total=0.1;
m=m_total/(n_points2*n_points)*ones(n_points,n_points2);
m_uav=0.032;

[X, Y] = ndgrid(1:1:n_points, 1:1:n_points2);
x_visible_points = [X(:), Y(:)];
[X, Y] = ndgrid(1:4:n_points, 1:4:n_points2);
x_actuators= [X(:), Y(:)];
n_actuators= size(x_actuators,1);
n_visible_points=size(x_visible_points,1);
x_actuators_2=zeros(n_actuators,3);
x_visible_points_2=zeros(n_visible_points,3);

for i=1:1:n_actuators
    m(x_actuators(i,1),x_actuators(i,2))=m_uav;
    x_actuators_2(i,1)=6*n_points*(x_actuators(i,2)-1)+6*(x_actuators(i,1)-1)+1;
    x_actuators_2(i,2)=6*n_points*(x_actuators(i,2)-1)+6*(x_actuators(i,1)-1)+2;
    x_actuators_2(i,3)=6*n_points*(x_actuators(i,2)-1)+6*(x_actuators(i,1)-1)+3;
end

for i=1:1:n_visible_points
    x_visible_points_2(i,1)=6*n_points*(x_visible_points(i,2)-1)+6*(x_visible_points(i,1)-1)+1;
    x_visible_points_2(i,2)=6*n_points*(x_visible_points(i,2)-1)+6*(x_visible_points(i,1)-1)+2;
    x_visible_points_2(i,3)=6*n_points*(x_visible_points(i,2)-1)+6*(x_visible_points(i,1)-1)+3;
end
%x_visible_points_2=x_actuators_2;

% SIMULATION PARAMETERS
delta = 0.02;
steps_change=30/0.02;
iter = steps_change*10;

% INTIAL STATE
x = zeros(n_points*n_points2*6,1);
x5=x;

x_z=0:0.1:n_points2*n_points;
for i=1:1:n_points
    for j=1:1:n_points2
        x(6*n_points*(j-1)+6*(i-1)+1)=l0*(i-1);
        x(6*n_points*(j-1)+6*(i-1)+2)=l0*(j-1);
        x5(6*n_points*(j-1)+6*(i-1)+1)=l0*(i-1);
        x5(6*n_points*(j-1)+6*(i-1)+2)=l0*(j-1);
        x5(6*n_points*(j-1)+6*(i-1)+3)=x_z(6*n_points*(j-1)+6*(i-1)+1);
    end
end
x(1:6:end)=x(1:6:end);
x(2:6:end)=x(2:6:end)-0.5;
x_estimated=x;


% LINEARIZED SYSTEM
T0=1;
theta0=0.0;
phi0=0.0;

%funciona

A= A_linearized(x5*1.02,k,k2,c1,c2,m,l0,n_points,n_points2);
B = zeros(n_points*n_points2*2*3,3*n_actuators); 
B2 = zeros(n_points*n_points2*2*3,3*n_actuators); 
B_matrix_3=[sin(theta0),0,T0*cos(theta0);-cos(theta0)*sin(phi0), -T0*cos(phi0)*cos(theta0),  T0*sin(phi0)*sin(theta0);cos(phi0)*cos(theta0), -T0*cos(theta0)*sin(phi0), -T0*cos(phi0)*sin(theta0)];
for i=1:1:n_actuators
    B(6*n_points*(x_actuators(i,2)-1)+6*(x_actuators(i,1)-1)+4:6*n_points*(x_actuators(i,2)-1)+6*(x_actuators(i,1)-1)+6,3*i-2:3*i) = B_matrix_3;
    B2(6*n_points*(x_actuators(i,2)-1)+6*(x_actuators(i,1)-1)+4:6*n_points*(x_actuators(i,2)-1)+6*(x_actuators(i,1)-1)+6,3*i-2:3*i) = eye(3);
end    
C = zeros(6*n_visible_points,n_points*n_points2*2*3);
for i=1:1:n_visible_points
    C(6*(i-1)+1:6*i,x_visible_points_2(i,1):x_visible_points_2(i,1)+5)=eye(6);
end
D=0;

sys_cont = ss(A, B, C, D);
sys_discrete = c2d(sys_cont, delta);
[A_d, B_d, C_d, D_d] = ssdata(sys_discrete);

% CONTROL PARAMETERS
Q = 1*eye(6*n_points*n_points2); 
for yu=1:1:n_points*n_points2
    Q(6*yu-3,6*yu-3) = 2500*eye(1);
    Q(6*yu-5:6*yu-4,6*yu-5:6*yu-4) = 700*eye(2);
    Q(6*yu-2:6*yu-1,6*yu-2:6*yu-1) = 10*eye(2);
end
for yu=1:1:n_actuators
    Q(x_actuators_2(yu,1)+5,x_actuators_2(yu,1)+5) = 1900*eye(1);
    Q(x_actuators_2(yu,1)+3:x_actuators_2(yu,1)+4,x_actuators_2(yu,1)+3:x_actuators_2(yu,1)+4) = 400*eye(2);
end

R = 100*eye(3*n_actuators);
for yi=1:1:n_actuators
    R(3*yi-1:3*yi,3*yi-1:3*yi)=20*eye(2);
end


Pe_d=zeros(6*n_points*n_points2);
L_d=zeros(6*n_points*n_points2,6*n_visible_points);
Qe_d = 50*eye(6*n_points*n_points2); 
for yu=1:1:n_visible_points
    Qe_d(x_visible_points_2(yu,1):x_visible_points_2(yu,1)+5,x_visible_points_2(yu,1):x_visible_points_2(yu,1)+5) = 60*eye(6);
end
Re_d = 60*eye(6*n_visible_points); 
for yu=1:1:n_visible_points
    Re_d(6*(yu-1)+4:6*(yu-1)+6,6*(yu-1)+4:6*(yu-1)+6) = 140*eye(3);
end

[K,P,E] = dlqr(A_d,B_d,Q,R);

% DATA MATRICES
u_save=zeros(3*n_actuators,iter);
u0_save=zeros(3*n_actuators,iter);
u1_save=zeros(3*n_actuators,iter);
u2_save=zeros(3*n_actuators,iter);
y_save=zeros(6*n_visible_points,iter);
y2_save=zeros(6*n_visible_points,iter);
x_save=zeros(6*n_points*n_points2,iter);
x_e_save=zeros(6*n_points*n_points2,iter);
x_d_save=zeros(6*n_points*n_points2,iter);
error_real=zeros(1,iter);
error_virtual_d=zeros(1,iter);
error_ed=zeros(1,iter);

xd=x;

x_g = linspace(-0.375, 0.375, (n_points));
y_g = linspace(-0.375, 0.375, (n_points2));

[Y_g,X_g] = meshgrid(x_g, y_g);

x_g_vector0 = reshape(X_g, [], 1);
y_g_vector0 = reshape(Y_g, [], 1);

% Define Gaussian parameters
Amp = 1.1;  % amplitude
x0 = 0.0; % center of Gaussian in x
y0 = 0.0; % center of Gaussian in y
sigma_x = 0.575; % standard deviation in x
sigma_y = 0.575; % standard deviation in y

% Calculate the 2D Gaussian
z_g_vector0 = Amp * exp(-((x_g_vector0 - x0).^2 / (2 * sigma_x^2) + (y_g_vector0 - y0).^2 / (2 * sigma_y^2)))-0.5;
x_g_vector0=x_g_vector0+0.5;

max_z_change0 = max(z_g_vector0);
z_g_vector1=-(z_g_vector0-max_z_change0);

[x_change0,y_change0,z_change0,x_dot_change0,y_dot_change0,z_dot_change0]=changes(xd(1:6:end),xd(2:6:end),xd(3:6:end),x_g_vector0,y_g_vector0,z_g_vector1,steps_change,delta,2/3,1/2,n_points,n_points2);

[x_change1,y_change1,z_change1,x_dot_change1,y_dot_change1,z_dot_change1]=changes(x_g_vector0,y_g_vector0,z_g_vector1,0.9*xd(1:6:end)+0.05,0.9*xd(2:6:end),xd(3:6:end)+0.25,steps_change/3,delta,1,1/8,n_points,n_points2);
[x_change1_5,y_change1_5,z_change1_5,x_dot_change1_5,y_dot_change1_5,z_dot_change1_5]=changes(0.9*xd(1:6:end)+0.05,0.9*xd(2:6:end),xd(3:6:end)+0.25,x_g_vector0,y_g_vector0,z_g_vector0,2*steps_change/3,delta,2/3,1/8,n_points,n_points2);

x_change1_5(:,int32(steps_change*2/3*2/3))=x_change1_5(:,int32(steps_change*2/3*2/3)+1);
y_change1_5(:,int32(steps_change*2/3*2/3))=y_change1_5(:,int32(steps_change*2/3*2/3)+1);
z_change1_5(:,int32(steps_change*2/3*2/3))=z_change1_5(:,int32(steps_change*2/3*2/3)+1);

x_change1=[x_change1,x_change1_5];
y_change1=[y_change1,y_change1_5];
z_change1=[z_change1,z_change1_5];
x_dot_change1=[x_dot_change1,x_dot_change1_5];
y_dot_change1=[y_dot_change1,y_dot_change1_5];
z_dot_change1=[z_dot_change1,z_dot_change1_5];

factor=11/12;
% Define the circle's radius
radius = 0.5;
factor2=4;
% Total time for the lap
T0 = factor2*steps_change*delta*factor;  % seconds
N=int32(factor2*steps_change*(factor));
t = linspace(0, T0, N); % Time parameter from 0 to T0
dt = t(2) - t(1); % Time step
% Smooth angular displacement for a complete lap (from 0 to 2*pi)
theta = pi * (1 - cos(pi * t / T0)); % Start/stop with zero velocity, full 2*pi displacement
% Preallocate arrays to store translated positions
Xt = zeros(n_points2*n_points, N);
Yt = zeros(n_points2*n_points, N);
Zt = repmat(z_g_vector0, 1, N); % No Z translation

x_circle = radius * cos(theta);
y_circle = radius * sin(theta);
z_circle = zeros(size(theta)); % Circle in XY plane

for i = 1:N
    % Compute the circular path translation using smooth angular theta
    x_circle = radius * cos(theta(i));
    y_circle = radius * sin(theta(i));
    
    % Apply the translation to the square points
    Xt(:, i) = x_g_vector0 + x_circle-0.5;
    Yt(:, i) = y_g_vector0 + y_circle;
end

% Compute the velocity (numerical derivative of the position for all points)
x_dot_change2 = diff(Xt, 1, 2) / dt; % Derivative along the time axis (columns)
y_dot_change2 = diff(Yt, 1, 2) / dt;
z_dot_change2 = diff(Zt, 1, 2) / dt; % Z remains zero

% Pad velocity arrays to match the dimensions of position arrays (4xN)
x_dot_change2 = [x_dot_change2, zeros(n_points2*n_points, 1)]; % Padding with zero to match dimensions
y_dot_change2 = [y_dot_change2, zeros(n_points2*n_points, 1)];
z_dot_change2 = [z_dot_change2, zeros(n_points2*n_points, 1)];

x_change2=[Xt,repmat(x_g_vector0,1,int32(factor2*steps_change*(1-factor)))];
y_change2=[Yt,repmat(y_g_vector0,1,int32(factor2*steps_change*(1-factor)))];
z_change2=[Zt,repmat(z_g_vector0,1,int32(factor2*steps_change*(1-factor)))];
x_dot_change2=[x_dot_change2,repmat(zeros(n_points*n_points2,1),1,int32(factor2*steps_change*(1-factor)))];
y_dot_change2=[y_dot_change2,repmat(zeros(n_points*n_points2,1),1,int32(factor2*steps_change*(1-factor)))];
z_dot_change2=[z_dot_change2,repmat(zeros(n_points*n_points2,1),1,int32(factor2*steps_change*(1-factor)))];



xd_3 = [xd(1:6:end),xd(2:6:end),xd(3:6:end)+0.5];
% Define the rotation angle in degrees
theta = 45;

% Convert the angle to radians
theta_rad = deg2rad(theta);

% Calculate the center of the square
center = mean(xd_3, 1);

% Shift the square to the origin
shifted_vertices = xd_3 - center;

% Calculate the center of the square

R_y = [cos(theta_rad), 0, sin(theta_rad);
       0,              1, 0;
      -sin(theta_rad), 0, cos(theta_rad)];

% Rotate the square vertices
rotated_vertices = (R_y * shifted_vertices')';

% Translate the rotated vertices back to the original center
xd_rotated= rotated_vertices + center;

factor3=4;
[x_change3a,y_change3a,z_change3a,x_dot_change3a,y_dot_change3a,z_dot_change3a]=changes(x_change2(:,end),y_change2(:,end),z_change2(:,end),xd_rotated(:,1),xd_rotated(:,2),xd_rotated(:,3),steps_change*factor3,delta,11/12,1/8,n_points,n_points2);

x_change3=zeros(n_points2*n_points,int32(factor3*steps_change));
y_change3=zeros(n_points2*n_points,int32(factor3*steps_change));
z_change3=zeros(n_points2*n_points,int32(factor3*steps_change));
x_dot_change3=x_dot_change2+x_dot_change3a;
y_dot_change3=y_dot_change2+y_dot_change3a;
z_dot_change3=z_dot_change2+z_dot_change3a;
x_change3(:,1)=x_change2(:,end);
y_change3(:,1)=y_change2(:,end);
z_change3(:,1)=z_change2(:,end);
for i=2:steps_change*4
    x_change3(:,i)= x_change3(:,i-1)+delta*(x_dot_change3(:,i-1)+x_dot_change3(:,i))/2;
    y_change3(:,i)= y_change3(:,i-1)+delta*(y_dot_change3(:,i-1)+y_dot_change3(:,i))/2;
    z_change3(:,i)= z_change3(:,i-1)+delta*(z_dot_change3(:,i-1)+z_dot_change3(:,i))/2;
end



% INTIIAL DESIRED STATE
xd0=xd;

% INITIAL CONTROL LAW
ud  = zeros(3*n_actuators,1);
ud(1:3:end) = sum(sum(m))*g/n_actuators*ones(n_actuators,1);

% NONLINEAR MODEL
M =  eye(6*n_points*n_points2);
G = zeros(6*n_points*n_points2,1);
for j=1:1:n_points2
    for i =1:1:n_points
        pos_v=6*n_points*(j-1)+6*(i-1)+4:6*n_points*(j-1)+6*(i-1)+6;
        M(pos_v,pos_v)=m(i,j)*eye(3);
        G(6*n_points*(j-1)+6*(i-1)+6)=m(i,j)*g;
    end
end

u_lim_T=1*[0,0.75];
u_lim_M=10*pi/180*[-1,1];

xe_d = x;
u0  = zeros(3*n_actuators,1);
u  = zeros(3*n_actuators,1);

%%

M_inv = M^-1;

% SIMULATION LOOP
for step=1:1:iter
    K_spring = K_matrix(x,k,k2,c1,c2,l0,n_points,n_points2);

    K_spring3 = K_matrix(xe_d,k,k2,c1,c2,l0,n_points,n_points2);

    y=C_d*x;

    y2=y+0.005*(rand(n_visible_points*6,1)-0.5);
   
    xe_d = xe_d + delta*(M_inv*(K_spring3*xe_d-G+B2*u));
    A2= A_linearized(xe_d,k,k2,c1,c2,m,l0,n_points,n_points2);
    A_d = eye(size(A2)) + A2 * delta + (A2 * delta)^2 / 2 + (A2 * delta)^3 / 6;
    Pe_d=A_d*Pe_d*A_d'+Qe_d;

    L_d=Pe_d*C_d'*(C_d*Pe_d*C_d'+Re_d)^-1;
    Pe_d = (eye(6*n_points*n_points2)-L_d*C_d)*Pe_d;
    xe_d = xe_d + L_d*(y2 - C_d*xe_d);
    
    u0 = ud+K*(xd-xe_d);

    for kv=1:1:n_actuators
        if u0(3*kv)<u_lim_M(1)
            u0(3*kv)=u_lim_M(1);
        elseif u0(3*kv)>u_lim_M(2)
            u0(3*kv)=u_lim_M(2);
        end
        if u0(3*kv-1)<u_lim_M(1)
            u0(3*kv-1)=u_lim_M(1);
        elseif u0(3*kv-1)>u_lim_M(2)
            u0(3*kv-1)=u_lim_M(2);
        end
        if u0(3*kv-2)<u_lim_T(1)
            u0(3*kv-2)=u_lim_T(1);
        elseif u0(3*kv-2)>u_lim_T(2)
            u0(3*kv-2)=u_lim_T(2);
        end
    end

    u0_save(:,step)=u0;
    
    
    u1=u0;
    u1(1:3:end)=u1(1:3:end)+0.05*(rand(n_actuators,1)-0.5);
    u1(2:3:end)=u1(2:3:end)+0.01*(rand(n_actuators,1)-0.5);
    u1(3:3:end)=u1(3:3:end)+0.01*(rand(n_actuators,1)-0.5);

    u1_save(:,step)=u1;
    
    for kv=1:1:n_actuators
        T = u0(3*kv-2);
        phi = u0(3*kv-1);
        theta = u0(3*kv);
        u(3*kv-2:3*kv,1) = [T*sin(theta); -T*cos(theta)*sin(phi);T*cos(theta)*cos(phi)];
        T = u1(3*kv-2);
        phi = u1(3*kv-1);
        theta = u1(3*kv);
        u2(3*kv-2:3*kv,1) = [T*sin(theta); -T*cos(theta)*sin(phi);T*cos(theta)*cos(phi)];
    end
    u_save(:,step)=u;
    u2_save(:,step)=u2;
    d_x = M_inv*(K_spring*x-G+B2*u2);
    x = x + d_x*delta;

    x_save(:,step)=x;
    x_e_save(:,step)=xe_d;
    x_d_save(:,step)=xd;
    y_save(:,step)=y;
    y2_save(:,step)=y2;
   
    figure(1)
    plot3(0,0,0)
    hold on
    xGrid = reshape(xe_d(1:6:end), n_points,n_points2);
    yGrid = reshape(xe_d(2:6:end), n_points,n_points2);
    zGrid = reshape(xe_d(3:6:end), n_points,n_points2);
    surf(-xGrid, yGrid, zGrid, 'FaceColor', 'k');
    plot3(-x(1:6:end),x(2:6:end),x(3:6:end),'k*',LineWidth=1.0)
    plot3(-xd(1:6:end),xd(2:6:end),xd(3:6:end),'r*',LineWidth=1.0)
    vec_forces=B2*u;
    quiver3(-x(x_actuators_2(:,1)), x(x_actuators_2(:,2)), x(x_actuators_2(:,3)), 1*vec_forces(x_actuators_2(:,1)+3)/2, 1*vec_forces(x_actuators_2(:,2)+3)/2,1*vec_forces(x_actuators_2(:,3)+3)/2, 0, 'k', 'LineWidth', 1,'MaxHeadSize', 1);


   hold off
    ylim([-1.5,1.5])
    xlim([-1.5,1.5])
    zlim([-0.5,1.0])
    grid


    pause(0.005)

    
   
    if step>steps_change*0 && step<=steps_change*1
        xd(1:6:end)=x_change0(:,step);
        xd(2:6:end)=y_change0(:,step);
        xd(3:6:end)=z_change0(:,step);
        xd(4:6:end)=x_dot_change0(:,step);
        xd(5:6:end)=y_dot_change0(:,step);
        xd(6:6:end)=z_dot_change0(:,step);
    end
    if step>steps_change*1 && step<steps_change*2
        xd(1:6:end)=x_change1(:,step-steps_change*1);
        xd(2:6:end)=y_change1(:,step-steps_change*1);
        xd(3:6:end)=z_change1(:,step-steps_change*1);
        xd(4:6:end)=x_dot_change1(:,step-steps_change*1);
        xd(5:6:end)=y_dot_change1(:,step-steps_change*1);
        xd(6:6:end)=z_dot_change1(:,step-steps_change*1);
    end
    if step>steps_change*2 && step<steps_change*6
        xd(1:6:end)=x_change2(:,step-steps_change*2);
        xd(2:6:end)=y_change2(:,step-steps_change*2);
        xd(3:6:end)=z_change2(:,step-steps_change*2);
        xd(4:6:end)=x_dot_change2(:,step-steps_change*2);
        xd(5:6:end)=y_dot_change2(:,step-steps_change*2);
        xd(6:6:end)=z_dot_change2(:,step-steps_change*2);
    end
    if step>steps_change*6 && step<steps_change*10
        xd(1:6:end)=x_change3(:,step-steps_change*6);
        xd(2:6:end)=y_change3(:,step-steps_change*6);
        xd(3:6:end)=z_change3(:,step-steps_change*6);
        xd(4:6:end)=x_dot_change3(:,step-steps_change*6);
        xd(5:6:end)=y_dot_change3(:,step-steps_change*6);
        xd(6:6:end)=z_dot_change3(:,step-steps_change*6);
    end
   
    if step==1
        pause(0.5)
    end
    pos=[x(1:6:end),x(2:6:end),x(3:6:end)];
    pos_d=[xd(1:6:end),xd(2:6:end),xd(3:6:end)];
    error_real(step)=mean(vecnorm(pos-pos_d, 2, 2));

    pos_ed=[xe_d(1:6:end),xe_d(2:6:end),xe_d(3:6:end)];
    error_virtual_d(step)=mean(vecnorm(pos_ed-pos_d, 2, 2));
    error_ed(step)=mean(vecnorm(pos-pos_ed, 2, 2));
end


%%

t=0:delta:(iter-1)*delta;
figure(2)
plot(t,error_real,LineWidth=1.5)
hold on
plot(t,error_virtual_d,LineWidth=1.5)
hold off
title('Error')
legend('Error real','Error virtual DEKF' )
xlabel('Time(s)')
grid on

figure(3)
plot(t,error_ed,LineWidth=1.5)
title('Error estimation')
legend( 'DEKF')
xlabel('Time(s)')
grid on


%%
function A=A_linearized(x,k,k2,c1,c2,m,l0,n_points,n_points2)
A=zeros(n_points*n_points2*6,n_points*n_points2*6);
for j=1:1:n_points2
    for i=1:1:n_points
        pos_x=6*n_points*(j-1)+6*(i-1)+1:6*n_points*(j-1)+6*(i-1)+3;
        pos_v=6*n_points*(j-1)+6*(i-1)+4:6*n_points*(j-1)+6*(i-1)+6;
        A(pos_x,pos_v)=eye(3);
        if (i==1 && j==1)||(i==1 && j==n_points2)||(i==n_points && j==1)||(i==n_points && j==n_points2)
            c=c2;
        else
            c=c1;
        end
        A(pos_v,pos_v)=-c/m(i,j)*eye(3);
        aux_same = 0;
        if i<n_points
            aux_same = aux_same - k*eye(3)+ component_A(i,j,i+1,j,x,n_points,k,l0);
            A(pos_v,pos_x+6)= k*eye(3) -component_A(i,j,i+1,j,x,n_points,k,l0);
        end
        if i<n_points-1
            aux_same = aux_same - k*eye(3)+ 2*component_A(i,j,i+2,j,x,n_points,k,l0); 
            A(pos_v,pos_x+6*2)= k*eye(3) - 2*component_A(i,j,i+2,j,x,n_points,k,l0); 
        end
        if i>1
            aux_same = aux_same - k*eye(3) + component_A(i,j,i-1,j,x,n_points,k,l0); 
            A(pos_v,pos_x-6)= k*eye(3) - component_A(i,j,i-1,j,x,n_points,k,l0); 
        end
        if i>2
            aux_same = aux_same - k*eye(3) + 2*component_A(i,j,i-2,j,x,n_points,k,l0);
            A(pos_v,pos_x-6*2)= k*eye(3) - 2*component_A(i,j,i-2,j,x,n_points,k,l0);
        end
        if j<n_points2
            aux_same = aux_same - k*eye(3)+ component_A(i,j,i,j+1,x,n_points,k,l0);
            A(pos_v,pos_x+6*n_points)= k*eye(3) - component_A(i,j,i,j+1,x,n_points,k,l0);
        end
        if j<n_points2-1
            aux_same = aux_same - k*eye(3)+ 2*component_A(i,j,i,j+2,x,n_points,k,l0);
            A(pos_v,pos_x+2*6*n_points)= k*eye(3) - 2*component_A(i,j,i,j+2,x,n_points,k,l0);
        end
        if j>1
            aux_same = aux_same - k*eye(3)+ component_A(i,j,i,j-1,x,n_points,k,l0);
            A(pos_v,pos_x-6*n_points)= k*eye(3) - component_A(i,j,i,j-1,x,n_points,k,l0);
        end
        if j>2
            aux_same = aux_same - k*eye(3)+ 2*component_A(i,j,i,j-2,x,n_points,k,l0);
            A(pos_v,pos_x-2*6*n_points)= k*eye(3) - 2*component_A(i,j,i,j-2,x,n_points,k,l0);
        end
        if i<n_points && j<n_points2
            aux_same = aux_same - k2*eye(3)+ component_A(i,j,i+1,j+1,x,n_points,k2,l0*sqrt(2));
            A(pos_v,pos_x+6*n_points+6)= k2*eye(3) - component_A(i,j,i+1,j+1,x,n_points,k2,l0*sqrt(2));
        end
        if i>1 && j>1
            aux_same = aux_same - k2*eye(3)+ component_A(i,j,i-1,j-1,x,n_points,k2,l0*sqrt(2));
            A(pos_v,pos_x-6*n_points-6)= k2*eye(3) - component_A(i,j,i-1,j-1,x,n_points,k2,l0*sqrt(2));
        end
        if i<n_points && j>1
            aux_same = aux_same - k2*eye(3)+component_A(i,j,i+1,j-1,x,n_points,k2,l0*sqrt(2));
            A(pos_v,pos_x-6*n_points+6)= k2*eye(3) -component_A(i,j,i+1,j-1,x,n_points,k2,l0*sqrt(2));
        end
        if i>1 && j<n_points2
            aux_same = aux_same - k2*eye(3)+ component_A(i,j,i-1,j+1,x,n_points,k2,l0*sqrt(2));
            A(pos_v,pos_x+6*n_points-6)= k2*eye(3) - component_A(i,j,i-1,j+1,x,n_points,k2,l0*sqrt(2));
        end
        A(pos_v,pos_x)= aux_same;
    end
end
end

function A=K_matrix(x,k,k2,c1,c2,l0,n_points,n_points2)
A=zeros(n_points*n_points2*6,n_points*n_points2*6);
for j=1:1:n_points2
    for i=1:1:n_points
        pos_x=6*n_points*(j-1)+6*(i-1)+1:6*n_points*(j-1)+6*(i-1)+3;
        pos_v=6*n_points*(j-1)+6*(i-1)+4:6*n_points*(j-1)+6*(i-1)+6;
        A(pos_x,pos_v)=eye(3);
        if (i==1 && j==1)||(i==1 && j==n_points2)||(i==n_points && j==1)||(i==n_points && j==n_points2)
            c=c2;
        else
            c=c1;
        end
        A(pos_v,pos_v)=-c*eye(3);
        aux_same = 0;
        if i<n_points
            aux_same = aux_same - k*eye(3)+ k*l0/norm_matrix(i,j,i+1,j,x,n_points)*eye(3);

            A(pos_v,pos_x+6)= k*eye(3) - k*l0/norm_matrix(i,j,i+1,j,x,n_points)*eye(3);
        end
        if i<n_points-1
            aux_same = aux_same - k*eye(3) + 2*k*l0/norm_matrix(i,j,i+2,j,x,n_points)*eye(3);
            A(pos_v,pos_x+6*2)= k*eye(3) - 2*k*l0/norm_matrix(i,j,i+2,j,x,n_points)*eye(3);
        end
        if i>1
            aux_same = aux_same - k*eye(3) + k*l0/norm_matrix(i,j,i-1,j,x,n_points)*eye(3);
            A(pos_v,pos_x-6)= k*eye(3) - k*l0/norm_matrix(i,j,i-1,j,x,n_points)*eye(3);
        end
        if i>2
            aux_same = aux_same - k*eye(3) + 2*k*l0/norm_matrix(i,j,i-2,j,x,n_points)*eye(3);
            A(pos_v,pos_x-6*2)= k*eye(3) - 2*k*l0/norm_matrix(i,j,i-2,j,x,n_points)*eye(3);
        end
        if j<n_points2
            aux_same = aux_same - k*eye(3)+ k*l0/norm_matrix(i,j,i,j+1,x,n_points)*eye(3);
            A(pos_v,pos_x+6*n_points)= k*eye(3) - k*l0/norm_matrix(i,j,i,j+1,x,n_points)*eye(3);
        end
        if j<n_points2-1
            aux_same = aux_same - k*eye(3)+ 2*k*l0/norm_matrix(i,j,i,j+2,x,n_points)*eye(3);
            A(pos_v,pos_x+2*6*n_points)= k*eye(3) - 2*k*l0/norm_matrix(i,j,i,j+2,x,n_points)*eye(3);
        end
        if j>1
            aux_same = aux_same - k*eye(3)+ k*l0/norm_matrix(i,j,i,j-1,x,n_points)*eye(3);
            A(pos_v,pos_x-6*n_points)= k*eye(3) - k*l0/norm_matrix(i,j,i,j-1,x,n_points)*eye(3);
        end
        if j>2
            aux_same = aux_same - k*eye(3)+ 2*k*l0/norm_matrix(i,j,i,j-2,x,n_points)*eye(3);
            A(pos_v,pos_x-2*6*n_points)= k*eye(3) - 2*k*l0/norm_matrix(i,j,i,j-2,x,n_points)*eye(3);
        end
        if i<n_points && j<n_points2
            aux_same = aux_same - k2*eye(3)+ k2*l0*sqrt(2)/norm_matrix(i,j,i+1,j+1,x,n_points)*eye(3);
            A(pos_v,pos_x+6*n_points+6)= k2*eye(3) - k2*l0*sqrt(2)/norm_matrix(i,j,i+1,j+1,x,n_points)*eye(3);
        end
        if i>1 && j>1
            aux_same = aux_same - k2*eye(3)+ k2*l0*sqrt(2)/norm_matrix(i,j,i-1,j-1,x,n_points)*eye(3);
            A(pos_v,pos_x-6*n_points-6)= k2*eye(3) - k2*l0*sqrt(2)/norm_matrix(i,j,i-1,j-1,x,n_points)*eye(3);
        end
        if i<n_points && j>1
            aux_same = aux_same - k2*eye(3)+ k2*l0*sqrt(2)/norm_matrix(i,j,i+1,j-1,x,n_points)*eye(3);
            A(pos_v,pos_x-6*n_points+6)= k2*eye(3) - k2*l0*sqrt(2)/norm_matrix(i,j,i+1,j-1,x,n_points)*eye(3);
        end
        if i>1 && j<n_points2
            aux_same = aux_same - k2*eye(3)+ k2*l0*sqrt(2)/norm_matrix(i,j,i-1,j+1,x,n_points)*eye(3);
            A(pos_v,pos_x+6*n_points-6)= k2*eye(3) - k2*l0*sqrt(2)/norm_matrix(i,j,i-1,j+1,x,n_points)*eye(3);
        end
        A(pos_v,pos_x)= aux_same;
    end
end
end

function value= norm_matrix(i1,j1,i2,j2,x,n_points)
    aux=x(6*n_points*(j1-1)+6*i1-5:6*n_points*(j1-1)+6*i1-3);
    aux2=x(6*n_points*(j2-1)+6*i2-5:6*n_points*(j2-1)+6*i2-3);
    value = norm(aux-aux2);
end

function value=x_dot(i1,j1,i2,j2,x,n_points)
    aux=x(6*n_points*(j1-1)+6*i1-5:6*n_points*(j1-1)+6*i1-3);
    aux2=x(6*n_points*(j2-1)+6*i2-5:6*n_points*(j2-1)+6*i2-3);
    value=(aux-aux2)*(aux-aux2)';
end

function value=component_A(i1,j1,i2,j2,x,n_points,k,l0)
    value=k*l0/norm_matrix(i1,j1,i2,j2,x,n_points)*eye(3)-k*l0*x_dot(i1,j1,i2,j2,x,n_points)/norm_matrix(i1,j1,i2,j2,x,n_points)^3;
end


function [x_change0,y_change0,z_change0,x_dot_change0,y_dot_change0,z_dot_change0]=changes(x_i,y_i,z_i,x_g_vector,y_g_vector,z_g_vector,steps_change,delta,factor,factor2,n_points,n_points2)


distance_total_x = x_g_vector-x_i;
distance_total_y = y_g_vector-y_i;
distance_total_z = z_g_vector-z_i;

t_acc=factor*steps_change*delta*factor2;
t_decel=t_acc;
t_const = factor*steps_change*delta - t_acc - t_decel;

v_max_x = distance_total_x / (t_const + t_acc);
v_max_y = distance_total_y / (t_const + t_acc);
v_max_z = distance_total_z / (t_const + t_acc);

x_change0=zeros(n_points2*n_points,int32(factor*steps_change));
y_change0=zeros(n_points2*n_points,int32(factor*steps_change));
z_change0=zeros(n_points2*n_points,int32(factor*steps_change));
x_dot_change0=zeros(n_points2*n_points,int32(factor*steps_change));
y_dot_change0=zeros(n_points2*n_points,int32(factor*steps_change));
z_dot_change0=zeros(n_points2*n_points,int32(factor*steps_change));

for j=1:n_points2*n_points
    for i=1:factor*steps_change
        current_time = i*delta;
        if current_time <= t_acc
            % Acceleration phase (linear increase in velocity)
            if  distance_total_x(j)~=0
                x_change0(j,i) = x_i(j) + 0.5 * (v_max_x(j) / t_acc) * current_time^2 * (x_g_vector(j) - x_i(j)) / distance_total_x(j);
                x_dot_change0(j,i) = (v_max_x(j) / t_acc) * current_time * (x_g_vector(j) - x_i(j)) / distance_total_x(j);
            else
                x_change0(j,i) = x_i(j) ;
            end
            if distance_total_y(j)~=0
                y_change0(j,i) = y_i(j) + 0.5 * (v_max_y(j) / t_acc) * current_time^2 * (y_g_vector(j) - y_i(j)) / distance_total_y(j);
                y_dot_change0(j,i) = (v_max_y(j) / t_acc) * current_time * (y_g_vector(j) - y_i(j)) / distance_total_y(j);
            else
                y_change0(j,i) = y_i(j) ;
            end
            if distance_total_z(j)~=0
                z_change0(j,i) = z_i(j) + 0.5 * (v_max_z(j) / t_acc) * current_time^2 * (z_g_vector(j) - z_i(j)) / distance_total_z(j);
                z_dot_change0(j,i) = (v_max_z(j) / t_acc) * current_time * (z_g_vector(j) - z_i(j)) / distance_total_z(j);
            else
                z_change0(j,i) = z_i(j) ;
            end

        elseif current_time <= (factor*steps_change*delta - t_decel)
            % Constant velocity phase
            if  distance_total_x(j)~=0
                x_change0(j,i) = x_i(j) + 0.5 * v_max_x(j) * t_acc * (x_g_vector(j) - x_i(j)) / distance_total_x(j) ...
                    + v_max_x(j) * (current_time - t_acc) * (x_g_vector(j) - x_i(j)) / distance_total_x(j);
                x_dot_change0(j,i) = v_max_x(j)*(x_g_vector(j) - x_i(j))/ distance_total_x(j);
            else
                x_change0(j,i) = x_i(j) ;
            end
            if  distance_total_y(j)~=0
                y_change0(j,i) = y_i(j) + 0.5 * v_max_y(j) * t_acc * (y_g_vector(j) - y_i(j)) / distance_total_y(j) ...
                    + v_max_y(j) * (current_time - t_acc) * (y_g_vector(j) - y_i(j)) / distance_total_y(j);
                y_dot_change0(j,i) = v_max_y(j)*(y_g_vector(j) - y_i(j))/ distance_total_y(j);
            else
                y_change0(j,i) = y_i(j) ;
            end
            if  distance_total_z(j)~=0
                z_change0(j,i) = z_i(j) + 0.5 * v_max_z(j) * t_acc * (z_g_vector(j) - z_i(j)) / distance_total_z(j) ...
                    + v_max_z(j) * (current_time - t_acc) * (z_g_vector(j) - z_i(j)) / distance_total_z(j);
                z_dot_change0(j,i) = v_max_z(j)*(z_g_vector(j) - z_i(j))/ distance_total_z(j);
            else
                z_change0(j,i) = z_i(j) ;
            end

        else
            % Deceleration phase (linear decrease in velocity)
            t_remaining = factor*steps_change*delta - current_time;
            if  distance_total_x(j)~=0
                x_change0(j,i) = x_g_vector(j)- 0.5 * (v_max_x(j) / t_decel) * t_remaining^2 * (x_g_vector(j) - x_i(j)) / distance_total_x(j);
                x_dot_change0(j,i) = v_max_x(j)* t_remaining / t_decel *(x_g_vector(j) - x_i(j))/ distance_total_x(j);
            else
                x_change0(j,i) = x_i(j) ;
            end
            if  distance_total_y(j)~=0
                y_change0(j,i) = y_g_vector(j)- 0.5 * (v_max_y(j) / t_decel) * t_remaining^2 * (y_g_vector(j) - y_i(j)) / distance_total_y(j);
                y_dot_change0(j,i) = v_max_y(j)* t_remaining / t_decel *(y_g_vector(j) - y_i(j))/ distance_total_y(j);
            else
                y_change0(j,i) = y_i(j) ;
            end
            if  distance_total_z(j)~=0
                z_change0(j,i) = z_g_vector(j)- 0.5 * (v_max_z(j) / t_decel) * t_remaining^2 * (z_g_vector(j) - z_i(j)) / distance_total_z(j);
                z_dot_change0(j,i) = v_max_z(j)* t_remaining / t_decel *(z_g_vector(j) - z_i(j))/ distance_total_z(j);
            else
                z_change0(j,i) = z_i(j) ;
            end

        end
    end
end
    x_change0=[x_change0,repmat(x_g_vector,1,int32(steps_change*(1-factor)))];
    y_change0=[y_change0,repmat(y_g_vector,1,int32(steps_change*(1-factor)))];
    z_change0=[z_change0,repmat(z_g_vector,1,int32(steps_change*(1-factor)))];
    x_dot_change0=[x_dot_change0,repmat(zeros(n_points*n_points2,1),1,int32(steps_change*(1-factor)))];
    y_dot_change0=[y_dot_change0,repmat(zeros(n_points*n_points2,1),1,int32(steps_change*(1-factor)))];
    z_dot_change0=[z_dot_change0,repmat(zeros(n_points*n_points2,1),1,int32(steps_change*(1-factor)))];
end