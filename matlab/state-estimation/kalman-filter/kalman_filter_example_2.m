% KALMAN_FILTER - updates a system state vector estimate based upon an
%                 observation, using a discrete Kalman filter.
%
% File inspired from "Learning the Kalman Filter" by Michael C. Kleder
% http://ch.mathworks.com/matlabcentral/fileexchange/5377-learning-the-kalman-filter
%
% FORM OF EQUATIONS USED IN THIS FILE
%
% STATE UPDATE
%
% x = Ax + Bu + w  
%
% MEASUREMENT UPDATE
% z = Hx + v       

% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is gaussian noise with covariance R
%
% VECTOR VARIABLES:
%
% s.x = state vector estimate. In the input struct, this is the
%       "a priori" state estimate. In the output struct,
%       this is the "a posteriori" state estimate.
% s.z = observation vector
% s.u = input control vector, optional (defaults to zero).
%
% MATRIX VARIABLES:
%
% state_.A = state transition matrix (defaults to identity).
% state_.P = covariance of the state vector estimate. In the input struct,
%       this is "a priori," and in the output it is "a posteriori."
%       (required unless autoinitializing as described below).
% state_.B = input matrix, optional (defaults to zero).
% state_.Q = process noise covariance (defaults to zero).
% state_.R = measurement noise covariance (required).
% state_.H = observation matrix (defaults to identity).
%
% NORMAL OPERATION:
%
% (1) define all state definition fields: A,B,H,Q,R
% (2) define intial state estimate: x,P
% (3) obtain observation and control vectors: z,u
% (4) call the filter to obtain updated state estimate: x,P
% (5) return to step (3) and repeat
%
% INITIALIZATION:
%
% If an initial state estimate is unavailable, then it can be obtained
% from the first observation as follows, provided that there are the
% same number of observable variables as state variables. This "auto-
% intitialization" is done automatically if state_.x is absent or NaN.
%
% x = inv(H)*z
% P = inv(H)*R*inv(H')
%
% This is mathematically equivalent to setting the initial state estimate
% covariance to infinity.
%

%% Problem Description %%
% This specific example considers a simple bi-wheeled robot model
% (each wheel can rotate independently to the other, and they both have
% the same radius while the robot Center-of-Mass lies at the center of
% their inter-axial distance)
% 
% The considered system states are:
% omega_1 : rotational speed of wheel 1 [rad/s]
% omega_2 : rotational speed of wheel 2 [rad/s]
% v_b   : velocity of the robot CoM [m/s] 
% theta : orientation angle of the robot [rad]
%
% The assembled state vector is:
% x = [omega_1 , omega_2 , v_b , theta]^T
% 
% Linear velocity of wheels i=1,2 : v_i = r * omega_i (where r the common radius)
% Linear velocity of robot : v_b = (v_1 + v_2)/2 
% Orientation angular speed of robot : omega_b = (v1 - v2)/d (where d the wheel inter-axle distance)
% Orientation angle of  robot : theta_{k} = theta_{k-1} + omega_b * Ts (where Ts the discrete sampling time) 
%
% Two measurement case setups are implemented:
% Case1: We only get measuements for states omega_i (rotational speed encoders on wheels)
% Case2: We get the measuremetns of Case1 as well as measurements for the states v_b (odometry) and theta (magnetometer)
%
%% Setup:
clear all;
clc;
close all;

% System (Model) Sampling Time
Ts = 0.1;
% Axle-to-axle Wheel distance
d  = 1.0;
% Wheel radius
r  = 0.25;

% Process uncertainty - noise stdev 
Q_rho = deg2rad(1); %[rad/s]
Q_v = 0.1; %[m/s]
Q_theta = deg2rad(1); %[rad]

% Measurement uncertainty - noise stdev 
R_rho = deg2rad(5); %wheel rotational speed encoders [rad/s]
R_v = 0.1; %robot velocity odometry [m/s] (only Case2)
R_theta = deg2rad(5); % robot orientation magnetometer [rad] (only Case2)

% Measurement case {1,2}
R_case = 2; 
% Note: try setting wheel_max=Inf for both measurement cases

%% Problem Definition:
clear state_
% Initial state
state_.x = [0;0;0;0];
% State Update Matrix (ss acquired from Problem Description)
state_.A = [1 0 r*Ts/d -r*Ts/d;
            0 0 r/2    r/2;
            0 0 1      0; 
            0 0 0      1];

% Define a process noise (stdev) % variance, hence stdev^2
state_.Q = [Q_theta^2 0     0       0;
            0         Q_v^2 0       0;
            0         0     Q_rho^2 0;
            0         0     0       Q_rho^2];
% Scale process noise model (if necessary)
state_.Q = state_.Q * 1;

% Case1: We only have measurements from rotational speed encoders (wheel rotational speed)
if (R_case==1)
% Define the measurement matrix:
state_.H = [0 0 1 0;  %encoder1 measures wheel speed (omega_1)
            0 0 0 1]; %encoder2 measures wheel speed (omega_2)
% Define a measurement error (stdev) % variance, hence stdev^2
state_.R = [R_rho^2 0;
            0       R_rho^2];
end
% Case2: We have measurements from rotational speed encoders (wheel rotational speed)
% & robot vecocity (odometry) & robot orientation (magnetometer) 
if (R_case==2)
% Define the measurement matrix:
state_.H = [0 0 1 0;  %encoder1 measures wheel speed (omega_1)
            0 0 0 1;  %encoder2 measures wheel speed (omega_2)
            1 0 0 0;  %magnetometer measures orientation (theta)
            0 1 0 0]; %odometry measures velocity (v_b)
% Define a measurement error (stdev) % variance, hence stdev^2
state_.R = [R_rho^2 0       0         0;
            0       R_rho^2 0         0;
            0       0       R_theta^2 0;
            0       0       0         R_v^2];  
end

% Scale sensor noise model (if necessary)
state_.R = state_.R * 1;

% Do not define any system input (control) functions:
state_.B = [0;0;0;0];
state_.u = 0;
%Large initial MSE covariace (uncertainty)
state_.P = 1000*ones(4); 

%% Generate signals:
% Run from 0[s] to tend[s] (final time)
tend = 20;
tvec = [0:Ts:tend];

% Real (True,Noise-Less) state profile (wheels 1,2 rotational speeds) 
wheel_1_slope = 1; %increase/decrease slope [(m/s)/s]
wheel_2_slope = 0.25 %increase/decrease slope [(m/s)/s]
wheel_max = 1.0/r; %[m/s]->[rad/s] max speed
%wheel_max = Inf; %no max

% Create slope profiles for states
x3 = tvec * wheel_1_slope;
x4 = tvec * wheel_2_slope;
% Constrain slope profiles for states
x3 = max(x3,-wheel_max*ones(size(x3)));
x3 = min(x3,wheel_max*ones(size(x3)));
x4 = max(x4,-wheel_max*ones(size(x4)));
x4 = min(x4,wheel_max*ones(size(x4)));
%Simulate remainder of system dynamics to obtain Real (True,Noise-less) state profile
x = [0;0;x3(1);x4(1)];
u = 0;
x1 = [];
x2 = [];
for i=1:length(tvec)
    x_new = state_.A*x + state_.B*u;
    x1 = [x1 x_new(1)];
    x2 = [x2 x_new(2)];
    x = [x_new(1);x_new(2);x3(i);x4(i)];   
end

% Artificially create measurements by superimposing measuement noise (possible scaled to more/less than sensor model)
z1 = x3 + randn([1 length(x3)]) * R_rho * 1;
z2 = x4 + randn([1 length(x4)]) * R_rho * 1;
z3 = x1 + randn([1 length(x1)]) * R_theta * 1;
z4 = x2 + randn([1 length(x2)]) * R_v * 1;

figure;
hold on;
plot(tvec,x1,'Color',[0 0 0.5],'LineWidth',1); %theta
plot(tvec,x2,'Color',[0 0.5 0],'LineWidth',1); %v_b 
plot(tvec,x3,'Color',[0.75 0 0],'LineWidth',1); %omega_1
plot(tvec,x4,'Color',[0.75 0 0.75],'LineWidth',1); %omega_2
plot(tvec,z1,'Color',[1 0 0],'LineWidth',2.5); %omega_1_meas
plot(tvec,z2,'Color',[1 0 1],'LineWidth',2.5); %omega_2_meas
if (R_case==1) %Measurements of Case1
  legend('\theta','v_b','\omega_1','\omega_2','Z\omega_1','Z\omega_2');
end
if (R_case==2) %Extra Measurements of Case2
  plot(tvec,z3,'Color',[0 0 1],'LineWidth',2.5); %theta_meas
  plot(tvec,z4,'Color',[0 1 0],'LineWidth',2.5); %v_b_meas
  legend('\theta','v_b','\omega_1','\omega_2','Z\omega_1','Z\omega_2','Z\theta','Zv_b');
end

%% Filter operation:
for i=1:length(tvec)-1
   if (R_case==1) %Case1: base measurements from wheel speed encoders
     state_(end).z = [z1(i);z2(i)];
   end
   if (R_case==2) %Case2: 2 extra measurements from magnetometer and odometry
     state_(end).z = [z1(i);z2(i);z3(i);z4(i)];
   end
   state_(end+1)=kalman_filter(state_(end)); % perform a Kalman filter iteration
end

figure
hold on
grid on

% Plot true (noise-less) state data:
plot(tvec,x1,'Color',[0 0 1],'LineWidth',2.0); %theta
plot(tvec,x2,'Color',[0 0.5 0],'LineWidth',2.0); %v_b
plot(tvec,x3,'Color',[0.75 0   0],'LineWidth',2.0); %omega_1
plot(tvec,x4,'Color',[0.75 0   0.75],'LineWidth',2.0); %omega_2

% Plot a-posteriori state estimates:
state_x = [];
for i=1:length(state_)
    state_x = [state_x state_(i).x];
end
plot(tvec,state_x(1,:),'Color',[0 0 0.75],'LineWidth',0.5);  %theta_estim
plot(tvec,state_x(2,:),'Color',[0 0.75 0],'LineWidth',0.5); %v_b_estim
plot(tvec,state_x(3,:),'Color',[0.75 0   0],'LineWidth',0.5); %omega_1_estim
plot(tvec,state_x(4,:),'Color',[0.75 0   0.75],'LineWidth',0.5); %omega_2_estim

% Plot measurement data:
state_z = [];
for i=1:length(state_)
    state_z = [state_z state_(i).z];
end
plot(tvec,state_z(1,:),'Color',[1 0 0],'LineStyle','none','Marker','.'); %omega_1_meas
plot(tvec,state_z(2,:),'Color',[1 0 1],'LineStyle','none','Marker','.'); %omega_2_meas
if (R_case==1) %Measurements of Case1
  legend('\theta','v_b','\omega_1','\omega_2','KF\theta','KFv_b','KF\omega_1','KF\omega_2','Z\omega_1','Z\omega_2');
end
if (R_case==2) %Extra Measurements of Case2
  plot(tvec,state_z(3,:),'Color',[0 0 1],'LineStyle','none','Marker','.'); %theta_meas
  plot(tvec,state_z(4,:),'Color',[0 1 0],'LineStyle','none','Marker','.'); %v_b_meas
  legend('\theta','v_b','\omega_1','\omega_2','KF\theta','KFv_b','KF\omega_1','KF\omega_2','Z\omega_1','Z\omega_2','Z\theta','Zv_b');
end
hold off

%% Extra plot of 2D-position in x/y axes
% Information is not part of system states, we derive it from integrating theta-rotated velocity v_b 

figure;
hold on;
axis equal;
x_true = [0];
y_true = [0];
x_estim = [0];
y_estim = [0];
for i=2:length(tvec)
    x_true = [x_true x_true(end)+x2(i)*cos(x1(i))];
    y_true = [y_true y_true(end)+x2(i)*sin(x1(i))];
    x_estim = [x_estim x_estim(end)+state_x(2,i)*cos(state_x(1,i))];
    y_estim = [y_estim y_estim(end)+state_x(2,i)*sin(state_x(1,i))];
end
plot(x_true,y_true,'b');
plot(x_estim,y_estim,'r');
xlabel('x'); ylabel('y');
legend('True','KF-estimated');


