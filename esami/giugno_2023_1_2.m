%Giugno 2023
%es1 ABSOLUTE ANGLES 4R
clear all 

n=4;
syms q m l dc [n 1] real
syms dq ddq [n 1] real
syms theta [n 1] real


%mapping absolute angles
% thetai = qi - qi-1
B = [1, 0, 0, 0;
    1, 1, 0, 0;
    1, 1, 1, 0;
    1, 1, 1, 1];

theta= inv(B)*[q1; q2; q3; q4];
disp(theta)
%%
% ----------Dynamic model-----------
%if the robot were planar  
syms Ic1 Ic2 Ic3 Ic4 real

%the robot is polar so
%syms Ic_xx_1 Ic_yy_1 Ic_zz_1 Ic_xx_2 Ic_yy_2 Ic_zz_2 Ic_xx_3 Ic_yy_3 Ic_zz_3
%Ic1 = diag([Ic_xx_1, Ic_yy_1, Ic_zz_1]);
%Ic2 = diag([Ic_xx_2, Ic_yy_2, Ic_zz_2]);
%Ic3 = diag([Ic_xx_3, Ic_yy_3, Ic_zz_3]);
 
q = [q1; q2; q3; q4];

%In this case all joints are revolute so kinetic energy is composed of
%translational and rotational.

% joint 1 
%position of the CoMi in RF0
x1 = dc1*cos(q1);
y1 = dc1*sin(q1);

%velocity of CoMi
vx1 = diff(x1,q1)*dq1 + diff(x1,q2)*dq2 + diff(x1,q3)*dq3 + diff(x1,q4)*dq4;
vy1 = diff(y1,q1)*dq1 + diff(y1,q2)*dq2 + diff(y1,q3)*dq3 + diff(y1,q4)*dq4;
vc1 = [vx1; vy1];
disp("Velocity CoM")
disp(vc1);

%angular velocity
w1 = [0; 0; dq1]; %absolute angles

disp("Joint 1")
T1_tr = (1/2)* m1 * vc1'*vc1;
T1_rot = (1/2)* w1'* Ic1 * w1;
T1 = simplify(T1_tr + T1_rot);
disp("Kinetic energy")
disp(T1);

% joint 2 
%position of the CoMi in RF0
x2 = l1*cos(q1)+dc2*cos(q2);
y2 = l1*sin(q1)+dc2*sin(q2);

%velocity of CoMi
vx2 = diff(x2,q1)*dq1 + diff(x2,q2)*dq2 + diff(x2,q3)*dq3 + diff(x2,q4)*dq4;
vy2 = diff(y2,q1)*dq1 + diff(y2,q2)*dq2 + diff(y2,q3)*dq3 + diff(y2,q4)*dq4;
vc2 = [vx2; vy2];
disp("Velocity CoM")
disp(vc2);

%angular velocity
w2 = [0; 0; dq2]; %absolute angles

disp("Joint 2")
T2_tr = (1/2)* m2 * vc2'*vc2;
T2_rot = (1/2)* w2'* Ic2 * w2;
T2 = simplify(T2_tr + T2_rot);
disp("Kinetic energy")
disp(T2);


% joint 3
%position of the CoMi in RF0
x3 = l1*cos(q1)+l2*cos(q2)+dc3*cos(q3);
y3 = l1*sin(q1)+l2*sin(q2)+dc3*sin(q3);

%velocity of CoMi
vx3 = diff(x3,q1)*dq1 + diff(x3,q2)*dq2 + diff(x3,q3)*dq3 + diff(x3,q4)*dq4;
vy3 = diff(y3,q1)*dq1 + diff(y3,q2)*dq2 + diff(y3,q3)*dq3 + diff(y3,q4)*dq4;
vc3 = [vx3; vy3];
disp("Velocity CoM")
disp(vc3);

%angular velocity
w3 = [0; 0; dq3]; %absolute angles

disp("Joint 3")
T3_tr = (1/2)* m3 * vc3'*vc3;
T3_rot = (1/2)* w3'* Ic3 * w3;
T3 = simplify(T3_tr + T3_rot);
disp("Kinetic energy")
disp(T3);


% joint 4
%position of the CoMi in RF0
x4 = l1*cos(q1)+l2*cos(q2)+l3*cos(q3)+dc4*cos(q4);
y4 = l1*sin(q1)+l2*sin(q2)+l3*sin(q3)+dc4*sin(q4);

%velocity of CoMi
vx4 = diff(x4,q1)*dq1 + diff(x4,q2)*dq2 + diff(x4,q3)*dq3 + diff(x4,q4)*dq4;
vy4 = diff(y4,q1)*dq1 + diff(y4,q2)*dq2 + diff(y4,q3)*dq3 + diff(y4,q4)*dq4;
vc4 = [vx4; vy4];
disp("Velocity CoM")
disp(vc4);

%angular velocity
w4 = [0; 0; dq4]; %absolute angles

disp("Joint 4")
T4_tr = (1/2)* m4 * vc4'*vc4;
T4_rot = (1/2)* w4'* Ic4 * w4;
T4 = simplify(T4_tr + T4_rot);
disp("Kinetic energy")
disp(T4);


% Total Energy
T = T1+T2+T3+T4;
disp("Total kinetic energy")
disp(T)

disp("Inertia matrix for the robot")

%-------------Inertia Matrix--------------

q_dot = [dq1; dq2; dq3; dq4];

M=simplify(hessian(T,q_dot));
disp(M)

%%

%Defining B for the absolute angles

M_theta = B'*M*B;
disp(M_theta)

%%
%es 2

pe = [cos(q1)+cos(q2)+cos(q3)+cos(q4);
      sin(q1)+sin(q2)+sin(q3)+sin(q4)];
p2 = [cos(q1)+cos(q2);
      sin(q1)+sin(q2)];

%jacobians for the tasks
J1 = jacobian(pe, [q1; q2; q3; q4]);
J2 = jacobian(p2, [q1; q2; q3; q4]);
disp("Jacobians for the tasks and the rank")
disp(J1)
disp(rank(J1))
disp(J2)
disp(rank(J2))

q0 = [0, pi/6, -pi/3, -pi/3];

%jacobian evaluated at q0

J1_q0 = eval(subs(J1, [q1, q2, q3, q4], q0));
J2_q0 = eval(subs(J2, [q1, q2, q3, q4], q0));
disp("Jacobians evaluated at q0")
disp(J1_q0)
disp(J2_q0)

% a) executing end effector task minimizing the norm

ve = [0.4330; -0.75]; %task velocity e-e
vt = [-0.5; 0.8660]; %task velocity tip 2 link

q_dot_a = pinv(J1_q0)*ve;
disp("A) joint velocities")
disp(q_dot_a)
disp("error task 1")
e_a_1 = ve - J1_q0*q_dot_a;
disp(e_a_1)
disp("error task 2")
e_a_2 = vt - J2_q0*q_dot_a;
disp(e_a_2)
%%
% b) executing joint 2 tip task minimizing the norm


q_dot_b = pinv(J2_q0)*vt;
disp("B)")
disp(q_dot_b)
disp("error task 1")
e_b_1 = ve - J1_q0*q_dot_b;
disp(e_b_1)
disp("error task 2")
e_b_2 = vt - J2_q0*q_dot_b;
disp(e_b_2)

%%
% c) executing task simoultaneously TASK AUGMENTATION

Ja = [J1_q0;
      J2_q0];

ra = [ve;
      vt];
q_dot_c  = pinv(Ja)*ra;

disp("C)")
disp(q_dot_c)
disp("error task 1")
e_c_1 = ve - J1_q0*q_dot_c;
disp(e_c_1)
disp("error task 2")
e_c_2 = vt - J2_q0*q_dot_c;
disp(e_c_2)

%%
% d) priority task e-e

%initialization
J = {J1_q0, J2_q0};
r_dot = {ve, vt};

n = 2; %number of tasks

q_dot_d = zeros(4, 1); %4 is the number of joints
P = eye(4);

for i = 1:n
    Ji = J{i};
    r_dot_i = r_dot{i};

    Ji_proj = Ji * P;
    Ji_pinv = pinv(Ji_proj, 0.001);  % Compute the DAMPED pseudoinverse

    q_dot_d = q_dot_d + Ji_pinv * (r_dot_i - Ji * q_dot_d);  % Task correction
    P = P - Ji_pinv * Ji_proj;                         % Update nullspace
end

%We use a damped pseudoinverse when we get large velocities so that we can
%handle numerical instability. In this case we were probably near a
%singularity that's why we were gettinhg really high values 1.0e+15

disp("D)")
disp(q_dot_d)
disp("error task 1")
e_d_1 = ve - J1_q0*q_dot_d;
disp(e_d_1)
disp("error task 2")
e_d_2 = vt - J2_q0*q_dot_d;
disp(e_d_2)

%%
% e) priority task tip 2 link
%initialization
J = {J2_q0, J1_q0};
r_dot = {vt, ve};

n = 2; %number of tasks

q_dot_e = zeros(4, 1); %4 is the number of joints
P = eye(4);

for i = 1:n
    Ji = J{i};
    r_dot_i = r_dot{i};

    Ji_proj = Ji * P;
    Ji_pinv = pinv(Ji_proj, 0.001);  % Compute the DAMPED pseudoinverse

    q_dot_e = q_dot_e + Ji_pinv * (r_dot_i - Ji * q_dot_e);  % Task correction
    P = P - Ji_pinv * Ji_proj;                         % Update nullspace
end

%We use a damped pseudoinverse when we get large velocities so that we can
%handle numerical instability. In this case we were probably near a
%singularity that's why we were gettinhg really high values 1.0e+15

disp("D)")
disp(q_dot_e)
disp("error task 1")
e_e_1 = ve - J1_q0*q_dot_e;
disp(e_e_1)
disp("error task 2")
e_e_2 = vt - J2_q0*q_dot_e;
disp(e_e_2)

