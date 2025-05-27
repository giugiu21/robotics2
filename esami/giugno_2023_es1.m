%Giugno 2023
%es1 ABSOLUTE ANGLES
clear all 

n=4;
syms q m l dc [n 1] real
syms dq ddq [n 1] real
syms theta [n 1] real

%position of the CoMi in RF0

pcm_1=  [(dc1-l1)*cos(q1); (dc1-l1)*sin(q1)];
pcm_2=  [l1*cos(q1)+(dc2-l2)*cos(q2); l1*sin(q1)+(dc2-l2)*sin(q2)];
pcm_3=  [l1*cos(q1)+l2*cos(q2)+(dc3-l3)*cos(q3); l1*sin(q1)+l2*sin(q2)+(dc3-l3)*sin(q3)];
pcm_4=  [l1*cos(q1)+l2*cos(q2)+l3*cos(q3)+(dc4-l4)*cos(q4); l1*sin(q1)+l2*sin(q2)+l3*sin(q3)+(dc4-l4)*sin(q4)];


B = [1, 0, 0, 0;
    1, 1, 0, 0;
    1, 1, 1, 0;
    1, 1, 1, 1];


theta= inv(B)*[q1; q2; q3; q4];
disp(theta)


%%

%DH parameters
alpha = [0, 0, 0, 0]; %IMPORTANT to make pi symbolic
a = [l1, l2, l3, l4];
d = [0, 0, 0, 0];


A = cell(1, n);
A_tot = eye(4);


for i = 1:n
    alphai = alpha(i);
    ai = a(i);
    di = d(i);
    thetai = theta(i);

    A{i} =  [cos(thetai), -sin(thetai)*cos(alphai), sin(thetai)*sin(alphai), ai*cos(thetai);
            sin(thetai), cos(thetai)*cos(alphai), -cos(thetai)*sin(alphai), ai*sin(thetai);
            0, sin(alphai), cos(alphai), di;
            0, 0, 0, 1];

    %complete DH matrix for the robot

    A_tot = simplify(A_tot * A{i}); 

end 

%extracting from the homogeneous matrices the rotation matrices and
%positions for each joint

R1 = A{1}(1:3,1:3);
p1 = A{1}(1:3,4);

R2 = A{2}(1:3,1:3); 
p2 = A{2}(1:3,4);

R3 = A{3}(1:3,1:3); 
p3 = A{3}(1:3,4);

R4 = A{4}(1:3,1:3); 
p4 = A{4}(1:3,4);


disp("Rotation matrices")
disp(R1)
disp(R2)
disp(R3)
disp(R4)

disp("Origin's position in RFi-1")
pause

disp(p1)
disp(p2)
disp(p3)
disp(p4)

%position vectors from origin i-1 to origin i in frame i
% ip = i-1Ri' * i-1p !!! i am doing this computation here
% i-1p = i-1Ri * ip

r01 = simplify(R1'*p1);
r12 = simplify(R2'*p2);
r23 = simplify(R3'*p3);
r34 = simplify(R4'*p4);

disp("Position vectors from origin i-1 to origin i in frame i")
pause 


disp(r01)
disp(r12)
disp(r23)
disp(r34)


%!!!CoM position vector is the vector from the origin of frame i to the CMi
% expressed in frame i
%in this case we assume it to be in the kinematic axis of each link
%if we define di as the distance between origin of frame i to CMi

%Understanding the link axis from DH parameters:
%if alpha_i = 0 the link extends on x-axis of frame i
%if alpha_i = pi/2 the link extends on y_axis of frame i
%if alpha_i = pi the link extends on -x-axis of frame i

% to find the link lengths look at the ai parameter:
%in this case a2 and a3 are our link lengths 

%we define d1 d2 d3 as the distance between origin i-1 and position of CoMi

rc_1 = [dc1 - l1; 0; 0];

rc_2= [dc2 - l2; 0; 0];

rc_3= [dc3 - l3; 0; 0];

rc_4= [dc4 - l4; 0; 0];


disp("Position vectors in frame i for CoMi for each link")
disp(rc_1)
disp(rc_2)
disp(rc_3)
disp(rc_4)

disp("CoM position vectors in RF0")
pause


% ip = i-1Ri' * i-1p 
% i-1p = i-1Ri * ip !!i am doing this computation here

rc_01 = simplify((R1*rc_1) + p1);
disp(rc_01)
rc_02 = simplify((R1* (R2*rc_2)) + R1*p2 + p1);
disp(rc_02)
rc_03 = simplify((R1*(R2*(R3*rc_3))) + R1*(R2*p3) + R1*p2 + p1);
disp(rc_03)
rc_04 = simplify((R1*(R2*(R3*(R4*rc_4)))) + R1*(R2*(R3*p4)) + R1*(R2*p3) + R1*p2 + p1);
disp(rc_04)

disp("MOVING FRAMES ALGORITHM")
pause

%-------MOVING FRAMES ALGORITHM-------


%Defining the sigma as: 0 if revolute joint, 1 if prismatic
s1 = 0;
s2 = 0;
s3 = 0;
s4 = 0;


%if the robot were planar  
syms Ic1 Ic2 Ic3 Ic4 real

%the robot is polar so
%syms Ic_xx_1 Ic_yy_1 Ic_zz_1 Ic_xx_2 Ic_yy_2 Ic_zz_2 Ic_xx_3 Ic_yy_3 Ic_zz_3
%Ic1 = diag([Ic_xx_1, Ic_yy_1, Ic_zz_1]);
%Ic2 = diag([Ic_xx_2, Ic_yy_2, Ic_zz_2]);
%Ic3 = diag([Ic_xx_3, Ic_yy_3, Ic_zz_3]);


T = cell(1, n);
%w = cell(1, n);
v = cell(1, n);
vc = cell(1, n);

%set up to change if we have more joints
s={s1, s2, s3, s4};
dq = {dq1, dq2, dq3, dq4};

R = {R1, R2, R3, R4};
r = {r01, r12, r23, r34};
rc = {rc_1, rc_2, rc_3, rc_4};
m = {m1, m2, m3, m4};
Ic = {Ic1, Ic2, Ic3, Ic4};


%Absolute angles
w1 =  [0; 0; dq1];
w2 =  [0; 0; dq2];
w3 =  [0; 0; dq3];
w4 =  [0; 0; dq4];

w = {w1, w2, w3, w4};


for i = 1:n
    Ri = R{i};
    si = s{i};
    dqi = dq{i};
    ri = r{i};
    rci  = rc{i};
    mi = m{i};
    Ici = Ic{i};
    wi  = w{i};

    %initialization if first frame = 0
    if i == 1
        %wi_1 = 0;
        vi_1 = 0;
    else
        %wi_1 = w{i-1};
        vi_1 = v(i-1);
    end

    %w{i} = simplify(transpose(Ri)*(wi_1+(1-si)*dqi*[0; 0; 1]));

    %wi = w{i};

    v{i} = simplify(transpose(Ri)*(vi_1+(si)*dqi*[0; 0; 1]) + cross(wi,ri));

    vi = v{i};

    vc{i} = simplify(vi + cross(wi,rci));

    vci = vc{i};

    T{i} = simplify((1/2)*mi*transpose(vci)*vci + (1/2)*transpose(wi)*Ici*wi);
    
end 


for i = 1:n
    disp("Joint")
    disp(i)
    disp("Angular velocity:")
    disp(w{i})
    disp("Linear velocity:")
    disp(v{i})
    disp("CoM linear velocity:")
    disp(vc{i})
    disp("Kinetic energy:")
    disp(T{i})
end

disp("Total kinetic energy")
pause


%Total kinetic Energy of the Robot 
T_tot = T{1} +T{2} +T{3} + T{4};
T_tot=simplify(T_tot);
T_tot=collect(T_tot,dq1^2);
T_tot=collect(T_tot,dq2^2);
T_tot=collect(T_tot,dq3^3);
T_tot=collect(T_tot,dq4^4);
disp(T_tot)


disp("Inertia matrix for the robot")
pause


%-------------Inertia Matrix--------------
dq = [dq1; dq2; dq3; dq4];

M=simplify(hessian(T_tot,dq));
disp(M)
%%

%Defining B for the absolute angles

M_theta = B'*M*B

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

