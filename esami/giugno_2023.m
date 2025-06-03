%es 1 
%4R planar robot inertia matrix ABSOLUTE ANGLES
clear all
n = 4; 
syms q l dc m [n 1] real

syms dq ddq [n 1] real


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
x1 = (dc1-l1)*cos(q1);
y1 = (dc1-l1)*sin(q1);
rc_01 = [x1; y1];

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
T1_rot = (1/2)* w1' * Ic1 *w1;
T1 = simplify(T1_tr + T1_rot);
disp("Kinetic energy")
disp(T1);

% joint 2 
%position of the CoMi in RF0
x2 = l1*cos(q1)+(dc2-l2)*cos(q2);
y2 = l1*sin(q1)+(dc2-l2)*sin(q2);
rc_02 = [x2; y2];

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
T2_rot = (1/2)* w2' * Ic2 *w2;
T2 = simplify(T2_tr + T2_rot);
disp("Kinetic energy")
disp(T2);

% joint 3 
%position of the CoMi in RF0
x3 = l1*cos(q1)+l2*cos(q2)+(dc3-l3)*cos(q3);
y3 = l1*sin(q1)+l2*sin(q2)+(dc3-l3)*sin(q3);
rc_03 = [x3; y3];

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
T3_rot = (1/2)* w3' * Ic3 *w3;
T3 = simplify(T3_tr + T3_rot);
disp("Kinetic energy")
disp(T3);

% joint 4
%position of the CoMi in RF0
x4 = l1*cos(q1)+l2*cos(q2)+l3*cos(q3)+(dc4-l4)*cos(q4);
y4 = l1*sin(q1)+l2*sin(q2)+l3*sin(q3)+(dc4-l4)*sin(q4);
rc_04 = [x4; y4];

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
T4_rot = (1/2)* w4' * Ic4 *w4;
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

%mapping absolute angles to theta DH angles and substituting to the M(q)

B = [1, 0, 0, 0;
    1, 1, 0, 0;
    1, 1, 1, 0;
    1, 1, 1, 1];

theta = inv(B) * q;

disp("theta angles mapping")
disp(theta)


M_theta = B'*M*B;
disp("M(theta)")
disp(M_theta)


%%
%es 3 
%PRR planar robot control

n = 3;

syms q m [n 1] real
syms dc2 l2 real 

%positions of CoMs in RF0
rc_01 = [0; q1];
rc_02 = [(dc2-l2)*cos(q2);
        q1 + (dc2-l2)*sin(q2)];
rc_03 = [l2*cos(q2); q1+l2*sin(q2)];

syms g0 real
% vector gravity acceleration this depends on the drawing
%see how the vector is represented in RF0

g=[0; -g0];%in this case it was opposite of z0-axis

U = cell(1, n);
%position of each CoM in RF0
rc_0 = {rc_01, rc_02, rc_03};
m = {m1, m2, m3};

for i=1:n
    mi = m{i};
    rc_0i = rc_0{i};

    U{i} = simplify(-mi * transpose(g) * rc_0i);

end

for i = 1:n
    disp("Joint")
    disp(i)
    disp(U{i})
end

disp("Total potential energy for the robot")

U_tot = simplify(U{1} + U{2} + U{3});

disp(U_tot)

G=transpose(jacobian(U_tot,q));

disp("Gravity term for the robot")
disp(G)

q = [q1; q2; q3];
dGdq = jacobian(G, q); %max eigenvalue is on the diagonal and it is = alpha

%to have an asymptotically stable system the min eigenvalue of Kp needs to
%be >alpha
disp(dGdq'*dGdq)

%%
%es 4 prima parte
clear all

syms a1 a2 a3 real
syms q1 q2 real
syms dq1 dq2 real

n=2;

q = [q1; q2];
dq = [dq1; dq2];

%inertia matrix
M = [a1 + 2*a2*cos(q2), a3 + a2*cos(q2);
    a3 + a2*cos(q2), a3];

%-------------Cristoffel terms--------------

C = cell(1, n); %Christoffel matrices
c = cell(1, n);

for i=1:n

    Mi = M(:,i);
    qi = q(i);

    Ci = (1/2)*(jacobian(Mi,q) + jacobian(Mi,q)' - diff(M,qi));

    C{i} = Ci;

    ci = dq' * Ci * dq;
    c{i} = ci;

end 

c_v=[c{1}; c{2}];
disp(c_v)

%-------------Gravity term--------------

disp("Potential Energy for each link")

syms g0 dc1 dc2 m1 m2 real
% vector gravity acceleration this depends on the drawing
%see how the vector is represented in RF0

g=[0; -g0];%in this case it was opposite of y0-axis

U = cell(1, n);
%position of each CoM in RF0
rc_01 = [dc1*sin(q1); dc1*cos(q1-pi/2)];
rc_02 = [sin(q1) + dc2*cos(q1+q2); cos(q1-pi/2)+ dc2*sin(q1+q2)];
rc_0 = {rc_01, rc_02};

m = {m1, m2};

for i=1:n
    rc_0i = rc_0{i};
    mi = m{i};

    U{i} = simplify(-mi * transpose(g) * rc_0i);

end

for i = 1:n
    disp("Joint")
    disp(i)
    disp(U{i})
end


syms mp real

disp("Potential energy of the payload")
%position of the payload
pp = [cos(q1)+cos(q1+q2); sin(q1)+sin(q1+q2)];
U3 = simplify(-mp * transpose(g) * pp);

disp(U3)
disp("Total potential energy for the robot")

U_tot = simplify(U{1} + U{2} +U3);

disp(U_tot)

G=transpose(jacobian(U_tot,q));

disp("Gravity term")
disp(G)


