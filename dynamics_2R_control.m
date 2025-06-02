%dynamic model 2R
clear all 

n=2;
syms q [n 1] real
syms m l real
syms dq ddq [n 1] real


%if the robot were planar  
syms Ic1 Ic2 real


%the robot is polar so
%syms Ic_xx_1 Ic_yy_1 Ic_zz_1 Ic_xx_2 Ic_yy_2 Ic_zz_2 Ic_xx_3 Ic_yy_3 Ic_zz_3
%Ic1 = diag([Ic_xx_1, Ic_yy_1, Ic_zz_1]);
%Ic2 = diag([Ic_xx_2, Ic_yy_2, Ic_zz_2]);
%Ic3 = diag([Ic_xx_3, Ic_yy_3, Ic_zz_3]);
 
q = [q1; q2];

%In this case all joints are revolute so kinetic energy is composed of
%translational and rotational.

% joint 1 
%position of the CoMi in RF0
x1 = (l/2)*cos(q1);
y1 = (l/2)*sin(q1);
rc_01 = [x1; y1];

%velocity of CoMi
vx1 = diff(x1,q1)*dq1 + diff(x1,q2)*dq2;
vy1 = diff(y1,q1)*dq1 + diff(y1,q2)*dq2;
vc1 = [vx1; vy1];
disp("Velocity CoM")
disp(vc1);

%angular velocity
w1 = [0; 0; dq1]; 

disp("Joint 1")
T1_tr = (1/2)* m * vc1'*vc1;
T1_rot = (1/2)* w1' * Ic1 *w1;
T1 = simplify(T1_tr + T1_rot);
disp("Kinetic energy")
disp(T1);

% joint 2 
%position of the CoMi in RF0
x2 = l*cos(q1)+(l/2)*cos(q1+q2);
y2 = l*sin(q1)+(l/2)*sin(q1+q2);
rc_02 = [x2; y2];

%velocity of CoMi
vx2 = diff(x2,q1)*dq1 + diff(x2,q2)*dq2;
vy2 = diff(y2,q1)*dq1 + diff(y2,q2)*dq2;
vc2 = [vx2; vy2];
disp("Velocity CoM")
disp(vc2);

%angular velocity
w2 = [0; 0; dq1+dq2]; 

disp("Joint 2")
T2_tr = (1/2)* m * vc2'*vc2;
T2_rot = (1/2)* w2' * Ic2 *w2;
T2 = simplify(T2_tr + T2_rot);
disp("Kinetic energy")
disp(T2);

%%

% Total Energy
T = T1+T2;
disp("Total kinetic energy")
disp(T)
disp("Inertia matrix for the robot")

%-------------Inertia Matrix--------------

q_dot = [dq1; dq2];

M=simplify(hessian(T,q_dot));
disp(M)


disp("Robot Coriolis Vector")

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

%%
%-------------Gravity term--------------

disp("Potential Energy for each link")

syms g0
% vector gravity acceleration this depends on the drawing
%see how the vector is represented in RF0

g=[g0; 0];%in this case it was opposite of z0-axis

U = cell(1, n);
%position of each CoM in RF0
rc_0 = {rc_01, rc_02};

for i=1:n
    rc_0i = rc_0{i};

    U{i} = simplify(-m * transpose(g) * rc_0i);

end

for i = 1:n
    disp("Joint")
    disp(i)
    disp(U{i})
end

disp("Total potential energy for the robot")
pause

U_tot = simplify(U{1} + U{2});

disp(U_tot)

G=transpose(jacobian(U_tot,q));

disp("Complete dynamic equation")

%-------------Final Dynamic Model--------------
syms u1 u2 real

disp(M)
disp("*")
disp(ddq)
disp("+")
disp(c_v)
disp("+")
disp(G)
disp("=")
disp([u1; u2])

%%

%position of the e-e

pe = [l*(cos(q1)+cos(q1+q2));
      l*(sin(q1)+sin(q1+q2))];

%jacobian for the robot
J = jacobian(pe, q);
disp("Jacobian for the robot")
disp(J)

%time derivative of the jacobian
%%

J_dot = sym(zeros(size(J))); % stessa dimensione di J

for i = 1:length(q)
    J_dot = J_dot + diff(J, q(i)) * dq(i);
end

J_dot = simplify(J_dot);

disp("Derivative w.r.t time of the Jacobian")
disp(J_dot)

