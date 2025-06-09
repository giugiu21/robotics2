%es 5 PPR robot planar dynamic model
clear all 

n=3;

syms q m dq ddq [n 1] real
syms l  dc3 real

q = [q1; q2; q3];


%if the robot were planar  
syms Ic3 real

% joint 1 
%position of the CoMi in RF0
x1 = q1;
y1 = 0;
rc_01 = [x1; y1];

%velocity of CoMi
vx1 = diff(x1,q1)*dq1 + diff(x1,q2)*dq2 + diff(x1,q3)*dq3;
vy1 = diff(y1,q1)*dq1 + diff(y1,q2)*dq2 + diff(y1,q3)*dq3;
vc1 = [vx1; vy1];
disp("Velocity CoM")
disp(vc1);

% no angular velocity since it is prismatic angular velocity

disp("Joint 1")
T1_tr = (1/2)* m1 * vc1'*vc1;
T1 = simplify(T1_tr);
disp("Kinetic energy")
disp(T1);

% joint 2 
%position of the CoMi in RF0
x2 = q1;
y2 = q2;
rc_02 = [x2; y2];

%velocity of CoMi
vx2 = diff(x2,q1)*dq1 + diff(x2,q2)*dq2 + diff(x2,q3)*dq3;
vy2 = diff(y2,q1)*dq1 + diff(y2,q2)*dq2 + diff(y2,q3)*dq3;
vc2 = [vx2; vy2];
disp("Velocity CoM")
disp(vc2);

% no angular velocity since it is prismatic angular velocity

disp("Joint 2")
T2_tr = (1/2)* m2 * vc2'*vc2;
T2 = simplify(T2_tr);
disp("Kinetic energy")
disp(T2);

% joint 3 
%position of the CoMi in RF0
x3 = q1 + dc3*cos(q3);
y3 = q2 + dc3*sin(q3);
rc_03 = [x3; y3];

%velocity of CoMi
vx3 = diff(x3,q1)*dq1 + diff(x3,q2)*dq2 + diff(x3,q3)*dq3;
vy3 = diff(y3,q1)*dq1 + diff(y3,q2)*dq2 + diff(y3,q3)*dq3;
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

% Total Energy
T = T1+T2+T3;
disp("Total kinetic energy")
disp(T)
disp("Inertia matrix for the robot")

%-------------Inertia Matrix--------------

q_dot = [dq1; dq2; dq3];

M=simplify(hessian(T,q_dot));
disp(M)

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

disp("Cristoffel terms")
disp(c_v)

disp("SKEW SYMMETRIC MATRIX Computation")

% Compute Mdot (time derivative of M)
%add if you have more joints in this case they were 3

dM = simplify(diff(M,q1)*dq1+diff(M,q2)*dq2 +diff(M,q3)*dq3);


%they should be n matrices
disp("Ci matrices")
disp(C{1})
disp(C{2})
disp(C{3})


S = [dq' * C{1};
      dq' * C{2};
      dq' * C{3}];

S = simplify(S);
disp("S(q, dq):")
disp(S)

disp("Skew Symmetric Matrix")
skew= simplify(dM-2*S);
disp(skew)
disp("The matrix is skew symmetric if these values are = 0")
disp(simplify([1; 1; 1]'*skew*[1; 1; 1]))%this should be = 0
disp(skew'+ skew)%this should be = 0






%%

% e-e position

pd = [l1*cos(q1)+l2*cos(q1+q2)+l3*cos(q1+q2+q3);
       l1*sin(q1)+l2*sin(q1+q2)+l3*sin(q1+q2+q3)];

J = jacobian(pd, q);
disp("Jacobian")
disp(J)

%potential energy for the robot

%positions of CoMs in RF0
rc_01 = [dc1*cos(q1); dc1*sin(q1)];

rc_02 = [l1*cos(q1) + dc2*cos(q1+q2);
        l1*sin(q1) + dc2*sin(q1+q2)];

rc_03 = [l1*cos(q1) + l2*cos(q1+q2) + dc3*cos(q1+q2+q3);
        l1*sin(q1) + l2*sin(q1+q2) + dc3*sin(q1+q2+q3)];

syms g0 real
% vector gravity acceleration this depends on the drawing
%see how the vector is represented in RF0

g=[0; -g0];%in this case it was opposite of y0-axis

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

q_dot = [dq1; dq2; dq3];


dGdq = jacobian(G, q);

disp(dGdq)

%derivata pd
px = pd(1)
py = pd(2)
vx = diff(px,q1)*dq1 + diff(px,q2)*dq2 + diff(px,q3)*dq3;
vy = diff(py,q1)*dq1 + diff(py,q2)*dq2 + diff(py,q3)*dq3;
r_dot = [vx; vy;];
disp(r_dot)

invG = inv(dGdq);

q_dot = simplify(invG*J'*inv(J*invG*J')*r_dot);
disp("Velocity of the joints for the task")
disp(q_dot)


