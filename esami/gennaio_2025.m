%es 1 dynamic model 2R polar divided in a, b, c, d, e points
clear all

% a)

n=2; %number of joints

syms q l dc m dq ddq [n 1] real


%DH parameters
alpha = [sym(pi)/2, 0]; %IMPORTANT to make pi symbolic
a = [0, l2];
d = [l1, 0];
q = [q1, q2];

A = cell(1, n);
T_tot = eye(4);


for i = 1:n
    alphai = alpha(i);
    ai = a(i);
    di = d(i);
    qi = q(i);

    A{i} =  [cos(qi), -sin(qi)*cos(alphai), sin(qi)*sin(alphai), ai*cos(qi);
            sin(qi), cos(qi)*cos(alphai), -cos(qi)*sin(alphai), ai*sin(qi);
            0, sin(alphai), cos(alphai), di;
            0, 0, 0, 1];

    %complete DH matrix for the robot

    T_tot = simplify(T_tot * A{i}); 

end 

%extracting from the homogeneous matrices the rotation matrices and
%positions for each joint

R1 = A{1}(1:3,1:3);
p1 = A{1}(1:3,4);

R2 = A{2}(1:3,1:3); 
p2 = A{2}(1:3,4);


disp("Rotation matrices")
disp(R1)
disp(R2)



%position vectors from origin i-1 to origin i in frame i
r01 = simplify(R1'*p1);
r12 = simplify(R2'*p2);


disp("Position vectors from origin i-1 to origin i in frame i")
disp(r01)
disp(r12)



%!!!CoM position vector is the vector from the origin of frame i to the CMi
% expressed in frame i?
%in this case we assume it to be in the kinematic axis of each link
%if we define di as the distance between origin of frame i to CMi

%Understanding the link axis from DH parameters:
%if alpha_i = 0 the link extends on x-axis of frame i
%if alpha_i = pi/2 the link extends on y_axis of frame i
%if alpha_i = pi the link extends on -x-axis of frame i

% to find the link lengths look at the ai parameter:
%in this case a2 and a3 are our link lengths 

%we define d1 d2 d3 as the distance between origin i-1 and position of CoMi

rc_1 = [0; dc1 - l1; 0];

rc_2= [dc2 - l2; 0; 0];


disp("Position vectors in frame i for CoMi for each link")
disp(rc_1)
disp(rc_2)


disp("CoM position vectors in RF0")
% ip = i-1Ri' * i-1p 
% i-1p = i-1Ri * ip !!i am doing this computation here

rc_01 = simplify((R1*rc_1) + p1);
disp(rc_01)
rc_02 = simplify((R1* (R2*rc_2)) + R1*p2 + p1);
disp(rc_02)


disp("MOVING FRAMES ALGORITHM")

%-------MOVING FRAMES ALGORITHM-------


%Defining the sigma as: 0 if revolute joint, 1 if prismatic
s1 = 0;
s2 = 0;


%if the robot were planar  
%syms Ic1 Ic2 real

%in this case it is polar so
syms Ic_xx_1 Ic_yy_1 Ic_zz_1 Ic_xx_2 Ic_yy_2 Ic_zz_2 real
Ic1 = diag([Ic_xx_1, Ic_yy_1, Ic_zz_1]);
Ic2 = diag([Ic_xx_2, Ic_yy_2, Ic_zz_2]);


T = cell(1, n);
w = cell(1, n);
v = cell(1, n);
vc = cell(1, n);

%set up to change if we have more joints
s={s1, s2};
dq = {dq1, dq2};

R = {R1, R2};
r = {r01, r12};
rc = {rc_1, rc_2};
m = {m1, m2};
Ic = {Ic1, Ic2};


for i = 1:n
    Ri = R{i};
    si = s{i};
    dqi = dq{i};
    ri = r{i};
    rci  = rc{i};
    mi = m{i};
    Ici = Ic{i};

    %initialization if first frame = 0
    if i == 1
        wi_1 = 0;
        vi_1 = 0;
    else
        wi_1 = w(i-1);
        vi_1 = v(i-1);
    end

    w{i} = simplify(transpose(Ri)*(wi_1+(1-si)*dqi*[0; 0; 1]));

    wi = w{i};

    v{i} = simplify(transpose(Ri)*(vi_1+(si)*dqi*[0; 0; 1]) + cross(wi,ri));

    vi = v{i};

    vc{i} = simplify(vi + cross(wi,rci));

    vci = vc{i};

    T{i} = (1/2)*mi*transpose(vci)*vci + (1/2)*transpose(wi)*Ici*wi;
    
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


%Total kinetic Energy of the Robot 
T_tot = T{1} +T{2};
T_tot=simplify(T_tot);
T_tot=collect(T_tot,dq1^2);
T_tot=collect(T_tot,dq2^2);
disp(T_tot)


disp("Inertia matrix for the robot")

%Robots inertia matrix
dq = [dq1; dq2];
M=simplify(hessian(T_tot,dq));
disp(M)

%%

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


%-------------Gravity term--------------

disp("Potential Energy for each link")

syms g0 real
% vector gravity acceleration this depends on the drawing
%see how the vector is represented in RF0

g=[0; 0; -g0];%in this case it was opposite of z0-axis

U = cell(1, n);
%position of each CoM in RF0
rc_0 = {rc_01, rc_02};

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

U_tot = simplify(U{1} + U{2});

disp(U_tot)

G=transpose(jacobian(U_tot,q));

disp("Gravity term for the robot")
disp(G)

disp("SKEW SYMMETRIC MATRIX Computation")

% Compute Mdot (time derivative of M)

dM = simplify(diff(M,q1)*dq1+diff(M,q2)*dq2);


%they should be n matrices
disp("Ci matrices")
disp(C{1})
disp(C{2})


S = [dq' * C{1};
      dq' * C{2}];

S = simplify(S);
disp("S(q, dq):")
disp(S)

disp("Skew Symmetric Matrix")
skew= simplify(dM-2*S);
disp(skew)
disp("The matrix is skew symmetric if these values are = 0")
disp(simplify([1; 1]'*skew*[1; 1]))%this should be = 0
disp(skew'+ skew)%this should be = 0

%%
% b)
disp("Linear Parametrization")

%-------------Minimal parametrization of the model--------------

ddq = [ddq1; ddq2];

%First you need to find it by hand and then you can substitute
syms a1 a2 a3 a4 real
a = [m2*dc2^2+Ic_yy_2-Ic_xx_2; Ic_xx_2+Ic_yy_1; m2*dc2^2+Ic_zz_2; dc2*m2];
a_symb = [a1; a2; a3; a4];

%since we changed M we rewrite it

M_subs = [a1*cos(q2)^2 + a2, 0;
       0, a3];

c_v_subs = [-2*dq1*dq2*cos(q2)*sin(q2)*a1;
            dq1^2*cos(q2)*sin(q2)*a1];
G_subs = [0;
          g0*cos(q2)*a4];

%tau = M*ddq+S2+G;
%tau_subs = subs(tau, a, a_symb);

tau_subs = M_subs*ddq+c_v_subs+G_subs;

Y = jacobian(tau_subs, a_symb);

Y=subs(Y, a_symb, a);
disp(Y)

%%
% c)
%-----------Regulation control law-----------

q=[q1; q2];

dgdq = jacobian(G, q);
disp("Partial derivative of the gravity term wrt q")
disp(dgdq)

%max eigenvalue is on the diagonal the max eigenvalue = alpha

disp(norm(dgdq))

%AtA = dgdq'*dgdq;
%disp(AtA)

%the min eigenvalue of Kp needs to be >alpha to have global asymptotical
%stability of the system with the controller u = Kp(qd -q)-Kd*dq +g(q)

%%
% d)
%----------Expression of the torque given a joint trajectory-------------

syms t w real
qd = [cos(3*t), 0];
dqd = diff(qd, t);
disp("desired joint velocity")
disp(dqd)
ddqd = diff(diff(qd, t), t);
disp("desired joint acceleration")
disp(ddqd)
disp("Expression of the torque")
tau = M*ddq+c_v+G;
disp(tau)
disp("Final expression of the torque")
tau_final = subs(tau, [q1, q2, dq1, dq2, ddq1, ddq2], [qd, dqd, ddqd]);
disp(tau_final)

% e)
%-------------Adaptive control law---------------
%was computed through the regressor matrix of the minimal linear
%parametrization

%%
%es 2 3P planar robot exercise in 3 parts a, b, c
clear all

n=3; %number of joints

% a) min norm solution q_dot
syms q m dq [n 1] real

%task end effector's position in a trajectory

pe = [q1+q3;
      q2];
%jacobian for the task
J = jacobian(pe, [q1; q2; q3]);
disp("Jacobian for the task e-e tracking")
disp(J)
disp("rank of J (if = 2 the jacobian is full rank")
disp(rank(J))

%task velocity r_dot = pe_dot
r_dot = [dq1 + dq3;
        dq2];


q_dot_a = pinv(J)*r_dot;
disp("Joint velocities q_dot_a:")
disp(q_dot_a)

%%
% b) minimizing the kinetic energy while computing the joint velocities
%finding the inertia matrix for the robot

q = [q1; q2; q3];

%In this case all joints are prismatic so kinetic energy is composed of
%only the translational element.

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

%angular velocity there's no angular velocity on prismatic joints

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

%angular velocity there's no angular velocity on prismatic joints

disp("Joint 2")
T2_tr = (1/2)* m2 * vc2'*vc2;
T2 = simplify(T2_tr);
disp("Kinetic energy")
disp(T2);

% joint 3 
%position of the CoMi in RF0
x3 = q1+q3;
y3 = q2;
rc_03 = [x3; y3];

%velocity of CoMi
vx3 = diff(x3,q1)*dq1 + diff(x3,q2)*dq2 + diff(x3,q3)*dq3;
vy3 = diff(y3,q1)*dq1 + diff(y3,q2)*dq2 + diff(y3,q3)*dq3;
vc3 = [vx3; vy3];
disp("Velocity CoM")
disp(vc3);

%angular velocity there's no angular velocity on prismatic joints

disp("Joint 3")
T3_tr = (1/2)* m3 * vc3'*vc3;
T3 = simplify(T3_tr);
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

%q_dot_b is computewd using the LQ optimization
invM = inv(M);

q_dot_b  = simplify(invM*J'*inv(J*invM*J')*r_dot);
disp("joint velocity that minimizes kinetic energy for the robot")
disp(q_dot_b)

%%
% c) cartesian inertia matrix = mi*I2

syms mi real

M2 = diag([m1 m2 m3]);


disp("I*mi")
disp(eye(2)*mi)
Mp = simplify(inv(J*inv(M2)*J'));
disp("Cartesian inertia matrix")
disp(Mp)


eq = Mp == mi * eye(2);
sol = solve(eq, [m1, m2, m3], 'Real', true, 'IgnoreAnalyticConstraints', true);
disp(sol.m1)
disp(sol.m2)
disp(sol.m3)

