%Midterm 2024 es 2
%2R robot polar
clear all

n = 2; %number of joints


syms l dc m q [n 1] real
syms dq [n 1] real
syms ddq [n 1] real


q = [q1; q2];
%dq = [dq1; dq2]; to define later
ddq = [ddq1; ddq2];

%DH parameters
alpha = [-(sym(pi)/2), 0]; %IMPORTANT to make pi symbolic
a = [0, l2];
d = [l1, 0];
theta = [q1, q2];

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


disp("Rotation matrices")
disp(R1)
disp(R2)


disp("Origin's position in RFi-1")
%pause

disp(p1)
disp(p2)


%position vectors from origin i-1 to origin i in frame i
% ip = i-1Ri' * i-1p !!! i am doing this computation here
% i-1p = i-1Ri * ip

r01 = simplify(R1'*p1);
r12 = simplify(R2'*p2);


disp("Position vectors from origin i-1 to origin i in frame i")
%pause 


disp(r01)
disp(r12)


%!!!CoM position vector is the vector from the origin of frame i to the CMi
% expressed in frame i
%in this case we assume it to be in the kinematic axis of each link
%if we define di as the distance between origin of frame i to CMi

%Understanding the link axis from DH parameters:
%if alpha_i = 0 the link extends on x-axis of frame i
%if alpha_i = pi/2 the link extends on y_axis of frame i
%if alpha_i = pi the link extends on -x-axis of frame i

%in this case they are given

syms rc_1_x rc_1_y rc_1_z rc_2_x real

rc_1 = [rc_1_x; rc_1_y; rc_1_z];

rc_2= [rc_2_x; 0; 0];



disp("Position vectors in frame i for CoMi for each link")
disp(rc_1)
disp(rc_2)


disp("CoM position vectors in RF0")
%pause


% ip = i-1Ri' * i-1p 
% i-1p = i-1Ri * ip !!i am doing this computation here

rc_01 = simplify((R1*rc_1) + p1);
disp(rc_01)
rc_02 = simplify((R1* (R2*rc_2)) + R1*p2 + p1);
disp(rc_02)



disp("MOVING FRAMES ALGORITHM")
%pause

%-------MOVING FRAMES ALGORITHM-------


%Defining the sigma as: 0 if revolute joint, 1 if prismatic
s1 = 0;
s2 = 0;



%if the robot were planar  
%syms Ic1 Ic2 Ic3 real

%the robot is polar so
syms Ic_xx_1 Ic_yy_1 Ic_zz_1 Ic_xx_2 Ic_yy_2 Ic_zz_2 real

syms Ic_yx_1 Ic_zx_1 Ic_xy_1 Ic_zy_1 Ic_xz_1 Ic_yz_1 real

Ic1 = [Ic_xx_1, Ic_xy_1, Ic_xz_1;
        Ic_xy_1, Ic_yy_1, Ic_yz_1;
        Ic_xz_1, Ic_yz_1, Ic_zz_1];

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
%pause


%Total kinetic Energy of the Robot 
T_tot = T{1} +T{2};
T_tot=simplify(T_tot);
T_tot=collect(T_tot,dq1^2);
T_tot=collect(T_tot,dq2^2);
disp(T_tot)


disp("Inertia matrix for the robot")
%pause


%-------------Inertia Matrix--------------
dq = [dq1; dq2];

M=simplify(hessian(T_tot,dq));
disp(M)

disp("Robot Centrifugal and Coriolis terms")
%pause

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
%pause

syms g0
% vector gravity acceleration this depends on the drawing
%see how the vector is represented in RF0

g=[-g0; 0; 0];%in this case it was opposite of x0-axis

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
%pause

U_tot = simplify(U{1} + U{2});

disp(U_tot)

G=transpose(jacobian(U_tot,q));

disp("Gravity term for the robot")
%pause

disp(G)


disp("Complete dynamic equation")
%pause

%-------------Final Dynamic Model--------------
syms u1 u2 u3 real

disp(M)
disp("*")
disp(ddq)
disp("+")
disp(c_v)
disp("+")
disp(G)
disp("=")
disp([u1; u2; u3])

disp("Linear Parametrization")
%pause

%%

%-------------Minimal parametrization of the model--------------

%First you need to find it by hand and then you can substitute
syms a1 a2 a3 a4 a5 a6 real
a = [Ic_yy_1 + Ic_xx_2 + m1*rc_1_x^2 + m1*rc_1_z^2; 
    Ic_yy_2 - Ic_xx_2 + m2*(l2 + rc_2_x)^2; 
    Ic_zz_2 + m2*(l2 + rc_2_x)^2; 
    -m2*(l2 +rc_2_x)^2;
    -m1*rc_1_z;
    -m1*rc_1_x];

%since we had to "change" M(q) with the factorization of:
%Ic_xx_2*sin(q2)^2 = Ic_xx_w*(1-cos(q2)^2
%let's redefine all the terms with the a values substituted

M_subs = [a1 + a2*cos(q2)^2, 0;
    0, a3];
c_v_subs = [-2*dq1*dq2*a2*cos(q2)*sin(q2);
            dq1^2*a2*cos(q2)*sin(q2)];
G_subs = [g0*(a5*cos(q1)+a6*sin(q1)+a4*cos(q2)*sin(q1));
            g0*a4*cos(q1)*sin(q2)];

a_symb = [a1; a2; a3; a4; a5; a6];

tau_subs = M_subs*ddq+c_v_subs+G_subs;
%tau_subs = subs(tau, a, a_symb);

disp(tau_subs)

Y = jacobian(tau_subs, a_symb);


Y=subs(Y, a_symb, a);
disp(Y)

%-----------------Skew Symmetric matrix------------
%pause 
disp("SKEW SYMMETRIC MATRIX Computation")

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


disp("Final tau with a joint trajectory given")
%pause

%-------------Final tau with a joint trajectory given--------------

syms t w real
qd = [2*t, pi/4];%joint trajectory given
dqd = diff(qd,t);%joint velocity
ddqd = diff(diff(qd,t), t);%joint acceleration

disp("Desired Joint Trajectory")
disp(qd)
disp("Desired Joint Velocity")
disp(dqd)
disp("Desired Joint Acceleration")
disp(ddqd)


disp("Expression of the Torque")
tau = M * ddq + c_v + G;
disp(tau)
disp("Final Expression of the Torque")
tau_final = subs(tau, [q1, q2, dq1, dq2, ddq1, ddq2], [qd, dqd, ddqd]);
disp(tau_final)

disp("Initial Configuration t=0")
disp("Initial joint trajectory")
disp(subs(qd, t, 0))
disp("Initial joint velocity")
disp(subs(dqd, t, 0))
disp("Initial joint acceleration")
disp(subs(ddqd, t, 0))
disp("Initial torque for the robot")
disp(subs(tau_final, t, 0))

disp("Final Configuartion t=pi")
disp("Final joint trajectory")
disp(subs(qd, t, pi))
disp("Final joint velocity")
disp(subs(dqd, t, pi))
disp("Final joint acceleration")
disp(subs(ddqd, t, pi))
disp("Final torque for the robot")
disp(subs(tau_final, t, pi))