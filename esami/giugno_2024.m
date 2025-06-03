%es 1
%Use moving frames algorithm to for a 3R robot 
%provide inertia matrix and a minimum linear parametrization
clear all 

n = 3;

syms q dc m a [n 1] real
syms dq ddq [n 1] real


q = [q1; q2; q3];
%dq = [dq1; dq2; dq3]; to define later
ddq = [ddq1; ddq2; ddq3];

%DH parameters
alpha = [sym(pi)/2, 0, 0]; %IMPORTANT to make pi symbolic
a = [0, a2, a3];
d = [0, 0, 0];
theta = [q1, q2, q3];

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

disp("Rotation matrices")
disp(R1)
disp(R2)
disp(R3)

disp("Origin's position in RFi-1")
pause

disp(p1)
disp(p2)
disp(p3)

%position vectors from origin i-1 to origin i in frame i
% ip = i-1Ri' * i-1p !!! i am doing this computation here
% i-1p = i-1Ri * ip

r01 = simplify(R1'*p1);
r12 = simplify(R2'*p2);
r23 = simplify(R3'*p3);

disp("Position vectors from origin i-1 to origin i in frame i")
pause 


disp(r01)
disp(r12)
disp(r23)

%%
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

%we define dc1 dc2 dc3 as the distance between origin i-1 and position of CoMi

rc_1 = [0; dc1; 0]; %coincident with the origin of its reference frame

rc_2= [dc2 - a2; 0; 0];

rc_3= [dc3 - a3; 0; 0];


disp("Position vectors in frame i for CoMi for each link")
disp(rc_1)
disp(rc_2)
disp(rc_3)

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


disp("MOVING FRAMES ALGORITHM")
pause


%-------MOVING FRAMES ALGORITHM-------


%Defining the sigma as: 0 if revolute joint, 1 if prismatic
s1 = 0;
s2 = 0;
s3 = 0;


%if the robot were planar  
%syms Ic1 Ic2 Ic3 real

%the robot is polar so
syms Ic_xx_1 Ic_yy_1 Ic_zz_1 Ic_xx_2 Ic_yy_2 Ic_zz_2 Ic_xx_3 Ic_yy_3 Ic_zz_3
Ic1 = diag([Ic_xx_1, Ic_yy_1, Ic_zz_1]);
Ic2 = diag([Ic_xx_2, Ic_yy_2, Ic_zz_2]);
Ic3 = diag([Ic_xx_3, Ic_yy_3, Ic_zz_3]);



T = cell(1, n);
w = cell(1, n);
v = cell(1, n);
vc = cell(1, n);

%set up to change if we have more joints
s={s1, s2, s3};
dq = {dq1, dq2, dq3};

R = {R1, R2, R3};
r = {r01, r12, r23};
rc = {rc_1, rc_2, rc_3};
m = {m1, m2, m3};
Ic = {Ic1, Ic2, Ic3};


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
pause


%Total kinetic Energy of the Robot 
T_tot = T{1} +T{2} +T{3};
T_tot=simplify(T_tot);
T_tot=collect(T_tot,dq1^2);
T_tot=collect(T_tot,dq2^2);
T_tot=collect(T_tot,dq3^3);
disp(T_tot)



disp("Inertia matrix for the robot")
pause


%-------------Inertia Matrix--------------
dq = [dq1; dq2; dq3];

M=simplify(hessian(T_tot,dq));
disp(M)


%%
%es 2 4R planar robot
%Primary task -> keeping the e-e at a position pd

clear all

n=4; % number of joints

syms q [n 1] real

q= [q1; q2; q3; q4];

%position of the e-e
pe = [cos(q1)+cos(q1+q2)+cos(q1+q2+q3)+cos(q1+q2+q3+q4);
       sin(q1)+sin(q1+q2)+sin(q1+q2+q3)+sin(q1+q2+q3+q4)];

%Jacobian for the task
J = jacobian(pe, [q1; q2; q3; q4]);
disp("Jacobian for the task")
disp(J)

%jacobian evaluated at values given q0
q0 = [0; pi/2; 0; -pi/4];

J_q0 = eval(subs(J, q, q0));
disp("Jacobian evaluated at q0")
disp(J_q0)

%rank of the jacobian
disp("Rank of the jacobian")
disp(rank(J_q0))

%%
%H(q) range function
q_min = [-pi/2; 0; -pi/4; -pi/4];
q_max = [pi/2; pi/2; pi/4; pi/4];

q_mid = (q_min + q_max)/2;
range = q_max - q_min;

H  = sum(((q-q_mid) ./ range).^2);

grad_H = gradient(H, q);
disp("gradient")
disp(grad_H)

%eval gradient at q0
grad_H_q0 = subs(grad_H, q, q0);
disp("Gradient evaluated at q0")
disp(grad_H_q0)

%velocity desired r_dot NOT given
%I WAS GIVEN THE POINT VELOCITY NOT THE TASK VELOCITY SO THIS IS WRONG
r_dot = [-1; -1];

%solution PG
q_dot = eval(pinv(J_q0)*r_dot + (eye(4) - pinv(J_q0)*J_q0)*grad_H_q0);

disp("solution in PG")
disp(q_dot)

