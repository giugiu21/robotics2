%PR planar robot in the vertical plane
%dynamic model
clear all

n=2;

syms q m l [n 1] real
syms dc2 real 
syms dq ddq [n 1] real


%if the robot were planar  
%syms Ic1 Ic2 real
syms Ic2 real %only second joint is revolute


%the robot is polar so
%syms Ic_xx_1 Ic_yy_1 Ic_zz_1 Ic_xx_2 Ic_yy_2 Ic_zz_2
%Ic1 = diag([Ic_xx_1, Ic_yy_1, Ic_zz_1]);
%Ic2 = diag([Ic_xx_2, Ic_yy_2, Ic_zz_2]);

q = [q1; q2];


%THE IMPORTANT THING IS THAT THE REFERENCE FRAMES ARE CONCORDANT

% joint 1 
%position of the CoMi in RF0

%for joint 1 the CoM is = to the origin of frame 1 = frame 0
x1 = q1;
y1 = 0;
rc_01 = [x1; y1];

%velocity of CoMi
vx1 = diff(x1,q1)*dq1 + diff(x1,q2)*dq2;
vy1 = diff(y1,q1)*dq1 + diff(y1,q2)*dq2;
vc1 = [vx1; vy1];
disp("Velocity CoM")
disp(vc1);

%angular velocity in frame 0
%first joint is prismatic so no rotational energy and no angular velocity

disp("Joint 1")
T1_tr = (1/2)* m1 * vc1'*vc1;
T1 = simplify(T1_tr);
disp("Kinetic energy")
disp(T1);


% joint 2 
%position of the CoMi in RFi
x2 = q1 + dc2*cos(q2);
y2 = dc2*sin(q2);
rc_02 = [x2; y2];

%velocity of CoMi
vx2 = diff(x2,q1)*dq1 + diff(x2,q2)*dq2;
vy2 = diff(y2,q1)*dq1 + diff(y2,q2)*dq2;
vc2 = [vx2; vy2];
disp("Velocity CoM")
disp(vc2);

%angular velocity
w2 = [0; 0; dq2]; 

disp("Joint 2")
T2_tr = (1/2)* m2 * vc2'* vc2;
T2_rot = (1/2)*  w2' * Ic2 * w2;
T2 = simplify(T2_tr + T2_rot);
disp("Kinetic energy")
disp(T2);


%%
%Total kinetic Energy of the Robot 
T_tot = T1 +T2;
T_tot=simplify(T_tot);
T_tot=collect(T_tot,dq1^2);
T_tot=collect(T_tot,dq2^2);
disp(T_tot)


disp("Inertia matrix for the robot")


%-------------Inertia Matrix--------------
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
pause

syms g0
% vector gravity acceleration this depends on the drawing
%see how the vector is represented in RF0

g=[g0; 0];%in this case it was concordant with x0-axis

U = cell(1, n);
%position of each CoM in RF0
rc_0 = {rc_01, rc_02};
m = {m1, m2};

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

disp("Gravity term for the robot g(q)")
disp(G)

%%
disp(simplify(M*[ddq1; 0]))

%%
clear all 

n=2;

syms q m l [n 1] real 
syms dq ddq [n 1] real

T1=  (1/2)*m1*dq1^2;
T2=  (1/2)*m2*(dq1^2 + dq2^2);

%Total kinetic Energy of the Robot 
T_tot = T1 +T2;
T_tot=simplify(T_tot);
T_tot=collect(T_tot,dq1^2);
T_tot=collect(T_tot,dq2^2);
disp(T_tot)


disp("Inertia matrix for the robot")


%-------------Inertia Matrix--------------
dq = [dq1; dq2];

M=simplify(hessian(T_tot,dq));
disp(M)


