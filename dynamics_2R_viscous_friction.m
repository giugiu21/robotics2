%2R dynamic model viscous friction
%planar robot

n=2;

syms q l dc m [n 1] real
syms dq ddq [n 1] real

syms Ic1 Ic2 real %robot is planar

q =  [q1; q2];

% joint 1
%position of the CoMi in RFi
x1 = dc1*cos(q1);
y1 = dc1*sin(q1);
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
T1_tr = (1/2)* m1 * vc1'* vc1;
T1_rot = (1/2)*  w1' * Ic1 * w1;
T1 = simplify(T1_tr + T1_rot);
disp("Kinetic energy")
disp(T1);


% joint 2 
%position of the CoMi in RFi
x2 = l1*cos(q1);
y2 = l1*sin(q1);
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
T2_tr = (1/2)* m2 * vc2'* vc2;
T2_rot = (1/2)*  w2' * Ic2 * w2;
T2 = simplify(T2_tr + T2_rot);
disp("Kinetic energy")
disp(T2);

%Total kinetic Energy of the Robot 
T_tot = T1 +T2;
T_tot=simplify(T_tot);
T_tot=collect(T_tot,dq1^2);
T_tot=collect(T_tot,dq2^2);
disp("Total kinetic energy for the robot")
disp(T_tot)


disp("Inertia matrix for the robot")


%-------------Inertia Matrix--------------
dq = [dq1; dq2];

M=simplify(hessian(T_tot,dq));
disp(M)
%%
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
disp("coriolis vector")
disp(c_v)

%-------------Gravity term--------------

disp("Potential Energy for each link")

syms g0 dc1 dc2 m1 m2 real
% vector gravity acceleration this depends on the drawing
%see how the vector is represented in RF0

g=[g0; 0];%in this case it was concordant with x0-axis

U = cell(1, n);
%position of each CoM in RF0
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


disp("Total potential energy for the robot")

U_tot = simplify(U{1} + U{2});

disp(U_tot)

G=transpose(jacobian(U_tot,q));

disp("Gravity term")
disp(G)

%%
%-------------Minimal parametrization of the model--------------

ddq = [ddq1; ddq2];

%viscous friction 
syms Fv1 Fv2 real
Fv_dq = [Fv1*dq1; Fv2*dq2];

%First you need to find it by hand and then you can substitute
syms a1 a2 a3 a4 a5 real
a = [m1*dc1^2+m2*l1^2+Ic1+Ic2; Ic2; g0*(m1*dc1+l1*m2); Fv1; Fv2];
a_symb = [a1; a2; a3; a4; a5];


tau = M*ddq+c_v+G+Fv_dq;
tau_subs = subs(tau, a, a_symb);

disp(tau_subs)

Y = jacobian(tau_subs, a_symb);

Y=subs(Y, a_symb, a);

disp("REGRESSOR MATRIX")
disp(Y)

