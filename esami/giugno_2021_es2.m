%giugno 2021
%es 2

clear all

n = 2; %number of joints


syms  m q [n 1]  real
syms d l Il real
syms dq  [n 1] real
syms ddq [n 1] real


q = [q1; q2];
%dq = [dq1; dq2]; to define later
ddq = [ddq1; ddq2];


pc2 = [d*cos(q2); q1+d*sin(q2);];


vcm2 = diff(pc2, q1)*dq1+diff(pc2, q2)*dq2;


disp(vcm2)

w2 =  [0; 0; dq2];

disp(w2);

T = (1/2)*m2*vcm2'*vcm2 + (1/2) * Il * w2'*w2;
T  = simplify(T);
%T=collect(T,dq1^2);
%T=collect(T,dq2^2);
disp(T)


%%
disp("Inertia matrix for the robot")


%-------------Inertia Matrix--------------
dq = [dq1; dq2];

M=simplify(hessian(T,dq));
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
syms k real 

disp("potential energy of the spring")
U =  (1/2)*k*q1^2;
disp(U)

%partial derivative of gradient Ue
dUdq  = jacobian(U, q);

disp(dUdq)


disp("Expression of the Torque")
tau = M * ddq + c_v + dUdq;
disp(tau)

%%
%------------Expression of the torque-----

M_subs = subs(M, [q1, q2, dq1, dq2], [0, 0, 0, 0]);
cv_subs = subs(c_v, [q1, q2, dq1, dq2], [0, 0, 0, 0]);
U_subs = subs(dUdq, [q1, q2, dq1, dq2], [0, 0, 0, 0]);

disp("Dynamc elements with q(0)=0 an dq(0)=0 substituted")
disp("inertia matrix")
disp(M_subs)
disp("centrifugal and coriolis terms")
disp(cv_subs)
disp("potential energy gradient")
disp(U_subs)

syms tau2 real
%since only the inertia matrix is different than 0 the dynamic model becomes:
%tau1 is = 0 since the spring has no torque
ddq = inv(M_subs)*[0; tau2];
disp("Acceleration of the joints with q(0)=0 an dq(0)=0")
disp(ddq)

%%
%----------Control law ?------------

syms Kd Kp real
syms qd1 qd2
dU  = [k*q1; 0];


disp("norm of the partial derivative")
disp(simplify(norm(dUdq)))

%%
%es 5
syms l1 l2 real 

pe = [l1*cos(q1)+l2*cos(q1+q2);
      l1*sin(q1)+l2*sin(q1+q2)];

J = jacobian(pe, [q1; q2]);
disp("jacobian for the robot")
disp(J)



