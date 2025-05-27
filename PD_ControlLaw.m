%Control law PD for a PR robot
clear all

syms q1 q2 real
syms l1 l2 real
syms dc2 real
syms dq1 dq2 real
syms ddq1 ddq2 real
syms m1 m2 real 

n = 2; %number of joints

q = [q1; q2];
dq = [dq1; dq2];
%ddq = [ddq1; ddq2];

%Position CoM for each link
%!!!CoM position vector is the vector from the origin of frame i to the CMi
% expressed in frame 0

rc_01 = [q1; 0; 0];
rc_02 = [q1 + dc2*cos(q2); dc2*sin(q2); 0];



%-------------Gravity term--------------

disp("Potential Energy for each link")
pause

syms g0
% vector gravity acceleration this depends on the drawing
%see how the vector is represented in RF0

g=[g0; 0; 0];%in this case it was coherent with x0-axis

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
pause

U_tot = simplify(U{1} + U{2});

disp(U_tot)

G=transpose(jacobian(U_tot,q));

disp("Gravity term for the robot")
pause

disp(G)

%--------------Control Law PD------------
syms u Kp Kd real
qd = [0; pi];


Gqd = subs(G, q, qd);
u = Kp*(qd-q)-Kd*dq+ Gqd;
%%
%partial derivative ∂g/∂q
dgdq = jacobian(G, q);

disp("Partial Derivative G")
disp(dgdq)

%--------------Eigenvalues of the partial derivative
syms lambda

% Form J^T * J
JTJ = simplify(transpose(dgdq) * dgdq);  


eq = det(lambda * eye(size(JTJ)) - JTJ);

% Solve for lambda
lambdas = solve(eq == 0, lambda);
disp('Eigenvalues (lambda):')
disp(lambdas)



%%
disp("Norm of the Partial Derivative G")
disp(simplify(norm(dgdq)))
