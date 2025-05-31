%Es1 cartesian control
clear all
%----------a------------
n = 3; %num joints
syms q m dc [n 1] real
syms l1 l3 real
syms g0 real

g = [0; -g0]; %gravity vector

%center of mass position for each link in RF0

r01 = [dc1*cos(q1); 
        dc1*sin(q1)];

r02 = [l1*cos(q1)-dc2*sin(q1); 
        l1*sin(q1)+dc2*cos(q1)];

r03 = [l1*cos(q1)-q2*sin(q1)+dc3*cos(q1+q3); 
        l1*sin(q1)+q2*cos(q1)+dc3*sin(q1+q3)];

r0 = {r01, r02, r03};

%potential energy for each link
U = cell(1, n);
q = [q1; q2; q3];
m = {m1, m2, m3};


for i=1:n
    mi = m{i};
    r0i = r0{i};

    U{i} = simplify(-mi * transpose(g) * r0i);

end

disp("potential energy per link")
for i=1:n
    disp("link")
    disp(i)
    disp(U{i})
end    
U_tot = simplify(U{1}+U{2}+U{3});

disp("total potential energy")
disp(U_tot)

disp("Gravity term dynamic model")
G = transpose(jacobian(U_tot,q));
disp(G)


%----------b-------

%e-e position
p = [l1*cos(q1)-q2*sin(q1)+l3*cos(q1+q3); 
    l1*sin(q1)+q2*cos(q1)+l3*sin(q1+q3)];

%jacobian for the robot
J = jacobian(p, [q1, q2, q3]);
disp("Jacobian for the robot")
disp(J)

%control law (cartesian control) -> u = J'*Kp(pd-p) + Kd*q_dot + G(q)

%----------c---------

disp("rank jacobian")
disp(rank(J))

%determinant of the minors

J_1 = [-sin(q1), -l3*sin(q1 + q3);
         cos(q1),  l3*cos(q1 + q3)];

J_2 = [- l3*sin(q1 + q3) - q2*cos(q1) - l1*sin(q1), -l3*sin(q1 + q3);
         l3*cos(q1 + q3) + l1*cos(q1) - q2*sin(q1), l3*cos(q1 + q3)];

J_3 = [- l3*sin(q1 + q3) - q2*cos(q1) - l1*sin(q1), -sin(q1);
        l3*cos(q1 + q3) + l1*cos(q1) - q2*sin(q1),  cos(q1)];

disp("Determinant of the minors")
disp("J-1")
disp(simplify(det(J_1)))
disp("J-2")
disp(simplify(det(J_2)))
disp("J-3")
disp(simplify(det(J_3)))


%these det are all 0 if q2=0 and q3=0 or q3=pi

qs = [q1; 0; 0];


%p in a singular configuration for the joints
ps = subs(p, q, qs);
disp("e-e pos in a singular configuration for the joints")
disp(ps)
Js = subs(J, q, qs);
disp("Jacobian a singular configuration for the joints")
disp(Js)


disp("null space of the singular jacobian transpose")
disp(null(Js'))

%!!!The cartesian controller fails if Kp*(pd-p)∈ N(J') and we are in a singular configuration 

