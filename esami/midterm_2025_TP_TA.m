%Midterm 2025 es 1

%---------------Task Priority Algorithm----------
syms q1 q2 q3 q4 q5 real

n = 3; %number of tasks
q = [q1; q2; q3; q4; q5];
q0 = [0; pi/3; -pi/4; pi/2; -pi/3];

%In this ex. the angles qi are all ABSOLUTE ANGLES

%Task 1 (e-e motion)
p1 = [cos(q1)+cos(q2)+cos(q3)+cos(q4)+cos(q5);
    sin(q1)+sin(q2)+sin(q3)+sin(q4)+sin(q5)] ; %e-e position

J1 = jacobian(p1, q); %Task 1 jacobian

%Task 2 (tip of link 3 motion)
p2 = [cos(q1)+cos(q2)+cos(q3);
    sin(q1)+sin(q2)+sin(q3)] ; %tip link 3 position

J2 = jacobian(p2, q); %Task 2 jacobian

%Task 2 (tip of link 2 motion)
p3 = [cos(q1)+cos(q2);
    sin(q1)+sin(q2)] ; %tip link 2 position

J3 = jacobian(p3, q); %Task 3 jacobian

%jacobians for the tasks
disp("Jacobians for the tasks in order")
disp(J1)
disp(J2)
disp(J3)

%task velocities
r_dot_1 = [2; 3];
r_dot_2 = [2; -0.5];
r_dot_3 = [-1; 0];

%jacobians evaluated at q0
J1 = double(subs(J1, q, q0));
J2 = double(subs(J2, q, q0));
J3 = double(subs(J3, q, q0));

J = {J1, J2, J3}; %jacobians for the tasks
r_dot = {r_dot_1, r_dot_2, r_dot_3}; %task velocities for each task

%--------------TASK PRIORITY ALGORITHM-------------------

%initialization
q_dot = zeros(5, 1); %5 is the number of joints
P = eye(5);

for i = 1:n
    Ji = J{i};
    r_dot_i = r_dot{i};

    Ji_proj = Ji * P;
    Ji_pinv = pinv(Ji_proj, 0.001);  % Compute the DAMPED pseudoinverse

    q_dot = q_dot + Ji_pinv * (r_dot_i - Ji * q_dot);  % Task correction
    P = P - Ji_pinv * Ji_proj;                         % Update nullspace
end

%We use a damped pseudoinverse when we get large velocities so that we can
%handle numerical instability. In this case we were probably near a
%singularity that's why we were gettinhg really high values 1.0e+15

q_dot_ts = q_dot;

disp('Joint velocities computed with Task Priority')
disp(q_dot_ts)

%computing the error for each task
e1 = r_dot_1  - J1*q_dot_ts;
e2 = r_dot_2  - J2*q_dot_ts;
e3 = r_dot_3  - J3*q_dot_ts;

disp('Error in norm for each tasks velocity')
disp(norm(e1))
disp(norm(e2))
disp(norm(e3))

%------------------TASK AUGMENTATION-------------------
%In this case all the tasks were to follow a certain velocity so we don't
%need to come up with an objective function to augment the jacobian

JA = [J1; J2; J3];
r_dot_A = [r_dot_1; r_dot_2; r_dot_3];

q_dot_ta = pinv(JA) * r_dot_A;

disp('Joint velocities computed with Task Augmentation')
disp(q_dot_ta)

%computing the error for each task
e1 = r_dot_1  - J1*q_dot_ta;
e2 = r_dot_2  - J2*q_dot_ta;
e3 = r_dot_3  - J3*q_dot_ta;

disp('Error in norm for each tasks velocity')
disp(norm(e1))
disp(norm(e2))
disp(norm(e3))

%---------------Mapping Absolute Angles with Theta Angles--------------

