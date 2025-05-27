%FIRST HOMEWORK FOR ALL STUDENTS (lambda not fixed)
clear all
% length of the links

l1 = 1; 
l2 = 1;

% initial value of the parameter lambda (see the lecture on task control with RCM constraint)

lambda0=0.5 

% initial manipulator configuration

q0=[pi/4; 0]

%positions of points p1 and p2 in the picture, respectively the position of the first joint and of the end-effector

p1=[l1*cos(q0(1));l1*sin(q0(1))];


p2=[l1*cos(q0(1))+l2*cos(q0(1)+q0(2));l1*sin(q0(1))+l2*sin(q0(1)+q0(2))];

 
% final desired orientation of the second link

theta_d=pi/8 ;

%position of the trocar

x_tr=p1(1)+lambda0*(p2(1)-p1(1));

y_tr=p1(2)+lambda0*(p2(2)-p1(2));

rcm = [x_tr; y_tr];

%NEWTON METHOD FOR INVERSE KINEMATICS

max_iter = 10;       % Maximum number of iterations
tol = 1e-6;           % Tolerance for convergence (acceptable residual norm)


% Initial configuration for [q1; q2; lambda]
x = [q0; lambda0];

for i = 1:max_iter
    q1 = x(1);
    q2 = x(2);
    lambda = x(3);
    
    % Forward kinematics
    p1 = [l1*cos(q1); l1*sin(q1)];
    p2 = [p1(1) + l2*cos(q1+q2); p1(2) + l2*sin(q1+q2)];
    p_rcm = p1 + lambda*(p2 - p1);
    
    % Constraint function
    f = [p_rcm(1) - x_tr; %x_rcm == x_trocar
         p_rcm(2) - y_tr; %y_rcm == y_trocar
         q1 + q2 - theta_d]; %the orientation of the e-e == theta_d


    %Analytic Jacobian
    J = [-l1*sin(q1) - lambda*l2*sin(q1+q2),     -lambda*l2*sin(q1+q2),     (p2(1) - p1(1));
     l1*cos(q1) + lambda*l2*cos(q1+q2),      lambda*l2*cos(q1+q2),     (p2(2) - p1(2));
     1,                                      1,                        0];


    % Update state
    x = x - pinv(J)*f;

    % Lambda belongs to [0, 1] interval
    x(3) = max(0, min(1, x(3)));

    disp(['Iter ', num2str(i), ' Lambda value: ', num2str(x(3)), ' Error norm: ', num2str(norm(f))]); %displaying lambda value and error value for each iteration

    % Convergence check
    if norm(f) < tol
        disp(['Converged at iteration ', num2str(i)]);
        break;
    end
end


% Output final solution
q_sol = x(1:2);
lambda_sol = x(3);

% Display final joint values
fprintf('Final joint angles:\nq1 = %.4f rad\nq2 = %.4f rad\n', q_sol(1), q_sol(2));
fprintf('Final lambda value:\nlambda= %.4f \n', lambda_sol);


%TRAJECTORY PLANNING (cubic trajectory rest-to-rest)

t_total = 5;           % total animation time (seconds)
dt = 0.05;             % time step
time = 0:dt:t_total;   % time vector

q_start = [q0(1); q0(2)];  % initial joint position
q_end   = [q_sol(1); q_sol(2)];      % final joint position computed with newton method
delta_q = q_end - q_start;


q_traj = zeros(2, length(time));

for i = 1:length(time)
    t = time(i);
    tau = t / t_total;
    q_traj(:, i) = q_start + (3*tau^2 - 2*tau^3) .* delta_q;
end


figure;
axis equal;
axis([-2.5 2.5 -2.5 2.5]);
grid on;

for i = 1:length(time)
    q1 = q_traj(1,i);
    q2 = q_traj(2,i);
    
    p0 = [0; 0];
    p1 = p0 + [l1*cos(q1); l1*sin(q1)];
    p2 = p1 + [l2*cos(q1 + q2); l2*sin(q1 + q2)];

    % Compute trocar point at this time step
    p_trocar = p1 + lambda_sol * (p2 - p1);
   
    
    % Plot links
    clf;
    plot([p0(1) p1(1)], [p0(2) p1(2)], 'b-o', 'LineWidth', 2); hold on;
    plot([p1(1) p2(1)], [p1(2) p2(2)], 'r-o', 'LineWidth', 2);
    
    % Plot desired orientation
    quiver(x_tr, y_tr, cos(theta_d), sin(theta_d), 'k--', 'LineWidth', 1.2, 'MaxHeadSize', 0.5, 'Tag', 'arm');
    
    % Draw fixed RCM point
    plot(x_tr, y_tr, 'ko', 'MarkerSize', 8, 'Color', [1, 0.5, 0], 'LineWidth', 2);

    % Plot current trocar point as 'x'
     plot(p_trocar(1), p_trocar(2), 'kx', 'MarkerSize', 8, 'Color', [0.4, 0.3, 1], 'LineWidth', 1.5);

    
    title(['Time: ', num2str(time(i), '%.2f'), 's']);
    axis equal;
    axis([-2.5 2.5 -2.5 2.5]);
    grid on;
    drawnow;
end

disp('Press any key to close the figure.');
pause;  
close(gcf); 
