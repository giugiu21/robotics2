close all
clear all
clc

%Find mapping between omega and phi dot
function omega = mapping_omega_phiDot(flag, sequence, angles)
    syms a b g dalfa dbeta dgamma 
    a = angles(1,1);
    b = angles(1,2);
    g = angles(1,3);
    if flag == 'r'
        R = rpy_rotation(sequence, angles);
    end
    if flag == 'e'
        R = euler_rotation(sequence, angles);
    end
    Rdot=diff(R,a)*dalfa+diff(R,b)*dbeta+diff(R,g)*dgamma
    S_omega=simplify(Rdot*transpose(R))
    omega=[S_omega(3,2);S_omega(1,3);S_omega(2,1)]
end

%Axis-angle inverse probelm
function [rfinale,rfinale1]=inverseProblemR(R)
   seno1=real((1/2)*sqrt((R(1,2)-R(2,1))^2+(R(1,3)-R(3,1))^2+(R(2,3)-R(3,2))^2))
   seno2=real(-(1/2)*sqrt((R(1,2)-R(2,1))^2+(R(1,3)-R(3,1))^2+(R(2,3)-R(3,2))^2))
   cose=(trace(R)-1)/2;
   teta=atan2(seno1,cose)
   teta1=atan2(seno2,cose)
   if(seno1==0)
       rx=sqrt((R(1,1)+1)/2);
       rx1=-sqrt((R(1,1)+1)/2);
       ry=sqrt((R(2,2)+1)/2);
       ry1=-sqrt((R(2,2)+1)/2);
       rz=sqrt((R(3,3)+1)/2);
       rz1=-sqrt((R(3,3)+1)/2);
       rfinale=[0; 0; 0];
       rfinale1=[0; 0; 0];
       if(2*rx*ry==R(1,2))
          disp("sono1")
          rfinale(1,1)=rx;
          rfinale(2,1)=ry;
       end
       if(2*rx*ry1==R(1,2))
          disp("sono2")
          rfinale(1,1)=rx;
          rfinale(2,1)=ry1;
       end
       if(2*rx1*ry==R(1,2))
          disp("sono3")
          rfinale(1,1)=rx1;
          rfinale(2,1)=ry;
       end
%        if(2*rx1*ry1==R(1,2))
%           disp("sono4")
%           rfinale(1,1)=rx1;
%           rfinale(2,1)=ry1;
%        end
       if(2*rz*rx==R(1,3))
          disp("sono5")
          rfinale(1,1)=rx;
          rfinale(3,1)=rz;
       end
       if(2*rz*rx1==R(1,3))
          disp("sono6")
          rfinale(1,1)=rx;
          rfinale(3,1)=rz;
       end
       if(2*rz1*rx==R(1,3))
          disp("sono7")
          rfinale(1,1)=rx;
          rfinale(3,1)=rz1;
       end
%        if(2*rz1*rx1==R(1,3))
%           disp("sono8")
%           rfinale(1,1)=rx1;
%           rfinale(3,1)=rz1;
%        end
       rfinale1=-rfinale;
       %R1=rfinale*rfinale'+(eye(3)-rfinale*rfinale')*cos(teta)+S*sin(teta);
       %R2=rfinale1*rfinale1'+(eye(3)-rfinale1*rfinale1')*cos(teta)+S*sin(teta);
   else 
       rfinale=(1/(2*sin(teta)))*[R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)];
       rfinale1=(1/(2*sin(teta1)))*[R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)];
   end
end

%Axis-angle matrix
function R=axis_angle_matrix(r,teta)
    S=[0 -r(3,1) r(2,1); r(3,1) 0 -r(1,1); -r(2,1) r(1,1) 0];
    R=r*r'+(eye(3)-r*r')*cos(teta)+S*sin(teta);
end

%Inverse problem of Axis-angle
function [theta, r] = get_theta_r(R)
% [theta, r] = get_theta_r(R) takes as inputs:
%   -R: A rotation matrix, tipycally a 3x3 matrix
% and outputs:
%   -theta: The radiants the rotation has been performed of
%   -r: The vector the rotation has been performed about
% Remember to use eval(symbolic function) to use this function
    risp = input("You can use either the positive or the negative sin(theta), this will correspond to 2 different solutions, which one you prefer? (pos, neg)\n", "s");
    
    if risp == "pos"
        sintheta= sqrt((R(1,2)-R(2,1))^2 + (R(1,3)-R(3,1))^2 + (R(2,3)-R(3,2))^2)/2;
        costheta= (R(1,1)+R(2,2)+R(3,3)-1)/2
        theta = atan2(sintheta,costheta);
        if theta == 0
            disp('Theta = 0, rotation axis is not defined, hence there is not a given solution');
            r = -1;
        else
            fprintf("Theta is %.2f (should be either +/- pi), hence we are in a singular case -> sin(theta) == 0\n", theta);
        
            r = (1/(2*sintheta))*[R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2);];
        end
    else
        sintheta= -sqrt((R(1,2)-R(2,1))^2 + (R(1,3)-R(3,1))^2 + (R(2,3)-R(3,2))^2)/2;
        costheta= (R(1,1)+R(2,2)+R(3,3)-1)/2
        theta = atan2(sintheta,costheta);
        
        if theta == 0
            disp('Theta = 0, rotation axis is not defined, hence there is not a given solution');
            r = -1;
        else
            fprintf("Theta is %.2f (should be either +/- pi), hence we are in a singular case -> sin(theta) == 0\n", theta);
        
            r = (1/(2*sintheta))*[R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2);];
        end
        
        
    end
end

%Construction of Rotation Matrix x,y,z
function R = elem_rot_mat(axis, s)
    switch axis
        case {"x", "X"}
            R = [1       0        0;
                 0    cos(s)   -sin(s);
                 0    sin(s)    cos(s)];
        case {"y", "Y"}
            R = [cos(s)     0   sin(s);
                   0	    1     0;
                -sin(s)     0   cos(s)];
        case {"z", "Z"}
            R = [cos(s)  -sin(s)    0;
                 sin(s)   cos(s)    0;
                   0        0       1];
        otherwise
            disp("First Parameter should be either 'x', 'y', 'z' or any of those capitalized")
    end
end

%Construction of euler rotation
function R = euler_rotation(sequence, angles)
    if strlength(sequence) ~= 3
        disp("Sequence not valid, must be of lenght three.")
        return;
    end
    
    sequence = lower(char(sequence));
    if (sequence(2) == sequence(1) || sequence(2) == sequence(3))
        disp("Two consecutive rotation along the same axis are not valid.")
        return
    end
    
    R = elem_rot_mat(sequence(1), angles(1)) * elem_rot_mat(sequence(2), angles(2)) * elem_rot_mat(sequence(3), angles(3));
end

%From euler rotation to angles
function [phi, theta, psi] = euler_rotation_inverse(sequence, R)
% [phi, theta, psi] = euler_rotation_inverse(sequence, R) takes as inputs:
%   -sequence: a string which specifies how the euler-rotation has been computed, e.g. "xyx"
%   -R: the rotation to be decomposed, should be a 3x3 matrix
% and outputs:
%   -phi: the radiants of the first rotation
%   -theta: the radiants of the second rotation
%   -psi: the radiants of the third rotation
% Upon execution the function will ask for an input which will determine
% which solution must be displayed (there is always a pair of solutions)
%
% Euler rotations work about moving-axes
    
    if isa(R, 'sym')
        disp("R must by non-symbolic matrix")
        phi = -1; theta = -1; psi = -1;
        return
    end
    
    if strlength(sequence) ~= 3
        disp("Sequence not valid, must be of length three.")
        return;
    end
    
    if (sequence(2) == sequence(1) || sequence(2) == sequence(3))
        disp("Two consecutive rotation along the same axis are not valid.")
        return
    end
    
    risp = input("There is always a pair of solutions to the inverse euler problem. Do you want to use the positive or negative sin(theta)? (pos, neg)\n", "s");
    cond = risp == "pos";
    
    switch lower(sequence)
        case "xyx"
            theta = atan2(sqrt(R(1, 2)^2 + R(1, 3)^2), R(1, 1))*cond + atan2(-sqrt(R(1, 2)^2 + R(1, 3)^2), R(1, 1))*(1-cond);
            if (abs(sin(theta)) <= 1e-6)
                disp('Singular case: sin(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(1, 2)/sin(theta), R(1, 3)/sin(theta));
            phi = atan2(R(2, 1)/sin(theta), -R(3, 1)/sin(theta));
            
        case "xyz"
            theta = atan2(R(1, 3), sqrt(R(1, 1)^2 + R(1, 2)^2))*cond + atan2(R(1, 3), -sqrt(R(1, 1)^2 + R(1, 2)^2))*(1-cond);
            if (abs(cos(theta)) <= 1e-6)
                disp('Singular case: cos(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(-R(1, 2)/cos(theta), R(1, 1)/cos(theta));
            phi = atan2(-R(2, 3)/cos(theta), R(3, 3)/cos(theta));
            
        case "xzx"
            theta = atan2(sqrt(R(1, 2)^2 + R(1, 3)^2), R(1, 1))*cond + atan2(-sqrt(R(1, 2)^2 + R(1, 3)^2), R(1, 1))*(1-cond);
            if (abs(sin(theta)) <= 1e-6)
                disp('Singular case: sin(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(1, 3)/sin(theta), -R(1, 2)/sin(theta));
            phi = atan2(R(3, 1)/sin(theta), R(2, 1)/sin(theta));
            
        case "xzy"
            theta = atan2(-R(1, 2), sqrt(R(1, 1)^2 + R(1, 3)^2))*cond + atan2(-R(1, 1), -sqrt(R(1, 3)^2 + R(2, 1)^2))*(1-cond);
            if (abs(cos(theta)) <= 1e-6)
                disp('Singular case: cos(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(1, 3)/cos(theta), R(1, 1)/cos(theta));
            phi = atan2(R(3, 2)/cos(theta), R(2, 2)/cos(theta));
            
        case "yxy"
            theta = atan2(sqrt(R(2, 3)^2 + R(2, 1)^2), R(2, 2))*cond + atan2(-sqrt(R(2, 3)^2 + R(2, 1)^2), R(2, 2))*(1-cond);
            if (abs(sin(theta)) <= 1e-6)
                disp('Singular case: sin(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(2, 1)/sin(theta), -R(2, 3)/sin(theta));
            phi = atan2(R(1, 2)/sin(theta), R(3, 2)/sin(theta));
            
        case "yxz"
            theta = atan2(-R(2, 3), sqrt(R(2, 2)^2 + R(2, 1)^2))*cond + atan2(-R(2, 3), -sqrt(R(2, 2)^2 + R(2, 1)^2))*(1-cond);
            if (abs(cos(theta)) <= 1e-6)
                disp('Singular case: cos(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(2, 1)/cos(theta), R(2, 2)/cos(theta));
            phi = atan2(R(1, 3)/cos(theta), R(3, 3)/cos(theta));
            
        case "yzx"
            theta = atan2(R(2, 1), sqrt(R(2, 2)^2 + R(2, 3)^2))*cond + atan2(R(2, 1), -sqrt(R(2, 2)^2 + R(2, 3)^2))*(1-cond);

            if (abs(cos(theta)) <= 1e-6)
                disp('Singular case: cos(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(-R(2, 3)/cos(theta), R(2, 2)/cos(theta));
            phi = atan2(-R(3, 1)/cos(theta), R(1, 1)/cos(theta));
            
        case "yzy"
            theta = atan2(sqrt(R(2, 1)^2 + R(2,3)^2), R(2, 2))*cond + atan2(-sqrt(R(2, 1)^2 + R(2,3)^2), R(2, 2))*(1-cond);
            if (abs(sin(theta)) <= 1e-6)
                disp('Singular case: sin(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(2, 3)/sin(theta), R(2, 1)/sin(theta));
            phi = atan2(R(3, 2)/sin(theta), -R(1, 2)/sin(theta));
            
        case "zxy"
            theta = atan2(R(3, 2), sqrt(R(3, 1)^2 + R(3,3)^2))*cond + atan2(R(3, 2), -sqrt(R(3, 1)^2 + R(3,3)^2))*(1-cond);
            if (abs(cos(theta)) <= 1e-6)
                disp('Singular case: cos(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(-R(3, 1)/cos(theta), R(3, 3)/cos(theta));
            phi = atan2(-R(1, 2)/cos(theta), R(2, 2)/cos(theta));
            
        case "zxz"
            theta = atan2(sqrt(R(1, 3)^2 + R(2, 3)^2), R(3, 3))*cond + atan2(-sqrt(R(1, 3)^2 + R(2, 3)^2), R(3, 3))*(1-cond);
            if (abs(sin(theta)) <= 1e-6)
                disp('Singular case: sin(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(3, 1)/sin(theta), R(3, 2)/sin(theta));
            phi = atan2(R(1, 3)/sin(theta), -R(2, 3)/sin(theta));
            
        case "zyx"
            theta = atan2(-R(3, 1), sqrt(R(3, 2)^2+R(3, 3)^2))*cond + atan2(-R(3, 1), -sqrt(R(3, 2)^2+R(3, 3)^2))*(1-cond);
            if (abs(cos(theta)) <= 1e-6)
                disp('Singular case: cos(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(3, 2)/cos(theta), R(3, 3)/cos(theta));
            phi = atan2(R(2, 1)/cos(theta), R(1, 1)/cos(theta));
            
        case "zyz"
            theta = atan2(sqrt(R(3, 1)^2 + R(3, 2)^2), R(3, 3))*cond + atan2(-sqrt(R(3, 1)^2 + R(3, 2)^2), R(3, 3))*(1-cond);
            if (abs(sin(theta)) <= 1e-6)
                disp('Singular case: sin(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(3, 2)/sin(theta), -R(3, 1)/sin(theta));
            phi = atan2(R(2, 3)/sin(theta), R(1, 3)/sin(theta));
            
        otherwise
            disp("Invalid sequence")
    end
    
end

%Construction of RPY rotation
function R = rpy_rotation(sequence, angles)

% R = rpy_rotation(sequence, angles) takes as inputs:
%   -sequence: a string which specifies which elementary rotation matrices
%              must be multiplied together to obtain the desired rotation
%   -angles: the radiants of the three rotation angles
% and outputs:
%   -R: The desired rotation
% RPY rotations work about fixed-axes

    sequence = char(sequence);
    R = euler_rotation(flip(sequence), flip(angles));
end

%From RPY rotation to angles
function [phi, theta, psi] = rpy_rotation_inverse(sequence, R)
% [phi, theta, psi] = rpy_rotation_inverse(sequence, R) takes as inputs:
%   -sequence: a string which specifies how the RPY-rotation has been computed, e.g. "xyx"
%   -R: the rotation to be decomposed, should be a 3x3 matrix
% and outputs:
%   -phi: the radiants of the first rotation
%   -theta: the radiants of the second rotation
%   -psi: the radiants of the third rotation
% Upon execution the function will ask for an input which will determine
% which solution must be displayed (there is always a pair of solutions)
%
% RPY rotations work about fixed-axes

    sequence = char(sequence);
    [psi, theta, phi] = euler_rotation_inverse(flip(sequence), R);
end

%Check if R is orthonotmal
function isOrthonormal = isOrthonormalMatrix(A)
    % Verifica che la matrice A sia ortonormale
    
    % Verifica che la matrice sia quadrata
    [m, n] = size(A);
    if m ~= n
        error('La matrice deve essere quadrata per essere ortonormale.');
    end
    
    % Verifica l'ortogonalità
    orthoCheck = norm(A' * A - eye(n), 'fro');
    
    % Verifica che i vettori delle colonne siano unitari
    normCheck = norm(A * A' - eye(n), 'fro');
    
    % Imposta la tolleranza
    tolerance = 1e-4;
    
    % Verifica complessiva
    isOrthonormal = (orthoCheck < tolerance) && (normCheck < tolerance);
end

%Inverse problem axis-angle SINGULAR CASE!
function [rfinale,rfinale1,R1,R2]=verificareR(R)
   seno1=(1/2)*sqrt((R(1,2)-R(2,1))^2+(R(1,3)-R(3,1))^2+(R(2,3)-R(3,2))^2)
   seno2=real(-(1/2)*sqrt((R(1,2)-R(2,1))^2+(R(1,3)-R(3,1))^2+(R(2,3)-R(3,2))^2));
   cose=(trace(R)-1)/2;
   teta=atan2(seno1,cose)
   teta1=atan2(seno2,cose)
   if(seno1==0)
       rx=sqrt((R(1,1)+1)/2);
       rx1=-sqrt((R(1,1)+1)/2);
       ry=sqrt((R(2,2)+1)/2);
       ry1=-sqrt((R(2,2)+1)/2);
       rz=sqrt((R(3,3)+1)/2);
       rz1=-sqrt((R(3,3)+1)/2);
       rfinale=[0; 0; 0];
       rfinale1=[0; 0; 0];
       if(2*rx*ry==R(1,2))
          disp("sono1")
          rfinale(1,1)=rx;
          rfinale(2,1)=ry;
       end
       if(2*rx*ry1==R(1,2))
          disp("sono2")
          rfinale(1,1)=rx;
          rfinale(2,1)=ry1;
       end
       if(2*rx1*ry==R(1,2))
          disp("sono3")
          rfinale(1,1)=rx1;
          rfinale(2,1)=ry;
       end
%        if(2*rx1*ry1==R(1,2))
%           disp("sono4")
%           rfinale(1,1)=rx1;
%           rfinale(2,1)=ry1;
%        end
       if(2*rz*rx==R(1,3))
          rfinale(1,1)=rx;
          rfinale(3,1)=rz;
       end
       if(2*rz*rx1==R(1,3))
          rfinale(1,1)=rx;
          rfinale(3,1)=rz;
       end
       if(2*rz1*rx==R(1,3))
          rfinale(1,1)=rx;
          rfinale(3,1)=rz1;
       end
%        if(2*rz1*rx1==R(1,3))
%           disp("sono8")
%           rfinale(1,1)=rx1;
%           rfinale(3,1)=rz1;
%        end
       rfinale1=-rfinale;
       R1=2*rfinale*rfinale'-eye(3)
       R2=2*rfinale1*rfinale1'-eye(3)
   else 
       rfinale=(1/(2*sin(teta)))*[R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)];
       rfinale1=(1/(2*sin(teta1)))*[R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)];
       S1=[0 -rfinale(3,1) rfinale(2,1); rfinale(3,1) 0 -rfinale(1,1); -rfinale(2,1) rfinale(1,1) 0];
       S2=[0 -rfinale1(3,1) rfinale1(2,1); rfinale1(3,1) 0 -rfinale1(1,1); -rfinale1(2,1) rfinale1(1,1) 0];
       R1=rfinale*rfinale'+(eye(3)-rfinale*rfinale')*cos(teta)+S1*sin(teta);
       R2=rfinale1*rfinale1'+(eye(3)-rfinale1*rfinale1')*cos(teta1)+S2*sin(teta1);
   end
end

%Newton Method
function [q_out, guesses, cartesian_errors] = newton_method(q_in, desired_point, f_r, initial_guess, max_iterations, max_cartesian_error, min_joint_increment, max_closeness_singularity)
% [q_out, guesses, cartesian_errors] = newton(q_in, desired_point, f_r, 
%  initial_guess, max_iterations, max_cartesian_error, min_joint_increment, 
%  max_closeness_singularity) takes as inputs:
%   -q_in: The variables of the joints, e.g. [q1 q2 q3 q4]
%   -desired_point: the configuration we wish to reach
%   -f_r: The mapping from joints to points
%   -initial_guess: Initial configuration of joints
%   -max_iterations: The maximum number of iterations we can perform
%   -max_cartesian_error: The level of precision we wish to reach
%   -min_joint_increment: The minimum level of increment of accuracy 
%                          between successive iterations.
%   -max_closeness_singularity: How close to a singularity we can get
% and outputs:
%   -q_out: The best reached configuration
%   -guesses: The history of tested configurations
%   -cartesian_errors: The history of errors

    J = jacobian(f_r, q_in);
    
    n_dof = length(q_in);
    n_config = length(desired_point);
    
    %simple_case = n_config == n_dof;
    
    if n_config < n_dof
        disp("Redundant case, let's use the pseudo inverse");
        J_inv = pinv(J);
    else
        disp("Non-redundant case, let's use the inverse");
        J_inv = inv(J);
    end
    
    guesses = zeros(max_iterations + 1, n_dof);
    cartesian_errors = zeros(1, max_iterations + 1);
    
    guess = initial_guess;
    
    for i = 1:max_iterations
        guesses(i, :) = guess;
        
        error = norm(desired_point - subs(f_r, q_in, guess));
        cartesian_errors(i) = error;
        if error < max_cartesian_error
            fprintf("Finshed at iteration %d because the error was lower than the specified amount.\n", i);
            break
        end
        
        % if simple_case
        %     if abs(det(subs(J, q_in, guess))) <= max_closeness_singularity
        %         fprintf("Finished at iteration %d because too close to a singularity.\n", i);
        %         break
        %     end
        % else
        %     [~, D, ~] = svd(subs(J, q_in, guess));
        %     if min(min(D)) >= D(1, 1) - max_closeness_singularity
        %         fprintf("Finished at iteration %d because too close to a singularity.\n", i);
        %         break
        %     end
        % end
        
        new_guess = guess + subs(J_inv, q_in, guess) * (desired_point - subs(f_r, q_in, guess));
        
        % if norm(new_guess - guess) <= min_joint_increment
        %     fprintf("Finished at iteration %d because of a too little joints values increment.\n", i);
        %     break
        % end
        
        guess = eval(new_guess);
    end
    
    guesses = guesses(1:i, :);
    cartesian_errors = cartesian_errors(1:i);
    q_out = guess;
end

%Gradient Method
function [q_out, guesses, cartesian_errors] = gradient_method(q_in, desired_point, f_r, initial_guess, alpha, max_iterations, max_cartesian_error, min_joint_increment, max_closeness_singularity)
% [q_out, guesses, cartesian_errors] = gradient(q_in, desired_point, f_r, 
%  initial_guess, alpha, max_iterations, max_cartesian_error, 
%  min_joint_increment, max_closeness_singularity) takes as inputs:
%   -q_in: The variables of the joints, e.g. [q1 q2 q3 q4]
%   -desired_point: the configuration we wish to reach
%   -f_r: The mapping from joints to points
%   -initial_guess: Initial configuration of joints
%   -alpha: The "learning rate" of the algorithm
%   -max_iterations: The maximum number of iterations we can perform
%   -max_cartesian_error: The level of precision we wish to reach
%   -min_joint_increment: The minimum level of increment of accuracy 
%                          between successive iterations.
%   -max_closeness_singularity: How close to a singularity we can get
% and outputs:
%   -q_out: The best reached configuration
%   -guesses: The history of tested configurations
%   -cartesian_errors: The history of errors

    J = jacobian(f_r, q_in);
    
    n_dof = length(q_in);
    n_config = length(desired_point);
    
    simple_case = n_config == n_dof;
    
    J_trs = J';
    
    guesses = zeros(max_iterations + 1, n_dof);
    cartesian_errors = zeros(1, max_iterations + 1);
    
    guess = initial_guess;
    
    for i = 1:max_iterations
        guesses(i, :) = guess;
        
        error = norm(desired_point - subs(f_r, q_in, guess));
        cartesian_errors(i) = error;
        
        if error <= max_cartesian_error
            fprintf("Finished at iteration %d\nBecause the error was lower than the specified amount.\n", i);
            break
        end
        
        if simple_case
            if abs(det(subs(J, q_in, guess))) <= max_closeness_singularity
                fprintf("Finished at iteration %d\nBecause we ended up too close to a singularity.\n", i);
                break
            end
        else
            [~, D, ~] = svd(subs(J, q_in, guess));
            if min(min(D)) >= D(1, 1) - max_closeness_singularity
                fprintf("Finished at iteration %d\nBecause we ended up too close to a singularity.\n", i);
                break
            end
        end
        
        p_k = gradient(eval(subs(f_r, q_in, guess)));
        while -p_k'*gradient(eval(subs(f_r, q_in, guess + alpha*p_k))) > -0.9*p_k'*gradient(eval(subs(f_r, q_in, guess)))%| subs(f_r, q_in, guess + alpha*desired_point) > subs(f_r, q_in, guess) + 1e-4*alpha*desired_point'*gradient(eval(subs(f_r, q_in, guess)))
            alpha = 0.99*alpha;
        end
        
        new_guess = guess + alpha * subs(J_trs, q_in, guess) * (desired_point - subs(f_r, q_in, guess));
        
        if norm(new_guess - guess) <= min_joint_increment
            fprintf("Finished at iteration %d\nBecause of a too little increase between iterations.\n", i);
            break
        end
        
        guess = eval(new_guess);
    end

    guesses = guesses(1:i, :);
    cartesian_errors = cartesian_errors(1:i);
    q_out = guess;
end

%Geometric Jacobian
function [Jl, Ja] = geometric_jacobian(f_r, sequence, q_in, params)
% geoJ = geometric_jacobian() takes as inputs:
%   -f_r: mapping from joint space to cartesian space
%   -sequence: A string containing the sequence of 'r's and 'p's
%   -q_in: The list of the variables we want to derive wrt
%   -params: List of params needed to get the DHMatrix
% and outputs:
%   -geoJ: The resulting geometric jacobian

    [~, A] = DHMatrix(params);
    sequence = char(lower(sequence));
    n_dof = strlength(sequence);
    
    f_r = f_r(1:3, :);
    
    %to calculate Jl, we use analytic jacobian
    Jl = jacobian(f_r, q_in);
    
    Ja = sym(zeros(3, n_dof));
    
    cells = cell(1, n_dof);
    cells{1} = get_rotation_mat(A{1});
    for i = 2:n_dof
        cells{i} = cells{i-1}*get_rotation_mat(A{i});
    end
    
    for i = 1:n_dof
        c = sequence(i);
        
        if c ~= 'r' && c ~= 'p'
            disp("Sequence must be composed of only 'r's and 'p's")
            geoJ = -1;
            return
        end
        
        %if the joint is prismatic then the element in Ja is equals to zero
        %otherwise:
        %   Ja(:, i) = z_(i-1) = R0_1 * ... * R(i-2)_(i-1)*[0; 0; 1]
        if c == 'r'
            if i == 1
                Ja(:, i) = [0; 0; 1];
            else
                Ja(:, i) = cells{i-1}*[0; 0; 1];
            end
        end
    end
    Jl = simplify(Jl);
    Ja = simplify(Ja);
    %geoJ = simplify([Jl; Ja]);

end

%DH-Matrix
function [T, A] = DHMatrix(arrays)
% T = DHMatrix(arrays) takes as inputs:
%   -arrays: a n-vector of vectors composed like this: [alpha a d theta]
% and outputs:
%   -T: the product of all the matrices corresponding to each vector of arrays
% Remember that:
% cos(q1 + q2) = cos(q1)*cos(q2) - sin(q1)*sin(q2)
% sin(q1 + q2) = cos(q1)*sin(q2) + cos(q2)*sin(q1)
% making use of the simplify function these are converted automatically
    T = eye(4);
    nums = size(arrays);
    
    A = cell(1,nums(1));
    
    for i = 1:nums(1)
        line = arrays(i, :);
        R = [cos(line(4)) -cos(line(1))*sin(line(4)) sin(line(1))*sin(line(4)) line(2)*cos(line(4));
             sin(line(4)) cos(line(1))*cos(line(4)) -sin(line(1))*cos(line(4)) line(2)*sin(line(4));
             0 sin(line(1)) cos(line(1)) line(3);
             0 0 0 1;];
        A{i} = R;
        T = T * R;
    end

    if isa(T, 'sym')
        T = simplify(T);
    end
end

%From DH-Matricex get rotation matricex
function R = get_rotation_mat(A)
    R = A(1:3, 1:3);
end