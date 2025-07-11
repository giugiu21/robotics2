% Nome del file: Funzioni2.m
classdef funzioni2
    methods (Static)

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
            Rdot=diff(R,a)*dalfa+diff(R,b)*dbeta+diff(R,g)*dgamma;
            S_omega=simplify(Rdot*transpose(R));
            omega=[S_omega(3,2);S_omega(1,3);S_omega(2,1)];
        end

        %Newton Method
        function [q_out, guesses, cartesian_errors] = newton_method(q_in, desired_point, f_r, initial_guess, max_iterations, max_cartesian_error)
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
    end
end

    

