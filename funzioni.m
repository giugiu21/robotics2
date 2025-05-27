% Nome del file: Funzioni.m
classdef funzioni
    methods (Static)
        %Building DH matrix 
        % parameters for the array (alphai, ai, di, thetai)
        function T = DHMatrix(arrays)
            T = eye(4);
            nums = size(arrays);
            
            for i = 1:nums(1)
                line = arrays(i, :);
                R = [cos(line(4)) -cos(line(1))*sin(line(4)) sin(line(1))*sin(line(4)) line(2)*cos(line(4));
                     sin(line(4)) cos(line(1))*cos(line(4)) -sin(line(1))*cos(line(4)) line(2)*sin(line(4));
                     0 sin(line(1)) cos(line(1)) line(3);
                     0 0 0 1;];
                T = T * R;   
            end
        
            if isa(T, 'sym')
                T = simplify(T);
            end
        end

        %Axis-angle inverse probelm
        function [rfinale,rfinale1]=inverseProblemR(R)
           seno1=real((1/2)*sqrt((R(1,2)-R(2,1))^2+(R(1,3)-R(3,1))^2+(R(2,3)-R(3,2))^2));
           seno2=real(-(1/2)*sqrt((R(1,2)-R(2,1))^2+(R(1,3)-R(3,1))^2+(R(2,3)-R(3,2))^2));
           cose=(trace(R)-1)/2;
           teta=atan2(seno1,cose);
           teta1=atan2(seno2,cose);
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
                costheta= (R(1,1)+R(2,2)+R(3,3)-1)/2;
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
                costheta= (R(1,1)+R(2,2)+R(3,3)-1)/2;
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
         
            R = funzioni.elem_rot_mat(sequence(1), angles(1)) * funzioni.elem_rot_mat(sequence(2), angles(2)) * funzioni.elem_rot_mat(sequence(3), angles(3));
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
            R = funzioni.euler_rotation(flip(sequence), flip(angles));
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
            [psi, theta, phi] = funzioni.euler_rotation_inverse(flip(sequence), R);
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
    end
end


