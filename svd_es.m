%definining the jacobian matrix
J = [1 2 3; 2 4 6];

%SVD
[U, S, V] = svd(J);
disp('Singular value decomposition of the jacobian:');
disp(U);
disp(S);
disp(V);

%computation of the pseudoinverse using svd
J_pseudo = V * pinv(S) * U';
disp('Jacobian Pseudoinverse using SVD:');
disp(J_pseudo);

%computation of the DLS jacobian using SVD
mi = 0.01; %mi factor >=0
S_dls = diag(diag(S)./ (diag(S).^2 + mi^2));


% Resize S_dls to match J's shape (S is usually m×n)
S_dls_padded = zeros(size(J')); % Create zero matrix of appropriate shape
S_dls_padded(1:size(S_dls,1), 1:size(S_dls,2)) = S_dls; % Place S_dls in top-left


disp('Damped Least-Squares Jacobian Pseudoinverse using SVD:');
J_dls = V * S_dls_padded * U';

disp(J_dls);

