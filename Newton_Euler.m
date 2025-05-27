%------------set up----------------
clear all

n = 3; %number of joints

syms l dc m q [n 1] real
syms dq [n 1] real
syms ddq [n 1] real


q = [q1; q2; q3];


%DH parameters
alpha = [sym(pi)/2, 0, 0]; %IMPORTANT to make pi symbolic
a = [0, l2, l3];
d = [l1, 0, 0];
theta = [q1, q2, q3];

A = cell(1, n);
A_tot = eye(4);


for i = 1:n
    alphai = alpha(i);
    ai = a(i);
    di = d(i);
    thetai = theta(i);

    A{i} =  [cos(thetai), -sin(thetai)*cos(alphai), sin(thetai)*sin(alphai), ai*cos(thetai);
            sin(thetai), cos(thetai)*cos(alphai), -cos(thetai)*sin(alphai), ai*sin(thetai);
            0, sin(alphai), cos(alphai), di;
            0, 0, 0, 1];

    %complete DH matrix for the robot

    A_tot = simplify(A_tot * A{i}); 

end 

%extracting from the homogeneous matrices the rotation matrices and
%positions for each joint

R1 = A{1}(1:3,1:3);
p1 = A{1}(1:3,4);

R2 = A{2}(1:3,1:3); 
p2 = A{2}(1:3,4);

R3 = A{3}(1:3,1:3); 
p3 = A{3}(1:3,4);

disp("Rotation matrices")
disp(R1)
disp(R2)
disp(R3)

disp("Origin's position in RFi-1")
pause

disp(p1)
disp(p2)
disp(p3)

%position vectors from origin i-1 to origin i in frame i
% ip = i-1Ri' * i-1p !!! i am doing this computation here
% i-1p = i-1Ri * ip

r01 = simplify(R1'*p1);
r12 = simplify(R2'*p2);
r23 = simplify(R3'*p3);

disp("Position vectors from origin i-1 to origin i in frame i")
pause 


disp(r01)
disp(r12)
disp(r23)


%!!!CoM position vector is the vector from the origin of frame i to the CMi
% expressed in frame i
%in this case we assume it to be in the kinematic axis of each link
%if we define di as the distance between origin of frame i to CMi

%Understanding the link axis from DH parameters:
%if alpha_i = 0 the link extends on x-axis of frame i
%if alpha_i = pi/2 the link extends on y_axis of frame i
%if alpha_i = pi the link extends on -x-axis of frame i

% to find the link lengths look at the ai parameter:
%in this case a2 and a3 are our link lengths 

%we define d1 d2 d3 as the distance between origin i-1 and position of CoMi

rc_1 = [0; dc1 - l1; 0];

rc_2= [dc2 - l2; 0; 0];

rc_3= [dc3 - l3; 0; 0];


disp("Position vectors in frame i for CoMi for each link")
disp(rc_1)
disp(rc_2)
disp(rc_3)

disp("CoM position vectors in RF0")
pause

% ip = i-1Ri' * i-1p 
% i-1p = i-1Ri * ip !!i am doing this computation here

rc_01 = simplify((R1*rc_1) + p1);
disp(rc_01)
rc_02 = simplify((R1* (R2*rc_2)) + R1*p2 + p1);
disp(rc_02)
rc_03 = simplify((R1*(R2*(R3*rc_3))) + R1*(R2*p3) + R1*p2 + p1);
disp(rc_03)



%-------------------NEWTON EULER---------------------------

disp("--------Newton Euler Algorithm-------")
pause

%----------forward pass------------

%Defining the sigma as: 0 if revolute joint, 1 if prismatic
s1 = 0;
s2 = 0;
s3 = 0;


%the robot is planar  
syms Ic1 Ic2 Ic3 real

%if the robot werw polar so
%syms Ic_xx_1 Ic_yy_1 Ic_zz_1 Ic_xx_2 Ic_yy_2 Ic_zz_2 Ic_xx_3 Ic_yy_3 Ic_zz_3 real
%Ic1 = diag([Ic_xx_1, Ic_yy_1, Ic_zz_1]);
%Ic2 = diag([Ic_xx_2, Ic_yy_2, Ic_zz_2]);
%Ic3 = diag([Ic_xx_3, Ic_yy_3, Ic_zz_3]);


%set up to change if we have more joints
s={s1, s2, s3};
dq = {dq1, dq2, dq3};
ddq = {ddq1, ddq2, ddq3};

R = {R1, R2, R3};
r = {r01, r12, r23};
rc = {rc_1, rc_2, rc_3};
m = {m1, m2, m3};
Ic = {Ic1, Ic2, Ic3};



%intialization
w_prev = [0; 0; 0];
w_dot_prev = [0; 0; 0];
a_prev = [0; 0; 0];


w_vector = sym(zeros(3,n));
w_dot_vector = sym(zeros(3,n));
a_vector = sym(zeros(3,n));
a_c_vector = sym(zeros(3,n));


% --------- starting forward pass -------------
    for i=1:n
        Ri = R{i};
        dqi = dq{i};
        ddqi = ddq{i};
        ri = r{i};

        w = Ri' * (w_prev + dqi * [0;0;1]);
        w = vpa(simplify(w));
        fprintf("the w of link %d is:",i);
        disp(w)
        w_vector(:,i) = w;

        w_dot = Ri' * (w_dot_prev + ddqi * [0;0;1] + cross(dqi*w_prev,[0;0;1]));
        fprintf("the w_dot of link %d is:",i);
        disp(w_dot)
        w_dot_vector(:,i) = w_dot;

        a = Ri' * a_prev + cross(w_dot,ri) + cross(w,ri);
        a = vpa(simplify(a));
        fprintf("the a of link %d is:",i);
        disp(a);
        a_vector(:,i) = a;

        a_c = a + cross(w_dot,ri) + cross(w,cross(w,ri));
        fprintf("the a_c of link %d is:",i);
        disp(a_c);
        a_c_vector(:,i) = a_c;

        %update the prev
        w_prev = w;
        W_dot_prev = w_dot;
        a_prev = a;
    end

%%

%starting backward pass
    %---------------- initialization --------------------
    f_next=[0;0;0];
    tau_next=[0;0;0];

    i_g = sym(zeros(3,n));

    f_vector = sym(zeros(3,n));
    tau_vector = sym(zeros(3,n));

    R{n+1} = eye(3);

    %--------------- starting backward pass --------------
    for i=n:-1:1

        Ici = Ic{i};
        Rii = R{i+1};
        ri = r{i};
        rci = rc{i};
        mi = m{i};

        if i == n
            f = mi * (a_c_vector(:,i) - i_g(:,i));
        else
            f = Rii * f_next + mi * (a_c_vector(:,i) - i_g(:,i));
        end
        f = vpa(simplify(f));
        fprintf("the f of link %d is:",i);
        disp(f);
        f_vector(:,i) = f;

        if i == n
            tau = -cross(f,(ri + rci)) + Ici * w_dot_vector(:,i) + cross(w_vector(:,i),(Ici * w_vector(:,i)));
        else
            tau = Rii * tau_next + cross((Rii * f_next),rci)-cross(f,(ri+rci))+ Ici * w_dot_vector(:,i)+cross(w_vector(:,i),(Ici * w_vector(:,i)));

        end
        tau=vpa(simplify(tau));
        fprintf("the tau of link %d is:",i);
        disp(tau);
        tau_vector(:,i) = tau;
        
        f_next = f;
        tau_next = tau;
    end

    u_vector=sym(zeros(3,n));

    for i=1:n
        si = s{i};
        Ri = R{i};

        if si == 1 %if the joint is prismatic
            ui = f_vector(:,i)' * (Ri' * [0;0;1]);
        else
            ui = tau_vector(:,i)' * (Ri' * [0;0;1]);
        end

        ui = vpa(simplify(ui));

        fprintf("u of link %d is:",i);
        disp(ui)
        u_vector(:,i) = ui;
    end

