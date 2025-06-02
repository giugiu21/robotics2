%dynamic model 2R polar

%?????

clear all

n= 2;

syms q l m dc [n 1] real

syms dq ddq [n 1] real


%if the robot were planar  
%syms Ic1 Ic2 real


%the robot is polar so
syms Ic_xx_1 Ic_yy_1 Ic_zz_1 Ic_xx_2 Ic_yy_2 Ic_zz_2
%i know that Ic_yy_2 = Ic_zz_2
Ic1 = diag([Ic_xx_1, Ic_yy_1, Ic_zz_1]);
Ic2 = diag([Ic_xx_2, Ic_yy_2, Ic_zz_2]);

 
q = [q1; q2];

%%

%In this case all joints are revolute so kinetic energy is composed of
%translational and rotational.

%THE IMPORTANT THING IS THAT THE REFERENCE FRAMES ARE CONCORDANT

% joint 1 
%position of the CoMi in RF0

%for joint 1 the CoM is = to the origin of frame 1 = frame 0
x1 = 0;
y1 = 0;
z1 = 0;
rc_01 = [x1; y1; z1];

%velocity of CoMi
vx1 = diff(x1,q1)*dq1 + diff(x1,q2)*dq2;
vy1 = diff(y1,q1)*dq1 + diff(y1,q2)*dq2;
vz1 = diff(z1,q1)*dq1 + diff(z1,q2)*dq2;
vc1 = [vx1; vy1; vz1];
disp("Velocity CoM")
disp(vc1);

%angular velocity in frame 0
w1 = [0; dq1; 0]; 

disp("Joint 1")
T1_tr = (1/2)* m1 * vc1'*vc1;
T1_rot = (1/2)* w1'* Ic1 * w1;
T1 = simplify(T1_tr + T1_rot);
disp("Kinetic energy")
disp(T1);


% joint 2 
%position of the CoMi in RFi
x2 = dc2*cos(q2)*cos(q1);
y2 = dc2*sin(q2)*sin(q1);
z2 = dc2*sin(q2);
rc_02 = [x2; y2; z2];

%velocity of CoMi
vx2 = diff(x2,q1)*dq1 + diff(x2,q2)*dq2;
vy2 = diff(y2,q1)*dq1 + diff(y2,q2)*dq2;
vz2 = diff(z2,q1)*dq1 + diff(z2,q2)*dq2;
vc2 = [vx2; vy2; vz2];
disp("Velocity CoM")
disp(vc2);

%angular velocity
w2 = [sin(q2)*dq1; cos(q2)*dq1; dq2]; 

disp("Joint 2")
T2_tr = (1/2)* m2 * vc2'* vc2;
T2_rot = (1/2)*  w2' * Ic2 * w2;
T2 = simplify(T2_tr + T2_rot);
disp("Kinetic energy")
disp(T2);


%%
%Total kinetic Energy of the Robot 
T_tot = T1 +T2;
T_tot=simplify(T_tot);
T_tot=collect(T_tot,dq1^2);
T_tot=collect(T_tot,dq2^2);
disp(T_tot)


disp("Inertia matrix for the robot")


%-------------Inertia Matrix--------------
dq = [dq1; dq2];

M=simplify(hessian(T_tot,dq));
disp(M)
