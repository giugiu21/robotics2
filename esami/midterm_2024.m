
% Robotics 2 Midterm
% 24 April 2024
% Ex #1
%
% Projected Gradient (PG) method (with obstacle avoidance)
% and Task Priority (TP) method (with two variants)
% for a 3R planar robot
clear all, clc

disp('*** Projected Gradient and Task Priority (TP) methods for a 3R planar robot ***')
disp('' )
syms q1 q2 q3 real
q=[q1 q2 q3]';

disp('kinematics of the 3R planar robot (links of unitary length)')
pesimb=[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);sin(q1)+sin(q1+q2)+sin(q1+q2+q3)] 
pmsimb=[cos(q1)+cos(q1+q2);sin(q1)+sin(q1+q2)]
Jesimb=jacobian(pesimb,q) 
Jmsimb=jacobian(pmsimb,q)

disp('at the given configuration (as in the text figure)') 
q0=[0;pi/2;-pi/2]
disp('numerical kinematics')
pe=subs(pesimb,q,q0)
pm=subs(pmsimb,q,q0)
Je=subs(Jesimb,q,q0) 
Jm=subs(Jmsimb,q,q0)
disp('clearance')

C=[0;2];r=0.5;
H=norm(pm-C)-r

disp('desired end-effector velocity')
ve=[0 1]'

disp('projected gradient (PG) method')
nablaH=Jm'*(pm-C)/norm(pm-C)
nablaH=eval(nablaH)
alfa1=1 % works also for a different alfa1>0
Jepinv=pinv(Je)
dqPG=alfa1*nablaH+Jepinv*(ve-alfa1*Je*nablaH) 
dqPG=eval(dqPG)

disp('check PG solution')
vecheck=eval(Je*dqPG) 
vm=eval(Jm*dqPG)

%TPwithtwocases: A&B disp(‘Task Priority (TP) - case A’)
Pe=eye(3)-Jepinv*Je
vmA=vm
dqTPA=Jepinv*ve+pinv(Jm*Pe)*(vmA-Jm*Jepinv*ve) 
dqTPA=eval(dqTPA)

disp('check TP A solution')
vecheck=eval(Je*dqTPA) 
vmcheck=eval(Jm*dqTPA)

disp('Task Priority (TP) - case B')
alfa2=1 % minimum error = vm-Jm*dq TP B is obtained with alfa2=0.9 vm B=alfa2*(1-(r/norm(pm-C)))*(pm-C)
vmB=alfa2*(1-(r/norm(pm-C)))*(pm-C)
vmB=eval(vmB)
dqTPB=Jepinv*ve+pinv(Jm*Pe)*(vmB-Jm*Jepinv*ve)
eval(dqTPB)
