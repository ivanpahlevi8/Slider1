%=================%
%System Parameters%
%=================%
L2=0.150 ; % in meter
L3=0.350; % in meter
w2=150; % in rad/s
uPjO1=[0;0];
uPjO2=[-L2/2;0];
uPjA2=[L2/2;0];
uPjA3=[-L3/2;0];
uPjB3=[L3/2;0];
uPjB4=[0;0];

%==============%   
%Initial Values%
%==============%`
q=zeros(12,1);
q6init=pi/6;
%==============%

t=0:0.001:0.04;

NumData=size(t,2);
q_alltime=zeros(12,NumData);
v_alltime=zeros(12,NumData);
a_alltime=zeros(12,NumData);

for j=1:NumData
  
epsilon=1e-12;
delta_q_norm=1;
i_newton=1;

while abs(delta_q_norm)>epsilon

%Transformation Matrix%
A1=[cos(q(3)) -sin(q(3)); sin(q(3)) cos(q(3))];
dA1t1=[-sin(q(3)) -cos(q(3)); cos(q(3)) -sin(q(3))];
A2=[cos(q(6)) -sin(q(6)); sin(q(6)) cos(q(6))];
dA2t2=[-sin(q(6)) -cos(q(6)); cos(q(6)) -sin(q(6))];
A3=[cos(q(9)) -sin(q(9)); sin(q(9)) cos(q(9))];
dA3t3=[-sin(q(9)) -cos(q(9)); cos(q(9)) -sin(q(9))];
A4=[cos(q(12)) -sin(q(12)); sin(q(12)) cos(q(12))];
dA4t4=[-sin(q(12)) -cos(q(12)); cos(q(12)) -sin(q(12))];

%Constraint Equation%
C=zeros(12,1);
C45=[q(1);q(2)]+A1*uPjO1-[q(4);q(5)]-A2*uPjO2;
C67=[q(4);q(5)]+A2*uPjA2-[q(7);q(8)]-A3*uPjA3;
C89=[q(7);q(8)]+A3*uPjB3-[q(10);q(11)]-A4*uPjB4;
C1011=[q(11);q(12)];
C12=q(6)-q6init-w2*t(j);

C=[q(1);q(2);q(3);C45;C67;C89;C1011;C12];
C_norm=sqrt(sum(C(1:12).^2))/12;

%Constrain Jacobian Matrix%
unit22=[1 0; 0 1];

%Ground Constraint%
C0q=[1 zeros(1,11); 0 1 zeros(1,10); 0 0 1 zeros(1,9)];

%Pin Join%
CPjO=[unit22 dA1t1*uPjO1 -unit22 -dA2t2*uPjO2 zeros(2,6)];
CPjA=[zeros(2,3) unit22 dA2t2*uPjA2 -unit22 -dA3t3*uPjA3 zeros(2,3)];
CPjB=[zeros(2,6) unit22 dA3t3*uPjB3 -unit22 -dA4t4*uPjB4];

%Slider%
Cslider=[zeros(1,10) 1 zeros(1,1); zeros(1,11) 1];

%Driving Constraint%
Cdrive=[zeros(1,5) 1 zeros(1,6)];

Cq=[C0q;CPjO;CPjA;CPjB;Cslider;Cdrive];

delta_q=inv(Cq)*(-C);
delta_q_norm=sqrt(sum(delta_q(1:12).^2))/12;

q=q+delta_q;
  
  i_newton=i_newton+1;
 
  if  i_newton > 30 %limit newton rhapson iteration in case non-convergence%
    break
  end  
end

Ct=[zeros(11,1);-w2];
velocity=inv(Cq)*(-Ct);
QdPjO=A1*uPjO1*(velocity(3,1))^2-A2*uPjO2*(velocity(6,1))^2;
QdPjA=A2*uPjA2*(velocity(6,1))^2-A3*uPjA3*(velocity(9,1))^2;
QdPjB=A3*uPjB3*(velocity(9,1))^2-A4*uPjB4*(velocity(12,1))^2;
Qd=[zeros(3,1);QdPjO;QdPjA;QdPjB;zeros(3,1)];
accel=inv(Cq)*Qd;

q_alltime(:,j)=q;
v_alltime(:,j)=velocity;
a_alltime(:,j)=accel;

end

C_norm
delta_q_norm
i_newton

figure(1);
plot(t,q_alltime(10,:));
grid on;
xlabel('Time [s]');
ylabel('slider position [m]');
title('Position of slider vs time');

figure(2);
plot(q_alltime(6,:),q_alltime(10,:));
grid on;
xlabel('crank angle [rad]');
ylabel('slider position [m]');
title('Position of slider vs crank angle');

figure(3);
plot(q_alltime(6,:),v_alltime(10,:));
grid on;
xlabel('crank angle [rad]');
ylabel('slider velocity [m/s]');
title('velocity of slider vs crank angle');

figure(4);
plot(q_alltime(6,:),a_alltime(10,:));
grid on;
xlabel('crank angle [rad]');
ylabel('slider acceleration [m/s^2]');
title('acceleration of slider vs crank angle');

figure(5);
plot(q_alltime(6,:),q_alltime(9,:));
grid on;
xlabel('crank angle [rad]');
ylabel('CR angular disp [rad]');
title('CR angular disp vs crank angle');


