
%% Plot Two Orbits
%Initialize Parameters

t1=linspace(0,2*pi);
t2=linspace(0,2*pi);

A1=10;
P1=2;
phi1=pi/8;

A2=4;
P2=1;
phi2=-pi/7;

E1=(P1-A1)/2;
E2=((P1+A1)/2)*cos(t1);
E3=sqrt(P1*A1)*sin(t1);

F1=(P2-A2)/2;
F2=((P2+A2)/2)*cos(t2);
F3=sqrt(P2*A2)*sin(t2);

x1=(cos(phi1)*(E1+E2))+(sin(phi1)*E3);
y1=(-sin(phi1)*(E1+E2))+(cos(phi1)*E3);

x2=(cos(phi2)*(F1+F2))+(sin(phi2)*F3);
y2=(-sin(phi2)*(F1+F2))+(cos(phi2)*F3);

%Plot Trajectory of Orbits
figure(1)
plot(x1,y1)
hold on
plot(x2,y2)
title('Figure 1');
xlabel('x');
ylabel('y');

%% Plot Contour Lines
[t1,t2]=meshgrid(0:2*pi/100:2*pi, 0:2*pi/100:2*pi); 

A1=10;
P1=2;
phi1=pi/8;

A2=4;
P2=1;
phi2=-pi/7;

E1=(P1-A1)/2;
E2=((P1+A1)/2)*cos(t1);
E3=sqrt(P1*A1)*sin(t1);

F1=(P2-A2)/2;
F2=((P2+A2)/2)*cos(t2);
F3=sqrt(P2*A2)*sin(t2);

x1=(cos(phi1)*(E1+E2))+(sin(phi1)*E3);
y1=(-sin(phi1)*(E1+E2))+(cos(phi1)*E3);

x2=(cos(phi2)*(F1+F2))+(sin(phi2)*F3);
y2=(-sin(phi2)*(F1+F2))+(cos(phi2)*F3);

dis=1/2*((x1-x2).^2+(y1-y2).^2);

figure(2)
contour(t1,t2, dis)
title('Figure 2');
xlabel('t1');
ylabel('t2');

figure(3)
surf(t1,t2, dis)
title('Figure 3');
xlabel('t1');
ylabel('t2');
zlabel('dis');

%% Using the Method of Steepest Descent

epsilon=.005;
t1=4;
t1T=t1;
t2=2;
t2T=t2;
k=1;

x=t1;
y=t2;


E1=(P1-A1)/2;
E2=((P1+A1)/2)*cos(t1);
E3=sqrt(P1*A1)*sin(t1);

F1=(P2-A2)/2;
F2=((P2+A2)/2)*cos(t2);
F3=sqrt(P2*A2)*sin(t2);

x1=(cos(phi1)*(E1+E2))+(sin(phi1)*E3);
y1=(-sin(phi1)*(E1+E2))+(cos(phi1)*E3);

x2=(cos(phi2)*(F1+F2))+(sin(phi2)*F3);
y2=(-sin(phi2)*(F1+F2))+(cos(phi2)*F3);

dis=1/2*((x1-x2).^2+(y1-y2).^2);
disT=dis;

x1prime=(cos(phi1)*((P1+A1)/2)*(-sin(t1)))+(sin(phi1)*(sqrt(P1*A1))*cos(t1));
y1prime=(-sin(phi1)*((P1+A1)/2)*(-sin(t1)))+(cos(phi1)*(sqrt(P1*A1))*cos(t1));
x2prime=(cos(phi2)*((P2+A2)/2)*(-sin(t2)))+(sin(phi2)*(sqrt(P2*A2))*cos(t2));
y2prime=(-sin(phi2)*((P2+A2)/2)*(-sin(t2)))+(cos(phi2)*(sqrt(P2*A2))*cos(t2));

dt1=(x1-x2)*x1prime+(y1-y2)*y1prime;
dt2=(x1-x2)*(-x2prime)+(y1-y2)*(-y2prime);

gradnorm=sqrt(dt1.^2+dt2.^2);

lambda=linspace(max(.05,(t1-2*pi)/dt1), max(.05, t1/dt1));

while gradnorm>=epsilon
    while k<=100
    t1=t1-lambda(k)*dt1;
    t2=t2-lambda(k)*dt2;
    
    E2=((P1+A1)/2)*cos(t1);
    E3=sqrt(P1*A1)*sin(t1);
    
    F2=((P2+A2)/2)*cos(t2);
    F3=sqrt(P2*A2)*sin(t2);

    x1=(cos(phi1)*(E1+E2))+(sin(phi1)*E3);
    y1=(-sin(phi1)*(E1+E2))+(cos(phi1)*E3);

    x2=(cos(phi2)*(F1+F2))+(sin(phi2)*F3);
    y2=(-sin(phi2)*(F1+F2))+(cos(phi2)*F3);

    dis=1/2*((x1-x2).^2+(y1-y2).^2);
    
    if dis<disT
        disT=dis;
        lambdaT=lambda(k);
        
    end
    
    t1=t1T;
    t2=t2T;
    k=k+1;
    end
    
    t1T=t1-lambdaT*dt1;
    t2T=t2-lambdaT*dt2;
    
    t1=t1T;
    t2=t2T;
    
    E2=((P1+A1)/2)*cos(t1);
    E3=sqrt(P1*A1)*sin(t1);
    
    F2=((P2+A2)/2)*cos(t2);
    F3=sqrt(P2*A2)*sin(t2);

    x1=(cos(phi1)*(E1+E2))+(sin(phi1)*E3);
    y1=(-sin(phi1)*(E1+E2))+(cos(phi1)*E3);

    x2=(cos(phi2)*(F1+F2))+(sin(phi2)*F3);
    y2=(-sin(phi2)*(F1+F2))+(cos(phi2)*F3);
    
    x1prime=(cos(phi1)*((P1+A1)/2)*(-sin(t1)))+(sin(phi1)*(sqrt(P1*A1))*cos(t1));
    y1prime=(-sin(phi1)*((P1+A1)/2)*(-sin(t1)))+(cos(phi1)*(sqrt(P1*A1))*cos(t1));
    x2prime=(cos(phi2)*((P2+A2)/2)*(-sin(t2)))+(sin(phi2)*(sqrt(P2*A2))*cos(t2));
    y2prime=(-sin(phi2)*((P2+A2)/2)*(-sin(t2)))+(cos(phi2)*(sqrt(P2*A2))*cos(t2));

    dt1=(x1-x2)*x1prime+(y1-y2)*y1prime;
    dt2=(x1-x2)*(-x2prime)+(y1-y2)*(-y2prime);
    
    gradnorm=sqrt(dt1.^2+dt2.^2);
    
    lambda=linspace(max(.05,(t1-2*pi)/dt1), max(.05, t1/dt1));
    k=1;
   
    figure(4)
    tp1=mod(t1, 2*pi);
    tp2=mod(t2, 2*pi);
    
    plot(tp1, tp2, "o")
    hold on;
    title('Figure 4');
    xlabel('t1');
    ylabel('t2');
    
    x = [x; tp1];
    y = [y; tp2];
   
end

[t1,t2]=meshgrid(0:2*pi/100:2*pi, 0:2*pi/100:2*pi); 

A1=10;
P1=2;
phi1=pi/8;

A2=4;
P2=1;
phi2=-pi/7;

E1=(P1-A1)/2;
E2=((P1+A1)/2)*cos(t1);
E3=sqrt(P1*A1)*sin(t1);

F1=(P2-A2)/2;
F2=((P2+A2)/2)*cos(t2);
F3=sqrt(P2*A2)*sin(t2);

x1=(cos(phi1)*(E1+E2))+(sin(phi1)*E3);
y1=(-sin(phi1)*(E1+E2))+(cos(phi1)*E3);

x2=(cos(phi2)*(F1+F2))+(sin(phi2)*F3);
y2=(-sin(phi2)*(F1+F2))+(cos(phi2)*F3);

dis=1/2*((x1-x2).^2+(y1-y2).^2);

figure(5)
contour(t1,t2, dis)
title('Figure 5');
xlabel('t1');
ylabel('t2');
hold on;
plot(x,y)

figure (6)
t1=linspace(0,2*pi);
t2=linspace(0,2*pi);

A1=10;
P1=2;
phi1=pi/8;

A2=4;
P2=1;
phi2=-pi/7;

E1=(P1-A1)/2;
E2=((P1+A1)/2)*cos(t1);
E3=sqrt(P1*A1)*sin(t1);

F1=(P2-A2)/2;
F2=((P2+A2)/2)*cos(t2);
F3=sqrt(P2*A2)*sin(t2);

x1=(cos(phi1)*(E1+E2))+(sin(phi1)*E3);
y1=(-sin(phi1)*(E1+E2))+(cos(phi1)*E3);

x2=(cos(phi2)*(F1+F2))+(sin(phi2)*F3);
y2=(-sin(phi2)*(F1+F2))+(cos(phi2)*F3);

%% Plot Trajectory of Orbits

figure(6)
plot(x1,y1)
hold on
plot(x2,y2)
title('Figure 6');
xlabel('x');
ylabel('y');

A1=10;
P1=2;
phi1=pi/8;

A2=4;
P2=1;
phi2=-pi/7;

E1=(P1-A1)/2;
E2=((P1+A1)/2)*cos(x);
E3=sqrt(P1*A1)*sin(x);

F1=(P2-A2)/2;
F2=((P2+A2)/2)*cos(y);
F3=sqrt(P2*A2)*sin(y);

x1=(cos(phi1)*(E1+E2))+(sin(phi1)*E3);
y1=(-sin(phi1)*(E1+E2))+(cos(phi1)*E3);

x2=(cos(phi2)*(F1+F2))+(sin(phi2)*F3);
y2=(-sin(phi2)*(F1+F2))+(cos(phi2)*F3);

plot(x1,y1, "o")

plot(x2,y2, "o")

  






