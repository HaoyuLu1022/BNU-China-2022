P=[0.01,0.091,0,0,0,0.05,0.01,0.001,3];
[t,y]=ode45(@F,[0,300],P);
% plot(t,y(:,9),y(:,6), y(:, 7), y(:, 8), '-.');
plot(t, y(:, 6));
hold on;
plot(t, y(:, 7));
hold on;
plot(t, y(:, 8));
hold on;
plot(t, y(:, 9));
hold on;

function ydot=F(t,y)
D=2;v=200/1.1*10^(-9);dae=0.057;
ka=0.005 ;da=0.057 ;miu=0.011 ;k1=0.003;k_1=0.001 ;
ci=1;ki=1.03 ;pi=0.1 ;dmi=0.247 ;di=0.005 ;
cr=1;kr=0.79 ;pr=0.1 ;dmr=0.247 ;dr=0.2 ;
m1=4;km1=1;km2=1.1;km3=0.8;m2=1;m3=1.5;
cg=1;kg=0.5;pg=0.1;dmg=0.247;a1=0.3;b1=60;dg=0.01 ;
cb=1;kb=0.5;pb=0.1;dmb=0.247;a2=0.5;b2=35;db=0.01 ;
ct=1;kt=0.5;pt=0.1;dmt=0.247;a3=0.5;b3=35;dt=0.01 ;
N0=3;k=1/6;G0=0.05;n1=4;n2=1;n3=1.5;b0=001;t0=0.001;
ydot=zeros(9,1);
ydot(1)=D*y(9)*(-v*y(1)+y(2))-dae*y(1);
ydot(2)=D*(v*y(1)-y(2))+ka*y(3)-(da+miu)*y(2)-k1*y(4)*y(2)+k_1*y(5);
ydot(3)=ci*ki*pi/(dmi+miu)-(di+miu)*y(3);
ydot(4)=cr*kr*pr/(dmr+miu)-(dr+miu)*y(4)-k1*y(4)*y(2)+k_1*y(5);
ydot(5)=k1*y(4)*y(2)-k'*y(5)-y(5)*y(5)^m1/(km2+y(5)^m1)-y(5)*y(5)^m2/(km2+y(5)^m2)-y(5)*y(5)^m3/(km3+y(5)^m3);
ydot(6)=cg*kg*pg/(dmg+miu)*(a1-(b1-a1)*y(5)^m1)/(km1+y(5)^m1)-(dg+miu)*y(6);
ydot(7)=cb*kb*pb/(dmb+miu)*(a2+(b2-a2)*y(5)^m2/(km2+y(5)^m2))-(db+miu)*y(7);
ydot(8)=ct*kt*pt/(dmt+miu)*(a3+(b3-a3)*y(5)^m3/(km3+y(5)^m3))-(dt+miu)*y(8);
ydot(9)=miu*y(9)*(N0-y(9))+k*y(9)*(y(6)^n1/(y(6)^n1+G0^n1)+y(7)^n2/(y(7)^n2+b0^n2)+y(8)^n3/(y(8)^n3+t0^n3));
end