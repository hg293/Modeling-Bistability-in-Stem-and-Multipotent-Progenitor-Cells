clear all
clc;
clf;
gcf;

%defining model paraneters
a1 = 1;
a2 = 1;
b1 = 1;
b2 = 1;
th1 = 0.5;
th2 = 0.5;
k = 1;
n = 1;

%To find the solution to the differential equations, we have implemented
%ode45 to find the curve of the two genes. They show the "activity" of the
%two proteins
tspan = linspace(-100,100); % define the time span to plot the protein activity
f = @(t,x) [a1*(x(1)^n / (th1^n + x(1)^n)) + b1*(th2^n/ (th2^n + (x(2)^n))) - k*x(1);
    a2*(x(2)^n / (th1^n + x(2)^n)) + b2*(th2^n/ (th2^n + (x(1)^n))) - k*x(2)];
[t,xa] = ode45(f,tspan,[0 0]); %using ode45 to find the solution of the dynamic system with initial values=0

%%Set-up figure and axes
figure(1);
subplot(2,1,1);
plot(xa(:,1),'color','blue','marker','.');
xlabel ('t');
ylabel ('x = GATA1'); 
title(['Plot of the solution of ODE system for a1=',num2str(a1) ' a2= ',num2str(a2) ' b1=',num2str(b1) ' b2=',num2str(b2) ' \theta1=',num2str(th1) ' \theta2=',num2str(th2) ' k=',num2str(k) ' n=',num2str(n)]);
subplot(2,1,2);
plot(xa(:,2),'color','blue','marker','.');
xlabel ('t');
ylabel ('y = PU.1 ');

%Now we solve the two ODEs by equating it to 0 and finding the jacobian,
%eigenvalues and eigenvectors
syms t y1 y2
F = f(t,[y1;y2]) %initialise the matrix
[y1s,y2s] = solve(F(1),F(2),y1,y2)
figure(2);
plot(y1s,'color','blue','marker','o')
hold on
plot(y2s,'color','red','marker','o')

jac = jacobian(F,[y1;y2]) %finds the jacobian
A = subs(jac,{y1,y2},{y1s(2),y2s(2)}) %to find the eigen values

[eigvect,eigval] = eig(A) %result of eigen values and eigen vectors of the dynamic system
