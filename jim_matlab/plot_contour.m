clear all;

dir1 = 'output-09Mar2023-163802';
file1 = [dir1 '/final.dat'];
%file1 = [dir1 '/it.10N.dat'];

A = importdata(file1);
x = A(:,2);
y = A(:,3);
q = A(:,4);

%Triangulate data
tri = delaunay(x,y);

%Analytical Solution
nu = 0.01;
u = 0.8; v = 0.8;
t = 1;
xc = -0.5;
yc = -0.5;


pre1 = 1.0./(4.*t + 1.0);
num = (x - u.*t - xc).^2 + (y - v.*t - yc).^2;
fac1 = num.*pre1./nu;
q_a = pre1.*exp(-fac1);



levels = 0:.02:0.2;
figure(1)
tricontf(x,y,tri,q,levels)
colorbar
set(gca,'FontSize',22)
xlabel('x')
ylabel('y')
title('t = 1 numerical')


figure(2)
tricontf(x,y,tri,q_a,levels)
colorbar
set(gca,'FontSize',22)
xlabel('x')
ylabel('y')
title('t = 1 analytical (5.123) in Kopriva')
print -dpng analytical_t1.png

diff1 = q -q_a;
figure(3)
tricontf(x,y,tri,diff1)
colorbar
set(gca,'FontSize',22)
xlabel('x')
ylabel('y')
title('q - q_a')
print -dpng diff1.png


rel_err = norm(q-q_a)./norm(q_a)



%title('final time')