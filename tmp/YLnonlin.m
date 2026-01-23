function [A,hL,pts]=YLnonlin(rho,sigma,g,omegaRPM,a,R,thetacDEG)
%Solves the nonlinear Young-Laplace equation in a circular domain subject
%to a lateral body force in x-direction of magnitude F*g

% all units in cgs

%INPUTS:
% rho == fluid density
% sigma == surface tension
% g == gravity
% omegaRPM == orbital RPM
% a == orbital radius
% R == well radius
% thetacDEG == contact angle (degrees)

%OUTPUTS:
% A == interface area
% hL == interface height difference
% pts == matrix of [x,y,z] coordinates of interface

%example run:
% [A,hL,pts]=YLnonlin(1,72,981,300,1.9/2,2.78/2,100);

lc=sqrt(sigma/rho/g); %capillary length

omega=omegaRPM*2*pi/60; %converts RPM to rad/s
F=a*omega^2/g; %non-dimensional acceleration F

thetac=thetacDEG*pi/180; %converts degrees to radians

model = createpde;

%Geometry definition
gd = decsg([1;0;0;R]); 
%For a circle, the first row contains 1. 
%The second and third rows contain the x- and y-coordinates of the center. 
%The fourth row contains the radius of the circle.

geometryFromEdges(model,gd); %adds the 2-D geometry described in gd to the model container

%Define coefficients for PDE of form: -div(c*grad(u))+a*u=f
a = (1/lc)^2; 
pf=num2str(-F*a);
f = append(pf,'*x');
c = '1./sqrt(1+ux.^2+uy.^2)'; 

%Apply BC for contact angle
applyBoundaryCondition(model,'neumann','Edge',1:model.Geometry.NumEdges,'q',0,'g',tan(pi/2-thetac));

%Generate mesh
generateMesh(model,'GeometricOrder','linear','Hmax',0.02*R);

%Solve PDE
u = pdenonlin(model,c,a,f); 

%Plot Mesh
figure(2)
pdemesh(model); 
axis equal
xlabel('x (cm)')
ylabel('y (cm)')

%Plot Interface
figure(1)
pdeplot(model,'XYData',u,'ZData',u)
axis equal
xlabel('x (cm)')
ylabel('y (cm)')
zlabel('\eta (cm)')

%post-processing
pos=model.Mesh.Nodes; %get positions of all nodes

%interpolate data to regular grid for surfarea function
rv=linspace(0.001,R-0.001,50);
thetav=linspace(0,2*pi,50);
[r,theta]=meshgrid(rv,thetav);
[xq,yq,v]=pol2cart(theta,r,u);
vq = griddata(pos(1,:),pos(2,:),v,xq,yq); 

%compute surface area
A = surfarea(xq,yq,vq);

%compute hL
hL=griddata(pos(1,:),pos(2,:),u,-R+0.001,0)-griddata(pos(1,:),pos(2,:),u,R-0.001,0);

%assemble matrix of x,y,z positions of interface shape
pts=[pos(1,:)' pos(2,:)' u];
%surf2stl('interface',pos(1,:)',pos(2,:)',u) %output stl of interface