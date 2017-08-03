%% Part 1: Vector objects
% a is a scalar
a=1

% v is a 3x1 vector (also a tensor)
v=[1;2;3]

% D is a 3x3 tensor
D= [1,2,3; 4,6,8; 7,3,2]

%% Note that adding is incompatible due to dimension mismatch
v+D
%% But we can multiply them together because (3x3)(3x1) or (1x3)(3x3) match the inside dimension
D*v
v'*D
%% Two ways of defining a dot product with a vector
v'*v
dot(v,v)
%% Inner product is the same as a dot product for a vector
IP=v'*v
%% Dot product is defined differently in matlab for 3x3 vectors. 
% To get a double dot, we need to sum over the dot product
DD = sum(dot(D,D))
%% Outer product yields a tensor
OP=v*v'

%% Let's create a spatial grid for our fluids
m=10;
n=20;
[x,y]=meshgrid(linspace(0,10,m),linspace(0,10,n));
%% In the case of Poiseuille flow, there pressure decreases linearly in a single dimension
pressure=x;
figure(1);imagesc(pressure);
%% If we solve poiseulle flow, we get a parabolic profile for the x velocity
velocity=zeros(n,m,2); % set an empty vector n x m x 2, the size of our velocity tensor
velocity(:,:,1)=(10-y).*y; % the first component, the x-flow (or u) is parabolic
velocity(:,:,2)=zeros(n,m); % the second component, the y-flow (or v) is zero.
figure(2);quiver(x,y,velocity(:,:,1),velocity(:,:,2)); % a quiver plot is a common way to show a vector field
%% Note that the velocity field is a two component vector field, while the
% pressure field is a scalar field
size(velocity)
size(pressure)
%% Part 2: Fluid mechanics with Poiseuille flow as an example
clear all
%% Let's set up our poiseuille flow problem again.
h=1; % height of our channel
l=5; % length we want to calculate
umax=1; % max flow in the center of the channel
npts=20; % number of points in our grid
dx=l/(npts-1); % spacing of our points
dy=h/(npts-1);
x=NaN(npts+2,npts+2); %setting up an empty tensor of NaN values for easy visualization later one
y=x;
%filling in our x and y values, but leaving an NaN at the edges
[x(2:end-1,2:end-1),y(2:end-1,2:end-1)]=meshgrid(linspace(0,l,npts),linspace(0,h,npts));

u=4*umax/h^2*(h-y).*y;v=x*0; % filling in our x an y flow components, also labelled u, v
figure(3);quiver(x,y,u,v);
%% The stream function is another way of visualizing flow. It is denoted as psi
psi=4*umax/h^2*(h*y.^2/2-y.^3/3);
contour(psi)
%% The velocity field can be derived from the stream function by taking a 
% spatial derivative. Note that it is the same as our plot previously
[psix,psiy]=gradient(psi,dx,dy);
u1=psiy;
v1=-psix;
quiver(x,y,u1,v1);
%% Now let's take the laplacian of our vector field 
% it should yield our pressure gradient, based on the Stokes Equation
dpx=del2(u,dx,dy);
dpy=del2(v,dx,dy);
quiver(x,y,dpx,dpy);
%% When we take the laplacian of the laplacian of our Stream function,
% it should equal zero (within machine computational tolerances
% note that the biharmonic equation is another way to solve for flow in two
% dimensional flow systems
d4=del2(del2(psi,dx,dy),dx,dy);
imagesc(d4)
%% The gradient of each component, u, v, can be taken separately with respect to x and y
% The gradient operation yields 4 different components ux, uy, vx, and vy
[ux,uy]=gradient(u,dx,dy);
[vx,vy]=gradient(v,dx,dy);
%% In incompressible flow, our divergence should also be zero.
divergence = ux+vy;
imagesc(divergence);
%% Our rate of strain tensor is a tensor that incorporates information
% about spatial gradients of u and v. It essentially tells us how quickly
% to fluid is deforming at various locations in various directions
visc=1;
ros = zeros([size(ux),2,2]);
ros(:,:,1,1)=ux;
ros(:,:,2,1)=1/2*(uy+vx);
ros(:,:,1,2)=1/2*(uy+vx);
ros(:,:,2,2)=uy;
figure(1);imagesc(squeeze(ros(:,:,2,1)));
figure(2);imagesc(squeeze(ros(:,:,1,1)));
%% The dissipation function can be calculated from the rate of strain tensor
% The dissipation function is the amount of energy dissipated in viscous
% flow at each location
phi=2*visc*sum(dot(ros,ros,4),3);
%%
figure(1);imagesc(phi)
%% We can calculate the total energy dissipated by the system by summing up our 
% dissipation function over all of our points in space
Energy=nansum(phi(:))*dx*dy
