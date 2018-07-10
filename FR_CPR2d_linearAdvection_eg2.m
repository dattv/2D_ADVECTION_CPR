%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Solving 2-D Linear Advection with NDG or FR/CPR + Art Diffusion
%
%          du/dt + dF/dx + dG/dy = k*( d^2u/dx^2 + d^2u/dy^2 ),  
%                                   for (x,y) \in [a,b]x[c,d]
%                 where u = u(x,y)
%                       F = F(u) = Ux(x,y)*u : linear
%                       G = G(u) = Uy(x,y)*u : linear
%
%                   coded by Manuel Diaz, NTU, 2013.08.30
%                       last modified, NTU, 2015.06.05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] J.Hestaven & T. Warburton; Nodal Discontinuous Galerkin Methods,
%     Algorithms, Analysis, and Applications. Text in Applied Mathematics,
%     (2008)
% [2] D. Wirasaet, Et al.; A performance comparison of nodal
%     discontinuous galerkin methods on triangles and quadrilaterals. Int.
%     J. Numer. Fluids (2010) DOI: 10.1002/fld.2376 .
% [3] A flux reconstruction approach to high-order schemes including
%     Discontinuous Galerkin methods. H.T. Huynh, AIAA 2007.
% [4] A New Class of High-Order Energy Stable Flux Reconstruction Schemes.
%     P.E. Vincent, P. Castonguay, A. Jameson, JSC 2010.
% [5] High-Order Methods for Diffusion Equation with Energy Stable Flux
%     Reconstruction Scheme. K. Ou, P. Vincent, A Jameson, AIAA 2011-46.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Basic Implementation without RK intergration scheme over a
% square domain of size [0,1]x[0,1] using Quadrilateral elements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% Parameters
cfl = 0.5; % CFL condition (~1.6k iterations)
tEnd= 0.50; % final time
  P = 06;  	% degree of accuracy
nEx = 64;   % number of elements in x direction
nEy = 64;   % number of elements in y direction

%% PreProcess
% Advection coefs: Rigid body rotation (washing machine problem ;D)
Ux = @(x,y)  2*pi/tEnd*(y-0.5);
Uy = @(x,y) -2*pi/tEnd*(x-0.5);

% Add source Term ?
sourcefun = 'dont'; 
switch sourcefun
    case 'add'
        S = @(w) 0.01*w.^2;
    case 'dont'
        S = @(w) zeros(size(w));
end

% Add artificial viscosity?
ArtVisc = 'add';

% Build 1d mesh (LGL only!)
mesh = DGmeshQuad('Q4',[0,1],[0,1],nEx,nEy,P); x=mesh.x; y=mesh.y;
dx = max(max(x(:,2:end)-x(:,1:end-1)));
dy = max(max(y(:,2:end)-y(:,1:end-1))); h = max(dx,dy);
data.nE = mesh.nE; data.npf = mesh.eNodesPerFace*mesh.eFaces;
nE = mesh.nE; npf = mesh.eNodesPerFace*mesh.eFaces; nN = mesh.eNodes;

% Create comunication and boundary maps
map = DGmeshQuad.BuildMaps2Dv2(mesh);

% Load DG tools
NDG = DGtoolsQuad(mesh);
    V=NDG.V2D; invV=inv(V); J=NDG.J;
    Dr = NDG.D_kxi;   Ds = NDG.D_eta;
    rx = NDG.kxi_x;   sx = NDG.eta_x;
    ry = NDG.kxi_y;   sy = NDG.eta_y;

% Build surface normals
[nx,ny,sJ] = DGtoolsQuad.Normals(NDG); 
Fscale = sJ./J(mesh.FaceMask(:),:);

% Select Solution Method
method = 'FR'; disp(method);
switch method
    case 'NDG'
        % Build Lift Operator
        [LIFT,Emat] = DGtoolsQuad.lift(NDG);
    case {'FR','CPR','ESFR'}
        dg = CorrectionPolynomial('DGRight',P,mesh.solutionPoints);
        % Build Lift Operator
        LIFT = DGtoolsQuad.liftFR(mesh,dg);
end

% Set velocity Functions
a=Ux(x,y); b=Uy(x,y);

% Set initial time step
dt0=cfl*min([dx./abs(a(:));dy./abs(b(:))]);

% load IC
u0 = IC2d(x,y,2);

%% visualize IC mesh
figure; DGmeshQuad.plotVertices(mesh.EtoV,mesh.VX,mesh.VY);
figure; DGmeshQuad.plotmesh(y,x,u0);

%% Solver Loop

% Low storage Runge-Kutta coefficients
rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];

% Using a 4rd Order 5-stage SSPRK time integration
res_u = zeros(size(u0)); % Runge-Kutta residual storage

% Set initial time & load IC
t=0; u=u0; it=0; dt=dt0;

tic
while t < tEnd
    % Correction for final time step
    if t+dt>tEnd, dt=tEnd-t; end
    
    % iteration counter
    it=it+1;
    
    % RK45 scheme
    for RKs = 1:5
        t_local = t + rk4c(RKs)*dt;
        switch ArtVisc
            case 'add'
                RHS = AdvArtDiffRHS(u,a,b,Dr,Ds,nx,ny,Fscale,LIFT,map,nN,nE,npf,rx,ry,sx,sy,invV,V,h,P);
            case 'dont'
                RHS = AdvectionRHS(u,a,b,Dr,Ds,nx,ny,Fscale,LIFT,map,nN,nE,npf,rx,ry,sx,sy);
        end
        res_u = rk4a(RKs)*res_u + dt*RHS;
        u = u + rk4b(RKs)*res_u;
        % Apply Positivity limiter
        u = PPlimiter(u,invV,nN,nE);
    end
    
    % update time
    t=t+dt;
end
CPUtime=toc; disp(toc);
 
%% plot solution
figure; DGmeshQuad.plotmesh(y,x,u);
xlabel('x'); ylabel('y'); zlabel('u'); title([method,' with Artificial Diffusion']);