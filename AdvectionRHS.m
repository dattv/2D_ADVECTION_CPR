function RHS = AdvectionRHS(uc,a,b,Dr,Ds,nx,ny,Fscale,LIFT,map,nN,nE,npf,rx,ry,sx,sy)
% Kernel coded by Manuel Diaz. NTU, IAM, 2014.05.20.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] J.Hestaven & T. Warburton; Nodal Discontinuous Galerkin Methods,
%     Algorithms, Analysis, and Applications. Text in Applied Mathematics,
%     (2008)
% [2] D. Wirasaet, Et al.; A performance comparison of nodal
%     discontinuous galerkin methods on triangles and quadrilaterals. Int.
%     J. Numer. Fluids (2010) DOI: 10.1002/fld.2376 .
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(nN,nE); u(:,:)=uc; 

% Jumps at elements faces
du = zeros(npf,nE); du(:) = u(map.vM)-u(map.vP);

% Apply BCs
BCtype = 'Neumann';
switch BCtype
    case 'Periodic' % Apply Periodic BCs: (vM-vP)
        du(map.L) = u(map.vL)-u(map.vR);
        du(map.R) = u(map.vR)-u(map.vL);
        du(map.U) = u(map.vU)-u(map.vD);
        du(map.D) = u(map.vD)-u(map.vU);
    case 'Neumann'  % Apply Neumann BCs
        du(map.L) = 0;
        du(map.R) = 0;
        du(map.U) = 0;
        du(map.D) = 0;
    case 'Dirichlet'% Apply Dirichlet BCs
        uL = 0.1; du(map.L) = u(map.vL)-uL;
        uR = 0.1; du(map.R) = u(map.vR)-uR;
        uU = 0.1; du(map.U) = u(map.vU)-uU;
        uD = 0.1; du(map.D) = u(map.vD)-uD;
end

method = 'LF';
switch method
    case 'NDG' % by Hestaven and Warburton ref.[1]
        alpha = 0; % alpha = 0; upwind,  alpha = 1; central flux.
        du1 = (a*nx-(1-alpha)*abs(a*nx)).*du/2;
        du2 = (b*ny-(1-alpha)*abs(b*ny)).*du/2;
    case 'LF' % by Lax Friedrichs flux splitting
        du1 = nx.*(a*du)/2-abs(a).*du/2;
        du2 = ny.*(b*du)/2-abs(b).*du/2;
end

% compute fluxes f = F(u); g = G(u);
f=a*u; g=b*u;

% Compute Divergence
Div = rx.*(Dr*f) + sx.*(Ds*f) + ry.*(Dr*g) + sy.*(Ds*g);

% Compute RHS for the semi-discrete PDE
RHS = -Div + LIFT*( Fscale.*(du1+du2) );
