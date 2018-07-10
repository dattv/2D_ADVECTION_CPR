function dF = AdvArtDiffRHS(uc,a,b,Dr,Ds,nx,ny,Fscale,LIFT,map,nN,nE,npf,rx,ry,sx,sy,invV,V,h,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Compute the Advection RHS for using NDG or FR/CPR schemes
%
%                    RHS = -f(u)_x + e(x)*u_xx
%         where f = f(x) is our Corrected/Continuous Flux function
%
%              coded by Manuel Diaz, NTU, 2014.09.16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] Persson, Per-Olof, and Jaime Peraire. "Sub-cell shock capturing for
%     discontinuous Galerkin methods." AIAA paper 112 (2006): 2006. 
% [2] Kl?ckner, Andreas, Tim Warburton, and Jan S. Hesthaven. "Viscous
%     shock capturing in a time-explicit discontinuous Galerkin method."
%     Mathematical Modelling of Natural Phenomena 6.03 (2011): 57-83.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Artificial Diffusion per element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arrange input data
u=zeros(nN,nE); u(:,:)=uc; 

% ut: modal values
ut=invV*u; %ut=V\u; 

% total elements
n = (P+1)*(P+1);

% u_hat: truncated ut
ut(P+1:P+1:n,:)=0; ut((P+1)*P+1:n-1,:)=0; u_hat=V*ut;

% Smooth/Resolution Indicator
s = log10(dot((u-u_hat),(u-u_hat))./dot(u,u));

% clean 's' (NaN or -Inf values are considered smooth values)
s(s==-Inf)=-6;

% Parameters
k = 1.5;      % $\kappa$: empirically choosen to obtain a sharp profile
s0 = -4.0;    % $s_0$
eps0 = 0.05*(h/P)*max(abs([a(:);b(:)])); % $epsilon_0$

% Ranges
range1 = s<(s0-k);
range2 = (s0-k)<s & s<(s0+k);
range3 = s>(s0+k);

% epsilon
k = ones(n,1)*(0*range1 + eps0/2*(1+sin(pi*(s-s0)/(2*k))).*range2 + eps0.*range3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 1st Derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Jumps at elements faces 
du = zeros(npf,nE); du(:) = u(map.vM)-u(map.vP);
da = zeros(npf,nE); da(:) = a(map.vM);
db = zeros(npf,nE); db(:) = b(map.vM);

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

% Jump at elements faces by LDG central flux
dux = nx.*du/2;
duy = ny.*du/2;

% Compute Gradients of u,
dudx = rx.*(Dr*u) + sx.*(Ds*u); 
dudy = ry.*(Dr*u) + sy.*(Ds*u);

% Compute auxiliary variable q,
qx = k.*(dudx - LIFT*( Fscale.*dux ));
qy = k.*(dudy - LIFT*( Fscale.*duy ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2nd Derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Jumps at elements faces by central flux
dqx = zeros(npf,nE); dqx(:) = qx(map.vM)-qx(map.vP);
dqy = zeros(npf,nE); dqy(:) = qy(map.vM)-qy(map.vP);

% Apply BCs
switch BCtype
    case 'Periodic' % Apply Periodic BCs in 'q'
        dqx(map.L) = qx(map.vL)-qx(map.vR);
        dqx(map.R) = qx(map.vR)-qx(map.vL);
        dqx(map.U) = qx(map.vU)-qx(map.vD);
        dqx(map.D) = qx(map.vD)-qx(map.vU);
        dqy(map.L) = qy(map.vL)-qy(map.vR);
        dqy(map.R) = qy(map.vR)-qy(map.vL);
        dqy(map.U) = qy(map.vU)-qy(map.vD);
        dqy(map.D) = qy(map.vD)-qy(map.vU);
    case 'Neumann' % Apply Dirichlet BCs in 'q'
        qxL = 0.0; dqx(map.L) = qx(map.vL)-qxL;
        qxR = 0.0; dqx(map.R) = qx(map.vR)-qxR;
        qxD = 0.0; dqx(map.D) = qx(map.vD)-qxD;
        qxU = 0.0; dqx(map.U) = qx(map.vU)-qxU;
        qyL = 0.0; dqy(map.L) = qy(map.vL)-qyL;
        qyR = 0.0; dqy(map.R) = qy(map.vR)-qyR;
        qyD = 0.0; dqy(map.D) = qy(map.vD)-qyD;
        qyU = 0.0; dqy(map.U) = qy(map.vU)-qyU;
    case 'Dirichlet' % Apply Neumann BCs in 'q'
        qxL = qx(map.vL); dqx(map.L) = qx(map.vL)-qxL;
        qxR = qx(map.vR); dqx(map.R) = qx(map.vR)-qxR;
        qxD = qx(map.vD); dqx(map.D) = qx(map.vD)-qxD;
        qxU = qx(map.vU); dqx(map.U) = qx(map.vU)-qxU;
        qyL = qy(map.vL); dqy(map.L) = qy(map.vL)-qyL;
        qyR = qy(map.vR); dqy(map.R) = qy(map.vR)-qyR;
        qyD = qy(map.vD); dqy(map.D) = qy(map.vD)-qyD;
        qyU = qy(map.vU); dqy(map.U) = qy(map.vU)-qyU;
end

% Jumps at elements faces by Lax-Friedrichs flux
flux = nx.*(da.*du-dqx)/2-abs(da.*nx).*du/2 + ...
      ny.*(db.*du-dqy)/2-abs(db.*ny).*du/2;
   
% Compute Divergence of q
Divq = rx.*(Dr*(a.*u-qx)) + sx.*(Ds*(a.*u-qx)) + ...
       ry.*(Dr*(b.*u-qy)) + sy.*(Ds*(b.*u-qy));
   
% Compute RHS for the semi-discrete PDE
dF = -Divq + LIFT*( Fscale.*flux );
end