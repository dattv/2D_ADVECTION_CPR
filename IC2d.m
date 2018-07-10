function u0 = IC2d(x,y,option)
%% Initial Condition (IC)
% This subroutine creates an special IC for testing purposes.
% The present IC models a rectangular domain with four equally sized
% regions with diferrent initial state value, u, an adimensional value.
%
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.09.06

%% Select your IC

switch option
    case 1 % 4 quadrants in Domain [0,1]^2
        [nN,nE] = size(x); domain = 'square';
        
        switch domain
            case 'square' % for Quad elements in a square domain!
            
            % Calculate number of nodes and elements
            u0 = zeros(nN,nE);
            
            % Only valid for square domain
            zero = zeros(sqrt(nE)/2); one = ones(sqrt(nE)/2);
            
            % Initial Condition for our 2D domain
            u0(:,[zero,zero;zero,one]==1) = 0.50;
            u0(:,[zero,zero;one,zero]==1) = 0.70;
            u0(:,[one,zero;zero,zero]==1) = 0.10;
            u0(:,[zero,one;zero,zero]==1) = 0.90;
            
            case 'squareT' % for Triangular elements in a square domain!
            
            % Calculate number of nodes and elements
            u0 = zeros(nN,nE);
            
            % Only valid for square domain
            zero = zeros(sqrt(nE/2),sqrt(nE/2)/2); 
            one = ones(sqrt(nE/2),sqrt(nE/2)/2);
            
            % Initial Condition for our 2D domain
            u0(:,[zero,zero;zero,one]==1) = 0.50;
            u0(:,[zero,zero;one,zero]==1) = 0.70;
            u0(:,[one,zero;zero,zero]==1) = 0.10;
            u0(:,[zero,one;zero,zero]==1) = 0.90;
            
            case 'non-square' % domain is not square!
            
            % set center
            x0 = 0.5; y0 = 0.5;
            
            % Preallocate u0,
            u0 = zeros(nN,nE);
            
            % Initial Condition for our 2D domain
            u0(x >x0 & y >y0) = 0.50; % region 1
            u0(x<=x0 & y >y0) = 0.70; % region 2
            u0(x<=x0 & y<=y0) = 0.10; % region 3
            u0(x >x0 & y<=y0) = 0.90; % region 4 
        end
        
    case 2 % Square Jump
        % set center
        x0 = 0.5*(x(end)+x(1))+0.00; 
        y0 = 0.5*(y(end)+y(1))+0.25;
        
        % parameters
        Lx = x(end)-x(1);
        Ly = y(end)-y(1);
        
        % Preallocate u0
        u0 = ones(size(x));
        
        % Parameters of region
        x1 = x0+Lx/12; x2 = x0-Lx/12;
        y1 = y0+Ly/12; y2 = y0-Ly/12;
        
        % Build Jump
        u0(x>x1) = 0.001;     u0(x<x2) = 0.001;
        u0(y>y1) = 0.001;     u0(y<y2) = 0.001;
        
    case 3 % Sine*Cosine 2-D in Domain [-1,1]^2
        u0 = sin(pi*x).*cos(pi*y);
        
    case 4 % Gaussian Jump
        % set center
        x0 = (x(end)-x(1))/2; y0 = (y(end)-y(1))/2;
        
        % Gaussian
        u0 = 0.1 + 0.5*exp(-20*((x-x0).^2+(y-y0).^2));
        
    case 5 % Cylindrical Jump
        % set center
        x0 = 0.5*(x(end)+x(1)); 
        y0 = 0.5*(y(end)+y(1));
        y0 = y0 + 0.25
        
        % radious
        r = 0.15;
        
        % Gaussian
        u0 = 0.1*ones(size(x));
        u0(sqrt((x+x0).^2+(y+y0).^2)<r) = 1.0;
        
        u0 = zeros(size(x));
        u0(sqrt((x-x0).^2+(y-y0).^2)<r) = 1.0;
        
    case 6 % rectangle in x direction
        % set center
        x0 = 0.0;
        
        % parameters
        Lx = x(end)-x(1);
                
        % Preallocate u0
        u0 = ones(size(x));
        
        % Parameters of region
        x1 = x0+Lx/10; x2 = x0-Lx/10;
                
        % Build Jump
        u0(x>x1) = 0.1;     u0(x<x2) = 0.1;
                
	case 7 % rectangle in y direction
        % set center
        y0 = 0.0;
        
        % parameters
        Ly = y(end)-y(1);
        
        % Preallocate u0
        u0 = ones(size(y));
        
        % Parameters of region
        y1 = y0+Ly/10; y2 = y0-Ly/10;
        
        % Build Jump
        u0(y>y1) = 0.1;     u0(y<y2) = 0.1;
        
    case 8 % Riemann for range [-1,1]
        [nN,nE] = size(x);
        
        % Calculate number of nodes and elements
        u0 = zeros(nN,nE);
        
        % Only valid for square domain
        zero = zeros(sqrt(nE)/2); one = ones(sqrt(nE)/2);
        
        % Initial Condition for our 2D domain
        u0(:,[zero,one;zero,one]==1) = 0.10;
        u0(:,[one,zero;one,zero]==1) = 1.00;
        
    case 9 % Shu's 0.5 + Sin(x+y) for Burgers with periodic conditions
        % make sure is on a domina of size [0,2*pi]x[0,2*pi].
        u0 = 0.5 + sin(x+y);
        
    case 10 % Steady wave for a square domain of size [-1,1]x[-1,1]
        [nN,nE] = size(x);
        
        % Calculate number of nodes and elements
        u0 = zeros(nN,nE);
        
        % Only valid for square domain
        zero = zeros(sqrt(nE)/2); one = ones(sqrt(nE)/2);
        
        % Initial Riemann Condition for our 2D domain
        u0(:,[one,one;zero,zero]==1) = 0.10;
        u0(:,[zero,zero;one,one]==1) = 1.00;
        
    case 11 % Gaussian wave for [-3,3]
        % for diffusion test with Dirichlet BCs. See ref [1].
        mu = 0.01;
        xmid = (x(end) + x(1))/2;
        u0 = exp(-(x-xmid).^2/(4*mu));
        
    otherwise 
        error('case not listed :P')
end
