classdef DGtoolsQuad
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DGtoolsQuad class
    %   Build a Quadrilateral mesh for NDG polynomial reconstructions
    %
    %   by Manuel Diaz, NTU, 2014.06.15
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Refs.: 
    % [1] J.Hestaven & T. Warburton; Nodal Discontinuous Galerkin Methods,
    % Algorithms, Analysis, and Applications. Text in Applied Mathematics,
    % (2008)
    % [2] D. Wirasaet, Et al.; A performance comparison of nodal
    % discontinuous galerkin methods on triangles and quadrilaterals. Int.
    % J. Numer. Fluids (2010) DOI: 10.1002/fld.2376 .
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Implementation based on ideas of refs [1] and [2].
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        pDeg
        kxi
        eta
        x
        y
        nE
        Nnodes
        Nfaces
        nodesPerFace
        FaceMask
    end
    
    properties (Dependent = true, SetAccess = private)
        V1D
        V2D
        GradV2D_kxi
        GradV2D_eta
        M
        invM
        D_kxi
        D_eta
        legLeftEnd
        legRightEnd
        x_kxi
        x_eta
        y_kxi
        y_eta
        kxi_x
        eta_x
        kxi_y
        eta_y
        J
    end
    
    methods (Static)
        function legP = legendreP(kxi,l)
            % Construct Legendre Polynomials
            %**************************************************************
            % Compute the value of the legendre polinomial of degree 'pDeg'
            % for every column value 'kxi'
            %**************************************************************
            x = sym('x'); 
            temp = simplify(1/(2^l*factorial(l))*diff((x^2-1)^l,x,l));
            legP = subs(temp,kxi);
        end
        
        function dlegP = dlegendreP(kxi,l)
            % Construct derivatives of Legendre Polynomials
            %**************************************************************
            % Compute the derivative value of the legendre polinomial of
            % degree 'pDeg' for every column value 'kxi'
            %**************************************************************
            x = sym('x'); 
            legP = simplify(1/(2^l*factorial(l))*diff((x^2-1)^l,x,l));
            dlegP = subs(diff(legP,x),kxi);
        end
        
        function P = JacobiP(x,alpha,beta,N)
            % function [P] = JacobiP(x,alpha,beta,N)
            % Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
            %          (alpha+beta <> -1) at points x for order N and returns P[1:length(xp))]
            % Note   : They are normalized to be orthonormal.
            
            % Turn points into row if needed.
            xp = x; dims = size(xp);
            if (dims(2)==1); xp = xp'; end;
            PL = zeros(N+1,length(xp));
            
            % Initial values P_0(x) and P_1(x)
            gamma0 = 2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*...
                gamma(beta+1)/gamma(alpha+beta+1);
            PL(1,:) = 1.0/sqrt(gamma0);
            if (N==0); P=PL'; return; end;
            gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
            PL(2,:) = ((alpha+beta+2)*xp/2 + (alpha-beta)/2)/sqrt(gamma1);
            if (N==1); P=PL(N+1,:)'; return; end;
            
            % Repeat value in recurrence.
            aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));
            
            % Forward recurrence using the symmetry of the recurrence.
            for i=1:N-1
                h1 = 2*i+alpha+beta;
                anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*...
                    (i+1+beta)/(h1+1)/(h1+3));
                bnew = - (alpha^2-beta^2)/h1/(h1+2);
                PL(i+2,:) = 1/anew*( -aold*PL(i,:) + (xp-bnew).*PL(i+1,:));
                aold =anew;
            end;
            
            P = PL(N+1,:)';
        end
        
        function dP = dJacobiP(r, alpha, beta, N)
            % function [dP] = GradJacobiP(r, alpha, beta, N);
            % Purpose: Evaluate the derivative of the Jacobi polynomial of type (alpha,beta)>-1,
            %	       at points r for order N and returns dP[1:length(r))]
            dP = zeros(length(r), 1);
            if(N == 0)
                dP(:,:) = 0.0;
            else
                dP = sqrt(N*(N+alpha+beta+1))*DGtoolsQuad.JacobiP(r(:),alpha+1,beta+1, N-1);
            end;
        end
                
        function V = Vandermonde1D(pDeg,kxi) % Legendre Vandermonde Matrix
            % Construct Legendre Vandermonde Matrix, V,
            %**************************************************************
            % V is a matrix array in which each element is a Legendre
            % Polynomial of Degree n, n=0:pDeg evaluated at a column array
            % of points, kxi.
            %
            % Coded by Manuel Diaz 2012.12.05
            %**************************************************************
            V = zeros(pDeg);% Allocate
            for l = 0:pDeg 	% All Polynomial Degrees up to pDeg
                j = l + 1;      % Dummy index
                for i = 1:pDeg+1  % All Element Points in kxi
                    V(i,j) = DGtoolsQuad.legendreP(kxi(i),l);
                end
            end
        end
        
        function [P] = BilinearP(xi,eta,i,j)
            % function [P] = QuadP(xi,eta,i,j);
            % Purpose : Evaluate 2D orthonormal polynomial
            %           on quadrilateral element at (xi,eta) of order (i,j).
            h1 = DGtoolsQuad.JacobiP(xi,0,0,i); h2 = DGtoolsQuad.JacobiP(eta,0,0,j);
            P = h1.*h2;
        end
        
        function [dPdxi,dPdeta] = GradBilinearP(xi,eta,i,j)
            % function [P] = QuadP(xi,eta,i,j);
            % Purpose : Evaluate 2D orthonormal gradient of a polynomial
            %           on quadrilateral element at (xi,eta) of order (i,j)
            Pxi = DGtoolsQuad.JacobiP(xi,0,0,i);   dPxi = DGtoolsQuad.dJacobiP(xi,0,0,i);
            Peta = DGtoolsQuad.JacobiP(eta,0,0,j); dPeta = DGtoolsQuad.dJacobiP(eta,0,0,j);
            
            % xi-derivative
            dPdxi = dPxi.*Peta;
            
            % eta-derivative
            dPdeta = Pxi.*dPeta;
        end
        
        %%%%%%%%%%%
        % Normals %
        %%%%%%%%%%%
        
        function [nx,ny,sJ] = Normals(obj)
            % Compute Normals
                
            % interpolate geometric factors to face nodes
            fx_kxi = obj.x_kxi(obj.FaceMask(:),:); fx_eta = obj.x_eta(obj.FaceMask(:),:); 
            fy_kxi = obj.y_kxi(obj.FaceMask(:),:); fy_eta = obj.y_eta(obj.FaceMask(:),:);

            % build normals
            nx = zeros(obj.Nfaces*obj.nodesPerFace,obj.nE);
            ny = zeros(obj.Nfaces*obj.nodesPerFace,obj.nE);
            fid1 = (1:obj.nodesPerFace)'; 
            fid2 = (obj.nodesPerFace+1:2*obj.nodesPerFace)'; 
            fid3 = (2*obj.nodesPerFace+1:3*obj.nodesPerFace)'; 
            fid4 = (3*obj.nodesPerFace+1:4*obj.nodesPerFace)';

            % face 1
            nx(fid1,:) =  fy_kxi(fid1,:); 
            ny(fid1,:) = -fx_kxi(fid1,:);

            % face 2
            nx(fid2,:) =  fy_eta(fid2,:); 
            ny(fid2,:) = -fx_eta(fid2,:);

            % face 3
            nx(fid3,:) = -fy_kxi(fid3,:); 
            ny(fid3,:) =  fx_kxi(fid3,:);

            % face 4
            nx(fid4,:) = -fy_eta(fid4,:); 
            ny(fid4,:) =  fx_eta(fid4,:);

            % normalize
            sJ = sqrt(nx.*nx+ny.*ny); nx = nx./sJ; ny = ny./sJ;
        end
        
        %%%%%%%%%%%%%%%%%
        % LIFT OPERATOR %
        %%%%%%%%%%%%%%%%%
        
        function [LIFT,Emat] = lift(obj)
            % Compute LIFT Operator
            Emat = zeros(obj.Nnodes,obj.Nfaces*obj.nodesPerFace);

            % face 1
            face1 = obj.kxi(obj.FaceMask(:,1));
            V1D = DGtoolsQuad.Vandermonde1D(obj.pDeg,face1); 
            massEdge1 = inv(V1D*V1D');
            Emat(obj.FaceMask(:,1),1:obj.nodesPerFace) = massEdge1;

            % face 2
            face2 = obj.eta(obj.FaceMask(:,2));
            V1D = DGtoolsQuad.Vandermonde1D(obj.pDeg,face2);
            massEdge2 = inv(V1D*V1D');
            Emat(obj.FaceMask(:,2),obj.nodesPerFace+1:2*obj.nodesPerFace) = massEdge2;

            % face 3
            face3 = obj.kxi(obj.FaceMask(:,3));
            V1D = DGtoolsQuad.Vandermonde1D(obj.pDeg,face3); 
            massEdge3 = inv(V1D*V1D');
            Emat(obj.FaceMask(:,3),2*obj.nodesPerFace+1:3*obj.nodesPerFace) = massEdge3;

            % face 4
            face4 = obj.eta(obj.FaceMask(:,4));
            V1D = DGtoolsQuad.Vandermonde1D(obj.pDeg,face4); 
            massEdge4 = inv(V1D*V1D');
            Emat(obj.FaceMask(:,4),3*obj.nodesPerFace+1:4*obj.nodesPerFace) = massEdge4;

            % inv(mass matrix)*Surf_Integral_Sn (L_i,L_j)_{edge_n}
            LIFT = obj.V2D*(obj.V2D'*Emat);
        end
        
        function LIFT = liftFR(obj,dg)
            % Compute LIFT Operator
            vNodes = reshape(1:obj.eNodes,[obj.pDeg+1,obj.pDeg+1]);
            LIFT = zeros(obj.eNodes,obj.eFaces*obj.eNodesPerFace);
                       
            % nodes in face 1 
            for i = 1:obj.eNodesPerFace
                LIFT(vNodes(:,i),i) = dg.R;
            end
            
            % face 2
            for i = 1:obj.eNodesPerFace
                LIFT(vNodes(i,:),i+obj.eNodesPerFace) = -dg.L;
            end

            % face 3
            for i = 1:obj.eNodesPerFace
                LIFT(vNodes(:,i),i+2*obj.eNodesPerFace) = -dg.L;
            end

            % face 4
            for i = 1:obj.eNodesPerFace
                LIFT(vNodes(i,:),i+3*obj.eNodesPerFace) = dg.R;
            end
            
            % Adjust the sign of coefs to match the NDG's formulation.
            LIFT = -LIFT;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % INTERPOLATION OPERATOR %
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [IM,QUAD] = InterpMatrix2D(N,obj)
            % Build Local Interpolation Operator.
            % Here want to interpolate information form NDG nodes to FEM (uniformly
            % distributed) nodes in every element. Manuel Diaz. NTU, 2014.07.10
            
            % Build uniformly distributed nodes on reference quadrilateral
            kxi_U=linspace(-1,1,N); eta_U=linspace(-1,1,N); [kxi_U,eta_U]=meshgrid(kxi_U,eta_U);
            Np=N*N; kxi_FEM=reshape(kxi_U,Np,1); eta_FEM=reshape(eta_U,Np,1);  
            
            % Compute Vandermonde Matrix at (kxi_FEM,eta_FEM) points
            P = obj.pDeg;
            
            % Initialize the 2d Vandermonde Matrix, V_{ij} = phi_j(r_i, s_i);
            V_FEM = zeros(Np,(P+1)*(P+1));
            sk = 0;
            for i=0:P      % for every p-degree
                for j=0:P  % for every p-degree
                    sk = sk+1;      % step counter
                    V_FEM(:,sk) = DGtoolsQuad.BilinearP(kxi_FEM,eta_FEM,i,j);
                end
            end
            
            % Build interpolation matrix
            IM = V_FEM*inv(obj.V2D);
            
            % Build EtoV matrix
            quad=zeros((N-1)*(N-1),4); 
            ks=1;
            for i = 1:N-1;
                for j = 1:N-1;
                    k = j+(N)*(i-1);
                    quad(ks,:) = [k, (N)+k, (N)+k+1, k+1];
                    ks = ks+1;
                end
            end
            
            % build quadrangulation for all equally spaced nodes on all elements
            QUAD = [];
            for k=1:obj.nE
                QUAD = [QUAD; quad+(k-1)*Np];
            end
        end
        
    end %Static Methods
     
    methods
        function obj = DGtoolsQuad(mesh) % Constructor
            obj.pDeg = mesh.pDeg;
            obj.kxi = mesh.kxi;
            obj.eta = mesh.eta;
            % mesh data for metrics, normals and Lift OP.
            obj.x = mesh.x;
            obj.y = mesh.y;
            obj.FaceMask = mesh.FaceMask;
            obj.nodesPerFace = mesh.eNodesPerFace;
            obj.Nnodes = mesh.eNodes;
            obj.Nfaces = mesh.eFaces;
            obj.nE = mesh.nE;
        end

        %%%%%%%%%%%%%%%%
        %  MAIN TOOLS: %
        %%%%%%%%%%%%%%%%
        
        function V = get.V1D(obj) % Linear Jacobi Vandermonde Matrix
            % Construct Jacobi Vandermonde Matrix, V,
            %**************************************************************
            % V is a matrix array in which each element is a jacobi
            % Polynomial of Degree n, n=0:pDeg evaluated at a column array
            % of points, kxi & eta.
            %
            % Coded by Manuel Diaz 2012.12.05
            %**************************************************************\
            V = zeros(obj.pDeg);% Allocate
            for l = 0:obj.pDeg 	% All Polynomial Degrees up to pDeg
                j = l + 1;      % Dummy index
                for i = 1:obj.pDeg+1  % All Element Points in kxi
                    V(i,j) = obj.legendreP(obj.eta(i),l);
                end
            end
        end
        
        function V = get.V2D(obj) % Bilinear Jacobi Vandermonde Matrix
            % Construct Jacobi Vandermonde Matrix, V,
            %**************************************************************
            % V is a matrix array in which each element is a jacobi
            % Polynomial of Degree n, n=0:pDeg evaluated at a column array
            % of points, kxi & eta.
            %
            % Coded by Manuel Diaz 2012.12.05
            %**************************************************************
            V = zeros(obj.Nnodes,obj.Nnodes);
            s = 0;
            for i = 0:obj.pDeg      % for every p-degree
                for j = 0:obj.pDeg  % for every p-degree
                    s = s + 1;      % step counter
                    V(:,s) = obj.BilinearP(obj.kxi,obj.eta,i,j);
                end
            end
        end
                
        %%%%%%%%%%%%%%%
        % NODAL TOOLS %
        %%%%%%%%%%%%%%%
        
        function M = get.M(obj) % Mass Matrix
            M = inv(obj.V2D*obj.V2D');
        end
        
        function invM = get.invM(obj) % Inverse of Mass Matrix
            invM = obj.V2D*obj.V2D';
        end
        
        function V2Dkxi = get.GradV2D_kxi(obj)
            V2Dkxi = zeros(length(obj.kxi),obj.Nnodes);
            sk = 1;
            for i=0:obj.pDeg    % for every p-degree
                for j=0:obj.pDeg    % for every p-degree
                    [V2Dkxi(:,sk),~] = obj.GradBilinearP(obj.kxi,obj.eta,i,j);
                    sk = sk+1;
                end
            end
        end
        
        function V2Deta = get.GradV2D_eta(obj)
            V2Deta= zeros(length(obj.kxi),obj.Nnodes);
            sk = 1;
            for i=0:obj.pDeg    % for every p-degree
                for j=0:obj.pDeg    % for every p-degree
                    [~,V2Deta(:,sk)] = obj.GradBilinearP(obj.kxi,obj.eta,i,j);
                    sk = sk+1;
                end
            end
        end
        
        function nDkxi = get.D_kxi(obj)
            nDkxi = obj.GradV2D_kxi/obj.V2D;
        end
        
        function nDeta = get.D_eta(obj)
            nDeta = obj.GradV2D_eta/obj.V2D;
        end
        
        %%%%%%%%%%%
        % METRICS %
        %%%%%%%%%%%
        
        function output = get.x_kxi(obj)
            output = obj.D_kxi*obj.x;
        end
        
        function output = get.y_kxi(obj)
            output = obj.D_kxi*obj.y;
        end
        
        function output = get.x_eta(obj)
            output = obj.D_eta*obj.x;
        end
        
        function output = get.y_eta(obj)
            output = obj.D_eta*obj.y;
        end
        
        function Jacobian = get.J(obj)
            % x_kxi = D_kxi*x; % y_kxi = D_kxi*y;
            % x_eta = D_eta*x; % y_eta = D_eta*y;
            % J = x_kxi.*y_eta - x_eta.*y_kxi;
            Jacobian = obj.x_kxi.*obj.y_eta - obj.y_kxi.*obj.x_eta;
        end
        
         function output = get.kxi_x(obj)
            output = obj.y_eta ./ obj.J;
        end
        
        function output = get.eta_x(obj)
            output = -obj.y_kxi ./ obj.J;
        end
        
        function output = get.kxi_y(obj)
            output = -obj.x_eta ./ obj.J;
        end
        
        function output = get.eta_y(obj)
            output = obj.x_kxi ./ obj.J;
        end
        
    end % Methods
end % Class

