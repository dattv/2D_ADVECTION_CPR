classdef DGmeshQuad
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DGmeshQuad class
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
        meshingRoutine  % Meshing strategy
        nEx             % number Elements in x-direction
        nEy             % number Elements in y-direction
        nE              % Total number of Elements
        pDeg            % Polynomial degree
        VX              % elements x-coordinate of vertices
        VY              % elements y-coordinate of vertices
        nV              % number of vertices
        EtoV            % EtoV: elements to vertices (connectivity mat)
        EtoE            % EtoE: element# connects to element#
        EtoF            % EtoF: element# with face# to face#
        eNodes          % nodes per element
        eFaces          % faces per element
        eNodesPerFace   % nodes per face
        NODETOL         % tolerance to identify node with same coordinates
        solutionPoints  % Basic solution points
        x               % x-Global node coordinates
        y               % y-Global node coordinates
        kxi             % kxi-standard element coordiantes.
        eta             % eta-standard element coordinates.
        w_kxi           % quadrature weigths in kxi
        w_eta           % quadrature weigths in eta
        q               % quadrature weights (w_kxi)x(w_eta)
        idx
    end
    
    properties (Dependent = true, SetAccess = private)
        FaceMask    % mask with all nodes on the edge
        Fx          % Face nodes x-coordinates
        Fy          % Face nodes y-coordiantes
    end
    
    methods (Static)
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
        
        function [x,w,P]=GaussLobatto(kk)
        	N = kk-1; % Compute for the number of points kk
            % lglnodes.m
            %
            % Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde
            % matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
            % integration and spectral methods.
            %
            % Reference on LGL nodes and weights:
            %   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
            %   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
            %
            % Written by Greg von Winckel - 04/17/2004
            % Contact: gregvw@chtm.unm.edu
            
            % Truncation + 1
            N1=N+1;
            % Use the Chebyshev-Gauss-Lobatto nodes as the first guess
            x=-cos(pi*(0:N)/N)';
            % The Legendre Vandermonde Matrix
            P=zeros(N1,N1);
            % Compute P_(N) using the recursion relation
            % Compute its first and second derivatives and
            % update x using the Newton-Raphson method.
            xold=2;
            while max(abs(x-xold))>eps
                xold=x;
                P(:,1)=1;    P(:,2)=x;
                for k=2:N
                    P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
                end
                x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
            end
            w=2./(N*N1*P(:,N1).^2);
        end
        
        function [x] = JacobiGL(alpha,beta,N)
            % function [x] = JacobiGL(alpha,beta,N)
            % Purpose: Compute the N'th order Gauss Lobatto quadrature
            %          points, x, associated with the Jacobi polynomial,
            %          of type (alpha,beta) > -1 ( <> -0.5).
            x = zeros(N+1,1);
            if (N==1); x(1)=-1.0; x(2)=1.0; return; end;
            
            [xint,~] = DGmeshQuad.JacobiGQ(alpha+1,beta+1,N-2);
            x = [-1, xint', 1]';
        end
        
        function [x,w] = JacobiGQ(alpha,beta,N)
            
            % function [x,w] = JacobiGQ(alpha,beta,N)
            % Purpose: Compute the N'th order Gauss quadrature points, x,
            %          and weights, w, associated with the Jacobi
            %          polynomial, of type (alpha,beta) > -1 ( <> -0.5).
            if (N==0); x(1)= -(alpha-beta)/(alpha+beta+2); w(1) = 2; return; end;
            
            % Form symmetric matrix from recurrence.
            J = zeros(N+1); %#ok<NASGU>
            h1 = 2*(0:N)+alpha+beta;
            J = diag(-1/2*(alpha^2-beta^2)./(h1+2)./h1) + ...
                diag(2./(h1(1:N)+2).*sqrt((1:N).*((1:N)+alpha+beta).*...
                ((1:N)+alpha).*((1:N)+beta)./(h1(1:N)+1)./(h1(1:N)+3)),1);
            if (alpha+beta<10*eps); J(1,1)=0.0; end;
            J = J + J';
            
            %Compute quadrature by eigenvalue solve
            [V,D] = eig(J); x = diag(D);
            w = (V(1,:)').^2*2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*...
                gamma(beta+1)/gamma(alpha+beta+1);
        end
        
        function dP = dJacobiP(r, alpha, beta, N)
            % function [dP] = GradJacobiP(r, alpha, beta, N);
            % Purpose: Evaluate the derivative of the Jacobi polynomial of type (alpha,beta)>-1,
            %	       at points r for order N and returns dP[1:length(r))]
            
            dP = zeros(length(r), 1);
            if(N == 0)
                dP(:,:) = 0.0;
            else
                dP = sqrt(N*(N+alpha+beta+1))*DGmeshQuad.JacobiP(r(:),alpha+1,beta+1, N-1);
            end;
        end
        
        function [V1D] = Vandermonde1D(N,r)
            % function [V1D] = Vandermonde1D(N,r)
            % Purpose : Initialize the 1D Vandermonde Matrix, V_{ij} = phi_j(r_i);
            
            V1D = zeros(length(r),N+1);
            for j=1:N+1
                V1D(:,j) = DGmeshQuad.JacobiP(r(:), 0, 0, j-1);
            end;
        end
        
        function [mapM, mapP, vmapM, vmapP, vmapB, mapB] = BuildMaps2D(obj)
            % Purpose: compute elements faces-

            % number volume nodes consecutively
            nodeids = reshape(1:obj.nE*obj.eNodes, obj.eNodes, obj.nE);
            vmapM   = zeros(obj.eNodesPerFace, obj.eFaces, obj.nE); 
            mapM  = (1:obj.nE*obj.eNodesPerFace*obj.eFaces)';     
            mapP  = reshape(mapM,[obj.eNodesPerFace,obj.eFaces,obj.nE]);

            % indexes of face nodes on every element using volume node ordening
            for e1=1:obj.nE
                for f1=1:obj.eFaces
                    vmapM(:,f1,e1) = nodeids(obj.FaceMask(:,f1), e1);
                end
            end

            % Purpose: compute elements faces+
            vmapP = zeros(obj.eNodesPerFace, obj.eFaces, obj.nE);
            
            one = ones(1, obj.eNodesPerFace);
            for e1=1:obj.nE
                for f1=1:obj.eFaces
                    % find neighbor
                    k2 = obj.EtoE(e1,f1); f2 = obj.EtoF(e1,f1);

                    % reference length of edge
                    v1 = obj.EtoV(e1,f1); v2 = obj.EtoV(e1, 1+mod(f1,obj.eFaces));
                    refd = sqrt( (obj.VX(v1)-obj.VX(v2))^2 + (obj.VY(v1)-obj.VY(v2))^2 );

                    % find find volume node numbers of left and right nodes
                    vidM = vmapM(:,f1,e1); vidP = vmapM(:,f2,k2);
                    x1 = obj.x(vidM); y1 = obj.y(vidM); 
                    x2 = obj.x(vidP); y2 = obj.y(vidP);
                    x1 = x1*one;  y1 = y1*one;  x2 = x2*one;  y2 = y2*one;

                    % Compute distance matrix
                    D = (x1 -x2').^2 + (y1-y2').^2;
                    [idM, idP] = find(sqrt(abs(D))< obj.NODETOL*refd);
                    vmapP(idM,f1,e1) = vidP(idP); 
                    mapP(idM,f1,e1) = idP + (f2-1)*obj.eNodesPerFace+(k2-1)*obj.eFaces*obj.eNodesPerFace;
                end
            end

            % Purpose: reshape vmapM and vmapP to be vectors and create boundary node list
            vmapP = vmapP(:); vmapM = vmapM(:); mapP = mapP(:);
            mapB = find(vmapP==vmapM); vmapB = vmapM(mapB); % or by:
            %vmapB = vmapM(vmapP==vmapM);
        end
        
        function map = BuildMaps2Dv2(obj)
            % Purpose: Compute maps. Coded and Extended by Manuel Diaz 2014.06.17
            
            % 0.Purpose: load nodes coordinates
            xx = obj.x; yy=obj.y;
            
            % 1.Purpose: compute elements faces-
            % number volume nodes consecutively
            nodeids = reshape(1:obj.nE*obj.eNodes, obj.eNodes, obj.nE);
            vmapM   = zeros(obj.eNodesPerFace, obj.eFaces, obj.nE); 
            mapM  = (1:obj.nE*obj.eNodesPerFace*obj.eFaces)';     
            mapP  = reshape(mapM,[obj.eNodesPerFace,obj.eFaces,obj.nE]);
            
            % indexes of face nodes on every element using volume node ordening
            FMask = obj.FaceMask;
            for e1=1:obj.nE
                for f1=1:obj.eFaces
                    vmapM(:,f1,e1) = nodeids(FMask(:,f1), e1);
                end
            end
            
            % 2.Purpose: Build Boundary Maps using elements face- info!
            % Find a way to distinguish every side in our domain as a single boundary.
            % Let's Take advantage of EtoE table:
            Eid = (1:obj.nE)'; FaceB = zeros(obj.nE,obj.eFaces);
            for fi = 1:obj.eFaces % i-face
                FaceB(:,fi) = (obj.EtoE(:,fi)==Eid);
            end
            
            % find elements with common boundary face#
            ElemBF1 = find(FaceB(:,1)==1);
            ElemBF2 = find(FaceB(:,2)==1);
            ElemBF3 = find(FaceB(:,3)==1);
            ElemBF4 = find(FaceB(:,4)==1);
            
            % Lets re-use face Ids
            fid1 = (1:obj.eNodesPerFace)';
            fid2 = (obj.eNodesPerFace+1:2*obj.eNodesPerFace)';
            fid3 = (2*obj.eNodesPerFace+1:3*obj.eNodesPerFace)';
            fid4 = (3*obj.eNodesPerFace+1:4*obj.eNodesPerFace)';
            
            % We must take advantage of the volume nodes map for faces,
            vM = reshape(vmapM,[obj.eNodesPerFace*obj.eFaces,obj.nE]);
             M = reshape( mapM,[obj.eNodesPerFace*obj.eFaces,obj.nE]);
            
            % Build Boundary maps using faces's nodes Ids and the normal vectors Ids,
            
            % for faces-1
            map.vD = reshape(vM(fid1,ElemBF1),[],1);
             map.D = reshape( M(fid1,ElemBF1),[],1);
            
            % for faces-2
            map.vR = reshape(vM(fid2,ElemBF2),[],1);
             map.R = reshape( M(fid2,ElemBF2),[],1);
            
            % for faces-3
            map.vU = reshape(vM(fid3,ElemBF3),[],1);
             map.U = reshape( M(fid3,ElemBF3),[],1);
            
            % for faces-4
            map.vL = reshape(vM(fid4,ElemBF4),[],1);
             map.L = reshape( M(fid4,ElemBF4),[],1);

            % 3.Purpose: compute elements faces+
            vmapP = zeros(obj.eNodesPerFace, obj.eFaces, obj.nE);
            
            one = ones(1, obj.eNodesPerFace);
            for e1=1:obj.nE
                for f1=1:obj.eFaces
                    % find neighbor
                    k2 = obj.EtoE(e1,f1); f2 = obj.EtoF(e1,f1);

                    % reference length of edge
                    v1 = obj.EtoV(e1,f1); v2 = obj.EtoV(e1, 1+mod(f1,obj.eFaces));
                    refd = sqrt( (obj.VX(v1)-obj.VX(v2))^2 + (obj.VY(v1)-obj.VY(v2))^2 );

                    % find find volume node numbers of left and right nodes
                    vidM = vmapM(:,f1,e1); vidP = vmapM(:,f2,k2);
                    x1 = xx(vidM); y1 = yy(vidM); 
                    x2 = xx(vidP); y2 = yy(vidP);
                    x1 = x1*one;  y1 = y1*one;  x2 = x2*one;  y2 = y2*one;

                    % Compute distance matrix
                    D = (x1 -x2').^2 + (y1-y2').^2;
                    [idM, idP] = find(sqrt(abs(D))< obj.NODETOL*refd);
                    vmapP(idM,f1,e1) = vidP(idP); 
                    mapP(idM,f1,e1) = idP + (f2-1)*obj.eNodesPerFace+(k2-1)*obj.eFaces*obj.eNodesPerFace;
                end
            end

            % 4.Purpose: reshape vmapM and vmapP to be vectors and create boundary node list
            map.vP = vmapP(:); map.vM = vmapM(:); map.P = mapP(:); map.M = mapM(:);
            map.B = find(map.vP==map.vM); map.vB = vmapM(map.B); 
        end
        
        function plotmesh(x,y,z)
            % Take advantage of surf plot to plot our mesh
            nNodes = size(z,1); nElem = size(z,2);  
            nNx = sqrt(nNodes); nNy = nNx; % nodes per size
            id = reshape(1:nNodes*nElem,[nNy,nNx,nElem]);
            % plot a surface for every element
            hold on;
            for e = 1:nElem
                surf(y(id(:,:,e)),x(id(:,:,e)),z(id(:,:,e)));
            end
            hold off;
        end
        
        function plotcontourf(x,y,z)
            % Take advantage of surf plot to plot our mesh
            nNodes = size(z,1); nElem = size(z,2);  
            nNx = sqrt(nNodes); nNy = nNx; % nodes per size
            id = reshape(1:nNodes*nElem,[nNy,nNx,nElem]);
            % plot a surface for every element
            hold on;
            for e = 1:nElem
                contourf(y(id(:,:,e)),x(id(:,:,e)),z(id(:,:,e)));
            end
            hold off; axis([0,1,0,1])
        end
        
        function plotnodes(x,y)
            % Take advantage of scatter plot to plot every node in the mesh
            nNodes = size(x,1); nElem = size(x,2);
            hold on;
            % Plot a single point per node
            for e = 1:nElem
                scatter(x(:,e),y(:,e)); 
            end
            % Plot Nodes number
            for i = 1:nNodes*nElem
                text(x(i),y(i),int2str(i),'fontsize',8,....
                    'fontweight','bold','Color','r');
            end
            hold off;
        end
        
        function plotVertices(EtoV,VX,VY)
            nVert = size(VX,1); nElem = size(EtoV,1);
            patch('Faces',EtoV,'Vertices',[VX,VY],...
                'FaceColor','none','EdgeColor','k'); 
            hold on;
            % Plot Nodes number
            for i = 1:nVert
                text(VX(i),VY(i),int2str(i),'fontsize',14,....
                    'fontweight','bold','Color','k');
            end
            % Plot Element number
            for e = 1:nElem
                pos = [mean(VX(EtoV(e,:))),mean(VY(EtoV(e,:)))];
                text(pos(1),pos(2),int2str(e),'fontsize',14,...
                    'fontweight','bold','color','b');
            end
            hold off;
        end
        
    end % Static Methods
    
    methods
        function obj = DGmeshQuad(mesher,xRange,yRange,nEx,nEy,Pdeg) % The Constuctor
            
            switch mesher
                    case 'Q4' % Quadrilateral mesh
                        [vx,vy,E2V,nv,ne] = Q4meshGenerator(xRange,yRange,nEx,nEy);
                otherwise
                    error('Mesher name not in list')
            end
            
            % Re-order elements to ensure counter clockwise orientation
            ax = vx(E2V(:,1)); ay = vy(E2V(:,1));
            bx = vx(E2V(:,2)); by = vy(E2V(:,2));
            cx = vx(E2V(:,3)); cy = vy(E2V(:,3));
            dx = vx(E2V(:,4)); dy = vy(E2V(:,4));
            %Area = 0.5*(-ay.*bx+ax.*by-by.*cx+bx.*cy-cy.*dx+cx.*dy);
            Area = 0.5*(ax-dx).*(cy-dy)-(bx-dx).*(by-dy)-(cx-dx).*(ay-dy);
            i = find(Area<0); fprintf('problem elements: %1.0f\n',i);
            E2V(i,:) = E2V(i,[1 4 3 2]);
            
            % Main parameters of mesh
            obj.meshingRoutine = mesher;
            obj.nEx = nEx;
            obj.nEy = nEy;
            obj.nE  = nEx*nEy;
            obj.pDeg= Pdeg;
            obj.VX  = vx;
            obj.VY  = vy;
            obj.nV  = nv;
            obj.EtoV= E2V;
                       
            % Define quadrilateral mesh parameters
            Nfaces = 4; TotalFaces = Nfaces*ne;
            obj.eFaces = Nfaces; 
            obj.eNodes = (Pdeg+1)*(Pdeg+1);
            obj.eNodesPerFace = Pdeg+1; 
            obj.NODETOL= 1e-12;
            
            % Nodes in Bilinear Standard Element (BSE)
            [KXI,WXI] = DGmeshQuad.GaussLobatto(Pdeg+1);
            [nkxi,neta] = meshgrid(KXI,KXI); 
            [Wkxi,Weta] = meshgrid(WXI,WXI);
            obj.kxi = reshape(nkxi,[obj.eNodes,1]);
            obj.eta = reshape(neta,[obj.eNodes,1]);
            obj.w_kxi = reshape(Wkxi,[obj.eNodes,1]);
            obj.w_eta = reshape(Weta,[obj.eNodes,1]);
            obj.q = obj.w_kxi.*obj.w_eta;
            obj.solutionPoints = KXI;
            
            % Build EtoE and EtoV connectivities
            
            % List of local face to local vertex connections
            vn = [[1,2];[2,3];[3,4];[4,1]];
            
            % Build global face to node sparse array
            SpFToV = spalloc(TotalFaces, nv, 2*TotalFaces);
            sk = 1;
            for k=1:ne
                for face=1:Nfaces
                    SpFToV( sk, obj.EtoV(k, vn(face,:))) = 1; %#ok<SPRIX>
                    sk = sk+1;
                end
            end
            
            % Build global face to global face sparse array
            SpFToF = SpFToV*SpFToV' - 2*speye(TotalFaces);
            
            % Find complete face to face connections
            [faces1, faces2] = find(SpFToF==2);
            
            % Convert face global number to element and face numbers
            element1 = floor( (faces1-1)/Nfaces ) + 1; face1 = mod( (faces1-1), Nfaces ) + 1;
            element2 = floor( (faces2-1)/Nfaces ) + 1; face2 = mod( (faces2-1), Nfaces ) + 1;
            
            % Rearrange into Nelements x Nfaces sized arrays
            ind = sub2ind([ne, Nfaces], element1, face1);
            
            obj.EtoE = (1:ne)'*ones(1,Nfaces); obj.EtoF = ones(ne,1)*(1:Nfaces);
            obj.EtoE(ind) = element2; obj.EtoF(ind) = face2;
        end
                
        function X = get.x(obj)
            % build x-coordinates for every node
            va = obj.EtoV(:,1)'; vb = obj.EtoV(:,2)'; 
            vc = obj.EtoV(:,3)'; vd = obj.EtoV(:,4)';
            X = (1-obj.kxi).*(1-obj.eta)/4*obj.VX(va)' + ...
                (1+obj.kxi).*(1-obj.eta)/4*obj.VX(vb)' + ...
                (1+obj.kxi).*(1+obj.eta)/4*obj.VX(vc)' + ...
                (1-obj.kxi).*(1+obj.eta)/4*obj.VX(vd)';
        end
        
        function Y = get.y(obj)
            % build y-coordinates for every node
            va = obj.EtoV(:,1)'; vb = obj.EtoV(:,2)'; 
            vc = obj.EtoV(:,3)'; vd = obj.EtoV(:,4)';
            Y = (1-obj.kxi).*(1-obj.eta)/4*obj.VY(va)' + ...
                (1+obj.kxi).*(1-obj.eta)/4*obj.VY(vb)' + ...
                (1+obj.kxi).*(1+obj.eta)/4*obj.VY(vc)' + ...
                (1-obj.kxi).*(1+obj.eta)/4*obj.VY(vd)';
        end
        
        function FMask = get.FaceMask(obj)
            % find all the nodes that lie on each edge
            fmask1 = find( abs(obj.eta+1) < obj.NODETOL)';
            fmask2 = find( abs(obj.kxi-1) < obj.NODETOL)';
            fmask3 = find( abs(obj.eta-1) < obj.NODETOL)';
            fmask4 = find( abs(obj.kxi+1) < obj.NODETOL)';
            FMask = [fmask1;fmask2;fmask3;fmask4]';
        end
        
        function FaceXnodes = get.Fx(obj)
            % mask the x coordinates of the nodes that lie on each edge
            FaceXnodes = obj.x(obj.FaceMask(:),:); 
        end
        
        function FaceYnodes = get.Fy(obj)
            % mask the y coordinates of the nodes that lie on each edge
            FaceYnodes = obj.y(obj.FaceMask(:),:); 
        end
        
        function index = get.idx(obj)
           index = reshape(1:obj.nE*obj.eNodes,obj.eNodes,obj.nE); 
        end
        
    end % Private Methods
    
end % End Class