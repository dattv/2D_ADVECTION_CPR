classdef CorrectionPolynomial
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CORRECTIONPOLYNOMIAL class
    %   Build Correction Polynomials of degree K for CPR scheme.
    %
    %   by Manuel Diaz, NTU, 2013.10.12
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        pDeg
        pType
        kxi
        L
        R
    end
    
    properties (Dependent = true, SetAccess = private)
        xi % local coordiante
        P  % Correction Polynomail
        dP % derivate of correction Polynomial
    end
    
    methods (Static)
        function legP = LegendreP(l)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Legendre Polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : kDeg: Polynomial Degree requested
            % Output : legP: symbolic Legendre polynomial
            %
            x = sym('x'); 
            legP = simplify(1/(2^l*factorial(l))*diff((x^2-1)^l,x,l));
        end
                   
    end % Methods
           
    methods
        function obj = CorrectionPolynomial(type,Kdeg,solutionPoints)
            obj.pDeg = Kdeg;
            obj.pType = type;
            obj.kxi = solutionPoints;
        end
        
        function cpoly = get.P(obj)
            switch obj.pType
                case 'Legendre'
                    cpoly = obj.LegendreP(obj.pDeg);
                case 'LGL'
                    cpoly = obj.LobattoP(obj.pDeg);
                case 'RadauRight'
                    cpoly = obj.RadauRightP(obj.pDeg);
                case 'RadauLeft'
                    cpoly = obj.RadauLeftP(obj.pDeg);
                case 'DGRight'
                    cpoly = obj.NDGRightP(obj.pDeg);
                case 'DGLeft'
                    cpoly = obj.NDGLeftP(obj.pDeg);
                case 'SDRight'
                    cpoly = obj.SDRightP(obj.pDeg);
                case 'SDLeft'
                    cpoly = obj.SDLeftP(obj.pDeg);
                case 'HURight'
                    cpoly = obj.HURightP(obj.pDeg);
                case 'HULeft'
                    cpoly = obj.HULeftP(obj.pDeg);
                case 'CinftyRight'
                    cpoly = obj.RadauRightP(obj.pDeg);
                case 'CinftyLeft'
                    cpoly = obj.RadauLeftP(obj.pDeg);
                otherwise
                    error('correction polynomial not available')
            end
        end
        
        function RRP = RadauRightP(obj,kDeg)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load from table the coefs for Radau Polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : kDeg: Polynomial Degree requested
            % Output : RRP: Right Radau polynomial
            %
            RRP = (-1)^(kDeg)/2*(obj.LegendreP(kDeg) - obj.LegendreP(kDeg-1));
        end
            
        function RLP = RadauLeftP(obj,kDeg)
        	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load from table the coefs for Radau Polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : kDeg: Polynomial Degree requested
            % Output : RLP: Left Radau polynomial
            %
            RLP = (1/2)*(obj.LegendreP(kDeg) + obj.LegendreP(kDeg-1) );
        end
            
        function lobP = LobattoP(obj,kDeg)
        	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Lobatto Polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : kDeg: Polynomial Degree requested
            % Output : lobP: Symbolic Lobatto polynomial
            %   
            lobP = obj.LegendreP(kDeg) - obj.LegendreP(kDeg-2);
        end
        
        function DGRP = NDGRightP(obj,kDeg)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load from table the coefs for NDG correction polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : kDeg: Polynomial Degree requested
            % Output : DGRP: Right Radau polynomial
            %
            DGRP = (-1)^(kDeg)/2*(obj.LegendreP(kDeg) - obj.LegendreP(kDeg+1));
        end
            
        function DGLP = NDGLeftP(obj,kDeg)
        	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load from table the coefs for NDG correction polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : kDeg: Polynomial Degree requested
            % Output : DGLP: Left Radau polynomial
            %
            DGLP = (1/2)*(obj.LegendreP(kDeg) + obj.LegendreP(kDeg+1) );
        end
        
        function SDRP = SDRightP(obj,kDeg)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load from table the coefs for SD correction polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : kDeg: Polynomial Degree requested
            % Output : SDRP: Right Radau polynomial
            %
            x = sym('x');
            SDRP = (-1)^(kDeg)/2*(1-x)*obj.LegendreP(kDeg);
        end
            
        function SDLP = SDLeftP(obj,kDeg)
        	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load from table the coefs for SD correction polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : kDeg: Polynomial Degree requested
            % Output : SDLP: Left Radau polynomial
            %
            x = sym('x');
            SDLP = (1/2)*(1+x)*obj.LegendreP(kDeg);
        end
        
        function HURP = HURightP(obj,k)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load from table the coefs for Huyng correction polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : k: Polynomial Degree requested
            % Output : HURP: Right Radau polynomial
            %
            HURP = (-1)^(k)/2*(obj.LegendreP(k) - ...
               ((k+1)*obj.LegendreP(k-1) + k*obj.LegendreP(k+1))/(2*k+1));
        end
            
        function HULP = HULeftP(obj,k)
        	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load from table the coefs for Huyng correction polynomials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Input : kDeg: Polynomial Degree requested
            % Output : HULP: Left Radau polynomial
            %
            HULP = (1/2)*(obj.LegendreP(k) + ...
               ((k+1)*obj.LegendreP(k-1) + k*obj.LegendreP(k+1))/(2*k+1));
        end
        
        function dcorrection = get.dP(obj)
            x = sym('x'); dcorrection = diff(obj.P,x);
        end
        
        function dg = eval_dP(obj,solutionPoints)
            dg.Right = double(subs(obj.dP,solutionPoints));
            dg.Left = -flipud(dg.Right);
        end
        
        function dgR = get.R(obj)
            dgR = double(subs(obj.dP,obj.kxi));
        end
        
        function dgL = get.L(obj)
            dgL = -flipud(double(subs(obj.dP,obj.kxi)));
        end
        
    end % Methods
    
end % Class

