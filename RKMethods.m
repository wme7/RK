classdef RKMethods
    % RKMETHODS - Runge Kutta Methods
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       This is a collection of RK coefficients for my favorite
    %       Implicit & Implicit-Explicit (IMEX) Runge Kutta Methods
    %
    %               coded by Manuel Diaz, NTU, 2014.01.21
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Refs.
    % (1) Lorenzo Pareschi & Giovanni Russo; Journal of Scientific Computing,
    %     Vol.25, Nos.1/2, (2005) DOI: 10.1007/s10915-004-4636-4
    % (2) CA. Kennedy & MH. Carpenter; Applied Numerical Mathematics 44 (2003)
    %     139-181
    % (3) A. KVÆRNØ; Singly diagonal implicit runge-kutta methods with and
    %     explicit first stage, BIT Numerical Mathematics 44: 489–502, (2004).
    % (4) Mark H. Carpenter, Sally A. Viken, & Eric J. Nielsen; The
    %     efficiency of high order temporal schemes, AIAA 2003-0086. 
    %     Downloaded Profs secret page:
    %     http://fun3d.larc.nasa.gov/papers/carpenter.pdf 
    % (5) Sandra Pieraccini, Gabriella Puppo; Implicit-Explicit schemes for
    %     BGK kinetic equations (2007)
    % (6) Sigal Gottlieb, David I. Ketcheson, Chi-Wang Shu; Strong
    %     Stability Preserving Runge-Kutta and Multistep Time Discretizations,  
    %     World Scientific (2011), ISBN	9814289264, 9789814289269
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        method
        r       % ssp parameter
        s       % stages
        Ahat    % Explicit coefs
        bhat
        chat
        A       % Implicit coefs
        b
        btilde  % Extra parameter in ARKs
        c
    end
    
    methods
        function obj = RKMethods(name) % Constructor
            obj.method = name;
            % List of coeficients for EXplicit, Implicit, IM-EX and ARK
            % Runge Kutta schemes.
            switch obj.method
                case 'FE11'
                    %Forward Euler
                    obj.s=1; obj.r=1;
                    obj.Ahat=0;
                    obj.bhat=1; obj.chat=0;
                    
                case 'SSP22'
                    obj.s=2; obj.r=1;
                    obj.Ahat=[0 0; 1 0];
                    obj.bhat=[1/2 1/2]; obj.chat=sum(obj.Ahat,2)';
                    
                case 'SSP42'
                    obj.s=4; obj.r=3;
                    obj.Ahat=[0 0 0 0; 1/3 0 0 0; 1/3 1/3 0 0; 1/3 1/3 1/3 0];
                    obj.bhat=1/4*ones(4,1);	obj.chat=sum(obj.Ahat,2)';
                    
                case 'SSP33'
                    obj.s=3; obj.r=1;
                    obj.Ahat=[0 0 0; 1 0 0; 1/4 1/4 0];
                    obj.bhat=[1/6 1/6 2/3]; obj.chat=sum(obj.Ahat,2)';
                    
                case 'SSP43'
                    obj.s=4; obj.r=2;
                    obj.Ahat=[0 0 0 0; 1/2 0 0 0; 1/2 1/2 0 0; 1/6 1/6 1/6 0];
                    obj.bhat=[1/6 1/6 1/6 1/2]; obj.chat=sum(obj.Ahat,2)';
                    
                case 'SSP54' 
                    % Low storage Runge-Kutta coefficients
                    obj.s=5; obj.r=4;
                    obj.Ahat = [        0.0 ...
                        -567301805773.0/1357537059087.0 ...
                        -2404267990393.0/2016746695238.0 ...
                        -3550918686646.0/2091501179385.0  ...
                        -1275806237668.0/842570457699.0];
                    obj.bhat = [ 1432997174477.0/9575080441755.0 ...
                        5161836677717.0/13612068292357.0 ...
                        1720146321549.0/2090206949498.0  ...
                        3134564353537.0/4481467310338.0  ...
                        2277821191437.0/14882151754819.0];
                    obj.chat = [        0.0  ...
                        1432997174477.0/9575080441755.0 ...
                        2526269341429.0/6820363962896.0 ...
                        2006345519317.0/3224310063776.0 ...
                        2802321613138.0/2924317926251.0];
                    
                  %=================Implicit Methods=========================
                  
                case 'BE11'
                    %Backward Euler
                    obj.s=1; obj.r=1.e10;
                    obj.A=1;
                    obj.b=1; obj.c=1;
                    
                case 'SDIRK34' %3-stage, 4th order singly diagonally implicit (SSP)
                    obj.s=3; obj.r=1.7588;
                    g=0.5*(1-cos(pi/18)/sqrt(3)-sin(pi/18));
                    q=(0.5-g)^2;
                    obj.A=[g     0    0
                        0.5-g g    0
                        2*g  1-4*g g];
                    obj.b=[1/(24*q) 1-1/(12*q) 1/(24*q)];
                    obj.c=sum(obj.A,2)';
                                       
                case 'ESDIRK4' % From ref [4], Derivation can be found in ref [3]
                    obj.s=6; obj.r=4;
                    obj.A = [ ...
                        0, 0, 0, 0, 0, 0;
                        1/4,                       1/4, 0, 0, 0, 0;
                        8611/62500,               -1743/31250,         1/4, 0, 0, 0;
                        5012029/34652500,         -654441/2922500,     174375/388108,       1/4, 0, 0;
                        15267082809/155376265600, -71443401/120774400, 730878875/902184768, 2285395/8070912,	1/4, 0;
                        82889/524892,              0,                  15625/83664,         69875/102672,    -2260/8211, 1/4;];
                    obj.b = [...
                        82889.0/524892.0,...
                            0.0, ...
                        15625.0/83664.0,...
                        69875.0/102672.0,...
                        -2260.0/8211.0,...
                            1.0/4.0];
                    obj.bhat = [...
                        4586570599.0/29645900160.0,...
                                 0.0, ...
                         178811875.0/945068544.0,...
                         814220225.0/1159782912.0,...
                          -3700637.0/11593932.0,...
                             61727.0/225920.0];
                    obj.c=sum(obj.A,2);
                    
                    %=================Classical Methods=========================
                    
                    %Gauss-Legendre methods -- order 2s
                case 'GL1'
                    obj.r=2; obj.A=1/2; obj.b=1; obj.c=1/2;
                case 'GL2'
                    obj.r=0;
                    obj.A=[1/4 1/4-sqrt(3)/6
                           1/4+sqrt(3)/6 1/4];
                    obj.b=[1/2 1/2];
                    obj.c=[1/2-sqrt(3)/6 1/2+sqrt(3)/6];
                case 'GL3'
                    obj.r=0;
                    obj.A=[5/36 (80-24*sqrt(15))/360 (50-12*sqrt(15))/360
                           (50+15*sqrt(15))/360 2/9  (50-15*sqrt(15))/360
                           (50+12*sqrt(15))/360 (80+24*sqrt(15))/360 5/36];
                    obj.b=[5/18 4/9 5/18];
                    obj.c=[(5-sqrt(15))/10 1/2 (5+sqrt(15))/10];
                    
                    %Radau IA methods -- order 2s-1
                case 'RIA1'
                    obj.r=1; obj.A=1; obj.b=1; obj.c=0;
                case 'RIA2'
                    obj.r=0; 
                    obj.A=[1/4 -1/4
                           1/4 5/12];
                    obj.b=[1/4 3/4];
                    obj.c=[0 2/3];
                case 'RIA3'
                    obj.r=0;
                    obj.A=[1/9 (-1-sqrt(6))/18 (-1+sqrt(6))/18
                           1/9 (88+7*sqrt(6))/360 (88-43*sqrt(6))/360
                           1/9 (88+43*sqrt(6))/360 (88-7*sqrt(6))/360];
                    obj.b=[1/9 (16+sqrt(6))/36 (16-sqrt(6))/36];
                    obj.c=[0 (6-sqrt(6))/10 (6+sqrt(6))/10];
                    
                    %Radau IIA methods -- order 2s-1
                case 'RIIA1'
                    obj.r=1; obj.A=1; obj.b=1; obj.c=1;
                case 'RIIA2'
                    obj.r=0;
                    obj.A=[5/12 -1/12
                           3/4 1/4];
                    obj.b=[3/4 1/4];
                    obj.c=[1/3 1];
                case 'RIIA3'
                    obj.r=0;
                    obj.A=[(88-7*sqrt(6))/360 (296-169*sqrt(6))/1800 (-2+3*sqrt(6))/225
                           (296+169*sqrt(6))/1800 (88+7*sqrt(6))/360 (-2-3*sqrt(6))/225
                           (16-sqrt(6))/36 (16+sqrt(6))/36 1/9];
                    obj.b=[(16-sqrt(6))/36 (16+sqrt(6))/36 1/9 ];
                    obj.c=[(4-sqrt(6))/10 (4+sqrt(6))/10 1];
                    
                    %Lobatto IIIA methods -- order 2s-2
                case 'LIIIA2'
                    obj.r=0; 
                    obj.A=[0 0
                           1/2 1/2];
                    obj.b=[1/2 1/2];
                    obj.c=[0 1];
                case 'LIIIA3'
                    obj.r=0;
                    obj.A=[0 0 0
                           5/24 1/3 -1/24
                           1/6 2/3 1/6];
                    obj.b=[1/6 2/3 1/6];
                    obj.c=[0 12 1];
                    
                    %==============Implicit-Explicit Methods====================
                    
                case 'IMEX-SSP1' % From ref[5]: Not realy and RK scheme,                    
                    obj.s=1; obj.r=1;
                    obj.Ahat=0; 
                    obj.A=1;
                    obj.bhat=1;
                    obj.b=1;
                    
                case 'IMEX-SSP2' % From ref[1]: IMEX-SSP2(2,2,2)
                    gamma = 1-1/sqrt(2);
                    obj.s=2; obj.r=2;
                    obj.Ahat=[0,0;1,0];
                    obj.A=[gamma,0;1-2*gamma,gamma];
                    obj.bhat=[1/2 1/2];
                    obj.b=[1/2 1/2];
                    obj.chat=[0 1];
                    obj.c=[0 1-gamma];
                    
                case 'IMEX-SSP3' % From ref[1]: IMEX-SSP3(3,3,2)
                    gamma = 1-1/sqrt(2);
                    obj.s=3; obj.r=3;
                    obj.Ahat=[0,0,0;1,0,0;1/4,1/4,0];
                    obj.A=[gamma,0,0;1-2*gamma,gamma,0;1/2-gamma,0,gamma];
                    obj.bhat=[1/6 1/6 2/3];
                    obj.b=[1/6 1/6 2/3];
                    obj.chat=[0 1 1/2];
                    obj.c=[gamma, 1-gamma, 1/2];
                    
                case 'ARK3' % From ref[2]: ARK3(2)4L[2]SA-ERK & ARK3(2)4L[2]SA-ESDIRK
                    obj.s=4; obj.r=3;
                    obj.Ahat = [ ...
                        0, 0, 0, 0;
                        1767732205903/2027836641118, 0,0,0;
                        5535828885825/10492691773637, 788022342437/10882634858940, 0, 0;
                        6485989280629/16251701735622, -4246266847089/9704473918619, 10755448449292/10357097424841, 0];
                    obj.A = [ ...
                        0, 0, 0, 0;
                        1767732205903/4055673282236, 1767732205903/4055673282236,0,0;
                        2746238789719/10658868560708,-640167445237/6845629431997,1767732205903/4055673282236,0;
                        1471266399579/7840856788654,-4482444167858/7529755066697,11266239266428/11593286722821,1767732205903/4055673282236];
                    obj.bhat = [1471266399579/7840856788654,...
                        -4482444167858/7529755066697,...
                        11266239266428/11593286722821,...
                        1767732205903/4055673282236];
                    obj.b = obj.bhat;
                    obj.btilde = [2756255671327/12835298489170,...
                        -10771552573575/22201958757719,...
                        9247589265047/10645013368117,...
                        2193209047091/5459859503100];
                    obj.c = [1767732205903/2027836641118,3/5,1];

                case 'ARK4' % From ref[2]: ARK4(3)6L[2]SA-ERK & ARK4(3)6L[2]SA-ESDIRK
                    obj.s=6; obj.r=4;
                    obj.Ahat = [ ...
                        0, 0, 0, 0, 0, 0;
                        1/2, 0, 0, 0, 0, 0;
                        13861/62500,                 6889/62500, 0, 0, 0, 0;
                       -116923316275/2393684061468, -2731218467317/15368042101831,  9408046702089/11113171139209,   0, 0, 0;
                       -451086348788/2902428689909, -2682348792572/7519795681897,   12662868775082/11960479115383,  3355817975965/11060851509271, 0, 0;
                        647845179188/3216320057751,  73281519250/8382639484533,     552539513391/3454668386233,     3354512671639/8306763924573,  4040/17871, 0];
                    obj.A = [ ...
                        0, 0, 0, 0, 0, 0;
                        1/4,                       1/4, 0, 0, 0, 0;
                        8611/62500,               -1743/31250,         1/4, 0, 0, 0;
                        5012029/34652500,         -654441/2922500,     174375/388108,       1/4, 0, 0;
                        15267082809/155376265600, -71443401/120774400, 730878875/902184768, 2285395/8070912,	1/4, 0;
                        82889/524892,              0,                  15625/83664,         69875/102672,    -2260/8211, 1/4;];
                    obj.bhat = [...
                        82889.0/524892.0,...
                            0.0, ...
                        15625.0/83664.0,...
                        69875.0/102672.0,...
                        -2260.0/8211.0,...
                            1.0/4.0];
                    obj.b = obj.bhat;
                    obj.btilde = [...
                        4586570599.0/29645900160.0,...
                                 0.0, ...
                         178811875.0/945068544.0,...
                         814220225.0/1159782912.0,...
                          -3700637.0/11593932.0,...
                             61727.0/225920.0];
                    obj.c = [0, 1/2, 83/250, 31/50, 17/20, 1];
                    
                    %===================Miscellaneous Methods================
                    
                case 'Mid22'
                    %Midpoint 22 method
                    obj.s=2; obj.r=0.5;
                    obj.Ahat=[0  0
                              1/2 0];
                    obj.bhat=[0 1]; obj.chat=[0 1/2];
                    
                case 'MTE22'
                    %Minimal truncation error 22 method (Heun)
                    obj.s=2; obj.r=0.5;
                    obj.Ahat=[0  0
                              2/3 0];
                    obj.bhat=[1/4 3/4]; obj.chat=[0 2/3];
                    
                case 'CN22'
                    %Crank-Nicholson
                    obj.s=2; obj.r=2;
                    obj.A=[0  0
                           1/2 1/2];
                    obj.b=[1/2 1/2]; obj.c=[0 1];
                    
                case 'Heun33'
                    obj.s=3; obj.r=0;
                    obj.Ahat=[0 0 0; 1/3 0 0; 0 2/3 0];
                    obj.bhat=[1/4 0 3/4]; obj.chat=sum(obj.A,2)';
                    
                case 'RK44'  %Classical fourth order
                    obj.s=4; obj.r=0;
                    obj.Ahat=[0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
                    obj.bhat=[1/6 1/3 1/3 1/6]; obj.chat=sum(obj.Ahat,2)';
                otherwise
                    error('method not listed :( ')
            end
            disp('See refs for implementation and details. MD 2014 ;)');
            
        end % Constructor 
    end % Methods 
end % Class

