classdef Parameters
    % Simulation parameters which are independent of the grid 
    % are stored here. They are shared by all scripts.
    % Methods are the mode functions which are also shared
    
    properties
        ltyp = 100e-6; %typical length scale of the system
        Lz; %propagation direction
        d = 50e-6; %inter-plate distance
        A; %area
        t0 = 0.5e-12; %illumination start
        tau = 100e-15; %plasma scattering time (100fs)
        nSi = 3.418; %silicone index of refraction at THz frequencies
        ng = 3.874; %silicone group index at 1030 nm
        c = 299792458; %speed of light
        cSi = 299792458/3.418; %speed of light in silicone at THz frequencies
        vg = 299792458/3.874; %speed of light in silicone at 1030nm
        eps = 3.418^2*8.8541878e-12; % permittivity of silicone at THz frequencies
        hbar = 1.0545718e-34; % reduced Planck's constant
        e = 1.6021766e-19; %electron charge
        meff = 0.197*9.10938356e-31; %electron effective mass
        kB=1.3806e-23; %Boltzmann's constant
        T=294; %Temperature in K
        n; %max number of ny modes considered
        m; %number of kz modes for each ny
        M = []; %vector of m
        N = []; %vector of n, same size as p
        maxf; %maximum cutoff frequency
        N_TM; %total number of TM modes being considered
        dz; %mesh size in z direction
        dy; %mesh size in y direction
        z; %vector in z size (1, Lz/dz+1)
        y; %vector in y size (d/dy+1, 1
        wTM; %vector of TM frequencies [0nz 1nz 2nz 3nz]
    end %properties
    
    methods
        function obj = Parameters(max_f, Lz)
            obj.Lz = Lz;
            obj.A = obj.Lz*obj.d;
            obj.maxf = max_f;
            obj.n = 0:floor(2*obj.d*max_f/obj.cSi);
            obj.m = floor(real(obj.Lz/pi*sqrt((max_f*2*pi/obj.c*obj.nSi)^2 - (obj.n*pi/obj.d).^2)));
            obj.N_TM = sum(obj.m);
            obj.dy = obj.d/(sum(obj.m>0)-1)/16;
            obj.dz = obj.Lz/max(obj.m)/16;
            obj.z = 0:obj.dz:obj.Lz;
            if length(obj.n)>1 %if ny>0
                obj.y = (0:obj.dy:obj.d)';
            else %if ny=0, constant along y
                obj.y=0;
            end
            for i=1:length(obj.m)
                obj.M = [obj.M 1:obj.m(i)];
                obj.N = [obj.N obj.n(i)*ones(1,obj.m(i))];
            end
            obj.wTM = obj.w(obj.N,obj.M);
        end %Parameters

        function out = NBE(obj,w) %Bose Einstein distrib
            out = 1./(exp(obj.hbar/obj.kB/obj.T*w)-1);
        end
        
        function out = w(obj,n,m) %TM mode frequency
            out = obj.cSi*sqrt((n*pi/obj.d).^2 + (m*pi/obj.Lz).^2);
        end
        
        function out = kz(obj,n) %wave vector in z
            out = n*pi/obj.Lz;
        end
        
        function out = ky(obj,n) %wave vector in y
            out = n*pi/obj.d;
        end
        
        function out = k(obj,ny,nz) %total wave vector
            out = sqrt(obj.ky(ny).^2 + obj.kz(nz).^2);
        end
        
        function out = A_TMy(obj, ny, nz, y, z) %Modes in y-direction
            if ny==0 %TEM modes
                out = -sqrt(2/obj.A).*sin(obj.kz(nz).*z);
            else %TM modes
                out = -sqrt(4/obj.A)*obj.kz(nz)./obj.k(ny,nz).*sin(obj.kz(nz).*z).*cos(obj.ky(ny)*y);
            end
        end %A_TMz
        
        function out = A_TMz(obj, ny, nz, y, z) %Modes in z-direction
            if ny==0 %TEM modes return array of zeros with correct size
                out = zeros(size(y)).*zeros(size(z)); 
            else %TM modes
                out = -1j*sqrt(4/obj.A)*obj.ky(ny)./obj.k(ny,nz).*sin(obj.kz(nz).*z).*sin(obj.ky(ny)*y);
            end
        end %A_TMz
        
%%%%%%%% FDTD Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = A_y(obj,u)
            [Z,Y] = meshgrid(obj.z,obj.y(1:end-1)+obj.dy/2);
            Ay = complex(zeros(size(Z))); %allocate memory
            for i=1:length(obj.M) %sum over modes
               Ay = Ay + obj.A_TMy(obj.N(i),obj.M(i),Y,Z)*u(i); 
            end
            out = Ay;
        end %E_y
        
        function out = A_z(obj,u)
            [Z,Y] = meshgrid(obj.z(1:end-1)+obj.dz/2,obj.y);
            Az = complex(zeros(size(Z))); %allocate memory
            for i=1:length(obj.M) %sum over modes
               Az = Az + obj.A_TMz(obj.N(i),obj.M(i),Y,Z)*u(i); 
            end
            out = Az;
        end %E_z

        function out = E_y(obj,u)
            [Z,Y] = meshgrid(obj.z,obj.y(1:end-1)+obj.dy/2);
            Ey = complex(zeros(size(Z))); %allocate memory
            for i=1:length(obj.M) %sum over modes
               Ey = Ey - obj.A_TMy(obj.N(i),obj.M(i),Y,Z)*u(i); 
            end
            out = real(Ey);
        end %E_y
        
        function out = E_z(obj,u)
            [Z,Y] = meshgrid(obj.z(1:end-1)+obj.dz/2,obj.y);
            Ez = complex(zeros(size(Z))); %allocate memory
            for i=1:length(obj.M) %sum over modes
               Ez = Ez - obj.A_TMz(obj.N(i),obj.M(i),Y,Z)*u(i); 
            end
            out = real(Ez);
        end %E_z
        
        function out = H_x(obj,u)
            Ay = obj.A_y(u); Az = obj.A_z(u);
            Hx = (Az(2:end, :) - Az(1:end-1, :))/obj.dy - (Ay(:, 2:end) - Ay(:, 1:end-1))/obj.dz;
            out = real(Hx);
        end %Hx
        
    end %methods
    
end %classdef