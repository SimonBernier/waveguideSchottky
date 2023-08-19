function FDTD_tilted_Fig2(maxf, v, Es, Lz, Profile, tskip, zskip)
    % 2D FDTD on Staggered grid
    % Solving Maxwell's Equations with the plasma source term 
    % dJ/dt + gamma*J = ne(y,z,t)*e^2/me*E = eps*wp^2*E
    % Using experimental pulse as a source.

    % padding vectors for experimental source
    Npad_l = length(0:5.2e-6:0.5e-3); Npad_r = length(0.5e-3+1280*5.2e-6:5.2e-6:Lz+5.2e-6); 
    %The experimental profiles are numbered 1:13. 1:7 is varying start, 8:13 is varying stop
    start = [0.17 0.17 0.17 0.63 1.87 2.80 3,66]*1e-3+0.5e-3;
    stop  = [3.11 4.28 5.59 5.41 5.51 5.51 5.51]*1e-3+0.5e-3;
    a=start(Profile); b=stop(Profile); % get start and stop from profile #

    % parameters and physical units
    PA = Parameters(maxf, Lz);
    PU = PhysicalUnits(PA.ltyp);

    %Grid parameters
    dz = PU.ConvertSILengthToFDUnits(PA.dz); dy = PU.ConvertSILengthToFDUnits(PA.dy);
    y = PU.ConvertSILengthToFDUnits(PA.y); d=PU.ConvertSILengthToFDUnits(PA.d);
    z = PU.ConvertSILengthToFDUnits(PA.z);
    nz = length(z)-1; ny = length(y)-1; [Z,Y] = meshgrid(z,y);

    %time stepping
    S = 0.99; %stability factor
    dt= PA.nSi/sqrt(1/dz^2 + 1/dy^2)*S; %stability condition
    runTime = PA.t0 + (b-a)/min(v*PA.cSi, PA.cSi) + (7e-3-b)/PA.cSi; %adaptive runTime
    fprintf(sprintf('Runtime = %0.2f ps\n', runTime*1e12))
    T = PU.ConvertSITimeToFDUnits(runTime);
    n_time = ceil(T/dt); 

    %Display simulation and grid properties
    fprintf('nt = %d, ny = %d, nz = %d\n', n_time, ny, nz)
    fprintf('dt = %0.5e, dy = %0.5e, dz = %0.5e \n', dt, dy, dz)
    fprintf('Accuracy: time ~ %0.3e, y ~ %0.3e, z ~ %0.3e\n\n', dt^2, dy^2, dz^2)

    % Parameters
    LD = 2.8e-6; vcSi = v/PA.nSi;
    gamma = 1/PU.ConvertSITimeToFDUnits(PA.tau); %plasma scattering rate
    LD=PU.ConvertSILengthToFDUnits(LD); %dopant parameters
    tpulse = 190e-15; sigma= PU.ConvertSITimeToFDUnits( tpulse )/sqrt(8*log(2)); %plasma profile
    vg = 1/PA.ng; t0 = PU.ConvertSITimeToFDUnits(PA.t0);
    a=PU.ConvertSILengthToFDUnits(a);

    % Plasma functions to integrate with FDTD
    YA = Y(1:ny,:)+dy/2; ZA = Z(1:ny,:); YB = Y(:,1:nz); ZB = Z(:,1:nz)+dz/2; %preallocate memory
    H = @(x) 0.5*(1+erf(x/(sqrt(2)*sigma))); %plasma profile
    Es = PU.ConvertSIElectricFieldToFDUnits(Es);
    Elocal = Es*(-heaviside(LD-YA).*(LD-YA)/LD + heaviside(YA-d+LD).*(LD-d+YA)/LD); % Schottky barrier electric field

    P = readmatrix(sprintf('Fig2_P%d.dat',Profile)); %import exp intensity profile
    P = (2*pi*PU.ConvertSIFrequencyToFDUnits(P)).^2; %make this wp^2 = (2*pi*fp)^2;
    P = [zeros(1, Npad_l) P zeros(1, Npad_r)]; %pad with zeros on both sides (add 0.5mm on the left)
    Py = interp1(PU.ConvertSILengthToFDUnits(-5.2e-6:5.2e-6:Lz+5.2e-6),P,z);
    Pz = interp1(PU.ConvertSILengthToFDUnits(-5.2e-6:5.2e-6:Lz+5.2e-6),P,z(1:end-1)+dz/2);
    By = @(t,y,z) Py.*H((t-t0)-(d-y)/vg-(z-a)/vcSi); % gaussian excitation profile
    Bz = @(t,y,z) Pz.*H((t-t0)-(d-y)/vg-(z-a)/vcSi); % gaussian excitation profile

    %precompute constants for FDTD
    eps_r = PA.nSi^2; wp2dt=dt/(1+0.5*dt*gamma); %just for notation
    dt_eps=dt/eps_r; dt_dz = dt/dz; dt_dy = dt/dy; dt_dz_eps = dt/dz/eps_r;
    dt_dy_eps = dt/dy/eps_r; dtgamma=(1-0.5*dt*gamma)/(1+0.5*dt*gamma);

    %Preallocate memory
    EyTEM=zeros(ceil((n_time+1)/tskip),ceil(length(z)/zskip)); %used to store and plot data
    Ey = zeros(ny,nz+1); Ez = zeros(ny+1,nz); Hx = zeros(ny,nz); Jy=zeros(size(Ey)); Jz=zeros(size(Ez));
    EyTEM(1,:) = mean( Ey(:,1:zskip:end) ); %save initial conditions

    %FDTD time integration
    for tn=1:n_time 

        %eps*dE/dt = curl H - J
        Ey(:, 2:nz) = Ey(:, 2:nz) + dt_dz_eps*(Hx(:, 2:nz) - Hx(:, 1:nz-1)) - dt_eps*Jy(:, 2:nz);
        Ez(2:ny, :) = Ez(2:ny, :) - dt_dy_eps*(Hx(2:ny, :) - Hx(1:ny-1, :)) - dt_eps*Jz(2:ny, :);

        %dH/dt = -curl E
        Hx = Hx + dt_dz*(Ey(:, 2:nz+1) - Ey(:, 1:nz)) - dt_dy*(Ez(2:ny+1, :) - Ez(1:ny, :));

        %dJ/dt = -gamma J + eps*wp^2*E (second order)
        Jy = dtgamma*Jy + wp2dt*By(tn*dt, YA, ZA).*(Ey+Elocal);
        Jz = dtgamma*Jz + wp2dt*Bz(tn*dt, YB, ZB).*Ez;
        
        % store data
        if mod(tn,tskip)==0
            EyTEM(tn/tskip+1,:) = mean( Ey(:,1:zskip:end) ); %average over y
        end

    end

    writematrix(EyTEM, sprintf('Ey_Fig2_v%0.2f_P%d.dat',v,Profile))
