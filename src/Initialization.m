function [gas_v_rad, gas_v_theta, gas_label] =Initialization(gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label)
    global RESTART system;
    %ReadPrevDim() Esta funcion leera dims.dat y used_rad.dat, mas adelante
    %se vera eso
    
    [gas_v_rad, gas_v_theta] = InitEuler(gas_v_rad, gas_v_theta, gas_density, gas_energy);
    gas_label = InitLabel(gas_label);

    if (RESTART == 1)
       % checkrebin and others functions 
    end
    %WriteDim();   
    %initPotential();
end

function [gas_v_rad, gas_v_theta] = InitEuler(gas_v_rad, gas_v_theta, gas_density, gas_energy)
    global NRAD NSEC SoundSpeed Pressure Temperature Potential;
    %InitTransport()
    %RadMomP         = zeros([NRAD,NSEC]);
    %RadMomM         = zeros([NRAD,NSEC]);
    %ThetaMomP       = zeros([NRAD,NSEC]);
    %ThetaMomM       = zeros([NRAD,NSEC]);
    %Work            = zeros([NRAD,NSEC]);
    %QRStar          = zeros([NRAD,NSEC]);
    %ExtLabel        = zeros([NRAD,NSEC]);
    %VthetaRes       = zeros([NRAD,NSEC]);
    %Elongations     = zeros([NRAD,NSEC]);
    %TempShift       = zeros(1,NRAD*NSEC);
    %dq              = zeros(1,NRAD*NSEC);
    
    %InitViscosity()
    %DivergenceVelocity = zeros([NRAD,NSEC]);
    %DRR                = zeros([NRAD,NSEC]);
    %DRP                = zeros([NRAD,NSEC]);
    %DPP                = zeros([NRAD,NSEC]);
    %TAURR              = zeros([NRAD,NSEC]);
    %TAURP              = zeros([NRAD,NSEC]);
    %TAUPP              = zeros([NRAD,NSEC]);
    
    % 
    %RhoStar      = zeros([NRAD,NSEC]);
    %RhoInt       = zeros([NRAD,NSEC]);
    %VradNew      = zeros([NRAD,NSEC]);
    %VradInt      = zeros([NRAD,NSEC]);
    %VthetaNew    = zeros([NRAD,NSEC]);
    %VthetaInt    = zeros([NRAD,NSEC]);
    %EnergyNew    = zeros([NRAD,NSEC]);
    %EnergyInt    = zeros([NRAD,NSEC]);
    %TemperInt    = zeros([NRAD,NSEC]);
    Potential    = zeros([NRAD+1,NSEC]);
    Pressure     = zeros([NRAD+1,NSEC]);
    SoundSpeed   = zeros([NRAD+1,NSEC]);
    Temperature  = zeros([NRAD+1,NSEC]);
    %Qplus        = zeros([NRAD,NSEC]);
    InitComputeAccel();
    ComputeSoundSpeed(gas_density, gas_energy);
    ComputePressureField(gas_density, gas_energy);
    ComputeTemperatureField(gas_density, gas_energy);
    [gas_v_rad, gas_v_theta] = InitGasVelocities (gas_v_rad, gas_v_theta);
end

function InitComputeAccel()
    global CellAbscissa CellOrdinate NRAD NSEC Rmed;
    
    %for j=1:NSEC
    %    CellAbscissa(1:NRAD,j) = Rmed.* cos(2.0*pi*(j-1)/NSEC);
    %    CellOrdinate(1:NRAD,j) = Rmed.* sin(2.0*pi*(j-1)/NSEC);
    %end
    CellAbscissa = Rmed'* cos(2.0*pi*(0:NSEC-1)/NSEC);
    CellOrdinate = Rmed'* sin(2.0*pi*(0:NSEC-1)/NSEC);
    CellAbscissa(NRAD+1,:) = 0.0;
    CellOrdinate(NRAD+1,:) = 0.0;
    
end

function ComputeSoundSpeed(gas_density, gas_energy)
    global Rmed NSEC FLARINGINDEX G Adiabatic SoundSpeed ADIABATICINDEX NRAD; 
    
    aspectratio = zeros([1,NRAD]);
    aspectratio(1:NRAD) = AspectRatio(Rmed);

    if (strcmp(Adiabatic,'NO'))        
        for j=1:NSEC
            SoundSpeed(1:NRAD,j) = aspectratio.*sqrt(G*1.0./Rmed).*(Rmed.^FLARINGINDEX);
        end    
    else
       SoundSpeed = sqrt(ADIABATICINDEX*(ADIABATICINDEX-1.0).*gas_energy(1:NRAD)./gas_density(1:NRAD));
    end    

end

function aspectratio = AspectRatio(rad)
    global ASPECTRATIO PhysicalTime TRANSITIONRADIUS TRANSITIONWIDTH TRANSITIONRATIO;
    global LAMBDADOUBLING PhysicalTimeInitial;
    
    aspectratio = ASPECTRATIO;
    rmin = TRANSITIONRADIUS-TRANSITIONWIDTH*ASPECTRATIO;
    rmax = TRANSITIONRADIUS+TRANSITIONWIDTH*ASPECTRATIO;
    scale = 1.0+(PhysicalTime-PhysicalTimeInitial)*LAMBDADOUBLING;
    rmin = rmin * scale;
    rmax = rmax * scale;
    if (rad < rmin) 
        aspectratio = aspectratio * TRANSITIONRATIO;
    end
    
    if ((rad >= rmin) & (rad <= rmax))
        aspectratio = aspectratio *exp((rmax-rad)/(rmax-rmin)*log(TRANSITIONRATIO));
    end  
end

function ComputePressureField(gas_density, gas_energy)
    global ADIABATICINDEX Adiabatic SoundSpeed Pressure NRAD;    
    if (strcmp(Adiabatic,'NO'))
        Pressure(1:NRAD,:) = gas_density(1:NRAD,:).*SoundSpeed(1:NRAD,:).*SoundSpeed(1:NRAD,:);
    else
        Pressure(1:NRAD,:) = (ADIABATICINDEX-1.0).*gas_energy(1:NRAD,:);
    end
end

function ComputeTemperatureField(gas_density, gas_energy)
    global Temperature Adiabatic MU R ADIABATICINDEX Pressure NRAD;
    if (strcmp(Adiabatic,'NO'))
        Temperature(1:NRAD,:) = MU/R.*Pressure(1:NRAD,:)./gas_density(1:NRAD,:);
    else
        Temperature(1:NRAD,:) = MU/R*(ADIABATICINDEX-1.0).*gas_energy(1:NRAD,:)./gas_density(1:NRAD,:);
    end
end

function gas_label = InitLabel(gas_label)
    global NSEC NRAD Rmed system;
    
    x = zeros([NRAD,NSEC]);
    y = zeros([NRAD,NSEC]);
    xp = system{3,1}(1);
    yp = system{4,1}(1);
    rp = sqrt(xp*xp + yp*yp);
    rhill = rp * (system{2,1}(1)/3)^(1./3.);
    
    % Initialize label as you wish. In this example, label only takes
    % into account fluid elements inside the planet's Hill Sphere
    
    for i=1:NRAD
        x(i,:) = Rmed(i).* cos((0:NSEC-1)/2.0*pi*NSEC);
        y(i,:) = Rmed(i).* sin((0:NSEC-1)/2.0*pi*NSEC);        
    end
    
    distance = sqrt((x-xp).*(x-xp) + (y-yp).*(y-yp));
    for i=1:NRAD
        if (distance(i, :) < rhill)
            gas_label(i,:) = 1.0;
        else
            gas_label(i,:) = 0.0;
        end        
    end
end

function WriteDim()
    global RMAX NTOT NINTERM NRAD NSEC
    fileID = fopen('dims.dat','w');
    fprintf(fileID,'%d\t%d\t\t%d\t%d\t%f\t%d\t%d\t%d\n',0,0,0,0,RMAX, int16(NTOT/NINTERM), NRAD, NSEC);
    fclose(fileID);
end

function [gas_v_rad, gas_v_theta] = InitGasVelocities(gas_v_rad, gas_v_theta)
    global NRAD Rmed Rinf G SELFGRAVITY NSEC ASPECTRATIO FLARINGINDEX SIGMASLOPE OmegaFrame;
    global IMPOSEDDISKDRIFT SIGMA0 SigmaInf;
    
    r = Rmed(1:NRAD);
    ri = Rinf(1:NRAD);
    r(NRAD+1) = Rmed(NRAD);
    ri(NRAD+1) = Rinf(NRAD);
    
    viscosity(1:NRAD+1) = Fviscosity(r);
    
    if (strcmp(SELFGRAVITY,'NO'))  
      omega = sqrt(G*1.0./r./r./r);
      for j=1:NSEC
        gas_v_theta(:,j) = omega.*r.*sqrt(1.0-(ASPECTRATIO.^2.0).*r.^(2.0*FLARINGINDEX)*(1.0+SIGMASLOPE-2.0*FLARINGINDEX)) - r*OmegaFrame;
      end   
    end
    
    gas_v_rad(NRAD+1,:) = 0.0;
    
    for j=1:NSEC
        gas_v_rad(1:NRAD, j) = IMPOSEDDISKDRIFT*SIGMA0./SigmaInf(1:NRAD)./ri(1:NRAD) - 3.0.*viscosity(1:NRAD)./r(1:NRAD).*(-SIGMASLOPE+0.5);
    end
    
    gas_v_rad(1,:) = 0.0;
    
end

function viscosity = Fviscosity(rad)
    global VISCOSITY CAVITYRADIUS CAVITYWIDTH ASPECTRATIO LAMBDADOUBLING PhysicalTime;
    global PhysicalTimeInitial CAVITYRATIO;
    
    viscosity = VISCOSITY;
    
    rmin = CAVITYRADIUS-CAVITYWIDTH*ASPECTRATIO;
    rmax = CAVITYRADIUS+CAVITYWIDTH*ASPECTRATIO;
    scale = 1.0 +(PhysicalTime-PhysicalTimeInitial)*LAMBDADOUBLING;
    rmin = rmin*scale;
    rmax = rmax*scale;
    
    if (rad < rmin) 
        viscosity = viscosity * CAVITYRATIO;
    end
    if ((rad >= rmin) & (rad <= rmax))
        viscosity = viscosity * exp((rmax-rad)/(rmax-rmin)*log(CAVITYRATIO));   
    end
end