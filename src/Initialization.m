function Initialization(gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, system)

    global RESTART;
    %ReadPrevDim() Esta funcion leera dims.dat y used_rad.dat, mas adelante
    %se vera eso
    
    InitEuler(gas_v_rad, gas_v_theta, gas_density, gas_energy);
    InitLabel(gas_label, system);
    
    if (RESTART == 1)
       % checkrebin and others functions 
    end
    WriteDim();
    
end

function z = InitEuler(Vr, Vt, gas_density, gas_energy)
    global NRAD NSEC;
    %InitTransport()
    RadMomP         = zeros([NRAD,NSEC]);
    RadMomM         = zeros([NRAD,NSEC]);
    ThetaMomP       = zeros([NRAD,NSEC]);
    ThetaMomM       = zeros([NRAD,NSEC]);
    Work            = zeros([NRAD,NSEC]);
    QRStar          = zeros([NRAD,NSEC]);
    ExtLabel        = zeros([NRAD,NSEC]);
    VthetaRes       = zeros([NRAD,NSEC]);
    Elongations     = zeros([NRAD,NSEC]);
    TempShift       = zeros(1,NRAD*NSEC);
    dq              = zeros(1,NRAD*NSEC);
    
    %InitViscosity()
    DivergenceVelocity = zeros([NRAD,NSEC]);
    DRR                = zeros([NRAD,NSEC]);
    DRP                = zeros([NRAD,NSEC]);
    DPP                = zeros([NRAD,NSEC]);
    TAURR              = zeros([NRAD,NSEC]);
    TAURP              = zeros([NRAD,NSEC]);
    TAUPP              = zeros([NRAD,NSEC]);
    
    % 
    RhoStar      = zeros([NRAD,NSEC]);
    RhoInt       = zeros([NRAD,NSEC]);
    VradNew      = zeros([NRAD,NSEC]);
    VradInt      = zeros([NRAD,NSEC]);
    VthetaNew    = zeros([NRAD,NSEC]);
    VthetaInt    = zeros([NRAD,NSEC]);
    EnergyNew    = zeros([NRAD,NSEC]);
    EnergyInt    = zeros([NRAD,NSEC]);
    TemperInt    = zeros([NRAD,NSEC]);
    Potential    = zeros([NRAD,NSEC]);
    Pressure     = zeros([NRAD,NSEC]);
    SoundSpeed   = zeros([NRAD,NSEC]);
    Temperature  = zeros([NRAD,NSEC]);
    Qplus        = zeros([NRAD,NSEC]);
    gas_label2   = zeros([NRAD,NSEC]);
    
    InitComputeAccel();
    ComputeSoundSpeed(gas_density, gas_energy);
    ComputePressureField(gas_density, gas_energy);
    ComputeTemperatureField(gas_density, gas_energy);
    %InitGasVelocities (Vr, Vt)
    z=0;
    
end

function InitComputeAccel()

    global CellAbscissa CellOrdinate NRAD NSEC Rmed;
   
    CellAbscissa = zeros([NRAD,NSEC]);
    CellOrdinate = zeros([NRAD,NSEC]);
    
    for i=1:NRAD
        CellAbscissa(i,:) = Rmed(i).* cos(2.0*pi*(0:NSEC-1)/NSEC);
        CellOrdinate(i,:) = Rmed(i).* sin(2.0*pi*(0:NSEC-1)/NSEC);
    end
    
end

function ComputeSoundSpeed(gas_density, gas_energy)

    global Rmed FLARINGINDEX G NSEC Adiabatic SoundSpeed ADIABATICINDEX;
    
    if (~Adiabatic)
        for j=1:NSEC
            SoundSpeed(:,j) = AspectRatio(Rmed).*sqrt(G*1.0./Rmed).*Rmed.^FLARINGINDEX;
        end
    
    else
       SoundSpeed = sqrt(ADIABATICINDEX*(ADIABATICINDEX-1.0).*gas_energy./gas_density);
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
    
    if (rad >= rmin & rad <= rmax)
        aspectratio = aspectratio *exp((rmax-rad)/(rmax-rmin)*log(TRANSITIONRATIO));
    end  
end

function ComputePressureField(gas_density, gas_energy)
    global ADIABATICINDEX Adiabatic SoundSpeed Pressure;
    
    if (~Adiabatic)
        Pressure = gas_density.*SoundSpeed.*SoundSpeed;
    else
        Pressure = (ADIABATICINDEX-1.0).*gas_energy;
    end
end

function ComputeTemperatureField(gas_density, gas_energy)
    
    global Temperature Adiabatic MU R ADIABATICINDEX Pressure;
    if (~Adiabatic)
        Temperature = MU/R.*Pressure./gas_density;
    else
        Temperature = MU/R*(ADIABATICINDEX-1.0).*gas_energy./gas_density;
    end
end

function InitLabel(gas_label,system)
    global NSEC NRAD Rmed;
    
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
    
    gas_label;
end

function WriteDim()


  %char filename[200];
  %FILE 	*dim;
  %if (!CPU_Master) return;
  %sprintf (filename, "%sdims.dat", OUTPUTDIR);
  %if ((dim = fopen (filename, "w")) == NULL) {
  %  fprintf (stderr, "Unable to open %s. Program stopped\n", filename);
  %  prs_exit (1);
  %}
  %fprintf (dim,"%d\t%d\t\t%d\t%d\t%f\t%d\t%d\t%d\n",0,0,0,0,RMAX, NTOT/NINTERM, GLOBALNRAD, NSEC);
  %fclose (dim);
%}
end