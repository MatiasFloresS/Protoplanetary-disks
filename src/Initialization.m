function Initialization(gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, pla_sys)

    %ReadPrevDim() Esta funcion leera dims.dat y used_rad.dat, mas adelante
    %se vera eso
    
    InitEuler(gas_v_rad, gas_v_theta, gas_density, gas_energy);
    
    
end

function z = InitEuler(Vr, Vt, Rho, Energy)
    global NRAD NSEC Rmed CellAbscissa CellOrdenate;
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
    
    InitComputeAccel();
    ComputeSoundSpeed(Rho, Energy, SoundSpeed);
    z=0;
    
end

function InitComputeAccel()

    global CellAbscissa CellOrdinate NRAD NSEC Rmed;
   
    CellAbscissa = zeros([NRAD,NSEC]);
    CellOrdinate = zeros([NRAD,NSEC]);
    
    
    for i=1:NRAD
        for j=1:NSEC
            CellAbscissa(i,j) = Rmed(i).* cos(2.0*pi*(j)/NSEC);
            CellOrdinate(i,j) = Rmed(i).* sin(2.0*pi*(j)/NSEC);
        end
    end
end

function ComputeSoundSpeed(Rho, Energy, SoundSpeed)

    global Rmed;
    aspectratio = AspectRatio(Rmed);

end


function aspectratio = AspectRatio(rad)

    global ASPECTRATIO PhysicalTime TRANSITIONRADIUS TRANSITIONWIDTH TRANSITIONRATIO;
    global LAMBDADOUBLING;
    aspectratio = ASPECTRATIO;
    
    rmin = TRANSITIONRADIUS-TRANSITIONWIDTH*ASPECTRATIO;
    rmax = TRANSITIONRADIUS+TRANSITIONWIDTH*ASPECTRATIO;
    scale = 1.0+(PhysicalTime-PhysicalTimeInitial)*LAMBDADOUBLING;
    rmin = rmin * scale;
    rmax = rmax * scale;
    
    if (rad < rmin) 
        aspectratio = aspectratio * TRANSITIONRATIO;
    end
    
    if ((rad >= rmin) && (rad <= rmax))
        aspectratio = aspecrtatio *exp((rmax-rad)/(rmax-rmin)*log(TRANSITIONRATIO));
    end
    
end