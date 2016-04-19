function main()
    global force system filename SigmaMed SIGMASLOPE SIGMA0 ASPECTRATIO;
    global R MU FLARINGINDEX NRAD NSEC EnergyMed GuidingCenter;
    global PhysicalTime LAMBDADOUBLING TRANSITIONRADIUS TRANSITIONRATIO;
    global TRANSITIONWIDTH PhysicalTimeInitial ADIABATICINDEX;
    
    ADIABATICINDEX = 1.4;
    TRANSITIONRADIUS = 0.0;
    TRANSITIONRATIO = 1.0;
    TRANSITIONWIDTH = 1.0;
    LAMBDADOUBLING = 0.0;
    PhysicalTimeInitial = 0.0;
    PhysicalTime =0.0;
    FREQUENCY = 2;
    OMEGAFRAME = 1.0;
    Corotating = 1;
    GuidingCenter = 1; % YES
    
    R = 1.0; % Mean molecular weight
    MU = 1.0; % Universal Gas Constant in code units 
    ASPECTRATIO = 0.05;
    SIGMASLOPE = 0.0;
    SIGMA0 = 6.3661977237e-4;
    FLARINGINDEX = 0.0;
  
    filename = '/in/Jup.cfg';
    NRAD = 128;
    NSEC= 384;
    adiabatic = true;
    selfgravity = false;
    
    gas_density     = zeros([NRAD,NSEC]);
    gas_v_rad       = zeros([NRAD,NSEC]);
    gas_v_theta     = zeros([NRAD,NSEC]);
    gas_energy      = zeros([NRAD,NSEC]);
    gas_label       = zeros([NRAD,NSEC]); 
    
    dimfxy = 11;
    FillPolar1DArrays();
    AllocateForce(dimfxy);

    planets = InitPlanetarySystem(filename);
    sprintf('%d\n',planets);
    
    % InitGasDensity
    FillSigma();
    for j=1:NSEC
        gas_density(:,j) = SigmaMed;
    end
    
    if (adiabatic)
        % InitGasEnergy
        FillEnergy();
        for j=1:NSEC
            gas_energy(:,j) = EnergyMed;
        end
    end
    
    if (selfgravity)
        % If SelfGravity = YES or Z, planets are initialized feeling disk
        % potential. Only the surface density is required to calculate
        % the radial self-gravity acceleration. The disk radial and
        % azimutal velocities are not updated 
    end
    
    %system;
    %ListPlanets(system);
    
    OmegaFrame = OMEGAFRAME;
    
    if (Corotating == true) 
        OmegaFrame = GetsPsysInfo(system, FREQUENCY);
    end
    
    sprintf('%.18ld\n',OmegaFrame);
    
    Initialization (gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, system);
    force;
    %force{9,1}(2) especifico

    %for i=1:8
    %    system{i,1}
    %end
    
end