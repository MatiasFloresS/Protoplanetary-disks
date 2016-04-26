function main(ParameterFile,TellE)
    global force system  SigmaMed R MU  NRAD NSEC EnergyMed Corotating G ;
    global PhysicalTime PhysicalTimeInitial RESTART Adiabatic SELFGRAVITY OMEGAFRAME;
    
    ReadVariables(ParameterFile);
    if (strcmp(TellE,'YES'))
        TellEverything();
    end
    
    RESTART = 0;
    PhysicalTimeInitial = 0.0;
    PhysicalTime =0.0;
    FREQUENCY = 2;
    dimfxy = 11;
    
    G = 1.0;
    R = 1.0; % Mean molecular weight
    MU = 1.0; % Universal Gas Constant in code units 
    
    gas_density     = zeros([NRAD,NSEC]);
    gas_v_rad       = zeros([NRAD,NSEC]);
    gas_v_theta     = zeros([NRAD,NSEC]);
    gas_energy      = zeros([NRAD,NSEC]);
    gas_label       = zeros([NRAD,NSEC]); 
    
    FillPolar1DArrays();
    AllocateForce(dimfxy);
    InitPlanetarySystem();

    % InitGasDensity
    FillSigma();
    for j=1:NSEC
        gas_density(:,j) = SigmaMed;
    end
    
    if (Adiabatic)
        % InitGasEnergy
        FillEnergy();
        for j=1:NSEC
            gas_energy(:,j) = EnergyMed;
        end
    end
    
    if (SELFGRAVITY)
        % If SelfGravity = YES or Z, planets are initialized feeling disk
        % potential. Only the surface density is required to calculate
        % the radial self-gravity acceleration. The disk radial and
        % azimutal velocities are not updated 
    end
    %system;
    ListPlanets(system);
    OmegaFrame = OMEGAFRAME;
    
    if (strcmp(Corotating,'YES'))
        OmegaFrame = GetsPsysInfo(system, FREQUENCY);
    end
    
    Initialization (gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, system);
    force;
end

function ReadVariables(ParameterFile)
    cant_inputs = 500; % Se da memoria para 500 inputs, esto se puede modificar;
    A = cell(2,1);
    A{1,1} = [];
    A{2,1} = [];
    cero = zeros(1,cant_inputs);
    A{1,1} = cero;
    A{2,1} = cero;
    fprintf('%s\n',ParameterFile);
    c = 1;
    fid = fopen(ParameterFile, 'r');
    if fid == -1
        fprintf('Cannot open file: %s', ParameterFile)
        quit cancel
    else
        tline = fgetl(fid);
        while ischar(tline)
            matches = strfind(tline, '###');
            if( matches ~= 0)
                %do nothing
            else
                a = fscanf(fid, '%s', 1);
                if (~strcmp(a, ''))
                    if (~strcmp(a(1),'#'))
                        b = fscanf(fid, '%s', 1);
                        if ((b ~= str2double(b)) & (~isnan(str2double(b))))
                            b = str2double(b);  
                        end
                        A{c,1} = a;
                        A{c,2} = b;
                        c = c + 1;
                    end
                end
            end
            tline = fgetl(fid);
        end
        fclose(fid);        
    end
    InitVariables(A);
end

function InitVariables(A)
    % variables leidas
    global DT SIGMA0 NINTERM NTOT OUTPUTDIR INNERBOUNDARY LABELADVECTION;
    global TRANSPORT PLANETCONFIG MASSTAPER RADIALSPACING NRAD NSEC RMIN RMAX;
    global THICKNESSSMOOTHING ROCHESMOOTHING ASPECTRATIO VISCOSITY ALPHAVISCOSITY;
    global SIGMASLOPE RELEASERADIUS RELEASEDATE OMEGAFRAME DISK FRAME OUTERSOURCEMASS;
    global WRITEDENSITY WRITEVELOCITY WRITEENERGY WRITETEMPERATURE WRITEDIVV WRITEQPLUS;
    global INDIRECTTERM EXCLUDEHILL IMPOSEDDISKDRIFT FLARINGINDEX ECCENTRICITY CAVITYRADIUS;
    global CAVITYRATIO CAVITYWIDTH TRANSITIONRADIUS TRANSITIONRATIO TRANSITIONWIDTH;
    global LAMBDADOUBLING SELFGRAVITY CICPLANET FORCEDCIRCULAR ZMPLUS ADIABATIC;
    global ADIABATICINDEX COOLING COOLINGTIME0;
    
    % otras variables
    global FastTransport OpenInner NonReflecting Evanescent LogGrid Corotating;
    global GuidingCenter SelfGravity ZMPlus Write_Temperature Adiabatic SGZeroMode;
    
    B ={'DT', 'SIGMA0', 'NINTERM', 'NTOT', 'OUTPUTDIR', 'INNERBOUNDARY', 'LABELADVECTION', 'TRANSPORT', ...
        'PLANETCONFIG', 'MASSTAPER', 'RADIALSPACING', 'NRAD', 'NSEC', 'RMIN', 'RMAX', 'THICKNESSSMOOTHING', ...
        'ROCHESMOOTHING','ASPECTRATIO', 'VISCOSITY', 'ALPHAVISCOSITY', 'SIGMASLOPE','RELEASERADIUS', ...
        'RELEASEDATE', 'OMEGAFRAME', 'DISK', 'FRAME', 'OUTERSOURCEMASS', 'WRITEDENSITY', 'WRITEVELOCITY', ...
        'WRITEENERGY', 'WRITETEMPERATURE', 'WRITEDIVV', 'WRITEQPLUS', 'INDIRECTTERM', 'EXCLUDEHILL', ...
        'IMPOSEDDISKDRIFT', 'FLARINGINDEX', 'ECCENTRICITY', 'CAVITYRADIUS', 'CAVITYRATIO', 'CAVITYWIDTH', ...
        'TRANSITIONRADIUS', 'TRANSITIONRATIO', 'TRANSITIONWIDTH', 'LAMBDADOUBLING', 'SELFGRAVITY', ...
        'CICPLANET', 'FORCEDCIRCULAR', 'ZMPLUS', 'ADIABATIC', 'ADIABATICINDEX', 'COOLING', 'COOLINGTIME0'};

    DT = 1.0;     
    SIGMA0 = 173.0;
    NINTERM = 10.0;
    NTOT = 1501.0;
    OUTPUTDIR = '/~masset';
    INNERBOUNDARY = 'WALL';
    LABELADVECTION = 'NO';
    TRANSPORT = 'FAST';
    PLANETCONFIG = '/Systems/SolarSystem.cfg';
    MASSTAPER = 0.0000001;
    RADIALSPACING = 'ARITHMETIC';
    NRAD = 64.0;
    NSEC = 64.0;
    RMIN = 1.0;
    RMAX = 1.0;
    THICKNESSSMOOTHING = 0.0;
    ROCHESMOOTHING = 0.0;
    ASPECTRATIO = 0.05;
    VISCOSITY = 0.0;
    ALPHAVISCOSITY = 0.0;
    SIGMASLOPE = 0.0;
    RELEASERADIUS = 0.0;
    RELEASEDATE = 0.0;
    OMEGAFRAME = 0.0;
    DISK = 'YES';
    FRAME = 'FIXED';
    OUTERSOURCEMASS = 'NO';
    WRITEDENSITY = 'YES';
    WRITEVELOCITY = 'YES';
    WRITEENERGY = 'NO';  
    WRITETEMPERATURE = 'NO';
    WRITEDIVV = 'NO';
    WRITEQPLUS = 'NO';
    INDIRECTTERM = 'YES';
    EXCLUDEHILL = 'NO';
    IMPOSEDDISKDRIFT = 0.0;
    FLARINGINDEX = 0.0;
    ECCENTRICITY = 0.0;
    CAVITYRADIUS = 0.0;
    CAVITYRATIO = 1.0;
    CAVITYWIDTH = 1.0;
    TRANSITIONRADIUS = 0.0;
    TRANSITIONRATIO = 1.0;
    TRANSITIONWIDTH = 1.0;
    LAMBDADOUBLING = 0.0;
    SELFGRAVITY = 'NO';
    CICPLANET = 'NO';
    FORCEDCIRCULAR = 'NO';
    ZMPLUS = 'NO';
    ADIABATIC = 'NO';
    ADIABATICINDEX = 1.4;
    COOLING = 'NO';
    COOLINGTIME0 = 6.28;
    
    % other variables
    FastTransport = 'YES';
    GuidingCenter = 'NO';
    NonReflecting = 'NO';
    Corotating = 'NO';
    Evanescent = 'NO';
    Write_Temperature = 'NO';
    SelfGravity = 'NO';
    SGZeroMode = 'NO';
    Adiabatic = 'NO';
    size(A)
    for i=1:size(A)
        varname = upper(matlab.lang.makeValidName(A{i,1}));
        eval([varname '= A{i,2};']);
    end
    
    C = setdiff(B,upper(A(1:size(A))));
    fprintf('Secondary variables omitted :\n');
    a = size(C);

    for i=1:a(1,2)
        if (isa(eval(C{1,i}),'double'))
            fprintf('%s ;\t Default value : %g\n', C{i}, eval(C{1,i}));
        else
            fprintf('%s ;\t Default Value : %s\n', C{i}, eval(C{1,i}));
        end
    end
    
  if (strcmp(TRANSPORT,'S'))  % falta definir que caso va acá
      FastTransport = 'NO';
  end
  
  if (strcmp(INNERBOUNDARY,'OPEN'))
      OpenInner = 'YES';
  end
  
  if (strcmp(INNERBOUNDARY,'NONREFLECTING'))
      NonReflecting = 'YES';
  end
  
  if (strcmp(INNERBOUNDARY,'EVANESCENT'))
      Evanescent = 'YES';
  end
  
  if (strcmp(RADIALSPACING,'LOGARITHM'))
      LogGrid = 'YES';
  end
  
  if(strcmp(FRAME,'CENTER'))
      Corotating = 'YES';
  end
  
  if(strcmp(FRAME,'GUIDING-CENTER'))
      Corotating = 'YES';
      GuidingCenter = 'YES';
  end
 
  if(strcmp(SELFGRAVITY,'Z'))
      SelfGravity = 'YES';
      SGZeroMode = 'YES';
  end
  
  if(strcmp(ZMPLUS,'YES') && (~strcmp(SGZeroMode,'YES')))
      fprintf('This is not very meaningfull to involve the anisotropic pressure model ');
      fprintf('(ZMPlus=Yes) without taking into account the axisymmetric component of the ');
      fprintf('disk self-gravity. I decided to put ZMPlus = No. Please check again!\n');
      ZMPlus = 'NO';
  end
  
  if (strcmp(ADIABATIC,'YES'))
      Adiabatic = 'YES';
      Write_Temperature = 'YES';
  end
  
  if ((strcmp(Adiabatic,'YES')) && (ADIABATICINDEX == 1)) 
      fprintf('You cannot have Adiabatic = YES and AdiabatcIndex = 1. I decided to put ');
      fprintf('Adiabatic = No, to simulate a locally isothermal equation of state. Please ');
      fprintf('check that it what you really wanted to do!\n');
      Adiabatic = 'NO';
  end
      
  if ((ALPHAVISCOSITY ~= 0.0) && (VISCOSITY ~= 0.0)) 
      fprintf('You cannot use at the same time\n');
      fprintf('VISCOSITY and ALPHAVISCOSITY.\n');
      fprintf('Edit the parameter file so as to remove\n');
      fprintf('one of these variables and run again.\n');
      quit cancel;
  end
  
  if (ALPHAVISCOSITY ~= 0.0) 
      fprintf('Viscosity is of alpha type\n');
  end
  
  if ((THICKNESSSMOOTHING ~= 0.0) && (ROCHESMOOTHING ~= 0.0)) 
      fprintf('You cannot use at the same time\n');
      fprintf('ThicknessSmoothing and RocheSmoothing.\n');
      fprintf('Edit the parameter file so as to remove\n');
      fprintf('one of these variables and run again.\n');
      quit cancel;
  end
  
  if ((THICKNESSSMOOTHING <= 0.0) && (ROCHESMOOTHING <= 0.0))
      fprintf('A non-vanishing potential smoothing length is required.\n');
      fprintf('Please use either of the following variables:\n');
      fprintf ('ThicknessSmoothing *or* RocheSmoothing.\n');
      fprintf('before launching the run again.\n');
      quit cancel;
  end
  
  if (ROCHESMOOTHING ~= 0.0)
      fprintf('Planet potencial smoothing scales with their Hill sphere.\n');
  end
  
  if (~strcmp(OUTPUTDIR,'/~masset'))
      fprintf('new output %s',OUTPUTDIR);
  end
  
  %if (*(OUTPUTDIR+strlen(OUTPUTDIR)-1) != '/')
  %  strcat (OUTPUTDIR, "/");   
end

function TellEverything()
    global RMIN RMAX ASPECTRATIO G SIGMA0 SIGMASLOPE NRAD NSEC NINTERM DT;
    global Adiabatic LABELADVECTION NTOT;
    
    total_size = NRAD*NSEC*8;
    fprintf('\nDisc properties:\n');
    fprintf('----------------\n');
    fprintf('Inner Radius          : %g\n', RMIN);
    fprintf('Outer Radius          : %g\n', RMAX);
    fprintf('G                     : %g\n', G)
    fprintf('Aspect Ratio          : %g\n', ASPECTRATIO)
    fprintf('VKep at inner edge    : %.3g\n', sqrt(G*1.0*(1.-0.0)/RMIN))
    fprintf('VKep at outer edge    : %.3g\n', sqrt(G*1.0/RMAX))
    temp=2.0*pi*SIGMA0/(2.0-SIGMASLOPE)*((RMAX^(2.0-SIGMASLOPE)) -RMIN^(2.0-SIGMASLOPE));
    %correct this and what follows...
    fprintf('Initial Disk Mass             : %g\n', temp)
    temp=2.0*pi*SIGMA0/(2.0-SIGMASLOPE)*(1.0 -RMIN^(2.0-SIGMASLOPE));
    fprintf('Initial Mass inner to r=1.0  : %g \n', temp)
    temp=2.0*pi*SIGMA0/(2.0-SIGMASLOPE)*(RMAX^(2.0-SIGMASLOPE)- 1.0);
    fprintf('Initial Mass outer to r=1.0  : %g \n', temp)
    fprintf('Travelling time for acoustic density waves :\n')
    temp = 2.0/3.0/ASPECTRATIO*((RMAX^1.5)-(RMIN^1.5));
    fprintf('* From Rmin to Rmax  : %.2g = %.2f orbits ~ %.1f outputs\n', temp, TellNbOrbits(temp), TellNbOutputs(temp))
    temp = 2.0/3.0/ASPECTRATIO*((RMAX^1.5)-(1.0^1.5));
    fprintf(' * From r=1.0 to Rmax: %.2g = %.2f orbits ~ %.1f outputs\n', temp, TellNbOrbits(temp), TellNbOutputs(temp))
    temp = 2.0/3.0/ASPECTRATIO*((1.0^1.5)-(RMIN^1.5));
    fprintf(' * From r=1.0 to Rmin: %.2g = %.2f orbits ~ %.1f outputs\n', temp, TellNbOrbits(temp), TellNbOutputs(temp))
    temp = 2.0*pi*sqrt(RMIN*RMIN*RMIN/G/1.0);
    fprintf('Orbital time at Rmin  : %.3g ~ %.2f outputs\n', temp, TellNbOutputs(temp))
    temp = 2.0*pi*sqrt(RMAX*RMAX*RMAX/G/1.0);
    fprintf('Orbital time at Rmax  : %.3g ~ %.2f outputs\n', temp, TellNbOutputs(temp))
    fprintf('Sound speed :\n')
    fprintf(' * At unit radius     : %.3g\n', ASPECTRATIO*sqrt(G*1.0))
    fprintf(' * At outer edge      : %.3g\n', ASPECTRATIO*sqrt(G*1.0/RMAX))
    fprintf(' * At inner edge      : %.3g\n', ASPECTRATIO*sqrt(G*1.0/RMIN))
    fprintf('\nGrid properties:\n')
    fprintf('----------------\n')
    fprintf('Number of rings       : %d\n', NRAD)
    fprintf('Number of sectors     : %d\n', NSEC)
    fprintf('Total cells           : %d\n', NRAD*NSEC)
    fprintf('\nOutputs properties:\n')
    fprintf('-------------------\n')
    fprintf('Time increment between outputs : %.3f = %.3f orbits\n', NINTERM*DT, TellNbOrbits(NINTERM*DT))
    fprintf('At each output #i, the following files are written:\n')
    fprintf('gasdens[i].dat : %d bytes\n',total_size) % sizeof(double)
    fprintf('gasvrad[i].dat : %d bytes\n',total_size)
    fprintf('gasvtheta[i].dat : %d bytes\n',total_size)
    if (strcmp(Adiabatic,'YES'))
        fprintf('gasTemperature[i].dat : %d bytes\n',total_size)
    end
    
    if (strcmp(LABELADVECTION, 'YES'))
        fprintf('gaslabel[i].dat : %d bytes\n',total_size)
    end
    fprintf('There will be in total %d outputs\n', int16(NTOT/NINTERM))
    fprintf('(which correspond to an elapsed time = %.3f or to %.2f orbits)\n', NTOT*DT, TellNbOrbits(NTOT*DT))
    nbfileoutput = 3.0;
    if (strcmp(Adiabatic,'YES'))
        nbfileoutput = nbfileoutput + 1.0;
    end
    
    if (strcmp(LABELADVECTION,'YES'))
        nbfileoutput = nbfileoutput + 1.0;
    end
    
    temp =nbfileoutput*total_size;
    temp = temp * NTOT/NINTERM;
    temp = temp /(1024.0*1024.0);
    fprintf('So the code will produce ~%.2f Mbytes of data\n', temp)
    fprintf('Check (eg by issuing a df command) that you have enough disk space,\n')
    fprintf('otherwise you will get a system full and the code will stop.\n')
    
end

function z = TellNbOrbits(temp)
    global G;
    z = temp/2.0/pi*sqrt(G*1.0/1.0/1.0);
end

function z = TellNbOutputs(temp)
    global DT NINTERM;
    z = temp/DT/NINTERM;
end

