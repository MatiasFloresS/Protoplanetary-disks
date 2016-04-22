function main(ParameterFile)
    global force system  SigmaMed R MU  NRAD NSEC EnergyMed GuidingCenter;
    global PhysicalTime PhysicalTimeInitial RESTART ADIABATIC SELFGRAVITY OMEGAFRAME;
    
    ReadVariables(ParameterFile);
    
    RESTART = 0;
    PhysicalTimeInitial = 0.0;
    PhysicalTime =0.0;
    FREQUENCY = 2;
    Corotating = 1;
    GuidingCenter = 1; % YES
    
    R = 1.0; % Mean molecular weight
    MU = 1.0; % Universal Gas Constant in code units 
    
    gas_density     = zeros([NRAD,NSEC]);
    gas_v_rad       = zeros([NRAD,NSEC]);
    gas_v_theta     = zeros([NRAD,NSEC]);
    gas_energy      = zeros([NRAD,NSEC]);
    gas_label       = zeros([NRAD,NSEC]); 
    
    dimfxy = 11;
    FillPolar1DArrays();
    AllocateForce(dimfxy);

    planets = InitPlanetarySystem();
    sprintf('%d\n',planets);
    
    % InitGasDensity
    FillSigma();
    for j=1:NSEC
        gas_density(:,j) = SigmaMed;
    end
    
    if (ADIABATIC)
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

function ReadVariables(ParameterFile)
    sprintf('%s',ParameterFile);
    c = 1;
    fid = fopen(ParameterFile, 'r');
    if fid == -1
        sprintf('Cannot open file: %s', ParameterFile)
        quit cancel
    else
        tline = fgetl(fid);
        while ischar(tline)
            %disp(tline);
            matches = strfind(tline, '###');
            if( matches ~= 0)
                %do nothing
            
            else
                %disp(tline);
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
    return
end


function InitVariables(A)
    global DT SIGMA0 NINTERM NTOT OUTPUTDIR INNERBOUNDARY LABELADVECTION;
    global TRANSPORT PLANETCONFIG MASSTAPER RADIALSPACING NRAD NSEC RMIN RMAX;
    global THICKNESSSMOTHING ROCHESMOOTHING ASPECTRATIO VISCOSITY ALPHAVISCOSITY;
    global SIGMASLOPE RELEASERADIUS RELEASEDATE OMEGAFRAME DISK FRAME OUTERSOURCEMASS;
    global WRITEDENSITY WRITEVELOCITY WRITEENERGY WRITETEMPERATURE WRITEDIVV WRITEQPLUS;
    global INDIRECTTERM EXCLUDEHILL IMPOSEDDISKDRIFT FLARINGINDEX ECCENTRICITY CAVITYRADIUS;
    global CAVITYRATIO CAVITYWIDTH TRANSITIONRADIUS TRANSITIONRATIO TRANSITIONWIDTH;
    global LAMBDADOUBLING SELFGRAVITY CICPLANET FORCEDCIRCULAR ZMPLUS ADIABATIC;
    global ADIABATICINDEX COOLING COOLINGTIME0;
    
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
    THICKNESSSMOTHING = 0.0;
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
    
    for i=1:size(A)
        varname = upper(genvarname(A{i,1}));
        eval([varname '= A{i,2};']);
    end
    
end