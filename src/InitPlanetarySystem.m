function cant_planets = InitPlanetarySystem()
    global system ;
    cant_planets = FindNumberOfPlanets();
    fprintf ('%d planet(s) found.\n', cant_planets);
    AllocPlanetsSystem(cant_planets);
    system{1} = cant_planets;
    information_planets(cant_planets);
    return
end

function cant_planets = FindNumberOfPlanets()
    global PLANETCONFIG;
    %fprintf('%s',PLANETCONFIG)
    num = 1;
    cant_planets = 0;
    filename = PLANETCONFIG;
    fid = fopen(filename, 'r');
    if fid == -1
        fprintf('Cannot open file: %s', PLANETCONFIG)
        quit cancel
    else
        tline = fgetl(fid);
        while ischar(tline)
            if num > 1
                cant_planets = cant_planets  + num -1;
            end
            %disp(tline)
            matches = strfind(tline, 'Others');
            if (num == length(matches))
                num=num+1;
            end
            tline = fgetl(fid);
        end
        fclose(fid);
        %fprintf('%d\n',cant_planets);
        return
    end
end

function AllocPlanetsSystem (cant_planets)
    global system ;
    system = cell(10,1);
    system{10,1} = [];

    % {
    %  int nb;
    %  real *mass;
    %  real *x;
    %  real *y;
    %  real *vx;
    %  real *vy;
    %  real *acc;
    %  char **name; no se ocupa, asi que se le da espacio de un array de 1x3
    %  boolean *FeelDisk;
    %  boolean *FeelOthers;
    %  }

    cero = zeros(1,cant_planets);
    cero2 = zeros(1);
    for i =1:10
        if (i == 1)
            system{i,1} = cero2;
        else
            system{i,1} = cero;
        end
    end

    for i = 1:10
        for j = 1:cant_planets
            if (i == 9 || i == 10)
                system{i,1}(j) = true;
            elseif (i == 1)
                % do nothing
            else
                system{i,1}(j) = 0.0;
            end
        end
    end   
end

function information_planets(cant_planets)
    global GlobalRmed Radii system G PLANETCONFIG;
    ECCENTRICITY = 0.0;
    G = 1;
    CICPlanet = false;
    num=0;
    num_planet=1;
    fid = fopen(PLANETCONFIG, 'r');
    if fid == -1
        fprintf('Cannot open file: %s', PLANETCONFIG)
        quit cancel
    else
        tline = fgetl(fid);
        while ischar(tline)
            matches = strfind(tline, 'Others');
            if (length(matches)==1)
                num = num + 1;
            end
            if (num == 1 && cant_planets ~=0)
                fscanf(fid, '%s', 1); % nombre de los planetas
                
                dist = fscanf(fid, '%d', 1);
                %distf = str2double(dist);
                if (CICPlanet == true) %initialization puts centered-in-cell planets (with excentricity = 0 only)
                    j = 0;
                    while ( GlobalRmed(j+1) < dist ) 
                        j=j+1;
                    end
                    dist = Radii(j+2);    
                end
                mass = fscanf(fid, '%f', 1);
                system{2,1}(num_planet) = mass;
                accret = fscanf(fid, '%f', 1);
                system{7,1}(num_planet) = accret;
                feeldis = fscanf(fid, '%s', 1);
                if (strcmp('NO',feeldis))
                    system{9,1}(num_planet) = false;
                else
                    system{9,1}(num_planet) = true;
                end
                
                feelothers = fscanf(fid, '%s', 1);
                
                if (strcmp('NO', feelothers))
                    system{10,1}(num_planet) = false;
                else
                    system{10,1}(num_planet) = true;
                end
                
                system{3,1}(num_planet) = real(dist)*(1.0+ECCENTRICITY);
                system{4,1}(num_planet) = 0.0;
                system{6,1}(num_planet) = real(sqrt(G*(1.0+mass)/dist))* sqrt((1.0-ECCENTRICITY)/(1.0+ECCENTRICITY));
                system{5,1}(num_planet) = system{6,1}(num_planet)* -0.0000000001;
                
                cant_planets = cant_planets - 1;
                num_planet = num_planet + 1;
            end
            tline = fgetl(fid);    
        end
        fclose(fid);
        
        %hillradius = system{3,1}(1) * (system{2,1}(1)/3)^(1./3.);
        return
    end
end