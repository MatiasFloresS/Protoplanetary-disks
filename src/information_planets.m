function information_planets(filename,cant_planets)
    global GlobalRmed Radii system G;
    ECCENTRICITY = 0.0;
    G = 1;
    sprintf('%s',filename);
    CICPlanet = false;
    num=0;
    num_planet=1;
    fid = fopen(filename, 'r');
    if fid == -1
        sprintf('Cannot open file: %s', filename)
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
                
                dist = fscanf(fid, '%s', 1);
                distf = str2double(dist);
                if (CICPlanet == true) %initialization puts centered-in-cell planets (with excentricity = 0 only)
                    j = 0;
                    while ( GlobalRmed(j+1) < distf ) 
                        j=j+1;
                    end
                    distf = Radii(j+2);    
                end
                mass = fscanf(fid, '%s', 1);
                massf = str2double(mass);
                system{2,1}(num_planet) = massf;
                accret = fscanf(fid, '%s', 1);
                accretf = str2double(accret);
                system{7,1}(num_planet) = accretf;
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
                
                system{3,1}(num_planet) = real(distf)*(1.0+ECCENTRICITY);
                system{4,1}(num_planet) = 0.0;
                system{6,1}(num_planet) = real(sqrt(G*(1.0+massf)/distf))* sqrt((1.0-ECCENTRICITY)/(1.0+ECCENTRICITY));
                system{5,1}(num_planet) = system{6,1}(num_planet)* -0.0000000001;
                
                
                cant_planets = cant_planets - 1;
                num_planet = num_planet + 1;
            end
            tline = fgetl(fid);    
        end
        fclose(fid);
        
        hillradius = system{3,1}(1) * (system{2,1}(1)/3)^(1./3.);
        return
    end
end