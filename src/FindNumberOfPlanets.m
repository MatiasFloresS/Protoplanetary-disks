function cant_planets = FindNumberOfPlanets(filename)
    sprintf('%s',filename);
    num = 1;
    cant_planets = 0;
    fid = fopen(filename, 'r');
    if fid == -1
        sprintf('Cannot open file: %s', filename)
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
        sprintf('%d\n',cant_planets);
        return
    end
end

