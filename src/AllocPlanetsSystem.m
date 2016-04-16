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