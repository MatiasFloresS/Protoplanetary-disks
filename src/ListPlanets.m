function ListPlanets (system)
    for i=1:system{1}
        fprintf('Planet number %d\n', i)
        fprintf('---------------\n')
        fprintf('x = %f\ty=%f\n',system{3,1}(i),system{4,1}(i))
        fprintf('vx = %f\tvy=%f\n',system{5,1}(i),system{6,1}(i))
        %[system{5,1}(i) system{6,1}(i)]
        if (system{7,1}(i) == 0.0)
            fprintf('Non-accreting\n')
        else
            fprintf('accretion time = %.20d\n', 1.0/system{7,1}(i))
        end
        
        if (system{9,1}(i) == 1)
            fprintf('Feels the disk potential\n')
        else
            fprintf('Does not feel the disk potential\n')
        end
        
        if (system{10,1}(i) == 1)
            fprintf('Feels the other planets potential\n')
        else
            fprintf('Does not feel the other planets potential \n')
        end
        fprintf('\n')
    end
end 
   