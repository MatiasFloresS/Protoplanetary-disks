function ListPlanets (system)
    
    for i=1:system{1}
        sprintf('Planet number %d\n', i)
        sprintf('---------------\n')
        sprintf('x = %f\ty=%f\n',system{3,1}(i),system{4,1}(i))
        sprintf('vx = %f\tvy=%f\n',system{5,1}(i),system{6,1}(i))
        %[system{5,1}(i) system{6,1}(i)]
        if (system{7,1}(i) == 0.0)
            sprintf('Non-accreting\n')
        
        else
            sprintf('accretion tiem = %.20d\n', 1.0/system{7,1}(i))
        end
        
        if (system{9,1}(i) == 1)
            sprintf('Feels the disk potential\n')
        else
            sprintf('Does not feel the disk potential\n')
        end
        
        if (system{10,1}(i) == 1)
            sprintf('Feels the other planets potential\n')
        else
            sprintf('Does not feel the other planets potential \n')
        end
        
        sprintf('\n')
    end
end 
   