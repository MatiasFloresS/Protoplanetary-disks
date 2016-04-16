function x = FillEnergy()
    global NRAD Rmed EnergyMed;
    EnergyMed = Energy(Rmed(1:NRAD));

end


function w = Energy(r)
    global SIGMASLOPE ASPECTRATIO R MU SIGMA0 FLARINGINDEX;
    ADIABATICINDEX = 1.4;
    double energy0;

    if (ADIABATICINDEX == 1.0)
        sprintf('The adiabatic index must differ from unity to initialize the gas internal energy. I must exit.')
        quit cancel
    
    else
        w = R/MU/(ADIABATICINDEX-1.0)*SIGMA0*(ASPECTRATIO.^2.0)*r.^(-SIGMASLOPE-1.0+2.0*FLARINGINDEX); %Energy0
    end
    


end