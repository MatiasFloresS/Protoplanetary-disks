function FillSigma()
    global Rmed Rinf SigmaMed SigmaInf;
    SigmaMed= Sigma(Rmed);    
    SigmaInf = Sigma(Rinf);
end 

function z = Sigma(r)
    global SIGMASLOPE SIGMA0;
    cavity = 1.0;
    CAVITYRADIUS= 0.0;
    ScalingFactor = 1;
    CAVITYRATIO = 1.0;

    if (r< CAVITYRADIUS)
        cavity = 1.0/CAVITYRATIO;
    end

    z = cavity*ScalingFactor*SIGMA0*r.^(-SIGMASLOPE);

end