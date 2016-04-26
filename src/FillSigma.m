function FillSigma()
    global Rmed Rinf SigmaMed SigmaInf;
    SigmaMed= Sigma(Rmed);    
    SigmaInf = Sigma(Rinf);
end 

function z = Sigma(r)
    global SIGMASLOPE SIGMA0 CAVITYRADIUS CAVITYRATIO;
    cavity = 1.0;
    ScalingFactor = 1;

    if (r< CAVITYRADIUS)
        cavity = 1.0/CAVITYRATIO;
    end
    z = cavity*ScalingFactor*SIGMA0*r.^(-SIGMASLOPE);
end