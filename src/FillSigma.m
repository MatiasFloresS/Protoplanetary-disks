function y = FillSigma()
    global Rmed Rinf SigmaMed NRAD;
    size(Rmed);
    SigmaMed = Sigma(Rmed(1:NRAD));
    SigmaInf = Sigma(Rinf(1:NRAD));

end 

function z = Sigma(r)
    global SIGMASLOPE SIGMA0;
    cavity = 1.0;
    CAVITYRADIUS= 0.0;
    ScalingFactor = 1;
    CAVITYRATIO = 1.0;

    if (r< CAVITYRADIUS)
        cavity = 1.0/CAVITIRATIO;
    end

    z = cavity*ScalingFactor*SIGMA0*r.^(-SIGMASLOPE);

    return 

end