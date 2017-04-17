function [gas_density, gas_v_theta] = transport(gas_v_rad, gas_v_theta, gas_density, dt)
    ComputeLRMomenta(gas_v_rad, gas_v_theta, gas_density);
    [gas_density] = OneWindRad(gas_density, gas_v_rad, dt);
    [gas_v_theta] = OneWindTheta(gas_density, gas_v_theta, dt);
    
end

%lista
function ComputeLRMomenta(gas_v_rad, gas_v_theta, gas_density)
    global RadMomP RadMomM ThetaMomP ThetaMomM NRAD OmegaFrame Rmed NSEC;
    
    RadMomP(1:NRAD, :) = gas_density(1:NRAD, :).*gas_v_rad(2:NRAD+1,:);
    RadMomM(1:NRAD, :) = gas_density(1:NRAD, :).*gas_v_rad(1:NRAD,:);
    
    thetatemp = bsxfun(@plus, gas_v_theta(1:NRAD, mod(1:NSEC, NSEC)+1), Rmed(1:NRAD)'.*OmegaFrame);
    thetatemp2 = bsxfun(@plus, gas_v_theta(1:NRAD, :), Rmed(1:NRAD)'.*OmegaFrame);
    
    ThetaMomP(1:NRAD,:) = bsxfun(@times, gas_density(1:NRAD, :).*thetatemp, Rmed(1:NRAD)');
    ThetaMomM(1:NRAD,:) = bsxfun(@times, gas_density(1:NRAD, :).*thetatemp2, Rmed(1:NRAD)');
end

function [gas_density] = OneWindRad(gas_density, gas_v_rad, dt)
    global densInt;
    
    ComputeStarRad(gas_v_rad, dt, 1, gas_density);
    densInt = gas_density;
    
    VanLeerRadial(gas_v_rad, dt, 1, []); %RadMomP
    VanLeerRadial(gas_v_rad, dt, 2, []); %RadMomM
    VanLeerRadial(gas_v_rad, dt, 3, []); %ThetaMomP
    VanLeerRadial(gas_v_rad, dt, 4, []); %ThetaMomM
    [gas_density] = VanLeerRadial(gas_v_rad, dt, 5, gas_density); %gas_density
end

function [gas_density] = VanLeerRadial(Vrad, dt, option, gas_density)
    global RadMomP NSEC densStar QRStar Rsup Rinf NRAD InvSurf RadMomM ThetaMomP ThetaMomM;
    DivisePolarGrid(option, gas_density);
    ComputeStarRad(Vrad,dt,2, []);
    
    dtheta = 2.0*pi/NSEC;
    
    varq = bsxfun(@times, dt*dtheta.*Rinf(1:NRAD)',QRStar(1:NRAD,:).*densStar(1:NRAD,:).*Vrad(1:NRAD,:)) - ...
               bsxfun(@times, dt*dtheta.*Rsup(1:NRAD)',QRStar(2:NRAD+1,:).*densStar(2:NRAD+1,:).*Vrad(2:NRAD+1,:));
           
    switch option
        case 1
            RadMomP(1:NRAD,:) = RadMomP(1:NRAD,:)+ bsxfun(@times, varq, InvSurf(1:NRAD)');
        case 2
            RadMomM(1:NRAD, :) = RadMomM(1:NRAD,:)+bsxfun(@times, varq, InvSurf(1:NRAD)');
        case 3
            ThetaMomP(1:NRAD, :) = ThetaMomP(1:NRAD,:)+bsxfun(@times, varq, InvSurf(1:NRAD)');
        case 4
            ThetaMomM(1:NRAD, :) = ThetaMomM(1:NRAD,:)+bsxfun(@times, varq, InvSurf(1:NRAD)');
        case 5
            gas_density(1:NRAD, :) = gas_density(1:NRAD,:)+bsxfun(@times, varq, InvSurf(1:NRAD)');
    end
    
    
end

function [gas_density] = DivisePolarGrid(option, gas_density)
    global RadMomP densInt Work RadMomM ThetaMomP ThetaMomM;
    
    switch option
        case 1
            Work = RadMomP./(densInt+1e-20);
        case 2
            Work = RadMomM./(densInt+1e-20);
        case 3
            Work = ThetaMomP./(densInt+1e-20);
        case 4
            Work = ThetaMomM./(densInt+1e-20);
        case 5
            Work = gas_density./(densInt+1e-20);
    end
end

function [gas_v_theta] = OneWindTheta(gas_density, gas_v_theta, dt)
    ComputeThetaElongations(gas_v_theta, dt);
    ComputeAverageThetaVelocities(dt);
    ComputeResiduals(dt);
    [gas_v_theta] = ComputeConstantResidual(gas_v_theta, dt);
    UniformTransport = 0;
    QuantitiesAdvection (gas_density, gas_v_theta, dt);
end

function ComputeThetaElongations(gas_v_theta, dt)
    global elong NRAD;
    elong = gas_v_theta(1:NRAD,:)*dt;
end

function ComputeAverageThetaVelocities(dt)
    global elong NRAD NSEC Vmed;
    
    invdt = 1.0/dt;
    moy = sum(elong(1:NRAD,:)*invdt, 2);
    
    Vmed(1:NRAD) = moy/NSEC;
end

function ComputeResiduals(dt)
    global VthetaRes NRAD Vmed elong;
    
    invdt = 1.0/dt;
    VthetaRes = bsxfun(@minus, elong(1:NRAD,:)*invdt,Vmed(1:NRAD)');
end

function [gas_v_theta] = ComputeConstantResidual(gas_v_theta, dt)
    global NRAD NSEC Vmed InvRmed Rmed NoSplitAdvection;
    
    invdt = 1.0/dt;
    dpinvns = 2.0*pi/NSEC;
    
    %maxfrac = 1.0;
    
    Ntilde = Vmed(1:NRAD).*InvRmed(1:NRAD)*dt*NSEC/2.0/pi;
    Nshift = floor(Ntilde+0.5);

    gas_v_theta(1:NRAD,:) = repmat((Ntilde(1:NRAD)-Nshift(1:NRAD)).*Rmed(1:NRAD)*invdt*dpinvns,NSEC,1)';
    NoSplitAdvection(1:NRAD) = 0;
    
end

function QuantitiesAdvection (gas_density, Vtheta, dt)
    global densInt;
    ComputeStarTheta(Vtheta, dt, 1, gas_density);
    densInt = gas_density;
    VanLeerTheta(Vtheta, dt, 1, []); %RadMomP
end

function VanLeerTheta(Vtheta, dt, option, gas_density)
    global NRAD Surf Rinf Rsup QRStar densStar;
    DivisePolarGrid(option, gas_density);
    ComputeStarTheta(Vtheta, dt, 2, []);
    
    dxrad = Rsup(1:NRAD)-Rinf(1:NRAD)*dt;
    invsurf = 1.0./Surf(1:NRAD);
    
    varq = bsxfun(@times, dxrad(1:NRAD)', QRStar(1:NRAD,:).*densStar(1:NRAD,:).*Vtheta(1:NRAD,:));
    varq = varq - dxrad(1:NRAD)'
    %{
    varq = bsxfun(@times, dt*dtheta.*Rinf(1:NRAD)',QRStar(1:NRAD,:).*densStar(1:NRAD,:).*Vrad(1:NRAD,:)) - ...
               bsxfun(@times, dt*dtheta.*Rsup(1:NRAD)',QRStar(2:NRAD+1,:).*densStar(2:NRAD+1,:).*Vrad(2:NRAD+1,:));
           
    switch option
        case 1
            RadMomP(1:NRAD,:) = RadMomP(1:NRAD,:)+ bsxfun(@times, varq, InvSurf(1:NRAD)');
        case 2
            RadMomM(1:NRAD, :) = RadMomM(1:NRAD,:)+bsxfun(@times, varq, InvSurf(1:NRAD)');
        case 3
            ThetaMomP(1:NRAD, :) = ThetaMomP(1:NRAD,:)+bsxfun(@times, varq, InvSurf(1:NRAD)');
        case 4
            ThetaMomM(1:NRAD, :) = ThetaMomM(1:NRAD,:)+bsxfun(@times, varq, InvSurf(1:NRAD)');
        case 5
            gas_density(1:NRAD, :) = gas_density(1:NRAD,:)+bsxfun(@times, varq, InvSurf(1:NRAD)');
    end
    %}
end

function ComputeStarTheta(Vtheta, dt, option, Qbase)
    global NSEC Rmed NRAD dqm2 dqp2 VthetaRes densStar QRStar Work;
    
    dxtheta = 2.0*pi/NSEC*Rmed(1:NRAD);
    invdxtheta = 1.0./dxtheta;
    
    if(option == 1)
        dqm2(1:NRAD, :) = (Qbase(1:NRAD, :) - Qbase(1:NRAD, mod((0:NSEC-1)+NSEC-1,NSEC)+1));
        dqp2(1:NRAD, :) = (Qbase(1:NRAD, mod(1:NSEC, NSEC)+1) - Qbase(1:NRAD, :));
    else
        dqm2(1:NRAD, :) = (Work(1:NRAD, :) - Work(1:NRAD, mod((0:NSEC-1)+NSEC-1,NSEC)+1));
        dqp2(1:NRAD, :) = (Work(1:NRAD, mod(1:NSEC, NSEC)+1) - Work(1:NRAD, :));
    end
    
    dqtemp = dqm2.*dqp2;
    
    
    I = find(dqtemp < 0.0);
    temp = dqtemp(:);
    temp(I) = 0.0;
    dqtemp = reshape(temp, NRAD, NSEC);
    
    I = find(dqtemp > 0.0);
    temp = dqtemp(:);
    dqptemp = dqp2(:);
    dqmtemp = dqm2(:);
    temp(I) = dqptemp(I).*dqmtemp(I)./(dqptemp(I)+dqmtemp(I));   
    dq = reshape(temp, NRAD, NSEC);
    dq(1:NRAD, :) = bsxfun(@times, dq(1:NRAD, :), invdxtheta(1:NRAD)');
    
    if(option == 1)
        ksi = VthetaRes(1:NRAD,:)*dt;
    else
        ksi = Vtheta(1:NRAD,:)*dt;
    end
    
    dxtemp = bsxfun(@plus, dxtheta(1:NRAD)', ksi(1:NRAD,:));
    dxtempm = bsxfun(@minus, dxtheta(1:NRAD)', ksi(1:NRAD, :));
    
    I = find(ksi <= 0.0);
    temp = ksi(:);
    if (option == 1)
        qbtemp = Qbase(:);
    else
        qbtemp = Work(:);
    end
    dqtemp = dq(:);
    dxtemp2 = dxtemp(:);
    temp2(I) = qbtemp(I) - (dxtemp2(I)).*dqtemp(I);
    
    I = find(ksi > 0);
    temp = ksi(:);
    if (option == 1)
        qbtemp = Qbase(:);
    else
        qbtemp = Work(:);
    end
    dqtemp = dq(:);
    dxtempm2 = dxtempm(:);
    temp2(I) = qbtemp(I) + (dxtempm2(I)).*dqtemp(mod((I-1)+NSEC-1,NSEC)+1);
    
    if(option == 1)
        densStar = reshape(temp2, NRAD, NSEC);
        densStar(1,:) = 0.0;
        densStar(NRAD+1,:) = 0.0;
    else
        QRStar = reshape(temp2, NRAD, NSEC);
        QRStar(1,:) = 0.0;
        QRStar(NRAD+1, :) = 0.0;
    end
    
end

function ComputeStarRad(Vrad, dt, option, Qbase)
    global densStar NRAD dq InvDiffRmed NSEC dqm dqp Rmed QRStar Work temp2;
    
    if(option == 1)
        dqm(2:NRAD-1, :) = bsxfun(@times, Qbase(2:NRAD-1,:) -Qbase(1:NRAD-2,:), InvDiffRmed(1:NRAD-2)');
        dqp(2:NRAD-1, :) = bsxfun(@times, Qbase(3:NRAD,:)-Qbase(2:NRAD-1,:),InvDiffRmed(2:NRAD-1)');
    else
        dqm(2:NRAD-1, :) = bsxfun(@times, Work(2:NRAD-1,:) -Work(1:NRAD-2,:), InvDiffRmed(1:NRAD-2)');
        dqp(2:NRAD-1, :) = bsxfun(@times, Work(3:NRAD,:)-Work(2:NRAD-1,:),InvDiffRmed(2:NRAD-1)');
    end
    dqtemp = dqm.*dqp;

    I = find(dqtemp < 0.0);
    temp = dqtemp(:);
    temp(I) = 0.0;
    dqtemp = reshape(temp, NRAD, NSEC);
    
    I = find(dqtemp > 0.0);
    temp = dqtemp(:);
    dqptemp = dqp(:);
    dqmtemp = dqm(:);
    temp(I) = 2.0.*dqptemp(I).*dqmtemp(I)./(dqptemp(I)+dqmtemp(I));   
    dq = reshape(temp, NRAD, NSEC);
 
    if(option == 1)
        vradtemp(2:NRAD,:) = Qbase(1:NRAD-1, :) + bsxfun(@minus, Rmed(2:NRAD)'-Rmed(1:NRAD-1)', Vrad(2:NRAD,:)*dt).*0.5.*dq(1:NRAD-1,:);   
        rmedtemp2 = bsxfun(@plus, Rmed(2:NRAD)'-Rmed(1:NRAD-1)', Vrad(1:NRAD-1,:)*dt);
        rmedtemp2(NRAD,:) = -Rmed(NRAD)'*Vrad(NRAD,:)*dt;
        vradtemp2(1:NRAD,:) = Qbase(1:NRAD,:) - rmedtemp2.*0.5.*dq(1:NRAD,:);
    else
        vradtemp(2:NRAD,:) = Work(1:NRAD-1, :) + bsxfun(@minus, Rmed(2:NRAD)'-Rmed(1:NRAD-1)', Vrad(2:NRAD,:)*dt).*0.5.*dq(1:NRAD-1,:);   
        rmedtemp2 = bsxfun(@plus, Rmed(2:NRAD)'-Rmed(1:NRAD-1)', Vrad(1:NRAD-1,:)*dt);
        rmedtemp2(NRAD,:) = -Rmed(NRAD)'*Vrad(NRAD,:)*dt;
        vradtemp2(1:NRAD,:) = Work(1:NRAD,:) - rmedtemp2.*0.5.*dq(1:NRAD,:);
    end
    
    I = find(Vrad(1:NRAD, :) > 0.0);
    
    
    if (~isempty(I))
        denstemp = vradtemp(:);
        temp2(I) = denstemp(I);
        temp2 = reshape(temp2, NRAD, NSEC);
    end
    
    I = find(Vrad(1:NRAD, :) <= 0.0);
    denstemp2 = vradtemp2(:);
    temp2(I) = denstemp2(I);
    
    if(option == 1)
        densStar = reshape(temp2, NRAD, NSEC);
        densStar(1,:) = 0.0;
        densStar(NRAD+1,:) = 0.0;
    else
        QRStar = reshape(temp2, NRAD, NSEC);
        QRStar(1,:) = 0.0;
        QRStar(NRAD+1, :) = 0.0;
    end
end

