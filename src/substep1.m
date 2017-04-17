function substep1(gas_v_rad, gas_v_theta, gas_density, dt)
    global Pressure NRAD InvDiffRmed NSEC Potential OmegaFrame Rinf InvRinf Rmed;
    global vradint vthetaint gradp gradphi gradp2 gradphi2 vt2;
    
    gradp(2:NRAD,:) = (Pressure(2:NRAD,:) - Pressure(1:NRAD-1,:)).*2./(gas_density(2:NRAD,:) + gas_density(1:NRAD-1,:));
    gradp = bsxfun(@times, gradp(2:NRAD,:), InvDiffRmed(1:NRAD-1)');
    
    gradphi(2:NRAD,:) = (Potential(2:NRAD,:) - Potential(1:NRAD-1,:));
    gradphi = bsxfun(@times, gradphi(2:NRAD,:),InvDiffRmed(1:NRAD-1)');
    
    vt2(2:NRAD, :) = gas_v_theta(2:NRAD, :) + gas_v_theta(2:NRAD, mod(1:NSEC, NSEC)+1) + gas_v_theta(1:NRAD-1, :)+ gas_v_theta(1:NRAD-1, mod(1:NSEC,NSEC)+1);
    vt2(2:NRAD, :) = vt2(2:NRAD,:)/4.0;
    
    vt2 = bsxfun(@plus, vt2(2:NRAD, :), Rinf(2:NRAD)'*OmegaFrame);
    vt2 = vt2.*vt2;
    gradp(2:NRAD, :) = gradp(1:NRAD-1,:);
    gradphi(2:NRAD, :) = gradphi(1:NRAD-1,:);
    vt2(2:NRAD, :) = vt2(1:NRAD-1, :);
    
    vt2(1,:) = 0.0;
    gradp(1,:) = 0.0;
    gradphi(1,:) = 0.0;
    dxtheta = 2.0*pi/double(NSEC)*Rmed(1:NRAD);
    invdxtheta = 1.0./dxtheta;

    vradtemp = bsxfun(@times, vt2(2:NRAD,:), InvRinf(2:NRAD)');
    vradtemp(2:NRAD,:) = vradtemp(1:NRAD-1,:);
    vradtemp(1,:) = 0.0;
    vradint(2:NRAD,:) = gas_v_rad(2:NRAD,:)+ dt*(-gradp(2:NRAD,:)-gradphi(2:NRAD,:) + vradtemp(2:NRAD,:));
    
    gradptemp(1:NRAD,:) = (Pressure(1:NRAD,:)-Pressure(1:NRAD, mod((0:NSEC-1)+NSEC-1,NSEC)+1))*2.0./(gas_density(1:NRAD,:) +...
        gas_density(1:NRAD, mod((0:NSEC-1)+NSEC-1,NSEC)+1));
    
    gradp2 = bsxfun(@times, gradptemp, invdxtheta');
    gradphitemp(1:NRAD,:) =  (Potential(1:NRAD,:)-Potential(1:NRAD, mod((0:NSEC-1)+NSEC-1,NSEC)+1));
    gradphi2 = bsxfun(@times, gradphitemp, invdxtheta');
    
    vthetaint(1:NRAD,:) = gas_v_theta(1:NRAD,:)-dt.*(gradp2(1:NRAD,:)+gradphi2(1:NRAD,:));
    %{
    for j=1:NSEC
        invdxtheta = 1.0./2.0*pi./Rmed(1:NRAD);
        gradp(2:NRAD,j) = (Pressure(2:NRAD,j) - Pressure(1:NRAD-1,j)).*2./(gas_density(2:NRAD,j) + gas_density(1:NRAD-1,j)).*InvDiffRmed(1:NRAD-1)';
        gradphi(2:NRAD,j) = (Potential(2:NRAD,j) - Potential(1:NRAD-1,j)).*InvDiffRmed(1:NRAD-1)';
        vt2(2:NRAD,j) = gas_v_theta(2:NRAD,j) + gas_v_theta(1:NRAD-1,mod(j,NSEC)+1) + gas_v_theta(2:NRAD, j)+ gas_v_theta(1:NRAD-1, mod(j,NSEC)+1);
        vt2(2:NRAD,j) = vt2(2:NRAD,j)/4.0+Rinf(2:NRAD)'*OmegaFrame;
        vt2(2:NRAD,j) = vt2(2:NRAD,j).*vt2(2:NRAD,j);
        vradint(2:NRAD,j) = gas_v_rad(2:NRAD,j)+dt.*(-gradp(2:NRAD,j)-gradphi(2:NRAD,j)+vt2(2:NRAD,j).*InvRinf(2:NRAD)');
        
        gradp2(1:NRAD,j) = (Pressure(1:NRAD,j) - Pressure(1:NRAD, mod((j-2)+NSEC,NSEC)+1)) *2.0./(gas_density(1:NRAD,j) + ...
            gas_density(1:NRAD, mod((j-2)+NSEC,NSEC)+1)).*invdxtheta';
        gradphi2(1:NRAD,j) = (Potential(1:NRAD,j) - Potential(1:NRAD, mod((j-2)+NSEC,NSEC)+1)).*invdxtheta';
        
        vthetaint(1:NRAD,j) = gas_v_theta(1:NRAD,j)-dt.*(gradp2(1:NRAD,j)+gradphi2(1:NRAD,j));
    end
    %}
    
    ComputeViscousTerms(gas_density);
    UpdateVelocitiesWithViscosity(gas_density, dt);
    ApplySubKeplerianBoundary();    
end

function ApplySubKeplerianBoundary()
    global G Rmed SIGMASLOPE FLARINGINDEX NRAD vthetaint OmegaFrame;
    aspectratio0 = AspectRatio(Rmed(1));
    aspectratioNr = AspectRatio(Rmed(NRAD));
    VKepIn = sqrt (G*1.0/Rmed(1)*(1.0 -(1.0+SIGMASLOPE-2.0*FLARINGINDEX)* aspectratio0^2.0)*Rmed(1)^(2.0*FLARINGINDEX));
    VKepOut = sqrt (G*1.0/Rmed(NRAD)*(1.0 -(1.0+SIGMASLOPE-2.0*FLARINGINDEX)* aspectratioNr^2.0)*Rmed(NRAD)^(2.0*FLARINGINDEX));
    
    vthetaint(1,:) = VKepIn-Rmed(1)*OmegaFrame;
    vthetaint(NRAD,:) = VKepOut-Rmed(NRAD)*OmegaFrame;
    
end

function UpdateVelocitiesWithViscosity(gas_density, DeltaT)
    global NRAD NSEC InvRmed InvDiffRsup Rinf Rsup Tpp Trp Trr InvRinf Rmed InvDiffRmed;
    global vthetaint vradint;
    phi = 2.0*pi/double(NSEC);
    invdphi = 1.0/phi;
    
    Trptemp = bsxfun(@times, Rsup(2:NRAD-1)', Trp(3:NRAD,:)) - bsxfun(@times, Rinf(2:NRAD-1)', Trp(2:NRAD-1, :));
    Trptemp2 = bsxfun(@times, Trptemp, InvDiffRsup(2:NRAD-1)');
    
    vthetatemp = (Trptemp2 + (Tpp(2:NRAD-1, :)- Tpp(2:NRAD-1, mod(((1:NSEC)-2)+NSEC,NSEC)+1))*invdphi + ...
            0.5.*(Trp(2:NRAD-1, :)+Trp(3:NRAD, :)))./(0.5.*(gas_density(2:NRAD-1, :)+gas_density(2:NRAD-1, mod(((1:NSEC)-2)+NSEC,NSEC)+1)));
  
    vthetaint(2:NRAD-1,:) = vthetaint(2:NRAD-1, :) + bsxfun(@times, DeltaT.*InvRmed(2:NRAD-1)', vthetatemp);
    
    Trrtemp = bsxfun(@times, Rmed(2:NRAD)', Trr(2:NRAD,:)) - bsxfun(@times, Rmed(1:NRAD-1)', Trr(1:NRAD-1, :));
    Trrtemp2 = bsxfun(@times, Trrtemp, InvDiffRmed(1:NRAD-1)');
    
    vradtemp = (Trrtemp2 + (Trp(2:NRAD, mod(1:NSEC,NSEC)+1)-Trp(2:NRAD,:))*invdphi - ...
            0.5*(Tpp(2:NRAD,:)+Tpp(1:NRAD-1,:)))./(0.5.*(gas_density(2:NRAD,:)+gas_density(1:NRAD-1, :)));
    
    vradint(2:NRAD,:) = vradint(2:NRAD,:) + bsxfun(@times, DeltaT.*InvRinf(2:NRAD)', vradtemp);
    
    %{
    for j=1:NSEC
        vthetaint(2:NRAD-1,j) = vthetaint(2:NRAD-1, j) + DeltaT.*InvRmed(2:NRAD-1)'.*((Rsup(2:NRAD-1)'.*Trp(3:NRAD,j)- ...
            Rinf(2:NRAD-1)'.*Trp(2:NRAD-1,j)).*InvDiffRsup(2:NRAD-1)' + (Tpp(2:NRAD-1, j)- Tpp(2:NRAD-1, mod((j-2)+NSEC,NSEC)+1))*invdphi + ...
            0.5.*(Trp(2:NRAD-1, j)+Trp(3:NRAD, j)))./(0.5.*(gas_density(2:NRAD-1, j)+gas_density(2:NRAD-1, mod((j-2)+NSEC,NSEC)+1)));
        
        vradint(2:NRAD,j) = vradint(2:NRAD,j) + DeltaT.*InvRinf(2:NRAD)'.*((Rmed(2:NRAD)'.*Trr(2:NRAD,j) - ...
            Rmed(1:NRAD-1)'.*Trr(1:NRAD-1,j)).*InvDiffRmed(1:NRAD-1)' + (Trp(2:NRAD, mod(j,NSEC)+1)-Trp(2:NRAD,j))*invdphi - ...
            0.5*(Tpp(2:NRAD,j)+Tpp(1:NRAD-1,j)))./(0.5.*(gas_density(2:NRAD,j)+gas_density(1:NRAD-1, j)));
    end
    %}
    
    %fprintf('%g\n',vradint(2:NRAD,1)')
    
end

function ComputeViscousTerms(gas_density)
    global NRAD vradint vthetaint InvDiffRsup NSEC InvRmed Rsup Rinf InvDiffRmed InvRinf;
    global Rmed Trp Tpp Trr Drr Dpp Drp divergence;
    
    % InvDiffRmed parte de 1 (en fargo parte de 1 el primer valor, en vez
    % de 0 
    
    dphi = 2.0*pi/double(NSEC);
    invdphi = 1.0/dphi;
    onethird = 1.0/3.0;
    
    Drrtemp(1:NRAD,:) = (vradint(2:NRAD+1,:)-vradint(1:NRAD,:));
    Drr = bsxfun(@times, Drrtemp, InvDiffRsup(1:NRAD)');
    
    Dpptemp(1:NRAD,:) = (vthetaint(1:NRAD, mod(1:NSEC,NSEC)+1)-vthetaint(1:NRAD,:))*invdphi;
    Dppt = bsxfun(@times, Dpptemp, InvRmed(1:NRAD)');
    
    Dpptemp2 = 0.5.*(vradint(2:NRAD+1,:)+vradint(1:NRAD,:));
    Dppt2 = bsxfun(@times, Dpptemp2, InvRmed(1:NRAD)');
    
    Dpp = Dppt + Dppt2;
    
    divtemp = bsxfun(@times, vradint(2:NRAD+1,:), Rsup(1:NRAD)');
    divtemp2 = bsxfun(@times, vradint(1:NRAD,:), Rinf(1:NRAD)');
    divergence = bsxfun(@times, divtemp - divtemp2, InvDiffRsup(1:NRAD)'.*InvRmed(1:NRAD)');
    
    
    divtemp = bsxfun(@times, (vthetaint(1:NRAD, mod(1:NSEC,NSEC)+1)-vthetaint(1:NRAD,:))*invdphi, InvRmed(1:NRAD)');
    divergence(1:NRAD,:) = divergence(1:NRAD, :) + divtemp(1:NRAD,:);
    
    Drptemp = bsxfun(@times, vthetaint(2:NRAD,:), InvRmed(2:NRAD)');
    Drptemp2 = bsxfun(@times, vthetaint(1:NRAD-1,:), InvRmed(1:NRAD-1)');
    
    Drptemp3 = bsxfun(@times, Rinf(2:NRAD)', Drptemp - Drptemp2);
    Drptemp4 = bsxfun(@times, Drptemp3, InvDiffRmed(1:NRAD-1)');
    Drptemp5 = bsxfun(@times, (vradint(2:NRAD,:)-vradint(2:NRAD,mod(((1:NSEC) -2)+NSEC,NSEC)+1))*invdphi, InvRinf(2:NRAD)');
    
    Drp(2:NRAD,:) = 0.5.*(Drptemp4 + Drptemp5);
    
    %{
    for j=1:NSEC
        Drr(1:NRAD,j) = (vradint(2:NRAD+1,j)-vradint(1:NRAD,j)).*InvDiffRsup(1:NRAD)';
        Dpp(1:NRAD,j) = (vthetaint(1:NRAD, mod(j,NSEC)+1)-vthetaint(1:NRAD,j))*invdphi.*InvRmed(1:NRAD)'+ ...
            0.5.*(vradint(2:NRAD+1,j)+vradint(1:NRAD,j)).*InvRmed(1:NRAD)';
        divergence(1:NRAD,j) = (vradint(2:NRAD+1,j).*Rsup(1:NRAD)'-vradint(1:NRAD,j).*Rinf(1:NRAD)').*InvDiffRsup(1:NRAD)'.*InvRmed(1:NRAD)';
        divergence(1:NRAD,j) = divergence(1:NRAD,j) + (vthetaint(1:NRAD, mod(j,NSEC)+1)-vthetaint(1:NRAD,j))*invdphi.*InvRmed(1:NRAD)';
	    Drp(2:NRAD,j) = 0.5.*(Rinf(2:NRAD)'.*(vthetaint(2:NRAD,j).*InvRmed(2:NRAD)'-vthetaint(1:NRAD-1,j).*InvRmed(1:NRAD-1)').*InvDiffRmed(1:NRAD-1)' + ...
            (vradint(2:NRAD,j)-vradint(2:NRAD, mod((j-2)+NSEC,NSEC)+1))*invdphi.*InvRinf(2:NRAD)');
    end
    %}
    viscosity(1:NRAD) = Fviscosity(Rmed);
    
    densVisc = bsxfun(@times, gas_density(1:NRAD,:), viscosity(1:NRAD)');
    
    Trr(1:NRAD,:) = 2.0.*densVisc.*(Drr(1:NRAD,:)-onethird.*divergence(1:NRAD,:));
    Tpp(1:NRAD,:) = 2.0.*densVisc.*(Dpp(1:NRAD,:)-onethird.*divergence(1:NRAD,:));
    Trptmp(2:NRAD,:) = 2.0*0.25.*(gas_density(2:NRAD,:)+gas_density(1:NRAD-1,:)+gas_density(2:NRAD, mod((-1:NSEC-2)+NSEC,NSEC)+1) + ...
            gas_density(1:NRAD-1, mod((-1:NSEC-2)+NSEC,NSEC)+1));
        
    Trp = bsxfun(@times, Trptmp(2:NRAD,:), viscosity(2:NRAD)');
    Trp(2:NRAD,:) = Trp(1:NRAD-1,:);
    
    Trp(2:NRAD,:) = Trp(2:NRAD,:).*Drp(2:NRAD,:);
    Trp(1,:) = 0.0;
    
    %{
    for j=1:NSEC
        Trr(1:NRAD,j) = 2.0.*gas_density(1:NRAD,j).*viscosity(1:NRAD)'.*(Drr(1:NRAD,j)-onethird.*divergence(1:NRAD,j)); 
        Tpp(1:NRAD,j) = 2.0.*gas_density(1:NRAD,j).*viscosity(1:NRAD)'.*(Dpp(1:NRAD,j)-onethird.*divergence(1:NRAD,j)); 
        Trp(2:NRAD,j) = 2.0*0.25.*(gas_density(2:NRAD,j)+gas_density(1:NRAD-1,j)+gas_density(2:NRAD, mod((j-2)+NSEC,NSEC)+1) + ...
            gas_density(1:NRAD-1, mod((j-2)+NSEC,NSEC)+1)).*viscosity(2:NRAD)'.*Drp(2:NRAD,j);
    end   
    %}
    
end

function viscosity = Fviscosity(rad)
    global VISCOSITY CAVITYRADIUS CAVITYWIDTH ASPECTRATIO LAMBDADOUBLING PhysicalTime;
    global PhysicalTimeInitial CAVITYRATIO;
    
    viscosity = VISCOSITY;
    
    rmin = CAVITYRADIUS-CAVITYWIDTH*ASPECTRATIO;
    rmax = CAVITYRADIUS+CAVITYWIDTH*ASPECTRATIO;
    scale = 1.0 +(PhysicalTime-PhysicalTimeInitial)*LAMBDADOUBLING;
    rmin = rmin*scale;
    rmax = rmax*scale;
    
    if (rad < rmin) 
        viscosity = viscosity * CAVITYRATIO;
    end
    if ((rad >= rmin) & (rad <= rmax))
        viscosity = viscosity * exp((rmax-rad)/(rmax-rmin)*log(CAVITYRATIO));   
    end
end

function aspectratio = AspectRatio(rad)
    global ASPECTRATIO PhysicalTime TRANSITIONRADIUS TRANSITIONWIDTH TRANSITIONRATIO;
    global LAMBDADOUBLING PhysicalTimeInitial;
    
    aspectratio = ASPECTRATIO;
    rmin = TRANSITIONRADIUS-TRANSITIONWIDTH*ASPECTRATIO;
    rmax = TRANSITIONRADIUS+TRANSITIONWIDTH*ASPECTRATIO;
    scale = 1.0+(PhysicalTime-PhysicalTimeInitial)*LAMBDADOUBLING;
    rmin = rmin * scale;
    rmax = rmax * scale;
    if (rad < rmin) 
        aspectratio = aspectratio * TRANSITIONRATIO;
    end
    
    if ((rad >= rmin) & (rad <= rmax))
        aspectratio = aspectratio *exp((rmax-rad)/(rmax-rmin)*log(TRANSITIONRATIO));
    end  
end