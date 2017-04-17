function FillPolar1DArrays()
    global GlobalRmed Radii Rmed Rinf NRAD NSEC RMAX RMIN Surf InvDiffRmed RADIALSPACING Rsup ...
        InvRmed InvSurf InvDiffRsup InvRinf;

    deltatheta = (RMAX-RMIN)/NRAD;
    
    if (strcmp(RADIALSPACING,'L')) 
        Radii = RMIN*exp((0:NRAD)./NRAD*log(RMAX/RMIN));

    elseif (strcmp(RADIALSPACING,'NO'))
        Radii = RMIN+deltatheta*(0:NRAD); 
    else 
        quit cancel;
    end

    %GlobalRmed = 2.0/3.0 * (Radii(2:NRAD+1).*Radii(2:NRAD+1).*Radii(2:NRAD+1)-...
    %    Radii(1:NRAD).*Radii(1:NRAD).*Radii(1:NRAD));
    %GlobalRmed = GlobalRmed./(Radii(2:NRAD+1).*Radii(2:NRAD+1) - Radii(1:NRAD).*Radii(1:NRAD));
    Rinf = Radii(1:NRAD);
    Rsup = Radii(2:NRAD+1);

    Rmed = 2.0/3.0 * (Rsup.*Rsup.*Rsup - Rinf.*Rinf.*Rinf);
    Rmed = Rmed ./(Rsup.*Rsup - Rinf.*Rinf);
    
    Surf = pi * (Rsup.*Rsup - Rinf.*Rinf)/NSEC;

    InvRmed = 1.0./Rmed;
    InvSurf = 1.0./Surf;
    InvDiffRsup = 1.0./(Rsup-Rinf);
    InvRinf = 1.0./Rinf;
    InvDiffRmed = 1.0./ (Rmed(2:NRAD) - Rmed(1:NRAD-1));   
    
    %{
    Rinf(NRAD+1) = Radii(NRAD+1);
    
    fileID = fopen('used_rad.dat','w');
    fwrite(fileID, Radii, 'double');
    fclose(fileID);
    %}
end

   
   