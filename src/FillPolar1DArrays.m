function FillPolar1DArrays()
    global GlobalRmed Radii Rmed Rinf NRAD NSEC;
    Radii = zeros(1,NRAD+1);
    RMIN = 0.4;
    RMAX = 2.5;

    prompt = 'LogGrid? ';
    logGrid = input(prompt,'s');
    deltatheta = (RMAX-RMIN)/NRAD;
    if (strcmp(logGrid,'YES')) 
        Radii = RMIN*exp([0:NRAD]./NRAD*log(RMAX/RMIN));
        sprintf('%.18g ',Radii);

    elseif (strcmp(logGrid,'NO'))
        Radii = RMIN+deltatheta*[0:NRAD]; 
        %format long
        %Radii'
        %sprintf('%.18g ', Radii)

    else 
        sprintf('end')
    end

    %for i=0:GLOBALNRAD
    %    ang=0:0.01:2*pi; 
    %    xp=Radii(i+1)*cos(ang);
    %    yp=Radii(i+1)*sin(ang);
    %    plot(xp,yp);
    %    hold on;
    %end

    GlobalRmed = 2.0/3.0 * (Radii(2:NRAD+1).*Radii(2:NRAD+1).*Radii(2:NRAD+1)-Radii(1:NRAD).*Radii(1:NRAD).*Radii(1:NRAD));
    GlobalRmed = GlobalRmed./(Radii(2:NRAD+1).*Radii(2:NRAD+1) - Radii(1:NRAD).*Radii(1:NRAD));
    
    %format long
    %GlobalRmed'
    %sprintf('%.18g ', GlobalRmed)
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

    %format long
    %Rmed'
    %sprintf('%.18g ', Rmed)
    
    %format long
    %Rinf'
    %sprintf('%.18g ', Rinf)
    
    fileID = fopen('used_rad.dat','w');
    fprintf(fileID,'%.18g \n',Radii);
    fclose(fileID);

end