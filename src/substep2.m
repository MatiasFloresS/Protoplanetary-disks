function substep2(gas_density, dt)
   global vradint CVNR NRAD NSEC vthetaint vradnew vthetanew Rmed InvDiffRmed;
   
   dv = vradint(2:NRAD+1,:) - vradint(1:NRAD,:);
   dvt = vthetaint(1:NRAD, mod(1:NSEC, NSEC)+1) - vthetaint(1:NRAD, :);
   
   I = find(dv >= 0.0);
   temp = dv(:);
   temp(I) = 0.0;
   dv = reshape(temp, NRAD, NSEC);
   
   I = find(dv < 0.0);
   temp = dv(:);
   temp2 = gas_density(1:NRAD, :);
   temp3 = temp2(:);
   
   temp(I) = CVNR*CVNR.*temp3(I).*temp(I).*temp(I);
   densint = reshape(temp, NRAD, NSEC);
   
   I = find(dvt >= 0.0);
   temp = dvt(:);
   temp(I) = 0.0;
   dvt = reshape(temp, NRAD, NSEC);
   
   I = find(dvt < 0.0);
   temp = dvt(:);
   temp2 = gas_density(1:NRAD, :);
   temp3 = temp2(:);
   
   temp(I) = CVNR*CVNR.*temp3(I).*temp(I).*temp(I);
   tempint = reshape(temp, NRAD, NSEC);
   
   dxtheta = 2.0*pi/double(NSEC).*Rmed(1:NRAD);
   invdxtheta = 1.0./dxtheta;
   
   densintdiff = bsxfun(@times, densint(2:NRAD,:)-densint(1:NRAD-1,:), InvDiffRmed(1:NRAD-1)');
   tempintdiff = bsxfun(@times, tempint(1:NRAD,:)-tempint(1:NRAD, mod((0:NSEC-1)+NSEC-1,NSEC)+1), invdxtheta');
   
   vradnew(2:NRAD,:) = vradint(2:NRAD,:)-dt*2.0./(gas_density(2:NRAD,:)+gas_density(1:NRAD-1,:)).*(densintdiff(1:NRAD-1,:));
   vthetanew(1:NRAD,:) = vthetaint(1:NRAD,:) - dt*2.0./(gas_density(1:NRAD,:)+gas_density(1:NRAD, mod((0:NSEC-1)+NSEC-1,NSEC)+1)).* ...
       (tempintdiff(1:NRAD,:));
end