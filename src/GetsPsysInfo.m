function p = GetsPsysInfo(system, action)
    global G GuidingCenter;

    Xplanet = 0;
    Yplanet = 0;
    x = system{3,1}(1);
    xc = x;
    y = system{4,1}(1);
    yc = y;
    vx = system{5,1}(1);
    vxc = vx;
    vy = system{6,1}(1);
    vyc = vy;
    m = system{2,1}(1)+1;
    h = x*vy-y*vx;
    d = sqrt(x*x+y*y);
    Ax = x*vy*vy-y*vx*vy-G*m*x/d;
    Ay = y*vx*vx-x*vx*vy-G*m*y/d;
    e = sqrt(Ax*Ax+Ay*Ay)/m;
    a = h*h/G/m/(1.-e*e);    
    
    if (e == 0.0)
        arg = 1.0;
    else
        arg = (1.0-d/a)/e;
    end
    
    if (abs(arg) >= 1.0)
        E = PI *(1.-arg/abs(arg))/2.;
    else
        E = acos((1.0-d/a)/e);
    end    
    
    if ((x*y*(vy*vy-vx*vx)+vx*vy*(x*x-y*y)) < 0) 
        E= -E;
    end
    
    M = E-e*sin(E);
    omega = sqrt(m/a/a/a);
    PerihelionPA=atan2(Ay,Ax);
    if (GuidingCenter == 1) % true or false 
        xc = a*cos(M+PerihelionPA);
        yc = a*sin(M+PerihelionPA);
        vxc = -a*omega*sin(M+PerihelionPA);
        vyc =  a*omega*cos(M+PerihelionPA);
    end
    
    if (e < 1e-8) 
        xc = x;
        yc = y;
        vxc = vx;
        vyc = vy;
    end    
    switch action
        
        case 0
            Xplanet = xc;
            Yplanet = yc;
            p = 0.0;

        case 1
            x = xc;
            y = yc;
            vx = vxc;
            vy = vyc;
            d2 = sqrt(x*x+y*y);
            d1 = sqrt(Xplanet*Xplanet+Yplanet*Yplanet);
            cross = Xplanet*y-x*Yplanet;
            Xplanet = x;
            Yplanet = y;
            p = asin(cross/(d1*d2));    
 
        case 2
            p = omega;
    end  
end