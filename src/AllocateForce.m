function AllocateForce(dimfxy)
    global force;
    force = cell(9,1);
    force{9,1} = []; 

    %{
     real fx_inner;
     real fy_inner;
     real fx_ex_inner;
     real fy_ex_inner;
     real fx_outer;
     real fy_outer;
     real fx_ex_outer;
     real fy_ex_outer;
     real *GlobalForce;
    %}

    cero = zeros(1,4*dimfxy);
    force{9,1}= cero;
end