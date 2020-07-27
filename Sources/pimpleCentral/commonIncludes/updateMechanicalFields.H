    K = 0.5*magSqr(U);
    
    //add convection part
    
   KChange = 
                fvc::div(phi_own,K)
                +
                fvc::div(phi_nei,K);

    KSource = 
            fvc::div
            (
                (
                    linearInterpolate
                    (
                        (-turbulence->devRhoReff()) & U
                    ) & mesh.Sf()
                )()
            );
    
    //add time derivative part
    if (LTS)
    {
        KChange +=
            (rho*K - rho.oldTime()*K.oldTime())*trDeltaT.ref();
        KSource +=
            (p - p.oldTime())*trDeltaT.ref();
    }
    else
    {
        KChange +=
            (rho*K - rho.oldTime()*K.oldTime())/runTime.deltaT();
        KSource +=
            (p - p.oldTime())/runTime.deltaT();
    }


