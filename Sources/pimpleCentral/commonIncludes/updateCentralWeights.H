    {
        if (noFluxWeighting)
        {
            phiv_own    = (fvc::interpolate(U, own, "reconstruct(U)")) & mesh.Sf();
            phiv_nei    = (fvc::interpolate(U, nei, "reconstruct(U)")) & mesh.Sf();
        }
        else
        {
            rho_own     = fvc::interpolate(rho, own, "reconstruct(rho)");
            rho_nei     = fvc::interpolate(rho, nei, "reconstruct(rho)");
            p_own     = fvc::interpolate(p, own, "reconstruct(p)");
            p_nei     = fvc::interpolate(p, nei, "reconstruct(p)");
            pS_own       = p_own*mesh.magSf()*mesh.magSf();
            pS_nei       = p_nei*mesh.magSf()*mesh.magSf();
            //correctCentralACMIInterpolation(rho_nei);
            phiv_own    = (fvc::interpolate(rho*U, own, "reconstruct(U)") / rho_own) & mesh.Sf();
            phiv_nei    = (fvc::interpolate(rho*U, nei, "reconstruct(U)") / rho_nei) & mesh.Sf();
        }
        
        if (mesh.changing())
        {
            phiv_own   -= mesh.phi();
            phiv_nei   -= mesh.phi();
        }

        c           = sqrt(thermo.Cp()/thermo.Cv() / psi);
        cf_own      = fvc::interpolate(c, own, "reconstruct(psi)");
        cf_nei      = fvc::interpolate(c, nei, "reconstruct(psi)");
        //correctCentralACMIInterpolation(cf_nei);
        cSf_own     = cf_own * mesh.magSf();
        cSf_own.setOriented(true);
        cSf_nei     = cf_nei * mesh.magSf();
        cSf_nei.setOriented(true);
        
        ap = max(max(phiv_own + cSf_own, phiv_nei + cSf_nei), v_zero);
        am = min(min(phiv_own - cSf_own, phiv_nei - cSf_nei), v_zero);
        aStar = (pS_nei - pS_own + rho_own*phiv_own*(am - phiv_own) - rho_nei*phiv_nei*(ap - phiv_nei))/(rho_own*(am - phiv_own) - rho_nei*(ap - phiv_nei));
		aStar_own = max(aStar, v_zero);
		aStar_nei = min(aStar, v_zero);
        if (fluxScheme == "HLL")
	    {
           alpha_own   = ap/(ap - am);
           aSf_own     = am*alpha_own;
           alpha_nei   = 1.0 - alpha_own;
           aSf_nei     = -aSf_own;
	    } else if (fluxScheme == "HLLC")
		{
			/*forAll(agrp_own, i)
			{
				if (aStar[i] >= 0.0)
				{
					agrp_own[i]   = 1.0;
					agrp_nei[i]   = 0.0; 
				} else 
				{
					agrp_own[i]   = 0.0;
					agrp_nei[i]   = 1.0; 
				} 
			}*/
			alpha_own   = aStar_own/(aStar_own - am);
			alpha_nei   = aStar_nei/(aStar_nei - ap);
			aSf_own     = am*alpha_own;
			aSf_nei     = ap*alpha_nei;			
			agrp_own    = aStar_own/max(aStar_own,v_eps);
			agrp_nei    = 1.0 - agrp_own;
			b_own       = (agrp_own - alpha_own)*(phiv_own - am)*(phiv_own - aStar_own)/mesh.magSf()/mesh.magSf();
			b_nei       = (agrp_nei - alpha_nei)*(phiv_nei - ap)*(phiv_nei - aStar_nei)/mesh.magSf()/mesh.magSf();
		}
    }
