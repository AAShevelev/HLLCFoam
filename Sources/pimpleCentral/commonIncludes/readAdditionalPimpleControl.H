{
	
    pimple.dict().lookup("updateEnergyInPISO") >> updateEnergyInPISO;
	
    pimple.dict().lookup("fluxScheme") >> fluxScheme;

    if ((fluxScheme == "HLL") || (fluxScheme == "HLLC"))
    {
    }
    else
    {
        FatalErrorIn
        (
        "readAdditionalPimpleControl.H "
        )   << "fluxScheme: " << fluxScheme
            << " is not a valid choice. "
            << "Options are: HLL, HLLC"
            << abort(FatalError);
    }
	
	if (pimple.dict().found("noFluxWeighting"))
    {
        pimple.dict().lookup("noFluxWeighting") >> noFluxWeighting;
    }
}
