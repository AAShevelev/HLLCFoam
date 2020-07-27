fvScalarMatrix hEqn
(
    fvm::ddt(rho,h)
    +
    fvm::div(phi_own,h)
    +
    fvm::div(phi_nei,h)
    -
    fvm::laplacian(turbulence->alphaEff(), h)
    +
	fvc::div(b_own*rho_own*aStar_own + b_nei*rho_nei*aStar_nei)
    +
    KChange
    ==
    KSource
    +
    fvOptions(rho,h)
);

hEqn.relax();

fvOptions.constrain(hEqn);

hEqn.solve();

fvOptions.correct(h);

thermo.correct();

Info << "max/min T: " << max(T).value() << "/" << min(T).value() << endl;
