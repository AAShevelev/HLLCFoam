MRF.correctBoundaryVelocity(U);
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho,U)
    +
    fvm::div(phi_own,U)
    +
    fvm::div(phi_nei,U)
    +
    MRF.DDt(rho,U)
    +
    turbulence->divDevRhoReff(U)
    +
    fvc::div((b_own*rho_own + b_nei*rho_nei)*mesh.Sf())
    ==
    fvOptions(rho, U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

#ifdef PIMPLE_RELAX
UEqn.relax();
#endif

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn == -gradp
    );

    fvOptions.correct(U);
}
