
//Smooth1D
/*
file.close();
scalar errL2 = 0;

forAll(err, i)
{
	err[i] = 1.0 + 1e-6*Foam::sin(2*M_PI*(mesh.C()[i][0]+ runTime.value())) - rho[i];
	errL2 += err[i]*err[i]*mesh.V()[i];
}
errL2 = Foam::sqrt(errL2);
Info << "ErrL2 = " << errL2 << endl;
*/

/*
//hErr

volScalarField rhoTrue = rho;
volScalarField heTrue = he;
forAll(rho, i)
	rhoTrue[i] = 1.0 + 1e-6*Foam::sin(2*M_PI*(mesh.C()[i][0]+ runTime.value())); 
he = (thermo.Cp()*pTrue*thermo.W())/(rhoTrue*1e3*Foam::constant::physicoChemical::R);

forAll(err, i)
{
	err[i] = heTrue[i] - he[i];
	errL2 += err[i]*err[i]*mesh.V()[i];
}
errL2 = Foam::sqrt(errL2);
Info << "ErrL2 = " << errL2 << endl;
*/

//Sod

file.close();
scalar EnFin = 0.0;

forAll(rho, i)
{
	
	EnFin += (rho[i]*(h[i] + K[i]) - p[i])*mesh.V()[i]; 
}
	
Info << "dE = " << En0 - EnFin << endl;


//Smooth2D
/*
file.close();
scalar errL2 = 0;

forAll(err, i)
{
	err[i] = 1.0 + 1e-6*Foam::sin(2*M_PI*(mesh.C()[i][0] + mesh.C()[i][1] - runTime.value()*(1.0 - 0.5))) - rho[i];
	errL2 += err[i]*err[i]*mesh.V()[i];
}
errL2 = Foam::sqrt(errL2);
Info << "ErrL2 = " << errL2 << endl;						   
*/
