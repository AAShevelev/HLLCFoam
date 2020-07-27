//Smooth1D
/*
dRho = 0;
forAll(rho, i)
{
	err[i] = 1.0 + 1e-6*Foam::sin(2*M_PI*(mesh.C()[i][0]+ runTime.value())) - rho[i];
	dRho += err[i]*err[i]*mesh.V()[i]; 
}
dRho = Foam::sqrt(dRho);
file <<  runTime.value() << " " << dRho << std::endl;
*/

//Sod

dEn = 0;
forAll(rho, i)
{
	dEn += (rho[i]*(h[i] + K[i]) - p[i])*mesh.V()[i]; 
}
dEn = En0 - dEn;
file <<  runTime.value() << " " << dEn.value() << std::endl;


//Smooth2D
/*
dRho = 0;
forAll(rho, i)
{
	err[i] = 1.0 + 1e-6*Foam::sin(2*M_PI*(mesh.C()[i][0] + mesh.C()[i][1] - runTime.value()*(1.0 - 0.5))) - rho[i]; 
	dRho += err[i]*err[i]*mesh.V()[i]; 
}
dRho = Foam::sqrt(dRho);
file <<  runTime.value() << " " << dRho << std::endl;
*/
