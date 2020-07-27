//Creating internal field

//Smooth1D test

//volScalarField pTrue = p;
/*
	forAll(rho, i)
		rho[i] = 1.0 + 1e-6*Foam::sin(2*M_PI*mesh.C()[i][0]); 
	he = (thermo.Cp()*p*thermo.W())/(rho*1e3*Foam::constant::physicoChemical::R);
	K = 0.5*magSqr(U);
	thermo.correct();
	
	volScalarField err(rho);
	scalar dRho = 0;	
	std::ofstream file("dRho.txt");
*/

//SOD test

	/*forAll(rho, i)
	{
		if (mesh.C()[i][0] < 0.5)
		{
			p[i] = 1.0;
			rho[i] = 1.0;
		}
		else
		{
			p[i] = 0.1;
			rho[i] = 0.125;
		}
		U[i][0] = 0.0;
		U[i][1] = 0.0;
		U[i][2] = 0.0;
	}
	he = (thermo.Cp()*p*thermo.W())/(rho*1e3*Foam::constant::physicoChemical::R);
	K = 0.5*magSqr(U);*/
	scalar En0 = 0;
	thermo.correct();
	forAll(rho, i)
	{	
		En0 += (rho[i]*h[i] - p[i])*mesh.V()[i]; 
	}
	
	dimensionedScalar dEn = En0; 
	std::ofstream file("Energy.txt");


//Smooth2D test
/*
	forAll(rho, i)
	{
		rho[i] = 1.0 + 1e-6*Foam::sin(2*M_PI*(mesh.C()[i][0] + mesh.C()[i][1])); 
		U[i][0] = 1.0;
		U[i][1] = -0.5;
		U[i][2] = 0.0;
		p[i] = 1.0;
	}
	he = (thermo.Cp()*p*thermo.W())/(rho*1e3*Foam::constant::physicoChemical::R);
	K = 0.5*magSqr(U);
	thermo.correct();
	
	volScalarField err(rho);
	scalar dRho = 0;	
	std::ofstream file("dRho.txt");
*/

//Shu-Osher test
/*
	forAll(rho, i)
	{
		if (mesh.C()[i][0] < -4.0)
		{
			U[i][0] = 2.629369;
			p[i] = 10.33333;
			rho[i] = 3.857143;
		}
		else
		{
			U[i][0] = 0.0;
			p[i] = 1.0;
			rho[i] = 1.0 + 0.2*Foam::sin(5.0*mesh.C()[i][0]); 
		}
		U[i][1] = 0.0;
		U[i][2] = 0.0;
	}

	he = (thermo.Cp()*p*thermo.W())/(rho*1e3*Foam::constant::physicoChemical::R);
	K = 0.5*magSqr(U);
	thermo.correct();
*/

//Acoustic test
/*
	forAll(rho, i)
	{
		U[i][0] = 0.0;
		U[i][1] = 0.0;
		U[i][2] = 0.0;
		rho[i] = 1.0 + 1e-6*Foam::exp(-2*(mesh.C()[i][0] - 0.5)*(mesh.C()[i][0] - 0.5) - 2*(mesh.C()[i][1] - 0.5)*(mesh.C()[i][1] - 0.5) ); 
		p[i] = 1.0;
	}
	he = (thermo.Cp()*p*thermo.W())/(rho*1e3*Foam::constant::physicoChemical::R);
	K = 0.5*magSqr(U);
	thermo.correct();
*/







