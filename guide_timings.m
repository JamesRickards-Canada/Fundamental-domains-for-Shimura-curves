P<x> := PolynomialRing(Rationals());

//We start with an easy one to make sure things are working
SetOutputFile("1_1_6.txt");
try
	F<y> := RationalsAsNumberField();
	A<i, j, k> := QuaternionAlgebra<F | 3, -1>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();


//Case 1
SetOutputFile("1_1_390.txt");
try
	F<y> := RationalsAsNumberField();
	A<i, j, k> := QuaternionAlgebra<F | -733, 775515>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 5
SetOutputFile("2_5_3724.txt");
try
	F<y> := NumberField(x^2-x-1);
	A<i, j, k> := QuaternionAlgebra<F | -16*y + 3, -6734*y + 5873>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 2
SetOutputFile("1_1_2145.txt");
try
	F<y> := RationalsAsNumberField();
	A<i, j, k> := QuaternionAlgebra<F | 1022, 9680385>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 7
SetOutputFile("3_81_1509.txt");
try
	F<y> := NumberField(x^3-3*x-1);
	A<i, j, k> := QuaternionAlgebra<F | -1, -4*y^2 - 80*y - 75>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 3
SetOutputFile("1_1_2021.txt");
try
	F<y> := RationalsAsNumberField();
	A<i, j, k> := QuaternionAlgebra<F | -1, 2021>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 8
SetOutputFile("3_257_423.txt");
try
	F<y> := NumberField(x^3-x^2-4*x+3);
	A<i, j, k> := QuaternionAlgebra<F | 580*y^2 + 656*y - 2851, -16*y^2 - 16*y + 7>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 9
SetOutputFile("4_14656_17.txt");
try
	F<y> := NumberField(x^4+2*x^3-4*x^2-6*x+1);
	A<i, j, k> := QuaternionAlgebra<F | -92*y^3 - 172*y^2 + 352*y + 221, 2*y^3 + 4*y^2 - 7*y - 18>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 13
SetOutputFile("6_371293_131.txt");
try
	F<y> := NumberField(x^6-x^5-5*x^4+4*x^3+6*x^2-3*x-1);
	A<i, j, k> := QuaternionAlgebra<F | 4*y^4 - 12*y^3 + 12*y - 47, -1>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 15
SetOutputFile("7_83934569_1.txt");
try
	F<y> := NumberField(x^7 - 8*x^5 + 19*x^3 - x^2 - 13*x + 1);
	A<i, j, k> := QuaternionAlgebra<F | -4*y^5 + y^4 + 20*y^3 - 2*y^2 - 16*y - 3, -10*y^6 - 2*y^5 + 58*y^4 + 21*y^3 - 74*y^2 - 20*y - 15>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 11
SetOutputFile("5_24217_415.txt");
try
	F<y> := NumberField(x^5-5*x^3-x^2+3*x+1);
	A<i, j, k> := QuaternionAlgebra<F | -104*y^4 + 116*y^3 + 672*y^2 - 660*y - 867, y^2 - y - 10>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 4
SetOutputFile("1_1_5111.txt");
try
	F<y> := RationalsAsNumberField();
	A<i, j, k> := QuaternionAlgebra<F | 2070, -45999>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 6
SetOutputFile("2_44_283.txt");
try
	F<y> := NumberField(x^2-11);
	A<i, j, k> := QuaternionAlgebra<F | -1, -7*y - 16>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 16
SetOutputFile("7_119292949_1.txt");
try
	F<y> := NumberField(x^7 - 2*x^6 - 6*x^5 + 8*x^4 + 12*x^3 - 5*x^2 - 6*x - 1);
	A<i, j, k> := QuaternionAlgebra<F | -4*y^6 + 18*y^5 - 11*y^4 - 37*y^3 + 25*y^2 + 22*y - 1, -21*y^6 + 17*y^5 + 176*y^4 - 20*y^3 - 400*y^2 - 183*y - 33>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 14
SetOutputFile("6_300125_491.txt");
try
	F<y> := NumberField(x^6-x^5-7*x^4+2*x^3+7*x^2-2*x-1);
	A<i, j, k> := QuaternionAlgebra<F | -1, -8*y^5 + 4*y^4 + 48*y^3 + 32*y^2 - 24*y - 59>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 12
SetOutputFile("5_240881_35.txt");
try
	F<y> := NumberField(x^5 - 2*x^4 - 5*x^3 + 9*x^2 + 5*x - 7);
	A<i, j, k> := QuaternionAlgebra<F | 8*y^4 - 104*y^3 + 164*y^2 + 188*y - 615, y - 11>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 17
SetOutputFile("8_309593125_19.txt");
try
	F<y> := NumberField(x^8 - 4*x^7 - x^6 + 17*x^5 - 5*x^4 - 23*x^3 + 6*x^2 + 9*x - 1);
	A<i, j, k> := QuaternionAlgebra<F | -1, -16*y^7 + 68*y^6 - 20*y^5 - 204*y^4 + 168*y^3 + 120*y^2 - 140*y - 3>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 18
SetOutputFile("8_456768125_11.txt");
try
	F<y> := NumberField(x^8 - 2*x^7 - 7*x^6 + 11*x^5 + 14*x^4 - 18*x^3 - 8*x^2 + 9*x - 1);
	A<i, j, k> := QuaternionAlgebra<F | -1, -12*y^7 + 16*y^6 + 76*y^5 - 48*y^4 - 96*y^3 + 28*y^2 + 16*y - 35>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();

//Case 10
SetOutputFile("4_1957_2511.txt");
try
	F<y> := NumberField(x^4 - 4*x^2 - x + 1);
	A<i, j, k> := QuaternionAlgebra<F | -1, -48*y^3 + 36*y^2 + 120*y - 111>;
	O := MaximalOrder(A);
	Discriminant(F);
	Norm(Discriminant(A));
	Gamma := FuchsianGroup(O);
	time FD := FundamentalDomain(Gamma);
	FD;
catch e
	print "Error:", e`Object;
end try;
UnsetOutputFile();