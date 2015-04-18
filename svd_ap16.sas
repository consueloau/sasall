

ODS TAGSETS.RTF STYLE= SEASIDE   FILE ="CHECKING_svd_apr172015.RTF";

TITLE "SVD";
proc iml;
U={-0.48 -0.02 -0.59  0.42 -0.31 0.8,
   -0.45  0.42  0.28  0.51  0.37 0.38,
   -0.26  0.29  0.46 -0.24  0.12 0.76,
   -0.55  0.13  0.05 -0.57 -0.46 -0.38,
   -0.41 -0.57 -0.16 -0.26  0.64  0,
   -0.14 -0.63  0.58  0.34 -0.36  0 } ;
   
   S={2.21    0      0      0     0,
         0  1.71     0      0     0,
		 0     0  1.31      0     0,
		 0     0     0   1.11     0,
		 0     0     0      0  0.48,  
		 0     0     0      0     0};
		
	VT= { -0.42	-0.57	-0.66	-0.25	-0.06,
	      0.23	0.49	-0.28	-0.70	-0.37,
		  -0.24	0.60	-0.53	0.32	0.44,
		  0.83	-0.27	-0.37	0.07	0.31,
		  0.12	0.05	-0.27	0.58	-0.75};

   
   print U;

   print S;
   
   print VT;
   
   chk1= (U1* S1 )*VT ;

  PRINT  'COMPARE WITH ORIGINAL MATRIX D';
   print chk1; *** should be same as D below;

*******; 
   TITLE ' REDUCED MATRICES SELECTED';
 U1 = U[, 1:2];
 S1 = S[1:2,1:2];
VT1 = VT[1:2,];

   print U1;
   
   print S1;
   
   print VT1;
   
   chk2 = (U1* S1 )*VT1;
PRINT 'COMPARE WITH REDUCED MATRIX FROM NOTES';
   print chk2;
****************************************************;

PRINT  'START SVD FROM ORIGINAL D MATRIX';
   D={1  0  1  0  0,
      1  1  0  0  0,
      0  1  0  0  0,
      0  1  1  0  0,
      0  0  1  1  0,
      0  0  0  1  1 } ;

*** GET SVD MATRICES;
call svd (Uv, qV, VV, d);
print uv qv vv;
CHK3 = (UV* diag(qv))*t(Vv);

PRINT 'COMPARE WITH PRIGINAL MATRIX D';
print CHK3;
**********;

*** GET REDUCED MATRICES;
 U1= Uv[, 1:2];
 S1= diag(qv[1:2]);
VT1= (VV)[,1:2 ];

   print U1;
   
   print S1;
   
   print VT1;
   
   CHK4 = (U1*S1)* t(VT1);

PRINT   'COMPARE WITH REDUCED MATRIX IN NOTES'  ;
   print CHK4 ;  *** COMPARE WITH RESULTS IN NOTES;
**********************************************************;
   *** SECOND EXAMPLE ************************************;
****  http://web.mit.edu/be.400/www/SVD/Singular_Value_Decomposition.htm  ***;
   PRINT 'EXAMPLE 2';
   D={2  4  ,
      1  3  ,
      0  0  ,
      0  0 } ;

*** GET D*T(D);
	  dd=D*t(D);

call svd (UD1, qD1, VD1, dd);
	  
   print dd ;

   print UD1;
   
   print qD1;
   
   print VD1;
   
   CHKD1= (UD1)*(diag(qD1))* t(VD1);
PRINT 'COMPARE WITH ORIGINAL MATRIX D*T(D) SECOND EX';
   print CHKD1;
   
*** GET T(D)*D;

	  dd2= t(D)*D;
	  call svd (UD2, qD2, VD2, dd2);

   print dd2;
   print UD2;
   
   print qD2;
   
   print VD2;
   
   CHKD2= (UD2)*(diag(qD2))* t(VD2);

PRINT 'COMPARE WITH ORIGINAL MATRIX  T(D0)*D SECOND EX';
   print CHKD2;

*****;
*** DO SVD IN MATRIX D;

call svd (UD, qD, VD, d);

print uD qD vD;

    CHKD = (UD* diag(qD))*t(VD);

PRINT 'COMPARE WITH ORIGINAL MATRIX D  SECOND EX';
   print CHKD;

*** REDUCE MATRICES  ***;

UDR= UD[, 1:2];
SDR=diag(qD[1:2]);
VDR=(VD)[,1:2 ];
   print UDR;
   
   print SDR;
   
   print VDR;
   
   CHKDR= (UDR*SDR)* t(VDR);

PRINT 'COMPARE WITH REDUCED MATRIX -  SECOND EX';
   print CHDR;
QUIT;

ODS TAGSETS.RTF CLOSE;
