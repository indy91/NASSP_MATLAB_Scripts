function GenerateSaturnIBLVDCLaunchAzimuthPolynomial

  %Generate the constants for the Saturn IB LVDC launch azimuth polynomial
  %The polynomial has the form: A_Z = A00 + A01*INCL + A10*NODE + A11*INCL*NODE

  %INPUTS
  LAT = 28.627151;  %Launchpad latitude
  INCL = 51.78;     %Desired orbit inclination
  ASC = 1;          %1 = ascending node launch, 0 = descending node launch
  V_ORBIT = 7835;   %Velocity at orbital insertion


  %CONSTANTS
  RAD = pi/180;
  DEG = 180/pi;
  LW_LIMIT = 3*60;                %Increment in time to generate the launch azimuth array, seconds
  D_INCL = 1.0;                   %Increment in inclination to generate the launch azimuth array, degrees
  R_EARTH = 6.373338e6;           %Radius of Earth
  EARTH_RATE = 7.29211514667e-5;  %Rotational rate of Earth
  DR = 3.1e-1;                    %Downrange angle for ascent, value is 17.762 degrees
  DTOPT = 6.0*60.0 + 33;

  %INTERNAL

  %Convert to internal unit
  LAT = LAT*RAD;
  INCL = INCL*RAD;
  D_INCL = D_INCL*RAD;

  %Generate the launch azimuths
  Incl_arr = [INCL INCL INCL+D_INCL INCL+D_INCL];
  DT_arr = [0 LW_LIMIT 0 LW_LIMIT];

  AZL_arr = zeros(4,1);

  for i=1:4
    [AZL_arr(i), Lambda_arr(i)] = OptimumLaunchAzimuth(LAT, Incl_arr(i), DT_arr(i), ASC, R_EARTH, EARTH_RATE, V_ORBIT, DR, DTOPT);
  endfor

  %AZL_arr/RAD
  %Lambda_arr/RAD

  %Calculate matrix
  for i=1:4
    A(i,1) = 1;
    A(i,2) = Incl_arr(i)*DEG;
    A(i,3) = Lambda_arr(i)*DEG;
    A(i,4) = Incl_arr(i)*Lambda_arr(i)*DEG*DEG;
  endfor

  %A
  COEF = A^-1*(AZL_arr*DEG);

  %Output prints
  printf("\nSATURN IB LVDC LAUNCH AZIMUTH POLYNOMIAL CALCULATOR 1.1\n");
   printf("\nLaunch latitude %.3f deg, Inclination %.3f deg\n", LAT*DEG, INCL*DEG);
  printf("\nSimulated launches:\n");
  for i=1:4
    printf("\nLaunch %d, inclination %.2f deg, delay of %.1f seconds from optimum\n", i, Incl_arr(i)*DEG, DT_arr(i));
    printf("Launch Azimuth %.4f deg, descending node angle %.4f deg\n", AZL_arr(i)*DEG, Lambda_arr(i)*DEG);
  endfor

  printf("\nLVDC Coefficients:\n");
  for i=1:4
    printf("LVDC_Ax[%d] %f\n",i-1,COEF(i));
  endfor
  printf("\nRTCC Constants:\n");
  printf("LAZCOE %f %f %f %f\n",COEF(1)*RAD, COEF(2), COEF(3), COEF(4)*DEG);

endfunction

  %To generate data for the LVDC launch azimuth polynomial
  function [AZL, DNA] = OptimumLaunchAzimuth(LAT_C, Incl, DT, ASC, R_EARTH, EARTH_RATE, v_orbit, DR, DTOPT)

    %Calculate speed at equator
    v_eqrot = EARTH_RATE*R_EARTH;

    %Inertial launch azimuth
    arg = cos(Incl)/cos(LAT_C);
    if arg > 1
      arg = 1;
    endif
    beta = asin(arg);
    if ASC == 0
      beta = pi - beta;
    endif

    %Compensate for rotating Earth
    v_rotx = v_orbit*sin(beta) - v_eqrot*cos(LAT_C);
    v_roty = v_orbit*cos(beta);
    AZP = atan2(v_rotx, v_roty);
    if AZP < 0
      AZP = AZP + 2*pi;
    endif

    %Now the Shuttle equations for accounting for launch off nominal time
    WE_DT = EARTH_RATE*DT;
    A = 2*asin(sin(WE_DT/2)*cos(LAT_C));
    if abs(A - DR) < 1e-4
      A = DR - 1e-4;
    endif

    PHI = atan2(1 + cos(WE_DT), sin(WE_DT)*sin(LAT_C));
    THET = asin(sin(A)*sin(PHI-beta)/sin(DR));
    ALP = 2*atan2(-sin((A - DR)/2)*cos((THET - PHI + beta)/2), -sin((A + DR)/2)*sin((THET - PHI + beta)/2));
    DELTA_PSI_TEMP = ALP - PHI - THET - beta;
    %DELTA_PSI_TEMP = MIDVAL(-DELTA_PSI_LIM, DELTA_PSI_TEMP, DELTA_PSI_LIM);
    AZL = AZP + DELTA_PSI_TEMP;

    %Descending node angle
    bias = EARTH_RATE*DTOPT;
    dlng = atan2(sin(LAT_C),cos(beta)/sin(beta));
    if dlng < 0
      dlng = dlng + 2*pi;
    endif
    if pi - beta < 0
      dlng = pi + dlng;
    endif
    h = pi - dlng + bias; %lng
    if h < 0
      h = h + 2*pi;
    endif
    DNA = h - WE_DT;

  endfunction

  function D = MIDVAL(A, B, C)
    if B < A
      D = A;
    elseif B > C
      D = C;
    else
      D = B;
    endif
  endfunction
