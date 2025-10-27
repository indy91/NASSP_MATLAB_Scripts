function LOI2PhantomVector

  %INPUTS
  %MJD of midnight of launch day, days
  GMTBASE = 40418.0;
  %Launch time (GMT), seconds
  t_launch = 13*3600 + 32*60;
  %Time (GET) of CSM pass over landing site for lunar landing, seconds
  t_land = 102*3600 + 47*60 + 11 - 4.2*60;
  %Landing site latitude, degrees
  lat_LS = 0.71388889;
  %Landing site longitude, degrees
  lng_LS = 23.70777778;
  %Landing site elevation, nautical miles
  alt_LS = -1.44;
  %Approach azimuth, degrees
  azi_land = -91;
  %Time (GET) of circular orbit, seconds
  t_circ = 125*3600;
  %Circular altitude, nautical miles
  h_circ = 60;
  %Time (GET) of LOI-2, seconds
  t_LOI2 = 79*3600;
  %Year defining the coordinate system (usually the year containing the January 1st closest to launch day)
  Epoch = 1970;

  %INTERNAL
  close all;
  Constants = ConstantsFunction();
  RAD = pi/180;
  U_X = [1 0 0]';
  U_Y = [0 1 0]';
  U_Z = [0 0 1]';
  r_M = 1738090.0;
  J2 = 0.207108e-3;
  J3 = -2.1e-5;
  J4 = 0.0;
  J22 = 0.20716e-4;
  C31 = 3.4e-5;
  mu = 0.4902778e13;
  R_E = 6378165.0;
  HRS = 3600.0;

  %Convert to internal units
  lat_LS = lat_LS*RAD;
  lng_LS = lng_LS*RAD;
  azi_land = azi_land*RAD;
  h_circ = h_circ*1852.0;
  alt_LS = alt_LS*1852.0;

  %Preliminary calculations
  GMT_land = t_launch + t_land;
  MJD_land = GMTBASE + GMT_land/24.0/3600.0;
  GMT_circ = t_launch + t_circ;
  GMT_LOI2 = t_launch + t_LOI2;
  MAT_J2000_BRCS = J2000EclToBRCS(Epoch);

  %Initial guess
  r_LS = r_M + alt_LS;
  r0 = r_LS + h_circ;
  gamma0 = 0.0;
  x0 = [r0/R_E,gamma0];

  %Enforce semi-major axis of desired height
  a0 = r0;

  %Test
  %phi(x);

  printf("Generating LOI-2 phantom state vector\n");
  printf("Converging on zero eccentricity:\n");
  seconds1 = time();

  %Optimize
  lb = [x0(1)-10.0*1852.0/R_E,-5.0*RAD];
  ub = [x0(1)+10.0*1852.0/R_E,5.0*RAD];
  maxiter = 100;
  tol = 1e-5;
  [xf, obj, info, iter, nf, lambda] = sqp (x0, @phi, [], [], lb, ub, maxiter, tol);
  %To test circular at input time
  %xf = x0;

  %Get solution state vector
  [R_land_ecl, V_land_ecl] = CalcInitialStateVector(xf);

  %Propagate to time of LOI-2
  [R_LOI2, V_LOI2] = Propagate(R_land_ecl, V_land_ecl, GMT_land, GMT_LOI2);

  %Apollo 11 Test
  %R_LOI2 = [0.22541097 0.16822345 0.071353489]';
  %V_LOI2 = [0.57390857 -0.65773378 -0.28288983]';
  %R_LOI2 = MAT_J2000_BRCS'*R_LOI2*R_E;
  %V_LOI2 = MAT_J2000_BRCS'*V_LOI2*R_E/HRS;

  %Propagate to time of circularization and store data
  [R_circ, V_circ, t_arr, y_arr] = Propagate(R_LOI2, V_LOI2, GMT_LOI2, GMT_circ);

  for i=1:length(y_arr);
    R1 = [y_arr(i,1) y_arr(i,2) y_arr(i,3)]';
    h_arr(i) = norm(R1) - r_LS;
  endfor

  plot((t_arr-t_launch)/3600, h_arr/1852);
  title("Orbit from LOI-2 to Circularization");
  xlabel("GET in hours");
  ylabel("Altitude in nautical miles");

  %Output state vector
  R_out = MAT_J2000_BRCS*R_LOI2;
  V_out = MAT_J2000_BRCS*V_LOI2;
  R_out_ER = R_out/R_E;
  V_out_ER = V_out/R_E*HRS;
  HRS_LOI2 = floor(GMT_LOI2/3600.0);
  MIN_LOI2 = floor(GMT_LOI2/60.0 - HRS_LOI2*60.0);
  SECS_LOI2 = GMT_LOI2 - MIN_LOI2*60.0 - HRS_LOI2*3600.0;


  seconds2 = time();
  printf("\nThe calculation took %f seconds to complete\n", seconds2 - seconds1);

  printf("\nOutput State Vector:\n");
  printf("Time: %f\n", t_LOI2);
  printf("%f %f %f %f %f %f\n\n", R_out(1), R_out(2), R_out(3), V_out(1), V_out(2), V_out(3));
  printf("MED Format:\n");
  printf("S84,LEM,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.0f:%.0f:%05.2f,ILHU001,MCI;\n\n", R_out_ER(1), R_out_ER(2), R_out_ER(3), V_out_ER(1), V_out_ER(2), V_out_ER(3), HRS_LOI2, MIN_LOI2, SECS_LOI2);

  function f = phi(x)
    %x(1) = r
    %x(2) = FPA

    %Calculate initial state vector
    [R_land_ecl, V_land_ecl] = CalcInitialStateVector(x);
    %Propagate to t_circ
    [R_circ, V_circ] = Propagate(R_land_ecl, V_land_ecl, GMT_land, GMT_circ);
    %Calculate dependent variables
    r = norm(R_circ);
    H = cross(R_circ, V_circ);
    E = cross(V_circ, H)/mu-R_circ/r;
    e = norm(E);

    f = e;

    printf("Radius %f Er, Flight Path Angle %f deg, Eccentricity: %f\n", x(1), x(2)/RAD, e);
  endfunction

  function [R_init, V_init] = CalcInitialStateVector(x)
    %Calculate radius in meters
    r = x(1)*R_E;
    %Calculate velocity
    v = sqrt(2.0*(mu/r - mu/(2.0*a0)));
    %Flight-path angle
    gamma = x(2);
    %Calculate selenographic state vector
    [R_land_sg, V_land_sg] = SphericalToCartesian(r, v, lat_LS, lng_LS, gamma, azi_land);
    %Convert to inertial coordinates
    Rot = oapiGetRotationMatrixEfficient(MJD_land);
    R_init = rhmul(Rot, R_land_sg);
    V_init = rhmul(Rot, V_land_sg);
  endfunction

  function [R1, V1, t_arr, y_arr] = Propagate(R0, V0, GMT0, GMT1)

    tspan = [GMT0 GMT1];
    ic = [R0(1) R0(2) R0(3) V0(1) V0(2) V0(3)];
    options = odeset('RelTol', 1e-8);
    [t_arr,y_arr] = ode45(@(t,y) GravityModel(t,y), tspan, ic, options);
    R1 = [y_arr(end,1) y_arr(end,2) y_arr(end,3)]';
    V1 = [y_arr(end,4) y_arr(end,5) y_arr(end,6)]';
  endfunction

  function dy = GravityModel(t, y)
    R = [y(1) y(2) y(3)]';
    V = [y(4) y(5) y(6)]';
    MJD = GMTBASE + t/24/3600;
    accel = R2Model(R, MJD);
    dy = [y(4) y(5) y(6) accel(1) accel(2) accel(3)];
  endfunction

  function accel = R2Model(R, MJD)
    accel = -mu*R/(norm(R)^3);

    Rot = oapiGetRotationMatrixEfficient(MJD);
    Rot = MatrixRH_LH(Rot);

    r = norm(R);
    U_R = unit(R);
    U_X = Rot*[1 0 0]';
    U_Y = Rot*[0 1 0]';
    U_Z = Rot*[0 0 1]';
    cos_lat = dot(U_R, U_Z);
    R_fixed = Rot'*R;
    x_M=R_fixed(1);
    y_M=R_fixed(2);
    z_M=R_fixed(3);

    P2 = 3*cos_lat;
    P3 = 0.5*(15*cos_lat^2-3);
    P4 = 1/3*(7*cos_lat*P3 - 4*P2);
    P5 = 1/4*(9*cos_lat*P4 - 5*P3);

    A_D = J2*(r_M/r)^2*(U_R*P3 - U_Z*P2);
    A_D = A_D + J3*(r_M/r)^3*(U_R*P4 - U_Z*P3);
    A_D = A_D + J4*(r_M/r)^4*(U_R*P5 - U_Z*P4);

    A_D += 3*J22*(r_M/r)^2*(-5*(x_M^2-y_M^2)/r/r*U_R+2*x_M/r*U_X-2*y_M/r*U_Y)...
                +3/2*C31*(r_M/r)^3*(5*x_M/r*(1-7*cos_lat^2)*U_R+(5*cos_lat^2-1)*U_X+10*x_M*z_M/r/r*U_Z);

    A_D *= mu/(r^2);
    accel = accel + A_D;
  endfunction

  function R = oapiGetRotationMatrixEfficient(t)
	  L_rel = Constants.L_0+2*pi*(t-Constants.t0)/Constants.T_p;
    SLREL = sin(L_rel);
    CLREL = cos(L_rel);
	  Rot3 = [CLREL 0 -SLREL;0 1 0;SLREL 0 CLREL];
	  R_rel = Rot3*Constants.Rot4;
	  phi = Constants.phi_0+2*pi*(t-Constants.t0)/Constants.T_s+(Constants.L_0-L_rel)*cos(Constants.e_rel);
    SPHI = sin(phi);
    CPHI = cos(phi);
	  R_rot = [CPHI 0 -SPHI;0 1 0;SPHI 0 CPHI];
	  Rot = R_rel*R_rot;
    R = Constants.R_ref*Rot;
  endfunction

  function [R, V] = SphericalToCartesian(r, v, lat, lng, gamma, azi)
    R = [cos(lat) * cos(lng) cos(lat) * sin(lng) sin(lat)]' * r;
    M = [cos(lat) * cos(lng) -sin(lng) -sin(lat) * cos(lng); cos(lat) * sin(lng) cos(lng) -sin(lat) * sin(lng); sin(lat) 0 cos(lat)];
	  V = M*([sin(gamma) cos(gamma) * sin(azi) cos(gamma) * cos(azi)]' * v);
  endfunction

  function B = MatrixRH_LH(A)
    %Function to convert a left handed matrix to a right handed matrix
	  B = [A(1,1) A(1,3) A(1,2); A(3,1) A(3,3) A(3,2); A(2,1) A(2,3) A(2,2)];
  endfunction

  function C = rhmul(A,b)
    C = zeros(3,1);
    C(1) = A(1,1)*b(1)+A(1,2)*b(3)+A(1,3)*b(2);
    C(2) = A(3,1)*b(1)+A(3,2)*b(3)+A(3,3)*b(2);
    C(3) = A(2,1)*b(1)+A(2,2)*b(3)+A(2,3)*b(2);
  endfunction

  function  C = rhtmul(A,b)
    C = zeros(3,1);
    C(1) = A(1,1)*b(1)+A(2,1)*b(3)+A(3,1)*b(2);
    C(2) = A(1,3)*b(1)+A(2,3)*b(3)+A(3,3)*b(2);
    C(3) = A(1,2)*b(1)+A(2,2)*b(3)+A(3,2)*b(2);
  endfunction

  function b = unit(a)
    b = a/norm(a);
  endfunction

  function [data] = ConstantsFunction()
    data.t0 = 51544.5;					      %LAN_MJD, MJD of the LAN in the "beginning"
    data.T_p = -6793.468728092782;    %Precession Period
    data.L_0 = 1.71817749;				    %LAN in the "beginning"
	  data.e_rel = 0.026699886264850;  %Obliquity/axial tilt of the earth in radians
	  data.phi_0 = 4.769465382;        %Sidereal Rotational Offset
	  data.T_s = 2360588.15/24/60/60;  %Sidereal Rotational Period
	  e_ref = 7.259562816e-005;   %Precession Obliquity
	  L_ref = 0.4643456618;       %Precession LAN

    data.R_ref = [cos(L_ref) 0 -sin(L_ref);0 1 0;sin(L_ref) 0 cos(L_ref)]*[1 0 0;0 cos(e_ref) -sin(e_ref); 0 sin(e_ref) cos(e_ref)];
    data.Rot4 = [1 0 0;0 cos(data.e_rel) -sin(data.e_rel);0 sin(data.e_rel) cos(data.e_rel)];
  endfunction

  function J_D = RTE_TJUDAT(Y, M, D)
    Y_apo = Y - 1900;
    TMM = [0 31 59 90 120 151 181 212 243 273 304 334];
    Z = floor(Y_apo/4);
    if Z == Y_apo/4
      Z = Z -1;
      for i=3:12
        TMM(i) = TMM(i) + 1;
      endfor
    endif
    J_D = 2415020.5 + 365*Y_apo + Z + TMM(M) + D - 1;
  endfunction

  function MJD = MJDOfNBYEpoch(epoch)
	  %Calculate MJD of Besselian epoch

	  A = 0.0929;
	  B = 8640184.542;
	  W1 = 1.720217954160054e-2;

	  E = epoch;
	  XN = floor((E - 1901) / 4);
	  C = -86400.0*(E - 1900) - 74.164;
	  T = 2.0 * C / (-B - sqrt(B*B - 4.0 * A*C));
	  DE = 36525.0*T - 365.0*(E - 1900) + 0.5 - XN;

	  JD = RTE_TJUDAT(epoch, 1, 0);
	  MJD = JD - 2400000.5 + DE;
  endfunction

  function b=MRz(a)
    ca = cos(a);
    sa = sin(a);
    b=[ca sa 0;-sa ca 0;0 0 1];
  endfunction

  function b=MRx(a)
    ca = cos(a);
    sa = sin(a);
    b=[1 0 0;0 ca sa;0 -sa ca];
  endfunction

  function Rot = J2000EclToBRCSMJD(mjd)
	  t1 = (mjd - 51544.5) / 36525.0;
	  t2 = t1 * t1;
	  t3 = t2 * t1;

	  t1 *= 4.848136811095359e-6;
	  t2 *= 4.848136811095359e-6;
	  t3 *= 4.848136811095359e-6;

	  i = 2004.3109*t1 - 0.42665*t2 - 0.041833*t3;
	  r = 2306.2181*t1 + 0.30188*t2 + 0.017998*t3;
	  L = 2306.2181*t1 + 1.09468*t2 + 0.018203*t3;

	  rot = -r - pi/2;
	  lan = pi/2 - L;
	  inc = i;
	  obl = 0.4090928023;

	  Rot = MRz(rot)*MRx(inc)*MRz(lan)*MRx(-obl);
  endfunction

  function M = J2000EclToBRCS(epoch)
	% Calculate the rotation matrix between J2000 and mean Besselian of epoch coordinate systems
	MJD = MJDOfNBYEpoch(epoch);
	M = J2000EclToBRCSMJD(MJD);
  endfunction

endfunction

