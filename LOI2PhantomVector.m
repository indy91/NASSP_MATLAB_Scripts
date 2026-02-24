function LOI2PhantomVector()

  % INPUTS
  % GMTBASE       = MJD of midnight of launch day, days
  % t_launch      = Launch time (GMT), seconds
  % t_land        = Time (GET) of CSM pass over landing site for lunar landing, seconds
  % lat_LS        = Landing site latitude, degrees
  % lng_LS        = Landing site longitude, degrees
  % alt_LS        = Landing site elevation, nautical miles
  % azi_land      = Approach azimuth, degrees
  % t_circ        = Time (GET) of circular orbit, seconds
  % h_circ        = Circular altitude, nautical miles
  % t_LOI2        = Time (GET) of LOI-2, seconds
  % Epoch         = Year defining the coordinate system (usually the year containing the January 1st closest to launch day)
  % GravityModel  = 0= R2 model, 1= L1 model
  % RelTolVal     = State vector propagation accuracy. Larger value is faster, but less accurate. Reasonable results starting with 1e-5
  % ConvergenceTol= Tolerance for the convergence of the eccentricity. Default value of 1e-5. Larger value = higher tolerance

  % Apollo 11 LOI-2
  GMTBASE       = 40418.0;
  t_launch      = 13*3600 + 32*60;
  t_land        = 102*3600 + 47*60 + 11 - 4.2*60;
  lat_LS        = 0.71388889;
  lng_LS        = 23.70777778;
  alt_LS        = -1.44;
  azi_land      = -91;
  t_circ        = 125*3600;
  h_circ        = 60;
  t_LOI2        = 79*3600;
  Epoch         = 1970;
  GravityModel  = 0;
  RelTolVal     = 1e-6;
  ConvergenceTol= 1e-5;

  CircPhantomVectorCalcs(GMTBASE, t_launch, t_land, lat_LS, lng_LS, alt_LS, azi_land, t_circ, h_circ, t_LOI2, Epoch, GravityModel, RelTolVal, ConvergenceTol);
endfunction
