function Generate_RTCC_CSM_CG_Table

  % INPUTS
  % 0 = Readable table, 1 = RTCC constants file, 2 = MED input, 3 = RTCC constructor code (for default values)
  OutputFormat = 0;

  % Minimum mass for output table, in pounds
  MinimumMass = 23200.0;
  % Maximum mass for output table, in pounds
  MaximumMass = 65000.0;
  % Number of data points (maximum 40 for RTCC table)
  Num = 20;

  % Mass of empty Command Module, in pounds
  Mass_CM = 12055.0;
  % CG location of empty Command Module (currently hardcoded in NASSP), in inches
  CG_CM = [1041.7 -0.4 5.6];
  % Mass of empty Service Module, in pounds
  Mass_SM = 9453.1;
  % CG location of empty Service Module ("EmptySMCG" variable in mission file), in inches
  CG_SM = [914.5916 -6.6712 12.2940];
  % CM RCS mass, in pounds
  Mass_CM_RCS = 244.71;
  % SM RCS mass, in pounds
  Mass_SM_RCS = 1344.8;
  % CG location of SM RCS propellant (hardcoded in NASSP), in inches
  CG_SM_RCS = [941.8 0 0];

  % INTERNAL

  % CONSTANTS
  % X-component of CG location of SPS propellant tanks, as a function of mass (pounds) in tank, in inches.
  oxid_store_tank_param = [-8.385141e-9 0.0061750118 838.7809363];
  fuel_store_tank_param = [-2.144163e-8 0.0098738581 838.7809363];
  oxid_sump_tank_param = [-2.599892e-9 0.0047770151 839.7803146];
  fuel_sump_tank_param = [-6.63916e-9 0.0076383695 839.7803146];

  % Y and Z components of CG location of SPS propellant tanks, in inches.
  oxid_store_tank_CG_YZ = [14.8 47.8];
  fuel_store_tank_CG_YZ = [-14.8 -47.8];
  oxid_sump_tank_CG_YZ = [48.3 6.6];
  fuel_sump_tank_CG_YZ = [-48.3 -6.6];

  % General constants

  % Conversion from pounds to kgs
  lbs = 0.453597;
  % Ratio of oxidizer to fuel in SPS tanks
  SPS_Propellant_Ratio_Loaded = 1.6;
  % Total oxidizer storage tank capacity, in pounds
  SPS_Oxid_Storage_Tank_Capacity = 11668.6;
  % Total fuel storage tank capacity, in pounds
  SPS_Fuel_Storage_Tank_Capacity = 7315.7;
  % Total oxidizer sump tank capacity, in pounds
  SPS_Oxid_Sump_Tank_Capacity = 14616.5;
  % Total fuel sump tank capacity, in pounds
  SPS_Fuel_Sump_Tank_Capacity = 9163.9;

  % Total sump tank capacity, in pounds
  SPS_Sump_Tanks_Capacity = SPS_Oxid_Sump_Tank_Capacity + SPS_Fuel_Sump_Tank_Capacity;%23068.1;

  % CALCULATIONS
  printf("RTCC CSM CG TABLE GENERATOR 1.0\n\n");
  % if BottomOfTanks == 1
  %   printf("[CSM propellant on bottom of tanks]\n");
  % else
  %   printf("[CSM propellant on top of tanks]\n");
  % endif

  if OutputFormat == 0
    printf(" WEIGHT      X        Y        Z     \n");
    printf("  LB       INCHES   INCHES   INCHES  \n");
  endif

  % Start with minimum weight
  Mass = MinimumMass;
  % Generate step length
  DeltaMass = (MaximumMass - MinimumMass)/(Num - 1);
  for i=1:Num
    % Calculate SPS propellant mass
    SPS_Propellant_Mass = Mass - Mass_CM - Mass_SM - Mass_CM_RCS - Mass_SM_RCS;
    % Error check
    if SPS_Propellant_Mass < 0
      printf("Error: Input minimum mass would result in no SPS propellant by %.2f lbs\n", -SPS_Propellant_Mass);
      return;
    endif

    % Split up in sump vs. storage tank weights
    if SPS_Propellant_Mass > SPS_Sump_Tanks_Capacity
      SPS_Sump_Mass = SPS_Sump_Tanks_Capacity;
      SPS_Storage_Mass = SPS_Propellant_Mass - SPS_Sump_Tanks_Capacity;
    else
      SPS_Sump_Mass = SPS_Propellant_Mass;
      SPS_Storage_Mass = 0;
    endif

    % Calculate oxidizer and fuel components
    oxstorem = SPS_Storage_Mass*SPS_Propellant_Ratio_Loaded/(1+SPS_Propellant_Ratio_Loaded);
    fuelstorem = SPS_Storage_Mass*1/(1+SPS_Propellant_Ratio_Loaded);
    oxsumpm = SPS_Sump_Mass*SPS_Propellant_Ratio_Loaded/(1+SPS_Propellant_Ratio_Loaded);
    fuelsumpm = SPS_Sump_Mass*1/(1+SPS_Propellant_Ratio_Loaded);

    % Calculate SPS tank CG locations
    SPS_Oxid_Storage_CG = GetSPSOxidStorageTankCGLocation(oxstorem);
    SPS_Fuel_Storage_CG = GetSPSFuelStorageTankCGLocation(fuelstorem);
    SPS_Oxid_Sump_CG = GetSPSOxidSumpTankCGLocation(oxsumpm);
    SPS_Fuel_Sump_CG = GetSPSFuelSumpTankCGLocation(fuelsumpm);

    % Calculate CG location
    CG = (CG_CM*(Mass_CM + Mass_CM_RCS) + CG_SM*Mass_SM + SPS_Oxid_Storage_CG*oxstorem + SPS_Oxid_Sump_CG*oxsumpm...
         + SPS_Fuel_Storage_CG*fuelstorem + SPS_Fuel_Sump_CG*fuelsumpm + CG_SM_RCS*Mass_SM_RCS)/Mass;

    if OutputFormat == 0
      printf("%.2f   %.2f    %+.2f    %+.2f\n", Mass,CG(1), CG(2), CG(3));
    elseif OutputFormat == 1
      printf("MHVCCG %d %.2f %f %f %f\n", i-1,Mass,CG(1), CG(2), CG(3));
    elseif OutputFormat == 2
      printf("M11,CSM,%d,%d,%.2f,%f,%f,%f;\n",Num, i, Mass, CG(1), CG(2), CG(3));
    else
      printf("MHVCCG.Weight[%d] = %.2f*0.453597;\nMHVCCG.CG[%d] = _V(%f, %f, %f)*0.0254;\n", i-1,Mass,i-1,CG(1), CG(2), CG(3));
    endif
    % Next weight
    Mass = Mass + DeltaMass;
  endfor


  function CG = GetSPSOxidStorageTankCGLocation(mass)
    CG = GetSPSTankCGLocation(mass, SPS_Oxid_Storage_Tank_Capacity, oxid_store_tank_param, oxid_store_tank_CG_YZ);
  endfunction

  function CG = GetSPSFuelStorageTankCGLocation(mass)
    CG = GetSPSTankCGLocation(mass, SPS_Fuel_Storage_Tank_Capacity, fuel_store_tank_param, fuel_store_tank_CG_YZ);
  endfunction

  function CG = GetSPSOxidSumpTankCGLocation(mass)
    CG = GetSPSTankCGLocation(mass, SPS_Oxid_Sump_Tank_Capacity, oxid_sump_tank_param, oxid_sump_tank_CG_YZ);
  endfunction

  function CG = GetSPSFuelSumpTankCGLocation(mass)
    CG = GetSPSTankCGLocation(mass, SPS_Fuel_Sump_Tank_Capacity, fuel_sump_tank_param, fuel_sump_tank_CG_YZ);
  endfunction

  function CG = GetSPSTankCGLocation(mass, capacity, x_coefs, yz_constants)
    % mass: mass in tank in pounds
    % capacity: total mass in pounds
    % x_coefs:
    % yz_constants:
    cgx = x_coefs(1)*mass^2+x_coefs(2)*mass+x_coefs(3);
    CG = [cgx yz_constants(1) yz_constants(2)];

  endfunction

endfunction

