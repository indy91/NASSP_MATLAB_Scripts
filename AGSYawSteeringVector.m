function AGSYawSteeringVector

  % Script for calculating the AGS yaw steering vector (WBX, WBY, WBZ) using an input yaw angle

  Angle = -20.0; % Yaw left = positive

  % INTERNAL
  ang = (180 - Angle)*pi/180;
  vec = [sin(ang) cos(ang) 0.0]';

  % To AEA
  WBX = DoubleToAEA(vec(1));
  WBY = DoubleToAEA(vec(2));
  WBZ = DoubleToAEA(vec(3));

  % To DEDA
  WBX_DEDA = AEAToDEDA(WBX);
  WBY_DEDA = AEAToDEDA(WBY);
  WBZ_DEDA = AEAToDEDA(WBZ);

  fprintf("\nInput angle: %f degrees\n\n", Angle);
  fprintf("NAME  AEA OCT.  DEDA     AEA VALUE\n");
  fprintf("WBX    %06o  %s   %+e\n", WBX, SignedOctal(WBX_DEDA), vec(1));
  fprintf("WBY    %06o  %s   %+e\n", WBY, SignedOctal(WBY_DEDA), vec(2));
  fprintf("WBZ    %06o  %s   %+e\n", WBZ, SignedOctal(WBZ_DEDA), vec(3));

endfunction

function txt = SignedOctal(val)
  if val >= 0
    txt = sprintf("+%05o", val);
  else
    txt = sprintf("-%05o", abs(val));
  endif
endfunction


function val2 = DoubleToAEA(val)
  val2 = int32(val*2^16);
  if val2 < 0
    val2 = 262144 - abs(val2);
  endif
endfunction

function val2 = AEAToDEDA(val)
  if val >= 131072
    val2 = -((val - 131072)/4);
  else
    val2 = (val/4);
  endif
endfunction

