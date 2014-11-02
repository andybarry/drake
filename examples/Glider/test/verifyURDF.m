function verifyURDF()
% Tests the URDF implementation of the glider model versus the GliderPlant
% that was adapted from Rick Cory's original model.  The wing areas may not
% exactly match the physical glider used--they were initially picked with
% a known working system so the dynamics would match

  % state:  
  %  x(1) - x position
  %  x(2) - z position
  %  x(3) - pitch (theta)
  %  x(4) - elevator (phi)
  %  x(5) - x velocity
  %  x(6) - z velocity
  %  x(7) - pitch velocity (thetadot)
  % input:
  %  u(1) - elevator velocity (phidot)
tmp = addpathTemporary(fullfile(pwd,'..'));

gp = GliderPlant();

%disp('constructing a Planar glider')
options.floating = true;
%p = TimeSteppingRigidBodyManipulator('Glider.URDF',.001, options);
p = PlanarRigidBodyManipulator('GliderBalanced.urdf', options);
%p2 = RigidBodyManipulator('GliderBalanced.urdf', options);
%v = p2.constructVisualizer();
%    [X Z Pitch El Vx Vz PitDot Velev]
% for i = 1:100
%   u0 = rand(1)-.5;
%   pitch = rand(1)-.5;
%   phi = rand(1)-.5;
%   xvel = rand(1)*3+4;
%   zvel = rand(1)*2-1;
%   pitchdot = rand(1)-.5;
  u0 = 2;
  pitch = 0;
  phi = 0;
  xvel = 4;
  zvel = 0;
  pitchdot = 0;
  x = [0 0   pitch   phi  xvel  zvel    pitchdot]';
  glider_xdot = gp.dynamics(0,x,u0)
  urdf_xdot = urdf_dynamics(p,x,u0)
  
  valuecheck(urdf_xdot,glider_xdot);
%end

  function xdot = urdf_dynamics(p,x,phidot)
    % URDF's are defined with positive pitch = down (y-axis rotation)
    %GliderPlant uses positive pitch = up.
    %because pitch axes are reversed between models--also for elevator(input)
    x = [x(1:2);-x(3);x(5:6);-x(7);-x(4)];
    u = -phidot;
    xdot = p.dynamics(0,x,u);
    xdot = [xdot(1:2);-xdot(3);-xdot(7);-xdot(4);xdot(5:6)];
  end

    function u = computeAccel(p, x, vel)
      qdd = 0;%10*(x(8)-vel);
      [H,C,B] = manipulatorDynamics(p,x(1:4),x(5:8),false);
      H11 = H(1:3,1:3);
      H12 = H(1:3,4);
      H21 = H(4,1:3);
      H22 = H(4,4);
      C1 = C(1:3);
      C2 = C(4);
      u = (H22-H21*inv(H11)*H12)*qdd + C2 - H21*inv(H11)*C1;
    end
end
