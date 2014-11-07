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

gp_lu = GliderPlantLinearizedU();

%disp('constructing a Planar glider')
options.floating = true;
%p = TimeSteppingRigidBodyManipulator('Glider.URDF',.001, options);
p = PlanarRigidBodyManipulator('GliderBalanced.urdf', options);

p2 = RigidBodyManipulator('GliderBalanced.urdf', options);
v = p2.constructVisualizer();

% first check with u = 0 since that should be basically perfect
% (no loss in the linearization from B*u)

%    [X Z Pitch El Vx Vz PitDot Velev]

%   u0 = -1;%rand(1)-.5;
%   pitch = .5;
%   phi = 0.5;
%   xvel = 3+4;
%   zvel = 0.2;
%   pitchdot = .5;
%   x = [0 0   pitch   phi  xvel  zvel    pitchdot]';
%   glider_xdot = gp.dynamics(0,x,u0)
%   
%   %glider_lu_xdot = gp_lu.dynamics(0,x,u0)
%   
%   %urdf_xdot = urdf_dynamics_planar(p,x,u0);
%   %valuecheck(urdf_xdot,glider_xdot, 1e-7);
%   
%   urdf_xdot_3d = urdf_dynamics_3d(p2,x,u0)
%   
%   
%   %glider_xdot - urdf_xdot_3d
%   
%   valuecheck(urdf_xdot_3d,glider_xdot, 1e-2);
%   
% return


fprintf('Testing');

for i = 1:80
  u0 = rand(1)-.5;
  pitch = rand(1)-.5;
  phi = rand(1)-.5;
  xvel = rand(1)*3+4;
  zvel = rand(1)*2-1;
  pitchdot = rand(1)-.5;
  x = [0 0   pitch   phi  xvel  zvel    pitchdot]';
  glider_xdot = gp.dynamics(0,x,u0);
  
  %urdf_xdot = urdf_dynamics_planar(p,x,u0);
  %valuecheck(urdf_xdot,glider_xdot, 1e-7);
  
  urdf_xdot_3d = urdf_dynamics_3d(p2,x,u0);
  valuecheck(urdf_xdot_3d,glider_xdot, 1e-2);
  fprintf('.');
end
disp(' passed.');


% now let u float as well, but we have to have decreased tolerance
for i = 1:100
  u0 = rand(1)-.5;
  pitch = rand(1)-.5;
  phi = rand(1)-.5;
  xvel = rand(1)*3+4;
  zvel = rand(1)*2-1;
  pitchdot = rand(1)-.5;
  x = [0 0   pitch   phi  xvel  zvel    pitchdot]';
  glider_xdot = gp.dynamics(0,x,u0);
  
  glider_xdot_lu = gp_lu.dynamics(0, x, u0);
  
  %urdf_xdot = urdf_dynamics_planar(p,x,u0);
  %valuecheck(urdf_xdot,glider_xdot, .5);
  
  urdf_xdot_3d = urdf_dynamics_3d(p2,x,u0);
  valuecheck(urdf_xdot_3d,glider_xdot, .5);
  
  
  valuecheck(urdf_xdot_3d, glider_xdot_lu, 1e-7);
  
  
  
end

  function xdot = urdf_dynamics_3d(p, x, phidot)
    % URDF's are defined with positive pitch = down (y-axis rotation)
    %GliderPlant uses positive pitch = up.
    %because pitch axes are reversed between models -- also for the
    %elevator
    
    x0 = x(1);
    y0 = 0;
    z0 = x(2);
    
    roll0 = 0;
    pitch0 = -x(3);
    yaw0 = 0;
    
    xdot0 = x(5);
    ydot0 = 0;
    zdot0 = x(6);
    
    rolldot0 = 0;
    pitchdot0 = -x(7);
    yawdot0 = 0;
    
    elv0 = -x(4);
    
    x0 = [ x0; y0; z0; roll0; pitch0; yaw0; xdot0; ydot0; zdot0; rolldot0;  pitchdot0; yawdot0; elv0 ];
    
    u = -phidot;
    xdot3d = p.dynamics(0,x0,u);
    
    xdot = xdot3d(1);
    ydot = xdot3d(2);
    zdot = xdot3d(3);
    
    rolldot = xdot3d(4);
    pitchdot_3d = xdot3d(5);
    yawdot = xdot3d(6);
    
    xddot = xdot3d(7);
    yddot = xdot3d(8);
    zddot = xdot3d(9);
    
    rollddot = xdot3d(10);
    pitchddot = xdot3d(11);
    yawddot = xdot3d(12);
    
    elevdot = xdot3d(13);
    
    
    xdot = [xdot; zdot; -pitchdot_3d; -elevdot; xddot; zddot; -pitchddot ];
    
    
  end


  function xdot = urdf_dynamics_planar(p,x,phidot)
    % URDF's are defined with positive pitch = down (y-axis rotation)
    %GliderPlant uses positive pitch = up.
    %because pitch axes are reversed between models
    
    % TODO
    x0 = x(1);
    y0 = 0;
    z0 = x(2);
    
    roll0 = 0;
    pitch0 = -x(3);
    yaw0 = 0;
    
    xdot0 = x(5);
    ydot0 = 0;
    zdot0 = x(6);
    
    rolldot0 = 0;
    pitchdot0 = -x(7);
    yawdot0 = 0;
    
    elv0 = x(4);
    
    x0 = [ x0; z0; pitch0; xdot0; zdot0; pitchdot0; elv0 ];
    
    u = -phidot;
    
    xdot2d = p.dynamics(0,x,u);
    
    xdot = xdot2d(1);
    zdot = xdot2d(2);
    
    pitchdot2d = xdot2d(3);
    
    xddot = xdot2d(4);
    zddot = xdot2d(5);
    
    pitchddot = xdot2d(6);
    
    elevdot = xdot2d(7);
    
    
    xdot = [xdot; zdot; -pitchdot2d; -elevdot; xddot; zddot; -pitchddot ];
    
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
