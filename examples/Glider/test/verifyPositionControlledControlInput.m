function verifyPositionControlledControlInput()
  % Tests the URDF implementation of the glider model versus the
  % GliderPlant.
  % Here we test a glider URDF that uses a position controlled control
  % input.

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

  gp_linearized = GliderPlantLinearizedElevAngle();

  %disp('constructing a Planar glider')
  options.floating = true;
  
  
  p = PlanarRigidBodyManipulator('GliderPositionControlled.urdf', options);

  p2 = RigidBodyManipulator('GliderPositionControlled.urdf', options);
  v = p2.constructVisualizer();

  % first check with u = 0 since that should be basically perfect
  % (no loss in the linearization from B*u)

  %    [X Z Pitch El Vx Vz PitDot Velev]

    u0 = 0;%rand(1)-.5;
    pitch = 0;
    phi = 1;
    xvel = 3+4;
    zvel = 0;
    pitchdot = 0;
    x = [0 0   pitch   phi  xvel  zvel    pitchdot]';
    %glider_xdot = gp.dynamics(0,x,u0)
    
    %glider_xdot_linearized = gp_linearized.dynamics(0,x,u0)
    
    glider_xdot = gp.dynamics(0, x, u0)
    
    %urdf_xdot = urdf_dynamics_planar(p,x,u0);
    %valuecheck(urdf_xdot,glider_xdot, 1e-7);
    
    urdf_xdot_3d = urdf_dynamics_3d(p2,x,u0)
    
    
    %glider_xdot - urdf_xdot_3d
    
    valuecheck(urdf_xdot_3d, glider_xdot, 1e-7);
    
  return


  real_dynamics_tol = 1e-2;
  linearized_dynamics_tol = 1e-7;

  disp(['Testing 3D with URDF against real dynamics, tolerance = ' num2str(real_dynamics_tol)]);

  
  for i = 1:80
    
    % for the position controlled control surface, the input to GliderPlant (velocity)
    % is always 0
    u0 = 0; 
    
    pitch = rand(1)-.5;
    phi = rand(1)-.5;
    xvel = rand(1)*3+4;
    zvel = rand(1)*2-1;
    pitchdot = rand(1)-.5;
    x = [0 0   pitch   phi  xvel  zvel    pitchdot]';

    glider_xdot = gp.dynamics(0,x,u0);

    urdf_xdot_3d = urdf_dynamics_3d(p2,x,u0);

    valuecheck(urdf_xdot_3d, glider_xdot, real_dynamics_tol);

    fprintf('.');
  end
  disp(' passed.');
  
% 
%   disp(['Testing 3D with URDF against linearized dynamics, tolerance = ' num2str(linearized_dynamics_tol)]);
% 
%   for i = 1:80
%     u0 = rand(1)-.5;
%     pitch = rand(1)-.5;
%     phi = rand(1)-.5;
%     xvel = rand(1)*3+4;
%     zvel = rand(1)*2-1;
%     pitchdot = rand(1)-.5;
%     x = [0 0   pitch   phi  xvel  zvel    pitchdot]';
% 
%     glider_xdot_lu = gp_lu.dynamics(0,x,u0);
% 
%     urdf_xdot_3d = urdf_dynamics_3d(p2,x,u0);
% 
%     valuecheck(urdf_xdot_3d,glider_xdot_lu, linearized_dynamics_tol);
%     
%     fprintf('.');
%   end
%   disp(' passed.');
%   
%   
%   disp(['Testing planar with URDF against real dynamics, tolerance = ' num2str(real_dynamics_tol)]);
%   for i = 1:80
%     u0 = rand(1)-.5;
%     pitch = rand(1)-.5;
%     phi = rand(1)-.5;
%     xvel = rand(1)*3+4;
%     zvel = rand(1)*2-1;
%     pitchdot = rand(1)-.5;
%     x = [0 0   pitch   phi  xvel  zvel    pitchdot]';
% 
%     glider_xdot = gp.dynamics(0,x,u0);
% 
%     urdf_xdot_2d = urdf_dynamics_planar(p,x,u0);
% 
%     valuecheck(urdf_xdot_2d,glider_xdot, real_dynamics_tol);
%     
%     fprintf('.');
%   end
%   disp(' passed.');
% 
% 
%   disp(['Testing planar with URDF against linearized dynamics, tolerance = ' num2str(linearized_dynamics_tol)]);
%   for i = 1:80
%     u0 = rand(1)-.5;
%     pitch = rand(1)-.5;
%     phi = rand(1)-.5;
%     xvel = rand(1)*3+4;
%     zvel = rand(1)*2-1;
%     pitchdot = rand(1)-.5;
%     x = [0 0   pitch   phi  xvel  zvel    pitchdot]';
% 
%     glider_xdot_lu = gp_lu.dynamics(0,x,u0);
% 
%     urdf_xdot_3d = urdf_dynamics_planar(p,x,u0);
% 
%     valuecheck(urdf_xdot_3d,glider_xdot_lu, linearized_dynamics_tol);
%     
%     fprintf('.');
%   end
%   disp(' passed.');

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
    
    if phidot ~= 0
      error('Non-zero phidot not allowed for position controlled surface.');
    end
    
    x0 = [ x0; y0; z0; roll0; pitch0; yaw0; xdot0; ydot0; zdot0; rolldot0;  pitchdot0; yawdot0 ];
    
    u = -x(4); % the position controlled input is the position of the glider's elevator
    
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
    
    elevdot = 0; % position controlled actuator means always 0 velocity
    
    
    xdot = [xdot; zdot; -pitchdot_3d; -elevdot; xddot; zddot; -pitchddot ];
    
    
  end


%   function xdot = urdf_dynamics_planar(p,x,phidot)
%     % URDF's are defined with positive pitch = down (y-axis rotation)
%     %GliderPlant uses positive pitch = up.
%     %because pitch axes are reversed between models
%     
%     x0 = x(1);
%     z0 = x(2);
%     
%     pitch0 = -x(3);
%     
%     xdot0 = x(5);
%     zdot0 = x(6);
%     
%     pitchdot0 = -x(7);
%     
%     elv0 = -x(4);
%     
%     x0 = [ x0; z0; pitch0; xdot0; zdot0; pitchdot0; elv0 ];
%     
%     u = -phidot;
%     
%     xdot2d = p.dynamics(0,x0,u);
%     
%     xdot = xdot2d(1);
%     zdot = xdot2d(2);
%     
%     pitchdot_2d = xdot2d(3);
%     
%     xddot = xdot2d(4);
%     zddot = xdot2d(5);
%     
%     pitchddot = xdot2d(6);
%     
%     elevdot = xdot2d(7);
%     
%     
%     xdot = [xdot; zdot; -pitchdot_2d; -elevdot; xddot; zddot; -pitchddot ];
%     
%   end

end
