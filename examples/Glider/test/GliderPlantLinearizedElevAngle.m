classdef GliderPlantLinearizedElevAngle < DrakeSystem
% Defines the dynamics for the perching glider.  Translated from Rick Cory's code in 
% https://svn.csail.mit.edu/russt/robots/perchingGlider/ 
  
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

  properties  % took parameters from Rick's "R1 = no dihedral" model
    S_w = 0.0885; % wing + fus + tail.
    S_e = 0.0147;
    l_w = 0;
    l_e = 0.022;
    l_h = 0.27;
    I = 0.0015;
    m = 0.08;
    g = 9.81; % m/s^2
    rho = 1.204; % kg/m^3
    xd = [0,0,pi/4,0,0,0,0]';% the goal
    phi_lo_limit = -0.9473; % lower elevator angle limit
    phi_up_limit = 0.4463; % upper elevator angle limit
%    talon_b = [0 -.043]'; % talon coordinates, about cg in body coords (meters).
  end
  
  methods
    function obj = GliderPlantLinearizedElevAngle()
      obj = obj@DrakeSystem(7,0,1,7,false,true);
      obj = setDirectFeedthrough(obj,0);
      ulimit = 13; % max servo velocity
      obj = setInputLimits(obj,-ulimit,ulimit);
      obj = setOutputFrame(obj,getStateFrame(obj));
    end
    
    function [xdot,df] = dynamics(obj,t,x,u)
      keyboard
      q=x(1:4); qdot=x(5:7); 
      m=obj.m; g=obj.g; rho=obj.rho; Sw=obj.S_w;
      Se=obj.S_e; I=obj.I; l=obj.l_h; lw=obj.l_w; le=obj.l_e;

      q1 = x(1,:); q2 = x(2,:);  q3 = x(3,:); q4 = x(4,:);
      qdot1 = x(5,:); qdot2 = x(6,:); qdot3 = x(7,:); qdot4=u(1);
      
      % here we linarized about elevator angle (q4 aka x(4) )
      
      if qdot4 ~= 0
        error('This plant does not support non-zero u because it is a test for the position controlled actuators and u controls the velocity of the surface.');
      end

      xwdot = qdot1 - lw*qdot3*sin(q3);
      zwdot=  qdot2 + lw*qdot3*cos(q3);
      alpha_w = q3 - atan2(zwdot,xwdot);
      Fw = rho*Sw*sin(alpha_w)*(zwdot^2+xwdot^2);
      
      diff_amount = 1e-5;
      q4_input = q4;
      
      % now change q4 to do numerical differentiation (linearize about 0)
      q4 = diff_amount;

      xedot = qdot1 + l*qdot3*sin(q3) + le*(qdot3+qdot4)*sin(q3+q4); 
      zedot = qdot2 - l*qdot3*cos(q3) - le*(qdot3+qdot4)*cos(q3+q4);
      alpha_e = q3+q4-atan2(zedot,xedot);
      Fe = rho*Se*sin(alpha_e)*(zedot.^2+xedot.^2);
      
      xdot_high = x;
      xdot_high(5,:)=-(Fw*sin(q3) + Fe*sin(q3+q4))/m;
      xdot_high(6,:)=(Fw*cos(q3) + Fe*cos(q3+q4))/m -g;
      xdot_high(7,:)=(Fw.*lw - Fe.*(l*cos(q4)+le))/I;

      
      q4 = -diff_amount;

      xedot = qdot1 + l*qdot3*sin(q3) + le*(qdot3+qdot4)*sin(q3+q4); 
      zedot = qdot2 - l*qdot3*cos(q3) - le*(qdot3+qdot4)*cos(q3+q4);
      alpha_e = q3+q4-atan2(zedot,xedot);
      Fe= rho*Se*sin(alpha_e)*(zedot.^2+xedot.^2);
      xdot_low = x;
      xdot_low(5,:)=-(Fw*sin(q3) + Fe*sin(q3+q4))/m;
      xdot_low(6,:)=(Fw*cos(q3) + Fe*cos(q3+q4))/m -g;
      xdot_low(7,:)=(Fw.*lw - Fe.*(l*cos(q4)+le))/I;
      
      q4 = 0;

      xedot = qdot1 + l*qdot3*sin(q3) + le*(qdot3+qdot4)*sin(q3+q4); 
      zedot = qdot2 - l*qdot3*cos(q3) - le*(qdot3+qdot4)*cos(q3+q4);
      alpha_e = q3+q4-atan2(zedot,xedot);
      Fe = rho*Se*sin(alpha_e)*(zedot.^2+xedot.^2);
      xdot_affine = x;
      xdot_affine(5,:)=-(Fw*sin(q3) + Fe*sin(q3+q4))/m;
      xdot_affine(6,:)=(Fw*cos(q3) + Fe*cos(q3+q4))/m -g;
      xdot_affine(7,:)=(Fw.*lw - Fe.*(l*cos(q4)+le))/I;

      % take a numerical derivative for elevator force
      xdot=x;
      
      dx5_du = (xdot_high(5,:) - xdot_low(5,:)) / (2*diff_amount);
      xdot(5,:) = dx5_du * q4_input + xdot_affine(5,:);
      
      dx6_du = (xdot_high(6,:) - xdot_low(6,:)) / (2*diff_amount);
      xdot(6,:) = dx6_du * q4_input + xdot_affine(6,:);
      
      dx7_du = (xdot_high(7,:) - xdot_low(7,:)) / (2*diff_amount);
      xdot(7,:) = dx7_du * q4_input + xdot_affine(7,:);
      
      
      xedot = qdot1 + l*qdot3*sin(q3) + le*(qdot3+qdot4)*sin(q3+q4); 
      zedot = qdot2 - l*qdot3*cos(q3) - le*(qdot3+qdot4)*cos(q3+q4);
      alpha_e = q3+q4-atan2(zedot,xedot);
      Fe = rho*Se*sin(alpha_e)*(zedot.^2+xedot.^2);

      xdot_true = x;
      xdot_true(5,:)=-(Fw*sin(q3) + Fe*sin(q3+q4))/m;
      xdot_true(6,:)=(Fw*cos(q3) + Fe*cos(q3+q4))/m -g;
      xdot_true(7,:)=(Fw.*lw - Fe.*(l*cos(q4)+le))/I;
      
      
      q4 = q4_input;

      xdot(4)=u(1);

      %xdot(5,:)=-(Fw*sin(q3) + Fe*sin(q3+q4))/m;

      %xdot(6,:)=(Fw*cos(q3) + Fe*cos(q3+q4))/m -g;

      %xdot(7,:)=(Fw.*lw - Fe.*(l*cos(q4)+le))/I;

      xdot(1,:)=x(5,:);

      xdot(2,:)=x(6,:);

      xdot(3,:)=x(7,:);
      
      % todo: enforce elevator constraints
      % phi_lo_limit < x(4) < phi_up_limit
      
      if (nargout>1)
        [df]= dynamicsGradients(obj,t,x,u,nargout-1);
      end
    end
    
    function [y,dy] = output(obj,t,x,u)
      y = x;
      if (nargout>1)
        dy=[zeros(obj.num_y,1),eye(obj.num_y),zeros(obj.num_y,obj.num_u)];
      end
    end
    
    function x = getInitialState(obj)
      x = [-3.5 .1 0 0 7.0 0 0]';
    end
    
  end  
  
end
