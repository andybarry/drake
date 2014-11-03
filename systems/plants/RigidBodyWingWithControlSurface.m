classdef RigidBodyWingWithControlSurface < RigidBodyWing & RigidBodyElementWithState
  % Implements functionality similar to RigidBodyWing but with a
  % control surface attached to the wing.
  %
  % URDF parsing is handled by RigidBodyWing.
      
  properties
    control_surface_chord; % chord of the control surface
    control_surface_min_deflection; % min deflection in radians
    control_surface_max_deflection; % max deflection in radians
    
    fCl_control_surface; % Interpolant or function_handle for the control surface, given aoa and u
    fCd_control_surface;
    fCm_control_surface; 
    dfCl_control_surface_du; % Interpolant or function_handle for the control surface, given aoa and u
    dfCd_control_surface_du;
    dfCm_control_surface_du; 
    
    control_surface_increment = 0.01; % resolution of control surface parameterization in radians

    control_surface_velocity_controlled = false; % set to true for the control surface to be velocity controlled instead of position controlled

  end
  
  methods
    
    function obj = RigidBodyWingWithControlSurface(frame_id, profile, chord, span, stall_angle, velocity, control_surface_chord, control_surface_min_deflection, control_surface_max_deflection, control_surface_velocity_controlled)
      % Constructor taking similar arguments to RigidBodyWing except
      % with the addition of control surface parameters
      %
      % @param control_surface_chord length of the control surface attached
      %   to the wing
      % @param control_surface_min_deflection minimum deflection of the
      %   control surface attached to the wing in radians (ex. -0.9)
      % @param control_surface_max_deflection maximum deflection of the
      %   control surface attached to the wing in radians (ex. 0.9)
      % @pararm control_surface_velocity_controlled set to true to make the
      %   control surface velocity controlled instead of position
      %   controlled (adds the control surface position as a state to the
      %   system)
      %   @default false
      
      
      % we need to be able to construct with no arguments per 
      % http://www.mathworks.com/help/matlab/matlab_oop/class-constructor-methods.html#btn2kiy
      
      obj = obj@RigidBodyWing(frame_id, profile, chord, span, stall_angle, velocity);

      if (nargin == 0)
        return;
      end
      
      obj.control_surface_chord = control_surface_chord;
      obj.control_surface_min_deflection = control_surface_min_deflection;
      obj.control_surface_max_deflection = control_surface_max_deflection;
      obj.control_surface_velocity_controlled = control_surface_velocity_controlled;
      
      obj.direct_feedthrough_flag = true;
      
      if control_surface_velocity_controlled
        %obj.joint_limits = [ obj.control_surface_min_deflection; obj.control_surface_max_deflection ]; 
        obj.input_limits = [ -inf; inf ];
      else
        obj.input_limits = [ obj.control_surface_min_deflection; obj.control_surface_max_deflection ];
      end
      
      
    end
    
    function n = getNumStates(obj)
      if (obj.control_surface_velocity_controlled)
        n = 1;
      else
        n = 0;
      end
    end
    
    function names = getCoordinateNames(obj)
      if (obj.control_surface_velocity_controlled)
        names={[obj.name,'_control_surface_angle']};
      else
        names = {};
      end
    end
    
    function q0 = getInitialState(obj)
      % will only be called if we're velocity controlled
      q0 = 0;
    end
    
    function [qdot,df] = dynamics(obj,manip,t,x,u)
      % will only be called if we're velocity controlled
      qdot = u(obj.input_num);
      df = 0*[t;x;u]'; df(1+manip.getNumStates()+obj.input_num) = 1;
    end
    
    function [obj, model] = onCompile(obj, model)
      % Called on compile.  Builds a new airfoil if needed
      % and also builds a new control surface parameterization if needed
      
      % call parent
      [obj, model] = onCompile@RigidBodyWing(obj, model);

      control_surface_area = obj.span .* obj.control_surface_chord;
      [fCl, fCd, fCm, dfCl, dfCd, dfCm ] = RigidBodyWing.flatplate(obj.rho,control_surface_area,obj.chord);

      obj.fCl_control_surface = @(aoa,u) fCl(aoa+u);
      obj.fCd_control_surface = @(aoa,u) fCd(aoa+u);
      obj.fCm_control_surface = @(aoa,u) fCm(aoa+u);

      obj.dfCl_control_surface_du = @(aoa,u) dfCl(aoa+u); % this is with respect to u, not aoa because u directly controls aoa
      obj.dfCd_control_surface_du = @(aoa,u) dfCd(aoa+u);
      obj.dfCm_control_surface_du = @(aoa,u) dfCm(aoa+u);
      
      % override fCM to accommodate shift in x-axis to aerodynamic center of
      % control surface (note: this is flat plate specific)
      % (lift and drag are uneffected)
      r = obj.chord / 2 + obj.control_surface_chord / 2; % TODO: check if this is correct
      obj.fCm_control_surface = @(aoa,u) obj.rho * r * sin(aoa+u) .* cos(aoa+u) * control_surface_area;
      obj.dfCm_control_surface_du = @(aoa,u) obj.rho * r * (cos(aoa+u).^2 - sin(aoa+u).^2) * control_surface_area;
    end
    
    function [force, B_force, dforce, dB_force] = computeSpatialForce(obj,manip,q,qd,xx)
      % Computes the forces from the wing including the control surface.
      % Returns the force from the wing along with the B matrix which the
      % matrix for a linearized input for the control surface.
      %
      % @param manip manipulator we are a part of
      % @param q state vector
      % @param qd q-dot, time derivative of the state vector
      %
      %
      % @retval force force from the wing that is independant of the
      %   control surface
      %
      % @retval B_force B matrix containing the linearized component of the
      %   force from the input (from the control surface's deflection)
      
      % first, call the parent class's computeSpatialForce to get the
      % u-invariant parts from the part of the wing that is not a control
      % surface

      [force, dforce] = computeSpatialForce@RigidBodyWing(obj, manip, q, qd); 
      
      kinsol = doKinematics(manip,q,true,true,qd);
      
      
      % get the coefficients for this point in state space
      
      wingvel_struct = RigidBodyWing.computeWingVelocity(obj.kinframe, manip, q, qd, kinsol);
      wingvel_rel = RigidBodyWing.computeWingVelocityRelative(obj.kinframe, manip, kinsol, wingvel_struct);
      
      wingvel_world_xz = wingvel_struct.wingvel_world_xz;
      wingYunit = wingvel_struct.wingYunit;
      
      airspeed = norm(wingvel_world_xz);
      
      aoa = -atan2(wingvel_rel(3),wingvel_rel(1));
      
      % lift is defined as the force perpendicular to the direction of
      % airflow, so the lift axis in the body frame is the axis
      % perpendicular to wingvel_world_xz
      % We can get this vector by rotating wingvel_world_xz by 90 deg:
      % cross(wingvel_world_xz, wingYunit)
      
      lift_axis_in_world_frame = cross(wingvel_world_xz, wingYunit);
      lift_axis_in_world_frame = lift_axis_in_world_frame / norm(lift_axis_in_world_frame);
      
      % drag axis is the opposite of the x axis of the wing velocity
      drag_axis_in_world_frame = -wingvel_world_xz / norm(wingvel_world_xz);
 
      
      % compute the u-invariant contribution from the control surface
      %
      % (if it is velocity controlled, that means this is the part coming
      % from the acutator at whatever aoa that state says it's at.
      %
      % if it's position controlled, this is the contribution of a flat
      % plate at u = 0.)

      if obj.control_surface_velocity_controlled
        control_surface_angle = xx(obj.extra_state_num);      
      else
        control_surface_angle = 0;
      end

      fCl_surface = obj.fCl_control_surface(aoa, control_surface_angle);
      fCd_surface = obj.fCd_control_surface(aoa, control_surface_angle);
      fCm_surface = obj.fCm_control_surface(aoa, control_surface_angle);


      lift_force = fCl_surface * airspeed * airspeed;
      drag_force = fCd_surface * airspeed * airspeed;

      torque_moment =  fCm_surface * airspeed * airspeed;

      lift_vector = lift_force * lift_axis_in_world_frame;
      drag_vector = drag_force * drag_axis_in_world_frame;

      % moment torque is around the wing's Y axis
      % build two forces (a couple) around the Y axis to express it
      % these forces should be on the X (forward) axis of the wing

      moment_location_in_body_frame1 = [1; 0; 0];
      moment_location_in_body_frame2 = [-1; 0; 0];

      moment_direction_in_body_frame1 = [0; 0; 1];
      moment_direction_in_body_frame2 = [0; 0; -1];

      [moment_location_in_world_frame1, J1] = forwardKin(manip, kinsol, obj.kinframe, moment_location_in_body_frame1);
      moment_axis_in_world_frame1 = forwardKin(manip, kinsol, obj.kinframe, moment_direction_in_body_frame1);
      moment_axis_in_world_frame1 = moment_axis_in_world_frame1 - moment_location_in_world_frame1;

      [moment_location_in_world_frame2, J2] = forwardKin(manip, kinsol, obj.kinframe, moment_location_in_body_frame2);
      moment_axis_in_world_frame2 = forwardKin(manip, kinsol, obj.kinframe, moment_direction_in_body_frame2);
      moment_axis_in_world_frame2 = moment_axis_in_world_frame2 - moment_location_in_world_frame2;

      moment_vector1 = torque_moment * moment_axis_in_world_frame1;
      moment_vector2 = torque_moment * moment_axis_in_world_frame2;
      
      x_unit = [1; 0; 0];
      y_unit = [0; 1; 0];
      z_unit = [0; 0; 1];

      force_x = dot(lift_vector, x_unit) + dot(drag_vector, x_unit) + dot(moment_vector1, x_unit) + dot(moment_vector2, x_unit);
      force_y = dot(lift_vector, y_unit) + dot(drag_vector, y_unit) + dot(moment_vector1, y_unit) + dot(moment_vector2, x_unit);
      force_z = dot(lift_vector, z_unit) + dot(drag_vector, z_unit) + dot(moment_vector1, z_unit) + dot(moment_vector2, x_unit);

      f_world_frame = [force_x; force_y; force_z];

      frame = getFrame(manip,obj.kinframe);
      f = manip.cartesianForceToSpatialForce(kinsol, frame.body_ind, zeros(3,1), f_world_frame);

      force(:, frame.body_ind) = force(:, frame.body_ind) + f; % add the forces with the rest of the wing


      if obj.control_surface_velocity_controlled
        warning('dforce does not include control surface stuff yet');
        
        
        
        % now compute the forces that depend on u
        
        r = obj.chord / 2 + obj.control_surface_chord / 2;
        v_air_x = wingvel_rel(1);
        v_air_z = wingvel_rel(3);
        
        % linearize lift about u = 0, where u is a velocity input
        
        theta = control_surface_angle;
        rho = obj.rho;
        S = obj.control_surface_chord * obj.span;
        cs_vel = 0; % linearize about control surface velocity = 0
        
        flatplate_dlift_du = (rho * r * S * (-3 * r^2 * cs_vel^2 * v_air_x * cos(3*theta) + 2*v_air_z * ...
          (3*r^2*cs_vel^2 + 2*v_air_z^2 + 3*r^2 * cs_vel^2 * cos(2*theta)) * sin(theta) ...
          + cos(theta) * (3 * r^2 * cs_vel^2 * v_air_x + 4 * v_air_x^3 + ...
          4 * r * cs_vel * (r^2 * cs_vel^2 + 3 * (v_air_x^2 + v_air_z^2)) * ...
          sin(theta))))/(4 * ((v_air_z + r * cs_vel * cos(theta))^2 + ...
          (v_air_x + r * cs_vel * sin(theta))^2)^(3/2));

        
        flatplate_ddrag_du = 2 * rho * S * r * sin(theta) * (r * cs_vel * sin(theta) + v_air_z);
        
        df_lift = flatplate_dlift_du;
        df_drag = flatplate_ddrag_du;
        dtorque_moment = 0; % todo?


      else

        % linearize about u = 0 where u is the position input
        u=0;
        % Cl = Cl_linear*u + Cl_affine  (but the affine terms are already
        % applied from above)
        Cl_linear = obj.dfCl_control_surface_du(aoa,u);
        Cd_linear = obj.dfCd_control_surface_du(aoa,u);
        Cm_linear = obj.dfCm_control_surface_du(aoa,u);

        % debug data to create plots

        %{
        % begin_debug
          control_surface_range = obj.getControlSurfaceRange();
          aoa_range = repmat(aoa, 1, length(control_surface_range));
          Cl = obj.fCl_control_surface(aoa_range, control_surface_range);
          Cd = obj.fCd_control_surface(aoa_range, control_surface_range);
          Cm = obj.fCm_control_surface(aoa_range, control_surface_range);

          figure(1)
          clf
          plot(rad2deg(control_surface_range), Cl)
          hold on
          plot(rad2deg(control_surface_range), Cl_linear * control_surface_range, 'r');
          xlabel('Control surface deflection (deg)');
          ylabel('Coefficient of lift');
          title(['aoa = ' num2str(rad2deg(aoa))]);

          figure(2)
          clf
          plot(rad2deg(control_surface_range), Cd)
          hold on
          plot(rad2deg(control_surface_range), Cd_linear * control_surface_range, 'g');
          xlabel('Control surface deflection (deg)');
          ylabel('Coefficient of drag');
          title(['aoa = ' num2str(rad2deg(aoa))]);

          figure(3)
          clf
          plot(rad2deg(control_surface_range), Cm);
          hold on
          plot(rad2deg(control_surface_range), Cm_linear * control_surface_range, 'k');
          xlabel('Control surface deflection (deg)');
          ylabel('Moment coefficient');
          title(['aoa = ' num2str(rad2deg(aoa))]);
        % end_debug
        %}

        df_lift = Cl_linear * airspeed*airspeed;
        df_drag = Cd_linear * airspeed*airspeed;
        dtorque_moment = Cm_linear * airspeed * airspeed;
        
      end
      
      % initalize B
      B_force = manip.B*0*q(1);
      

      % position of origin
      [~, J] = forwardKin(manip, kinsol, obj.kinframe, zeros(3,1));

      
      B_lift = df_lift * J' * lift_axis_in_world_frame;

      B_drag = df_drag * J' * drag_axis_in_world_frame;
      
      
      % use two forces in opposite directions one meter away to create a
      % torque (a couple)
      
      moment_location_in_body_frame1 = [1; 0; 0];
      moment_location_in_body_frame2 = [-1; 0; 0];
      
      moment_direction_in_body_frame1 = [0; 0; 1];
      moment_direction_in_body_frame2 = [0; 0; -1];
      
      [moment_location_in_world_frame1, J1] = forwardKin(manip, kinsol, obj.kinframe, moment_location_in_body_frame1);
      moment_axis_in_world_frame1 = forwardKin(manip, kinsol, obj.kinframe, moment_direction_in_body_frame1);
      moment_axis_in_world_frame1 = moment_axis_in_world_frame1 - moment_location_in_world_frame1;
      
      [moment_location_in_world_frame2, J2] = forwardKin(manip, kinsol, obj.kinframe, moment_location_in_body_frame2);
      moment_axis_in_world_frame2 = forwardKin(manip, kinsol, obj.kinframe, moment_direction_in_body_frame2);
      moment_axis_in_world_frame2 = moment_axis_in_world_frame2 - moment_location_in_world_frame2;
      
      B_moment = dtorque_moment * 0.5 * J1' * moment_axis_in_world_frame1 ...
       + dtorque_moment * 0.5 * J2' * moment_axis_in_world_frame2;
      
      B_force(:, obj.input_num) = B_lift + B_drag + B_moment; 
      
      dB_force = 0; %todo
    end
    
    function control_surface_range = getControlSurfaceRange(obj)
      % Returns a range of values between the minimum and maximum
      % deflection of the control surface
      %
      %
      % @retval control_surface_range the range as an array
      
      control_surface_range = obj.control_surface_min_deflection : obj.control_surface_increment : obj.control_surface_max_deflection;
    end
    
    function model = addWingVisualShapeToBody(obj, model, body)
      % Adds a visual shape of the wing to the model on the body given for
      % drawing the wing in a visualizer.
      %
      % @param model manipulator the wing is part of
      % @param body body to add the visual shape to
      %
      % @retval model updated model

      % call the parent wing's drawing system to add the main wing
      model = addWingVisualShapeToBody@RigidBodyWing(obj, model, body);
      
      
      % add another box for the control surface
      
      control_surface_height = 0.01;
      
      box_size = [ obj.control_surface_chord, obj.span, control_surface_height ];
      
      % get the xyz and rpy of the control surface
      origin = [-obj.chord/2 - obj.control_surface_chord/2; 0; 0];

      T = model.getFrame(obj.kinframe).T;
      R = T(1:3,1:3);
      
      pts = [origin; 1];
      
      xyz_rpy = [T(1:3,:)*pts; repmat(rotmat2rpy(R),1, 1)];
      
      xyz = xyz_rpy(1:3);
      
      rpy = xyz_rpy(4:6);
      
      shape = RigidBodyBox(box_size, xyz, rpy);
      
      shape = shape.setColor([1 .949 .211]);
      shape.name = [ obj.name '_urdf_shape' ];
      
      model = model.addVisualShapeToBody(body, shape);

    end
    
    function model = updateVisualShapes(obj, model)
      % Update visual shapes.  This should be called after a parameter
      % update that may change the automatic drawing of shapes (ie the area
      % of the drag force is changed).
      %
      % @param model RigidBodyManipulator this is a part of
      %
      % @retval model updated model
      
      if (obj.visual_geometry)
        
        % call parent class's version
        model = updateVisualShapes@RigidBodyWing(obj, model);
        
        % remove any existing shapes for this body
        model = model.removeShapeFromBody(obj.parent_id, [ obj.name '_urdf_shape' ]);
        
        % add new shapes
        model = addWingVisualShapeToBody(obj, model, obj.parent_id);
        
        
      end
      
    end
    
  end
  
  methods (Static)
    
    function [model, obj] = parseURDFNode(model,robotnum,node,options)
      % Parses URDF node for wing_with_control_surface
      %
      % @param model model we are buidling
      % @param robotnum
      % @param node URDF node to parse
      % @param options none right now
      %
      % @retval model updated model
      % @reval obj newly created RigidBodyWingWithControlSurface object
      
      
      name = char(node.getAttribute('name'));
      name = regexprep(name, '\.', '_', 'preservecase');

      elNode = node.getElementsByTagName('parent').item(0);
      parent = findLinkInd(model,char(elNode.getAttribute('link')),robotnum);

      xyz=zeros(3,1); rpy=zeros(3,1);
      elnode = node.getElementsByTagName('origin').item(0);
      if ~isempty(elnode)
        if elnode.hasAttribute('xyz')
          xyz = reshape(parseParamString(model,robotnum,char(elnode.getAttribute('xyz'))),3,1);
        end
        if elnode.hasAttribute('rpy')
          rpy = reshape(parseParamString(model,robotnum,char(elnode.getAttribute('rpy'))),3,1);
        end
      end
      
      
      
      wing_profile = char(node.getAttribute('profile'));

      if ~strcmpi(wing_profile, 'flat plate')
          error('wings with control surfaces that are not flat plates is unimplemneted at the current time.');
      end
      
      

      wing_chord = parseParamString(model,robotnum,char(node.getAttribute('chord')));
      wing_span = parseParamString(model,robotnum,char(node.getAttribute('span')));
      wing_stall_angle = parseParamString(model,robotnum,char(node.getAttribute('stall_angle')));
      nominal_speed = parseParamString(model,robotnum,char(node.getAttribute('nominal_speed')));
      control_surface_chord_urdf = parseParamString(model,robotnum,char(node.getAttribute('control_surface_chord')));
      control_surface_min_deflection_urdf = parseParamString(model,robotnum,char(node.getAttribute('control_surface_min_deflection')));
      control_surface_max_deflection_urdf = parseParamString(model,robotnum,char(node.getAttribute('control_surface_max_deflection')));
      control_surface_velocity_controlled_urdf = parseParamString(model,robotnum,char(node.getAttribute('control_surface_velocity_controlled')));
      visual_geometry_urdf = parseParamString(model,robotnum,char(node.getAttribute('visual_geometry')));
      
      if isempty(visual_geometry_urdf)
        visual_geometry_urdf = true;
      end
      
      [model,this_frame_id] = addFrame(model,RigidBodyFrame(parent,xyz,rpy,[name,'_frame']));
      
      obj = RigidBodyWingWithControlSurface(this_frame_id, wing_profile, wing_chord, ...
        wing_span, wing_stall_angle, nominal_speed, control_surface_chord_urdf, ...
        control_surface_min_deflection_urdf, control_surface_max_deflection_urdf, ...
        control_surface_velocity_controlled_urdf);
      
      obj.name = name;
      obj.visual_geometry = visual_geometry_urdf;
      obj.parent_id = parent;
      
      
      % bind parameters before drawing so we know what to draw
      obj = obj.bindParams(model, double(model.getParams()));
      bound_frame = model.getFrame(this_frame_id).bindParams(model, double(model.getParams()));
      model = model.setFrame(this_frame_id, bound_frame);

      if (visual_geometry_urdf)
        % add visual shapes for automatic drawing of the wing
       % model = obj.addWingVisualShapeToBody(model, parent);
      end

    end
  end
  
  
end
