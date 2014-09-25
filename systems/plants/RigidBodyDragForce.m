classdef RigidBodyDragForce < RigidBodyForceElement
  % Element that provides an aerodynamic drag force at a point
  
  properties
    
    coefficient_drag_x; % coefficient of drag in the x (forward) direction
    coefficient_drag_y; % coefficient of drag in the y (left/right) direction
    coefficient_drag_z; % coefficient of drag in the z (up/down) direction
    
    kinframe;
    
    %Air density for 20 degC dry air, at sea level
    rho = 1.204;
    
    area_x; % area of the drag surface in x
    area_y; % area of the drag surface in y
    area_z; % area of the drag surface in z
    
  end
  
  
  methods
    
    function obj = RigidBodyDragForce(frame_id, coefficient_drag_x, coefficient_drag_y, coefficient_drag_z, area_x, area_y, area_z)
      % Provides an aerodynamic drag force at a point
      %
      % Drag force is computed as:
      %  \f$ F_d = \frac{1}{2} \rho v^2 C_d A \f$
      %
      % where:
      % <pre>
      %  F_d: drag force
      %  rho: air density
      %  v: velocity
      %  C_d: coefficient of drag
      %  A: area
      % </pre>
      %  
      %
      % @param frame_id
      % @param coefficient_drag_x coefficient of drag in the x (forward) direction
      % @param coefficient_drag_y coefficient of drag in the y (left/right) direction
      % @param coefficient_drag_z coefficient of drag in the z (up/down) direction
      % @param area_x of the object in drag in the x direction
      % @param area_y of the object in drag in the y direction
      % @param area_z of the object in drag in the z direction
      %
      % @retval obj newly constructed object

      obj.coefficient_drag_x = coefficient_drag_x;
      obj.coefficient_drag_y = coefficient_drag_y;
      obj.coefficient_drag_z = coefficient_drag_z;
      
      obj.area_x = area_x;
      obj.area_y = area_y;
      obj.area_z = area_z;
      
      obj.kinframe = frame_id;
      
      
    end
    
    function [force] = computeSpatialForce(obj,manip,q,qd)
      % Computes the force and derivative of force for the aerodynamic
      % drag on this object.
      %
      % @param manip RigidBodyManipulator we are a part of
      % @param q state vector
      % @param qd time derivative of state vector
      %
      % @retval force force produced by the drag
      % @retval dforce gradient of force produced by the drag
      
      
      frame = getFrame(manip,obj.kinframe);
      kinsol = doKinematics(manip,q,true,true,qd);
      
      % Compute airspeed over this point in xyz
      
      [~, ~, ~, ~, ~, ~, wingvel_world_xyz ] = RigidBodySubWing.computeWingVelocity(obj.kinframe, manip, q, qd, kinsol);
      
      velocity_x = wingvel_world_xyz(1);
      velocity_y = wingvel_world_xyz(2);
      velocity_z = wingvel_world_xyz(3);
      
      force_x = -0.5 * obj.rho * velocity_x * velocity_x * obj.coefficient_drag_x * obj.area_x;
      force_y = -0.5 * obj.rho * velocity_y * velocity_y * obj.coefficient_drag_y * obj.area_y;
      force_z = -0.5 * obj.rho * velocity_z * velocity_z * obj.coefficient_drag_z * obj.area_z;
      
      f = sparse(6, getNumBodies(manip)) * q(1); % q(1) for taylorvar
      
      f(:, frame.body_ind) = [force_x; force_y; force_z; 0; 0; 0];
      
      force = manip.cartesianForceToSpatialForce(kinsol, frame.body_ind, zeros(3,1), f);
      
      
      full(force)
      % TODO
      dforce = 0;
      
      
    end
    
    
  end
  
  methods (Static)
    
    function [model,obj] = parseURDFNode(model,robotnum,node,options)
      % Parse URDF node for drag object.
      
      name = char(node.getAttribute('name'));
      name = regexprep(name, '\.', '_', 'preservecase');

      elNode = node.getElementsByTagName('parent').item(0);
      parent = findLinkInd(model,char(elNode.getAttribute('link')),robotnum);
      
      % get the xyz position of the element
      xyz=zeros(3,1);
      elnode = node.getElementsByTagName('origin').item(0);
      if ~isempty(elnode)
        if elnode.hasAttribute('xyz')
          xyz = reshape(parseParamString(model,robotnum,char(elnode.getAttribute('xyz'))),3,1);
        end
      end
      
      rpy = zeros(3,1);
      [model,this_frame_id] = addFrame(model,RigidBodyFrame(parent,xyz,rpy,[name,'_frame']));
      
      % get the drag coefficients of the element
      cdrag_x = parseParamString(model,robotnum,char(node.getAttribute('coefficient_drag_x')));
      cdrag_y = parseParamString(model,robotnum,char(node.getAttribute('coefficient_drag_y')));
      cdrag_z = parseParamString(model,robotnum,char(node.getAttribute('coefficient_drag_z')));
      
      this_area_x = parseParamString(model,robotnum,char(node.getAttribute('area_x')));
      this_area_y = parseParamString(model,robotnum,char(node.getAttribute('area_y')));
      this_area_z = parseParamString(model,robotnum,char(node.getAttribute('area_z')));
      
      
      obj = RigidBodyDragForce(this_frame_id, cdrag_x, cdrag_y, cdrag_z, this_area_x, this_area_y, this_area_z);
      obj.name = name;
      
    end
    
  end
  
  
end