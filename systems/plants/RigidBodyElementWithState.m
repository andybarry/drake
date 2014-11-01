classdef RigidBodyElementWithState < RigidBodyElement
  % Some rigid body elements will add to the state vector
  % at least for now, I can only see them adding positions
  
  methods (Abstract=true)
    n = getNumPositions(obj);
    names = getCoordinateNames(obj);
    x0 = getInitialPosition(obj,manip);
    [xdot,df] = dynamics(obj,manip,t,x,u);
  end
  
  properties
    position_num=[];
  end
end