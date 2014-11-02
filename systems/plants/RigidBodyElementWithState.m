classdef RigidBodyElementWithState < RigidBodyElement
  % Some rigid body elements will add to the state vector
  % at least for now, I can only see them adding positions
  
  methods (Abstract=true)
    names = getCoordinateNames(obj);
    n = getNumStates(obj);
    x0 = getInitialState(obj,manip);
    [xdot,df] = dynamics(obj,manip,t,x,u);
  end
  
  properties
    extra_state_num=[];
  end
end