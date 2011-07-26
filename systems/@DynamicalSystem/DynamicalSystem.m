classdef DynamicalSystem
% An interface class for a state-space dynamical system 
%  with a single (vector) input u, a single (vector) output y, and a single (vector) state x composed of a combination of continuous time and discrete time variables
                                                   
  methods (Abstract = true)
    % Returns the scalar number of continuous state variables, aka xc
    n = getNumContStates(obj);  

    % Returns the scalar number of discrete state variables, aka xd
    n = getNumDiscStates(obj);

    % Returns the scalar length of the single input vector, aka u
    n = getNumInputs(obj);
    
    % Returns the scalar length of the single output vector, aka y
    n = getNumOutputs(obj);
    
    % Returns a simulink compatible sample time for the system.  
    % See http://www.mathworks.com/help/toolbox/simulink/ug/brrdmmw-5.html
    % for a full description
    ts = getSampleTime(obj);

    % Returns a handle (the string) of a simulink model which implements the system
    mdl = getModel(obj);

    % Returns a (possibly random) state vector which is a feasible initial
    % condition for the system.  
    % @retval x0 initial state vector, containing discrete and continuous
    % states
    x0 = getInitialState(obj);

    % Implements the differential equation governing the dynamics of the
    % continuous state variables.  
    % @param t time (scalar)
    % @param x state vector, containing discrete and continuous states
    % @param u input vector
    % @retval xcdot derivative vector of ONLY the continuous states
    xcdot = dynamics(obj,t,x,u);
    
    % Implements the difference equation governing the dynamics of the
    % discrete state variables
    % @param t time (scalar)
    % @param x state vector, containing discrete and continuous states
    % @param u input vector
    % @retval xdn new value (e.g. xd[n+1])for ONLY the discrete states in the state vector
    xdn = update(obj,t,x,u);
    
    % Implements the output function
    % @param t time (scalar)
    % @param x state vector, containing discrete and continuous states
    % @param u input vector
    y = output(obj,t,x,u);
  end
  
  % construction methods
  methods
    function newsys = cascade(sys1,sys2,bsys1out)
      % Creates a new system with the output of system 1 connected to the input of system 2. 
      %
      % @param sys1 the first DynamicalSystem  
      % @param sys2 the second DynamicalSystem
      % @param bsys1out defines the output of system1 to be the output of the 
      %   new system.  Useful, for instance, when cascading a plant with a 
      %   visualizer (which has no outputs) so that you can still access the
      %   output of the plant.  @default false
      %   
      %
      % @retval newsys the new DynamicalSystem

      if (nargin<3) bsys1out=false; end
      mdl = ['Cascade_',datestr(now,'MMSSFFF')];  % use the class name + uid as the model name
      new_system(mdl,'Model');
      set_param(mdl,'SolverPrmCheckMsg','none');  % disables warning for automatic selection of default timestep
      
      load_system('simulink3');
      add_block('simulink3/Subsystems/Subsystem',[mdl,'/system1']);
      Simulink.SubSystem.deleteContents([mdl,'/system1']);
      Simulink.BlockDiagram.copyContentsToSubSystem(sys1.getModel(),[mdl,'/system1']);
      add_block('simulink3/Subsystems/Subsystem',[mdl,'/system2']);
      Simulink.SubSystem.deleteContents([mdl,'/system2']);
      Simulink.BlockDiagram.copyContentsToSubSystem(sys2.getModel(),[mdl,'/system2']);
      add_line(mdl,'system1/1','system2/1');

      if (getNumInputs(sys1)>0)
        add_block('simulink3/Sources/In1',[mdl,'/in']);
        add_line(mdl,'in/1','system1/1');
      end
      
      if (bsys1out)
        if (getNumOutputs(sys1)>0)
          add_block('simulink3/Sinks/Out1',[mdl,'/out']);
          add_line(mdl,'system1/1','out/1');
        end
      else
        if (getNumOutputs(sys2)>0)
          add_block('simulink3/Sinks/Out1',[mdl,'/out']);
          add_line(mdl,'system2/1','out/1');
        end
      end
      
      newsys = SimulinkModel(mdl);
      newsys.time_invariant_flag = sys1.time_invariant_flag && sys2.time_invariant_flag;
      newsys.simulink_params = catstruct(sys1.simulink_params,sys2.simulink_params);
    end
    function newsys = feedback(sys1,sys2)
      % Creates a new systems with sys1 and sys2 in a feedback interconnect 
      % (with sys1 on the forward path, and sys2 on the return path).  
      % the output of sys1 will be the output of the new system.
      if (getNumOutputs(sys1)<1 || getNumOutputs(sys2)<1 || getNumOutputs(sys1)~=getNumInputs(sys2) || getNumOutputs(sys2)~=getNumInputs(sys1))
        error('these systems can''t be combined in feedback because they don''t have a compatible number of inputs / outputs');
      end
      
      mdl = ['Feedback_',datestr(now,'MMSSFFF')];  % use the class name + uid as the model name
      new_system(mdl,'Model');
      set_param(mdl,'SolverPrmCheckMsg','none');  % disables warning for automatic selection of default timestep
      
      load_system('simulink3');
      add_block('simulink3/Subsystems/Subsystem',[mdl,'/system1']);
      Simulink.SubSystem.deleteContents([mdl,'/system1']);
      Simulink.BlockDiagram.copyContentsToSubSystem(sys1.getModel(),[mdl,'/system1']);
      add_block('simulink3/Subsystems/Subsystem',[mdl,'/system2']);
      Simulink.SubSystem.deleteContents([mdl,'/system2']);
      Simulink.BlockDiagram.copyContentsToSubSystem(sys2.getModel(),[mdl,'/system2']);

      add_line(mdl,'system1/1','system2/1');
      add_line(mdl,'system2/1','system1/1');

      add_block('simulink3/Sinks/Out1',[mdl,'/out']);
      add_line(mdl,'system1/1','out/1');
      
      newsys = SimulinkModel(mdl);      
    end
    function newsys = sampledData(sys,tsin,tsout)
      % Creates a new system which is a sampled data (e.g. discretized in time) version of the original system.  
      % This is accomplished by adding rate transition blocks to the inputs and outputs.
      % @param tsin the sample times of the input
      % @param tsout the sample times of the output
      % @retval newsys the new, sampled-data system
      % See getSampleTime for more details about sample times.
      if (nargin<3) tsout = tsin; end
      typecheck(tsin,'double'); 
      sizecheck(tsin,[1 1]);
      typecheck(tsout,'double');
      sizecheck(tsout,[1 1]);
      
      mdl = ['SampledData_',datestr(now,'MMSSFFF')];  % use the class name + uid as the model name
      new_system(mdl,'Model');
      set_param(mdl,'SolverPrmCheckMsg','none');  % disables warning for automatic selection of default timestep
      
      load_system('simulink3');
      add_block('simulink3/Subsystems/Subsystem',[mdl,'/system']);
      Simulink.SubSystem.deleteContents([mdl,'/system']);
      Simulink.BlockDiagram.copyContentsToSubSystem(sys.getModel(),[mdl,'/system']);

      if (getNumInputs(sys)>0)
        add_block('simulink3/Sources/In1',[mdl,'/in']);
        add_block('simulink/Signal Attributes/Rate Transition',[mdl,'/zoh1'],'OutPortSampleTime',num2str(tsin));
        add_line(mdl,'in/1','zoh1/1');
        add_line(mdl,'zoh1/1','system/1');
      end
      
      if (getNumOutputs(sys)>0)
        add_block('simulink3/Sinks/Out1',[mdl,'/out']);
        add_block('simulink/Signal Attributes/Rate Transition',[mdl,'/zoh2'],'OutPortSampleTime',num2str(tsout));
        add_line(mdl,'system/1','zoh2/1');
        add_line(mdl,'zoh2/1','out/1');
      end
      
      newsys = SimulinkModel(mdl);
      newsys.time_invariant_flag = sys.time_invariant_flag;
      newsys.simulink_params = sys.simulink_params;  
    end
  end
  
  % utility methods
  methods
    function n = getNumStates(obj)
      % Returns the total number of states (discrete + continuous) in the system
      n = getNumDiscStates(obj) + getNumContStates(obj);
    end
    function ts = getInputSampleTimes(obj)
      % Returns the sample time of the input
      % See getSampleTime for more details about sample times.
      mdl = getModel(obj);
      inport = find_system(mdl,'SearchDepth',1,'BlockType','Inport');
      ts = [];
      for i=1:length(inport)
        ts = [ts,Simulink.Block.getSampleTimes(inport{i}).Value'];
      end
    end
    function ts = getOutputSampleTimes(obj)
      % Returns the sample time of the output
      % See getSampleTime for more details about sample times.
      mdl = getModel(obj);
      outport = find_system(mdl,'SearchDepth',1,'BlockType','Outport');
      ts = [];
      for i=1:length(outport)
        ts = [ts,Simulink.Block.getSampleTimes(outport{i}).Value'];
      end
    end
    
    function tf = isDT(obj)
      % Returns true if the system has only one sample time [a b], with b>0
      ts = getSampleTime(obj);
      tf = (size(ts,2)==1 && ts(1)>0); % only one sample time and not continuous
    end
    function tf = isCT(obj)
      % Returns true if the system has only one sample time [a b], with b==0
      ts = getSampleTime(obj);
      tf = (size(ts,2)==1 && ts(1)==0); % only one sample time and continuous
    end
    
    function tf = isTI(obj)  
      % Returns true if the system is time-invariant
      tf = obj.time_invariant_flag;
    end
    function obj = setTIFlag(obj,bval)
      % Sets the time invariant flag.
      % Set TI=true if you know the system is time-invariant.  It simplifies many things.
      obj.time_invariant_flag = bval;
    end

    
    function u = wrapInput(obj,u)
      % Wraps all input signals with input_angle_flag==true to be inside [-pi,pi]
      i=find(obj.input_angle_flag);
      u(i) = mod(u(i)+pi,2*pi)-pi;
    end
    function x = wrapState(obj,x)
      % Wraps all state variables with state_angle_flag==true to be inside [-pi,pi]
      i=find(obj.state_angle_flag);
      x(i) = mod(x(i)-pi,2*pi)-pi;
    end
    function y = wrapOutput(obj,y)
      % Wraps all output signals with output_angle_flag==true to be inside [-pi,pi]
      i=find(obj.output_angle_flag);
      y(i) = mod(y(i)-pi,2*pi)-pi;
    end
    function obj = setAngleFlags(obj,in_angle_flag,state_angle_flag,out_angle_flag)
      % Set the input,state, and output angle flags
      % @param in_angle_flag a vector of length(getNumInputs), where
      % in_angle_flag(i)==true means that input i is an angle and should be
      % wrapped around 2*pi when appropriate
      % @param state_angle_flag the corresponding vector of length(getNumStates) 
      % @param out_angle_flag the corresponding vector of length(getNumOutputs) 
      % @retval obj the DynamicalSystem object with updated angle flags
      if (~isempty(in_angle_flag))
        typecheck(in_angle_flag,{'logical','double'});
        sizecheck(in_angle_flag,obj.getNumInputs());
        obj.input_angle_flag = in_angle_flag(:);
      end
      if (~isempty(state_angle_flag))
        typecheck(state_angle_flag,{'logical','double'});
        sizecheck(state_angle_flag,obj.getNumStates());
        obj.state_angle_flag = state_angle_flag(:);
      end
      if (~isempty(out_angle_flag))
        typecheck(out_angle_flag,{'logical','double'});
        sizecheck(out_angle_flag,obj.getNumOutputs());
        obj.output_angle_flag = out_angle_flag(:);
      end
    end
    
    function xs = stateVectorToStructure(obj,xv)
      % Converts the vector state xv to the structure xs for simulink state
      if (isempty(obj.structured_x))
        obj.structured_x = Simulink.BlockDiagram.getInitialState(getModel(obj));
      end
      xs = obj.structured_x;
      if (length(xs.signals)>1)
        l = {xs.signals(:).label};
        ind = [find(strcmp(l,'DSTATE')), find(strcmp(l,'CSTATE'))];
        c = 1;
        for i=ind
          d = xs.signals(i).dimensions;
          xs.signals(i).values = xv(c+[1:d]-1);
          c = c+d;
        end
      else
        xs.signals.values = xv;
      end
    end
    
    function xv = stateStructureToVector(obj,xs)
      % Converts the simulink state structure representation back to the vector state
      if (length(xs.signals)>1)
        l = {xs.signals(:).label};
        ind = [find(strcmp(l,'DSTATE')), find(strcmp(l,'CSTATE'))];
        c = 1;
        for i=ind
          d = xs.signals(i).dimensions;
          xv(c+[1:d]-1) = xs.signals(i).values;
          c = c+d;
        end
      else
        xv = xs.signals.values;
      end
      xv = xv';
    end
    
    function obj = setSimulinkParam(obj,varargin)
      % Sets parameters of the simulink model
      % Syntax
      % setSimulinkParam(obj,param_name,param_value[,param_name,param_value,...])
      % Se
      % @param param_name the string parameter name
      % @param param_value the *string* representing the parameter value
      % @retval obj the DynamicalSystem object with the updated parameters
      % 
      % Try 'doc set_param' in matlab to see the simulink model and block
      % parameters.
      % Note that this supports many of the same parameter values that you
      % would send to odeset.
      if (mod(length(varargin),2)) error('invalid input'); end        
      i=1;
      while(i<length(varargin))
        if (length(varargin{i+1})==0) 
          obj.simulink_params = rmfield(obj.simulink_params,varargin{i});
        else
          if (~ischar(varargin{i+1})) error('simulink params should all be strings (unfortunately)'); end
          obj.simulink_params = setfield(obj.simulink_params,varargin{i},varargin{i+1});
          set_param(obj.getModel(),varargin{i},varargin{i+1});
        end
        i=i+2;
      end
    end
    
    function [A,B,C,D,xdot0,y0] = linearize(obj,t0,x0,u0)
      % Linearize the system about an operating point (continuous time)
      % @param t0 nominal scalar time
      % @param x0 nomimal vector of discrete+continuous state
      % @param u0 nominal vector input
      % @retval A,B,C,D,xdot0,y0   linear dynamics such that xdot=Ax+Bu+xdot0, y=C*x+Du+y0

      % note: also worth looking into using linearize instead of linmod: http://www.mathworks.com/help/toolbox/slcontrol/ug/linearize.html
      mdl = getModel(obj);
      if (~strcmp(get_param(mdl,'SimulationStatus'),'stopped'))
        feval(mdl,[],[],[],'term');
      end
      [A,B,C,D] = linmod(mdl,x0,u0,[1e-5,t0]);
      
      if (nargout>4) 
        xdot0 = dynamics(obj,t0,x0,u0);
        if (length(xdot0)~=size(A,1)) 
          % linmod did something clever, like adding states to the model.  
          error('have to handle this case more carefully'); 
        end
        if (nargout>5)
          y0 = output(obj,t0,x0,u0);
        end
      end
    end
    function [A,B,C,D,xn0,y0] = dlinearize(obj,ts,t0,x0,u0)
      % Linearize the system about an operating point (discrete time)
      %
      % @param ts sample time for linearization
      % @param t0 nominal scalar time
      % @param x0 nomimal vector of discrete+continuous state
      % @param u0 nominal vector input
      %
      % @retval A,B,C,D,xn0,y0   linear dynamics such that x[n+1]=Ax[n]+Bu[n]+xn0, y[n]=C*x[n]+Du[n]+y0

      mdl = getModel(obj);
      if (~strcmp(get_param(mdl,'SimulationStatus'),'stopped'))
        feval(mdl,[],[],[],'term');
      end
      [A,B,C,D] = dlinmod(mdl,ts,x0,u0,[1e-5,t0]);
      if (nargout>4) 
        xn0 = update(obj,t0,x0,u0);
        if (length(xn0)~=size(A,1)) 
          % linmod did something clever, like adding states to the model.  
          error('have to handle this case more carefully'); 
        end
        if (nargout>5)
          y0 = output(obj,t0,x0,u0);
        end
      end
    end
    
    function runLCMPlant(obj,lcmCoder,options)
      % Runs the system as an LCM client which subscribes to u and publishes xhat
      % @param lcmCoder an LCMCoder object, which defines all the messages
      % @param options see the options for the runLCM() method
      if (nargin<3) options = struct(); end
      runLCM(obj,lcmCoder,'u','xhat',options);
    end
    function runLCMControl(obj,lcmCoder,options)
      % Runs the system as an LCM client which subscribes to xhat and publishes u
      % @param lcmCoder an LCMCoder object, which defines all the messages
      % @param options see the options for the runLCM() method
      if (nargin<3) options = struct(); end
      runLCM(obj,lcmCoder,'xhat','u',options);
    end
    function runLCMEstimator(obj,lcmCoder,options)
      % Runs the system as an LCM client which subscribes to y and publishes xhat
      % @param lcmCoder an LCMCoder object, which defines all the messages
      % @param options see the options for the runLCM() method
      if (nargin<3) options = struct(); end
      runLCM(obj,lcmCoder,'y','xhat',options);
    end
    function runLCMVisualizer(obj,lcmCoder,options)
      % Runs the system as an LCM client which subscribes to xhat and publishes nothing
      % @param lcmCoder an LCMCoder object, which defines all the messages
      % @param options see the options for the runLCM() method
      if (nargin<3) options = struct(); end
      runLCM(obj,lcmCoder,'xhat',[],options);
    end

  end
  
  properties (SetAccess=private,GetAccess=protected)
    time_invariant_flag = false;  % set to true if you know the system is time invariant
    simulink_params=struct();     % simulink model parameters
    structured_x;                 % simulink state structure (cached for efficiency)
  end
  properties (SetAccess=private,GetAccess=public)
    input_angle_flag = [];  % in_angle_flag(i)==true means that input i is an angle and should be wrapped around 2*pi when appropriate
    state_angle_flag = [];  % state_angle_flag(i)==true means that state i is an angle and should be wrapped around 2*pi when appropriate
    output_angle_flag = []; % out_angle_flag(i)==true means that output i is an angle and should be wrapped around 2*pi when appropriate
  end

end