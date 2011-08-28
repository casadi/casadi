casadi_main
casadi_primitive
casadi_primitive_tools
casadi_noncore

global casadi = struct();
names = {'casadi_main','casadi_primitive','casadi_primitive_tools', 'casadi_noncore'};
interfaces = { casadi_main casadi_primitive casadi_primitive_tools casadi_noncore};

for i=1:numel(names)
  name = names{i};
  s=completion_matches([name "."]);
  strlen = numel(name);
  for j=1:size(s,1)
    command = deblank(s(j,strlen+2:end));
    if ~isfield(casadi,command)
      casadi.(command) = subsref(interfaces{i},struct('type','.','subs',command));
    end
  end
end


casadi_helpers


