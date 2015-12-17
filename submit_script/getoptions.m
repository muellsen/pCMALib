% ---------------------------------------------------------------  
function opts=getoptions(inopts, defopts)
% OPTS = GETOPTIONS(INOPTS, DEFOPTS) handles an arbitrary number of
% optional arguments to a function. The given arguments are collected
% in the struct INOPTS.  GETOPTIONS matches INOPTS with a default
% options struct DEFOPTS and returns the merge OPTS.  Empty or missing
% fields in INOPTS invoke the default value.  Fieldnames in INOPTS can
% be abbreviated.
%
% The returned struct OPTS is first assigned to DEFOPTS. Then any
% field value in OPTS is replaced by the respective field value of
% INOPTS if (1) the field unambiguously (case-insensitive) matches
% with the fieldname in INOPTS (cut down to the length of the INOPTS
% fieldname) and (2) the field is not empty.
%
% Example:
%   In the source-code of the function that needs optional
%   arguments, the last argument is the struct of optional
%   arguments:
%
%   function results = myfunction(mandatory_arg, inopts)
%     % Define four default options
%     defopts.PopulationSize = 200;
%     defopts.ParentNumber = 50;
%     defopts.MaxIterations = 1e6;
%     defopts.MaxSigma = 1;
%  
%     % merge default options with input options
%     opts = getoptions(inopts, defopts);
%
%     % Thats it! From now on the values in opts can be used
%     for i = 1:opts.PopulationSize
%       % do whatever
%       if sigma > opts.MaxSigma
%         % do whatever
%       end
%     end
%   
%   For calling the function myfunction with default options:
%   myfunction(argument1, []);
%   For calling the function myfunction with modified options:
%   opt.pop = 100; % redefine PopulationSize option
%   opt.PAR = 10;  % redefine ParentNumber option
%   opt.maxiter = 2; % opt.max=2 is ambiguous and would result in an error 
%   myfunction(argument1, opt);

%
% 04/07/19: Entries can be structs itself leading to a recursive
%           call to getoptions. 
%

if nargin < 2 || isempty(defopts) % no default options available
  opts=inopts;
  return;
elseif isempty(inopts) % empty inopts invoke default options
  opts = defopts;
  return;
elseif ~isstruct(defopts) % handle a single option value
  if isempty(inopts) 
    opts = defopts;
  elseif ~isstruct(inopts)
    opts = inopts;
  else
    error('Input options are a struct, while default options are not');
  end
  return;
elseif ~isstruct(inopts) % no valid input options
  error('The options need to be a struct or empty');
end

  opts = defopts; % start from defopts 
  % if necessary overwrite opts fields by inopts values
  defnames = fieldnames(defopts);
  idxmatched = []; % indices of defopts that already matched
  for name = fieldnames(inopts)'
    name = name{1}; % name of i-th inopts-field
%     if isoctave
%       for i = 1:size(defnames, 1)
% 	idx(i) = strncmp(lower(defnames(i)), lower(name), length(name));
%       end
%     else
	idx = strncmpi(defnames, name, length(name));
%    end
    if sum(idx) > 1
      error(['option "' name '" is not an unambigous abbreviation. ' ...
	     'Use opts=RMFIELD(opts, ''' name, ...
	     ''') to remove the field from the struct.']);
    end
    if sum(idx) == 1
      defname  = defnames{find(idx)}; 
      if ismember(find(idx), idxmatched)
	error(['input options match more than ones with "' ...
	       defname '". ' ...
	       'Use opts=RMFIELD(opts, ''' name, ...
	       ''') to remove the field from the struct.']);
      end
      idxmatched = [idxmatched find(idx)];
      val = getfield(inopts, name);
      % next line can replace previous line from MATLAB version 6.5.0 on and in octave
      % val = inopts.(name);
      if isstruct(val) % valid syntax only from version 6.5.0
	opts = setfield(opts, defname, ...
	    getoptions(val, getfield(defopts, defname))); 
      elseif isstruct(getfield(defopts, defname)) 
      % next three lines can replace previous three lines from MATLAB 
      % version 6.5.0 on
      %   opts.(defname) = ...
      %      getoptions(val, defopts.(defname)); 
      % elseif isstruct(defopts.(defname)) 
	warning(['option "' name '" disregarded (must be struct)']); 
      elseif ~isempty(val) % empty value: do nothing, i.e. stick to default
	opts = setfield(opts, defnames{find(idx)}, val);
	% next line can replace previous line from MATLAB version 6.5.0 on
	% opts.(defname) = inopts.(name); 
      end
    else
      warning(['option "' name '" disregarded (unknown field name)']);
    end
  end

% ---------------------------------------------------------------  