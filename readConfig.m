
function [configOptionsStruct, configOptionsCell, allheaders, symbolicDefs] = readConfig(filename,varargin)
% Config parser for ASCII based files

%   Author: Tigran Mkhoyan
%   Delft University of Technology, 2017

% see below for description
headerDelim = [];
maxheaders = [];
commentDelim = [];
warnEnabled = [];
maxheaders = [];


structnamefieldfillelemn = [];
setOptargs;

fileID = fopen(filename);
charsizecomment = numel(commentDelim); 

fgets(fileID); %skip header
N = 1000; %lines

configOptionsCell = cell(N/10,2);
allheaders = cell(N/10,2);

headercount = 1;
loopcount = 1;
rowcount = 2; %to make sure first heder gets

noncharIdx = zeros(size(configOptionsCell));
%put back old
while (~feof(fileID) && loopcount < N && ((headercount+rowcount) <= (maxheaders+10)))
    
    %read
    nextline = strtrim(fgetl(fileID));
    numchar = numel(nextline);
    if (~isempty(nextline) && sum(nextline(1:min(numchar,charsizecomment)) == commentDelim) ~= charsizecomment)
        
        if (nextline([1 end]) == headerDelim)  % % new header
            
            rowcount = 2;
            configOptionsCell(1,headercount) = {nextline};    %header
            allheaders(headercount,1) = {nextline(2:end-1)}; %row two is header without the brakets
            headercount = headercount + 1;
            % optionList
        elseif (headercount>1)
            if (isempty(str2num(nextline))) %#ok<*ST2NM> %%&& ~strcmp(nextline,'[]'))
                %if (isempty(str2double(nextline))) %%&& ~strcmp(nextline,'[]'))
                configOptionsCell(rowcount,headercount-1) = {nextline}; %option
                noncharIdx(rowcount,headercount-1) = 0;
            else
                configOptionsCell(rowcount,headercount-1) = {str2num(nextline)}; %convert to numbers
                %configOptionsCell(rowcount,headercount-1) = {str2double(nextline)}; %convert to numbers
                noncharIdx(rowcount,headercount-1) = 1;
            end
            rowcount = rowcount + 1;
            allheaders(headercount-1,2) = {rowcount-2}; %update number of elements
        end
    end
    
    loopcount = loopcount + 1;
    
    %TODO fix the number of headers part
    if (headercount+rowcount) == (maxheaders+10)
        warning('Maxixmum number of header elements reached. Ignoring extra elemenst. Use read config with optional argumens to increase')
    end
end

%cleanup options
linenum = loopcount; %#ok<NASGU>
emptyRows = ~sum(~cellfun('isempty',configOptionsCell),2); %emty rows
configOptionsCell(emptyRows,:) = [];
noncharIdx(emptyRows,:) = [];
% checkheaders = configOptionsCell{~noncharIdx};
% checkvals    = configOptionsCell{2,k};

noncharIdx(1,:) = 1; % dont covert headers
noncharIdx = logical(noncharIdx);
allheaders = allheaders(1:size(configOptionsCell,2),:);

configOptionsCell(cellfun('isempty',configOptionsCell)) = {'~NaN'}; % fill all empty fields with nan

%TODO consider using function handle instead of parsing pwd. However then
%option must be accesed as opt.myopt();
%mydir = @() pwd;
%configOptionsCell(strcmp(configOptionsCell,'pwd')) = {mydir}; //to make
%it an aninymous function. However in usage option must be called with
%opt.mydir()
configOptionsCell(strcmp(configOptionsCell,'pwd')) = {pwd()};

%Update2.0: functionality header addition:
%c1 = configOptionsCell;

for k = 1:numel(configOptionsCell(1,:))
    checkheader = configOptionsCell{1,k};
    checkval    = configOptionsCell{2,k}; %TODO:takes only first value if the header has multiple values wont work
    if(ischar(checkval))
        configOptionsCell(~noncharIdx) = regexprep(configOptionsCell(~noncharIdx),checkheader,checkval);
    end
end
%check for double path //

configOptionsCell(~noncharIdx) = regexprep(configOptionsCell(~noncharIdx),'//','/');

% %in case [] is config option make sure its numeric not string
configOptionsCell(strcmp(configOptionsCell,'[]')) = {[]};
%create struct

%Update 2.2: recognize = sign and split header options
% for k = 1:numel(configOptionsCell(1,:))
%     %for m =1:numel(configOptionsCell())
%     checkval    = configOptionsCell{2,k}; %TODO:takes only first value if the header has multiple values wont work
%     if(ischar(checkval))
%         configOptionsCell(~noncharIdx) = regexprep(configOptionsCell(~noncharIdx),'=','split');
%     end
% end

%generate temp option cell and replace non char values
tmpconfigOptionsCell = configOptionsCell;
tmpconfigOptionsCell(noncharIdx) = {''};

%dynamicHeaderValIdx = find(~~cellfun(@sum,(regexp(tmpconfigOptionsCell,'=')))); % get the ids of cells to be replaced
%nonchaidx and configoptionscell are always same size
dynamicHeaderValIdx = (~~cellfun(@sum,(regexp(tmpconfigOptionsCell,'=')))); % get the ids of cells to be replaced


toSplitIdx = ~~(zeros(size(noncharIdx)));
toSplitIdx(dynamicHeaderValIdx) = true;


configOptionsCell(toSplitIdx) = strtrim(regexp(configOptionsCell(toSplitIdx),'=','split'));
toSplitHeaderIdx = sum(toSplitIdx)>=1; %indeces of headers with special = operator

%Update 3.0 functionality for simbolic addition
% find symbolic expressions and replace them with numeric values converted
% to strings. This works only with numeric values and variables under the headers preceeded by expression '= @'
[configOptionsCell(toSplitIdx), symbolicDefs]= evalSymbolicVar(configOptionsCell(toSplitIdx));
%TODO consider putting symbolicDefs in configOptionsStruct

%zz = regexp(configOptionsCell(~noncharIdx),'=','split');
%zz2 = regexp(configOptionsCell(:,9),'=','split');
% update: check if header contains same variable twice remove (only first one is used)
% [~,uniqueindeces] = unique(allheaders(:,1));
% uniqueindeces = sort(uniqueindeces); % keep order
% allheaders = allheaders(uniqueindeces,:);
% headercount = numel(uniqueindeces)+1;
% configOptionsCell = configOptionsCell(:,uniqueindeces);

%UPDATE for simulink parser
%replace headers names with valid fielnames for structure
charToreplace = {'/',' '};
configOptionsCell(1,:) = regexprep(configOptionsCell(1,:),charToreplace,structnamefieldfillelemn);
allheaders(:,1) = regexprep(allheaders(:,1),charToreplace,structnamefieldfillelemn);

[allheaders, uniqueindeces] = checkUniqness(allheaders);
headercount = numel(uniqueindeces)+1;
configOptionsCell = configOptionsCell(:,uniqueindeces);
%do the same for
%match options if given


configOptionsStruct = struct();
if isempty(optionsList)
    optionsList = allheaders(:,1); %use all headers in case no specific headers are provded to build the struct
else
    %replace headers names with valid fielnames for structure
    optionsList = checkUniqness(optionsList);
    optionsList= regexprep(optionsList,charToreplace,structnamefieldfillelemn);
end

%  configOptionsStruct = initStruct(optionsList);
optionsSelected = cell(max([allheaders{:,2}]+1),numel(optionsList));

if (headercount-1)>numel(optionsList)
    if (warnEnabled)
        warning(['Extra heading entries in config file, ',filename, ', ignored'])
    end
end

for i=1:numel(optionsList)
    
    columnmatched = strcmpi(allheaders(:,1),optionsList(i));
    if ~columnmatched
        error(['config option, ', optionsList{i}, ', not found'])
    else
        rowsmatched = 1:(1+allheaders{columnmatched,2});
    end
    optionsSelected(rowsmatched,i) = configOptionsCell(rowsmatched,columnmatched);
    
    %set struct fields
    fieldname  = optionsList{i};
    % fieldname(isspace(fieldname)) = structnamefieldfillelemn; already
    % done by regexprep
    if (max(rowsmatched)>2 || toSplitHeaderIdx(i)) % if we have special headers '=' put them al in cell container
        configOptionsStruct = setfield(configOptionsStruct,fieldname,configOptionsCell(rowsmatched(rowsmatched>=2),columnmatched)); %#ok<SFLD> %fieldvalue is cell
    else
        configOptionsStruct = setfield(configOptionsStruct,fieldname,configOptionsCell{rowsmatched(rowsmatched>=2),columnmatched}); %#ok<SFLD> %fieldvalue is just value
    end
    
end

configOptionsCell = optionsSelected;


%     function s = initStruct(fieldnames)
%     zz = fieldnames;
%     zz(:,2) = {[]};
%     zz = zz';
%     s = struct(zz{:});
%
%     end

    function [incell, uniqueIdx] = checkUniqness(incell)
        [~,uniqueIdx] = unique(incell(:,1)); % checks only first column
        uniqueIdx = sort(uniqueIdx); % keep order
        incell = incell(uniqueIdx,:); % returns full cell
    end

%TODO consider returning the struct s with sorted '=' variables
%irrespective of existence of symbolic definitions
    function  [original, sTempValContainer] = evalSymbolicVar(original) %returns
        %create a safe workspace to evaluate simbolic expressions
        vec = original; %temporary container
        vec = [vec{:}]'; %unfold array
        
        sTempValContainer = []; % if no symbolic expressions the struct will be returned as empty
        
        if(~isempty(vec))
            idx_var = ~~zeros(size(vec));
            idx_var(1:2:end) = true; % skip one from 1
            idx_val = ~idx_var; % skip one from 2 i.e shift by one
            
            vec_var = vec(idx_var);
            vec_val = vec(idx_val);
            
            %find symbolic expressions in Vals indicated by = @
            [~,idx_sym] = regexp(vec(idx_val),'@','once');
            idx_sym = ~cellfun(@isempty,idx_sym);
            
            if(sum(idx_sym)) % if we have symbolic expressions else dont bother
                
                %find in symbolic expressions patterns that contain variables and replace with appended variable sTempValContainer.* such that it can later be evaluated using eval
                % pattern \< 'var(k)' \> is appended to variable such that the search is constrained to whole word. This prevents long named variables partially beeing replaced by short variable names e.g. var a -> cls.alpha instead of s.clapha vaiables
                rvec_val = vec_val;
                
                rvec_val(idx_sym) = regexprep(vec_val(idx_sym), strcat('\<',vec_var,'\>'), strcat('sTempValContainer.',vec_var)); %accounts for nested relationships i.e. also replaces the definition with variables that are defined as symbolic
                
                idx_sym_matlab =  regexp(rvec_val,'@(','once'); % match definition @() then it is a matlab symbolic definition
                idx_sym_matlab = ~cellfun(@isempty,idx_sym_matlab);
                
                idx_sym_keepstr =  regexp(rvec_val,'@\[%s\]','once'); % match definition @() then it is must be kept as string without evaluating
                idx_sym_keepstr = ~cellfun(@isempty,idx_sym_keepstr);
                
                rvec_val = strtrim(regexprep(rvec_val,'@\[%s\]',''));
                vec_val = strtrim(regexprep(vec_val,'@\[%s\]',''));
                
                %remove @ from value field except @(..) match anything but '('
                %after @ '@[^(]' . This allows to evaluate matlab syntax
                %symbolic expressions e.g. @(x,y) (x^2 + y^2)
                vec_val = regexprep(vec_val,'@[^(]','');
                rvec_val = regexprep(rvec_val,'@[^(]','');
                
                %evaluate first non symbolic to obtain constants then sibolic
                %expressions
                sTempValContainer = struct();
                xval = [];  %#ok<NASGU>
                idx = find(~idx_sym);
                for l_count=1:numel(idx) %make sure l is not a var
                    try
                        %eval(combinedCell{idx(l)}); % we just need them in the workspace but eval wont work so put in structure
                        xval = str2num(vec_val{idx(l_count)});
                        sTempValContainer.(vec_var{idx(l_count)}) = xval; %will return empty value if not properly matched
                        
                        %TODO check if string flag works
                        if(isempty(xval) && ~idx_sym_keepstr(idx(l_count)))
                            warning('Check expression ''%s'' forgotten ''=@'' or wrong syntax. Value of variable %s is set to empty',...
                                [vec_var{idx(l_count)},' = ',vec_val{idx(l_count)}],vec_var{idx(l_count)});
                        end
                        %else it had keep string flag @[s]
                    catch
                        warning('Could not evaluate expression ''%s''. Check variable definition ''%s '' \n',...
                            [vec_var{idx(l_count)},' = ',vec_val{idx(l_count)}],vec_var{idx(l_count)});
                    end
                end
                %now we have all the variables irrespective of order in thestructure so we canevaluate expressions
                
                % First pass: see if can be evaluated otherwise store as failed attempt
                idx = find(idx_sym);
                idx_fail = idx_sym<0;
                for l_count=1:numel(idx)
                    try
                        %                     expression = regexprep(vec_val{idx(l)},)
                        if(idx_sym_keepstr(idx(l_count)))
                            xval = vec_val{idx(l_count)};
                            original{idx(l_count)}{2} = xval; %kept string
                        elseif(idx_sym_matlab(idx(l_count)))  %matlab symbolic definition
                            xval = eval(rvec_val{idx(l_count)});
                            original{idx(l_count)}{2} = xval; % we need to convert back to string the second entry that is the value of the variable
                            %      rvec_val{idx(l)}    = xval; % also update the value so that ecursive relationships are evaluated!
                        else
                            xval = eval(rvec_val{idx(l_count)});
                            original{idx(l_count)}{2} = num2str(xval); % we need to convert back to string the second entry that is the value of the variable
                        end
                        %else if matlab symbolic definition then dont convert to string just return the handle in array without evaluating
                        sTempValContainer.(vec_var{idx(l_count)}) = xval;
                        rvec_val{idx(l_count)}    = xval; % also update the value for symbolic expressions, these can then be evaluated
                    catch
                        idx_fail(idx(l_count)) = true;
                    end
                end
                % REMAINING PASSES: evaluating nesterd variables undefined symbolic variables. So the order wont matter!
                xval=[]; %#ok<NASGU>
                idx_fail_prev = idx_fail;
                while sum(idx_fail)
                    idx = find(idx_fail);
                    for l_count=1:numel(idx)
                        
                        try
                            xval = eval(rvec_val{idx(l_count)}); % try second time
                            original{idx(l_count)}{2} = num2str(xval);
                            
                            %store if succesfull
                            sTempValContainer.(vec_var{idx(l_count)}) = xval;
                            rvec_val{idx(l_count)}    = xval;
                            idx_fail(idx(l_count)) = false;
                        catch
                            idx_fail(idx(l_count)) = true;
                        end
                        
                    end
                    
                    if ~(sum(idx_fail)<sum(idx_fail_prev)) % if unchanged throw error warning and exit loop
                        
                        idx = find(idx_fail);
                        for m_count=1:numel(idx)
                            warning('Could not evaluate expression ''%s''. Check variable definition ''%s '' \n',...
                                [vec_var{idx(l_count)},' = ',vec_val{idx(l_count)}],vec_var{idx(l_count)});
                        end
                        break
                    end
                    
                    idx_fail_prev = idx_fail;
                    
                end
            end
        end
        %return original vector
    end

%UPDATE simulinkparser: now the runtimesettings can be set using two
%methods: 1. via providing ordered variable arguments. 2. Via providing
%name tag of runtime option followed by the value e.q.
%readconfig(...,'structnamefieldfillelemn','__',...) irrespective of the order of variable argument this usage of readConfig will pick the right configuration for the explicit setting
    function setOptargs
        numvarargs  = length(varargin);
        
        optargs = {cell(0), '//', '{}',true, '_', 100};
        %optargs{1:numvarargs} = varargin;
        %UPDATE simulink parser. Now one option is parced
        runtimeSettings = {'optionsList', 'commentDelim', 'headerDelim', 'warnEnabled', 'structnamefieldfillelemn', 'maxheaders'};
        idx = cell(size(runtimeSettings));
        idxval = idx;
        id = idx;
        %idxval = ~~(zeros(numvarargs,numel(runtimeSettings)));
        
        tmp = (zeros(size(runtimeSettings)));
        % check if explicit options are provided
        for n=1:numel(runtimeSettings)
            tmp(n) = 1;
            idx{n} = strcmpi(varargin,runtimeSettings{n})';
            idxval{n} = circshift(idx{n},1);
            id{n} = (tmp*sum(idx{n})==1)'; % (tmp*(sum(idx{n})==1))'; gives doubles array!
            tmp = tmp*0;
        end
        
        if sum(cellfun(@sum,idx)>1) %if duplicates
            str = strcat(runtimeSettings,',');
            error('runtime setting ''%s'' specified more than once',[str{sum(zz,1)>1}])
            
        elseif ~(sum(cellfun(@sum,idx))) % regular way of setting runtimesettings
            % set defaults for optional inputs
            if numvarargs > 10
                error('functions:randRange:TooManyInputs', ...
                    'requires atmost 10 optional input');
            end
            %TODO: make this work for all variable arguments
            [optargs{1:numvarargs}] = varargin{:};
            [optionsList, commentDelim, headerDelim, warnEnabled, structnamefieldfillelemn, maxheaders] = optargs{:};
        else % parce from strings order does not matter
            %set optional argumenst and apliy specific runtime settings for the remaining
            %idxval = circshift(idx,1); %value must be the enry afer
            for n=1:numel(runtimeSettings)
                [optargs{id{n}}] = varargin{idxval{n}};
            end
            [optionsList, commentDelim, headerDelim, warnEnabled, structnamefieldfillelemn, maxheaders] = optargs{:};
            
        end
    end

end

% ---------------------------------------------------------------------------------:

% Description/Usage:
% -A simle config parser for types of string,boolean, vector,cell and matrix values.
%   The parser uses simple headers with customizable header and comment delimiters.
%   In other words use a delimiter of your preference for comments ( e.g. // or %, or  ) and setting headers (e.g. {myoption1} or [myoption2]).
%   De delimiters can be provided to readConfig as second and third input for comments and headers respectively. First option is the path to the config file.
% -The parser will generate a stuct with settings that are accecible via same header names. Example:

%   in config file  --> set option with
%   {myoption1}
%   1 2 3
%   in matlab file  --> acces vector [1 2 3] opt struct with
%   opt.myoption1

%   The parser will automatically recognize the following types:

%  -(single) path/string    = string with pathname/caracter
%  -(multiple) path/        = 1xm cell with m pathnames listed after {header}
%  - boolean                = false or true logical (case sensitive!)
%  - vector                 = simple 1xm vector with m space separated values after header
%  - cell/matrix            = 2D matrix with m columns (space separated values) and n rows (lines after header). Note that the           function returns a cell. To convert to array simple do cell2mat(opt.mycellarray)


% -File type is also customizable (i.e. txt or conf or whatever)

% Usage:
% -create a txt file similar to config.txt
%   e.g.
% -put desired header names in a cell with strings (order is not relevant but caracters must match)
% -acces config with:
%     [options] = readConfig('path/to/config',headernames,'//','{}',false); % varargin: 1=commentdelim, 2=headerdelim 3: -
% -use additional options to customyze to needs
% loadOptionsList is a local function to create the header names. If this is not desired simply use:
% optionsList = {...
%      'path'
%      'paths'
%      'boolean'
%      'vector'
%      'cell'

% Usage of 2.0 functionality:
%
% - as of 2.0 the parser will recognize pwd command. so using:
% {mydir} --> {mydir}
% pwd     --> path/to/location/calling/master/function
% In other words if you call readConfig with mast.m located at ~/somedirectory/ the options.mydir = ~/somedirectory/
%
% Parser will now also allow variable additions. In other words if you define a variable {mydir} you can use this to construct a new variable composed partially or fully of {mydir}. Example:
%
% define:
%  {mydir}
% ~/somepath/
% use {mydir} to definecomposite path variable
%
% {myotherdir}
% {mydir}/../
%
% Now the {myotherdir} will point to a folder level higher than {mydir} (this is due to addition ../)
% This can also be use for addition of any string variables. Note that the functionality doesn't work yet for numeric values. In case there is interest i might consider adding that.