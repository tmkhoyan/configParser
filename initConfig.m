function initConfig(varargin)

mydir = pwd();

if isempty(varargin)
    pathConfigLib = '~/github/configParser';
else
    pathConfigLib= varargin{1};
end

a(1) = copyfile([pathConfigLib,'/master.m'],[mydir '/master.m']);
a(2) = copyfile([pathConfigLib,'/config/'], [mydir  '/config/']);

if sum(a)==2
    disp('Config files copied')
else
    warning('wrong path provided, files not copied.')
end
end