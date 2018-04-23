function initConfig()

mydir = pwd();

a(1) = copyfile('~/github/configParser/master.m',[mydir '/master.m']);
a(2) = copyfile('~/github/configParser/config/', [mydir  '/config/']);

if sum(a)==2
   disp('Config files copied') 
end
end