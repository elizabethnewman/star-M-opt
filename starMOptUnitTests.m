starMOptSetup;

myDir   = 'unitTests';
myFiles = dir(myDir);

for i = 1:length(myFiles)
    name = myFiles(i).name;
    if contains(name,'.m') && contains(name,'UnitTests')
        results = runtests(name); 
        assertSuccess(results);
    end
end
