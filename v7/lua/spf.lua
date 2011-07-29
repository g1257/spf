
require("myspf")

concurrency = myspf.ConcurrencyT(#arg,arg[3]);

print(concurrency:nprocs())

--print(c:rank())

filename = "input50.inp"
io = myspf.IoSimpleIn(filename);
print(io:filename())

engineParams = myspf.ParametersEngineT(io);

mp=myspf.ParametersModelT(io,engineParams);
print(mp)

license = "C";
if (concurrency:root()) then  print(license) end

len1 = engineParams.latticeLength
print(len1)
geometry = myspf.GeometryT(len1);

model = myspf.ModelT(engineParams,mp,geometry,concurrency)
algorithm = myspf.AlgorithmT(engineParams,model);
engine= myspf.EngineT(engineParams,model,algorithm,concurrency);

engine:main()
--print(g:length())


-- r = myspf.Random48D();
-- r:seed(3439483);
-- print(r:random())
-- print(r:random())


