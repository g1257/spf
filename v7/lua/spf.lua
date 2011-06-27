
require("myspf")
g = myspf.GeometrySquareD(4);

print(g:length())

c = myspf.ConcurrencySerialD(#arg,arg[3]);

print(c:nprocs())

print(c:rank())