import JuPyMake

JuPyMake.InitializePolymake()

print(JuPyMake.ExecuteCommand("print 1+1;"))
print(JuPyMake.ExecuteCommand('application "polytope";'))
print(JuPyMake.ExecuteCommand("$p = cube(4);"))
print(JuPyMake.ExecuteCommand("print $p->VOLUME;"))
