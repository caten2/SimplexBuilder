import sage.all
from sage.groups.perm_gps.permgroup_named import SymmetricGroup, AlternatingGroup, DihedralGroup
# from SimplexBuilder import GroupComplex

print("We build the collection of complexes for the symmetric group on 4 elements and print the associated documentation string.")
C = GroupComplex(SymmetricGroup(4))
print(C)
print

print("We can check to see what group belongs to C.")
print(C.group)
print

print("The object C won't be constructed if the group isn't abelian, but we can check anyway.")
print("This also illustrates a convenient way to perform arbitrary group methods on the group corresponding to our simplicial complex object.")
print(C.group.is_abelian())
print

print("We can look at the list of conjugacy class representatives C uses.")
print(C.representatives)
print

print("We can also get a dictionary showing the elements of each of those conjugacy classes, indexed by their representative elements.")
print(C.conjugacyClasses)
print

print("We examine the sheets about each element in a conjugacy class by giving the class's representative element. In this case we conisder the sheets about (2,4), which is listed as the first conjugacy class representative.")
print(C.Sheetbuilder(C.representatives[0]))
print

print("The detailedInfo dictionary contains information on each component of the construction.")
print("Components are assigned an integer label.")
print

print("A list of the sheet labels in the zeroth component of the construction follows.")
c = 0
print(C.detailedInfo[c]["raw data"])
print

print("One can also view more verbose information about the vertices in each polygonal sheet.")
print(C.detailedInfo[c]["parsed data"])
print

print("The numbers of vertices, edges, and polygonal faces in the component are given below.")
print(C.detailedInfo[c]["vertices"])
print(C.detailedInfo[c]["edges"])
print(C.detailedInfo[c]["faces"])
print
 
print("The number of edges in each polygonal face is given below.")
print(C.detailedInfo[c]["polygon"])
print
 
print("Finally, we can calulate the genus of the component.")
print(C.detailedInfo[c]["genus"])
print
 
C.GeneralInfo()