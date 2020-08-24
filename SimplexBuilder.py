# Simplex Builder
# Charlotte Aten 2015
# Based on a program by Mark Herman
# Constructs simplices that are described in "On a canonical construction of tesselated surfaces
# via finite group theory" by Mark Herman, Jonanthan Pakianathan, and Ergun Yalcin (2013)

import sage.all
from sage.matrix.constructor import zero_matrix
from sage.graphs.graph import Graph

class GroupComplex():
    """
    Examine a canonical simplicial complex construction using a given nonabelian group.
    """

    def __init__(self,group):
        # Require that group is abelian.
        assert group.is_abelian() == False
        # Return the underlying group.
        self.group = group
        # Return a list of noncentral elements in the group.
        self.oppCen = list(set(group).difference(set(group.center())))
        # Return a list of conjugacy class representatives for noncentral elements.
        self.representatives = self.Reps()
        # Return a dictionary of conjugacy classes of noncentral elements.
        self.conjugacyClasses = self.ConClasses()
        # Return a list of the polyhedral sheets in the construction.
        self.sheets = list(self.SheetIterator())
        # Return the detailed component information for the construction as a dictionary.
        self.detailedInfo = self.ComponentBuilder()
    def __repr__(self):
        return "The collection of simplicial complexes corresponding to the %s." % (self.group)

    # Return a list of conjugacy class representatives which are noncentral in the group.
    def Reps(self):
        lis = []
        for elem in self.group.conjugacy_classes_representatives():
            if elem in self.oppCen:
                lis.append(elem)
        return lis

    # Return a dictionary indexed by the conjugacy class representatives chosen, with entries containing lists of all elements in the given conjugacy class.
    def ConClasses(self):
        dic = {}
        for elem in self.representatives:
            dic[elem] = [elem]
        for candidate in set(self.oppCen)-set(self.representatives):
            for elem in self.oppCen:
                if elem**(-1)*candidate*elem in self.representatives:
                    dic[elem**(-1)*candidate*elem].append(candidate)
                    break
        return dic

    # Test whether two sheets given as lists have an edge in common.
    def SheetAdjacency(self,sheet0,sheet1):
        for elem0 in sheet0[:-1]:
            for elem1 in sheet1[:-1]:
                if [elem0,sheet0[sheet0.index(elem0)+1]] == [sheet1[sheet1.index(elem1)+1],elem1] or [elem0,sheet0[sheet0.index(elem0)+1]] == [elem1,sheet1[sheet1.index(elem1)+1]]:
                    return True
        return False

    # Return a dictionary of sheets indexed by the (type 2) elements of a conjugacy class given the chosen representative for that class.
    def Sheetbuilder(self,conclassrep):
        dic = {}
        for elem0 in self.conjugacyClasses[conclassrep]:
            counter = 0
            dic[elem0] = {}
            for elem1 in self.oppCen:
                if elem0*elem1 != elem1*elem0:
                    sheet_candidate = [elem1,elem1**(-1)*elem0]
                    while sheet_candidate[0] != sheet_candidate[-1]:
                        sheet_candidate.append(sheet_candidate[-1]**(-1)*elem0)
                    for key in dic[elem0].keys():
                        if set(sheet_candidate) == set(dic[elem0][key]):
                            sheet_candidate = 'already obtained'
                    if sheet_candidate != 'already obtained':
                        dic[elem0][counter] = sheet_candidate
                        counter = counter + 1
        return dic

    # Create an iterator for all the polygonal sheets in the construction.
    def SheetIterator(self):
        for rep in self.representatives:
            for elem in self.Sheetbuilder(rep):
                for sheetkey in self.Sheetbuilder(rep)[elem]:
                    yield [(elem,sheetkey),self.Sheetbuilder(rep)[elem][sheetkey]]

    # Return a dictionary with a numerical key for each component of the construction.
    def ComponentBuilder(self):
        """
        Create a dictionary containing information about the simplicial complex construction.
        The keys for the dictionary are nonnegative integers assigned to the components of the construction as they are created.
        Each of these keys has a dictionary as an entry, whose keys are as follows:
            "raw data" - Return the list created initially in constructing the complex. Sheets are labeled by numbers.
            "parsed data" - Return the list obtained by replacing each number in the "raw data" list with its corresponding sheet information.
            "vertices" - Return the number of vertices in the component.
            "edges" - Return the number of edges in the component.
            "faces" - Return the number of triangular faces in the component.
            "polygon" - Return the number of sides in each polygonal face of the component.
            "genus" - Return the genus of the component.
        """
        dic = {}
        size = len(self.sheets)
        # Create a sparse matrix whose ij entry is 1 if sheet number i and sheet number j are adjacent.
        mat = zero_matrix(size,size,sparse=True)
        for i in range(size):
            adjacentSheetsFound = sum(mat[i])
            while adjacentSheetsFound < len(self.sheets[i][1])+1:
                for j in range(i,size):
                    if self.SheetAdjacency(self.sheets[i][1],self.sheets[j][1]):
                        mat[i,j] = 1
                        mat[j,i] = 1
                        adjacentSheetsFound = adjacentSheetsFound + 1
        # Create a graph from the previously-created adjacency matrix.
        G = Graph(mat)
        # The components of this graph correspond exactly to the components of the construction.
        rawData = G.connected_components()
        # Create a dictionary entry with information about each component, as arbitrarily numbered by the graph method used previously.
        for component in range(len(rawData)):
            dic[component] = {}
            dic[component]["raw data"] = rawData[component]
            dic[component]["faces"] = len(dic[component]["raw data"])
            dic[component]["parsed data"] = [self.sheets[i] for i in rawData[component]]
            dic[component]["polygon"] = len(dic[component]["parsed data"][0][1]) - 1
            # We count the vertices in a component.
            vertexCount = 0
            # First find all the group elements appearing as type 1 vertices.
            putativeVertices = set()
            for sheet in dic[component]["parsed data"]:
                for vert in sheet[1]:
                    putativeVertices.add(vert)
            # Then for each element determine the number of distinct conjugacy classes it induces on the component.
            for vert in putativeVertices:
                containingSheetElements = set()
                for sheet in dic[component]["parsed data"]:
                    if vert in sheet[1]:
                        containingSheetElements.add(sheet[0][0])
                while containingSheetElements != set():
                    elem = containingSheetElements.pop()
                    conjugates = set([vert**(-i)*elem*vert**i for i in range(1,vert.order())])
                    containingSheetElements = containingSheetElements - conjugates
                    vertexCount = vertexCount + 1
            dic[component]["vertices"] = vertexCount
            # We count the edges in a component.
            edgeCount = 0
            size = dic[component]["faces"]
            for i in range(size-1):
                for j in range(i+1,size):
                    sheet0 = dic[component]["parsed data"][i][1]
                    sheet1 = dic[component]["parsed data"][j][1]
                    for s in range(dic[component]["polygon"]):
                        for t in range(dic[component]["polygon"]):
                            if [sheet0[s],sheet0[s+1]] == [sheet1[t],sheet1[t+1]] or [sheet0[s],sheet0[s+1]] == [sheet1[t+1],sheet1[t]]:
                                edgeCount = edgeCount + 1
            dic[component]["edges"] = edgeCount
            dic[component]["Euler characteristic"] = dic[component]["vertices"] - dic[component]["edges"] + dic[component]["faces"]
            dic[component]["genus"] = (dic[component]["Euler characteristic"]-2)/(-2)
        return dic
    
    def GeneralInfo(self):
        compNum = len(self.detailedInfo.keys())
        dic = {}
        for component in self.detailedInfo:
            if self.detailedInfo[component]["genus"] not in dic:
                dic[self.detailedInfo[component]["genus"]] = {self.detailedInfo[component]["polygon"]: 1}
            else:
                if self.detailedInfo[component]["polygon"] not in dic[self.detailedInfo[component]["genus"]]:
                    dic[self.detailedInfo[component]["genus"]][self.detailedInfo[component]["polygon"]] = 1
                else:
                    dic[self.detailedInfo[component]["genus"]][self.detailedInfo[component]["polygon"]] = dic[self.detailedInfo[component]["genus"]][self.detailedInfo[component]["polygon"]] + 1
        print("There are %s components in %s of %s different genus." % (compNum,self,len(dic.keys())))
        for genus in dic:
            print("Of those components of genus %s there are:" % (genus))
            for polygon in dic[genus]:
                print("    %s with %s-gonal faces" % (dic[genus][polygon],polygon))