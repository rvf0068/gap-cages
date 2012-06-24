# The point graph of the classical generalized quadrangle W(p^n)
InstallGlobalFunction( ClassicalGeneralizedQuadrangleW, function (p,n)
    local v,s1,s2,ti2,ady,H;
    H:=[[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]];
    v := FullRowSpace(GF(p^n),4);
    s1 := Elements(Subspaces(v,1));
    ady := function (u,v)
        local g1,g2;
        g1 := BasisVectors(Basis(u));
        g2 := BasisVectors(Basis(v));
        return g1<>g2 and IsZero(g1*H*TransposedMat(g2));
    end;
    return Graph(
                 Group(()),
                 s1,
                 TrivialAction,
                 ady,
                 true
                 );
end );

# terminology of J. A. Thas.
InstallGlobalFunction( SpanThas, function (g,x,y)
    local inters,a,b;
    a := ClosedNeighbourhood(g,x);
    b := ClosedNeighbourhood(g,y);
    inters := Intersection(a,b);
    return
      Filtered(Vertices(g),z->IsSubset(ClosedNeighbourhood(g,z),inters));
end );

# The point graph of the (q-1,q+1) quadrangle of Ahrens and Szekeres,
# with q=p^n
InstallGlobalFunction (AhrensSzkekeres, function (p,n)
    local w,points,typea,typeb,notnei,lines;
    w := ClassicalGeneralizedQuadrangleW(p,n);
    points := Filtered(Vertices(w),x->not(x in ClosedNeighbourhood(w,1)));
    typea := Filtered(CliqueComplex(w),x->not(1 in x));
    typeb := AsSet(List(points,y -> SpanThas(w,1,y)));
    lines := Union(typea,typeb);
    return BlockGraph(points,lines);
end );

InstallGlobalFunction (CageGirthSix, function (k)
    local d;
    d:=PGPointFlatBlockDesign(2,k-1,1);
    return BipartiteBlockGraph(Union(d.blocks),d.blocks);
end);

InstallGlobalFunction (CageGirthEight, function (k)
    local l,w;
    if IsPrimePowerInt(k-1) then
        l := Factors(k-1);
        w:= ClassicalGeneralizedQuadrangleW(l[1],Length(l));
        return VertexCliqueBipartite(w);
    else
        Print(k,"-1 is not a prime power\n");
        return fail;
    fi;
end);

# this seems not to be used
InstallGlobalFunction (IsTotallyIsotropicTwoDim, function(subs)
    local l,H;
    H:=[[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]];
    l := BasisVectors(Basis(subs));
    return IsZero(l[1]*H*TransposedMat([l[2]]));
end);
    
# a tree with one root of degree m, the other nonleaves of degree r
# useful for cages of girth 5 or 6
InstallGlobalFunction (TreeForCageGirthFive, function(r,m)
    local ord,h,i,j;
    ord := 1+m*r;
    h:=NullGraph(SymmetricGroup(ord));
    h:=NewGroupGraph(Group(()),h);
    for i in [1..m] do
        AddEdgeOrbit(h,[1,i+1]);
        for j in [1..r-1] do
            AddEdgeOrbit(h,[i+1,m+2+(i-1)*(r-1)+j-1]);
        od;
    od;
    return UnderlyingGraph(h);
end );

# a tree with one root of degree m, the other nonleaves of degree r
# useful for cages of girth 8
InstallGlobalFunction( TreeForCageGirthEight, function(r,m)
    local ord,h,i,j,k;
    ord := 1+m+m*(r-1)+m*(r-1)^2;
    h:=NullGraph(SymmetricGroup(ord));
    h:=NewGroupGraph(Group(()),h);
    for i in [1..m] do
        AddEdgeOrbit(h,[1,i+1]);
        for j in [1..r-1] do
            AddEdgeOrbit(h,[i+1,m+2+(i-1)*(r-1)+j-1]);
            for k in [1..r-1] do
                AddEdgeOrbit(h,[m+2+(i-1)*(r-1)+j-1,r*m+2+(i-1)*(r-1)^2+(j-1)*(r-1)+k-1]);
            od;
        od;
    od;
    return UnderlyingGraph(h);
end);

# 'graph adjunction'
# input is given by a graph g, and a list=pair of lists. The first
# list gives vertices of g, the second, and some edges of h
# output is a graph formed with vertices the disjoint union of V(g),
# V(h), in the resulting graph, g stays the same, the induced subgraph
# by h is a null graph, and if [x,[i,j]] is in the list, then vertex x
# in g is joined to both i and j
# for example
#g:=GraphAdjunction(TreeForCageGirthEight(3,3),
#[[11..22],[[1,5],[3,7],[2,6],[4,8],[1,2],[7,8],[3,4],[5,6],[1,4],[6,7],[2,3],[5,8]]])
# is the Tutte 8-cage
# Some lists that can be used with TreeForCageGirthEight(3,6):
# l1:= [ [ 20 .. 43 ], 
#  [ [ 2, 6 ], [ 4, 8 ], [ 10, 14 ], [ 12, 16 ], [ 1, 5 ], [ 3, 7 ], 
#      [ 9, 13 ], [ 11, 15 ], [ 1, 2 ], [ 7, 8 ], [ 9, 10 ], [ 15, 16 ], 
#      [ 2, 3 ], [ 5, 8 ], [ 10, 11 ], [ 13, 16 ], [ 1, 4 ], [ 6, 7 ], 
#      [ 9, 12 ], [ 14, 15 ], [ 3, 4 ], [ 5, 6 ], [ 11, 12 ], [ 13, 14] ] ];
# l2:=[ [ 20 .. 43 ], 
#  [ [ 2, 6 ], [ 10, 14 ], [ 4, 8 ], [ 12, 16 ], [ 1, 5 ], [ 9, 13 ], 
#      [ 3, 7 ], [ 11, 15 ], [ 1, 2 ], [ 9, 10 ], [ 7, 8 ], [ 15, 16 ], 
#      [ 2, 3 ], [ 10, 11 ], [ 5, 8 ], [ 13, 16 ], [ 1, 4 ], [ 9, 12 ], 
#      [ 6, 7 ], [ 14, 15 ], [ 3, 4 ], [ 11, 12 ], [ 5, 6 ], [ 13, 14 ] ] ];
# l3:=[ [ 20 .. 43 ], 
#  [ [ 1, 2 ], [ 4, 5 ], [ 6, 7 ], [ 9, 10 ], [ 2, 3 ], [ 5, 6 ], [ 7, 8 ], 
#      [ 10, 11 ], [ 1, 16 ], [ 3, 4 ], [ 11, 12 ], [ 14, 15 ], [ 8, 9 ], 
#      [ 12, 13 ], [ 3, 11 ], [ 15, 16 ], [ 1, 9 ], [ 13, 14 ], [ 4, 12 ], 
#      [ 8, 16 ], [ 5, 13 ], [ 7, 15 ], [ 2, 10 ], [ 6, 14 ] ] ];
# A list for TreeForCage(3,5), then one has to add an edge to make it
# a {3,5}-graph.
# l:=[ [ 17 .. 36 ],
#  [ [ 8, 9 ], [ 6, 11 ], [ 2, 3 ], [ 5, 14 ], [ 1, 2 ], [ 11, 12 ], [ 6, 7 ],
#      [ 13, 14 ], [ 2, 7 ], [ 9, 10 ], [ 3, 12 ], [ 5, 6 ], [ 4, 5 ],
#      [ 12, 13 ], [ 1, 10 ], [ 7, 8 ], [ 1, 14 ], [ 3, 4 ], [ 8, 13 ],
#      [ 10, 11 ] ] ]

InstallGlobalFunction( GraphAdjunction, function (g,l)
    local len,j,k,u,h;
    len := Length(l[1]);
    u := Union(l[2]);
    h := GraphUnion(g,NullGraph(SymmetricGroup(Length(u))));
    h:=NewGroupGraph(Group(()),h);
    for j in [1..len] do
        for k in [1..2] do
            AddEdgeOrbit(h,[l[1][j],OrderGraph(g)+l[2][j][k]]);
        od;
    od;
    return UnderlyingGraph(h);
end);

# A family of (D,g) graphs with somewhat few vertices

InstallGlobalFunction( GraphCloseToDgCage, function (n,k)
    local i,j,g;
    g:=GraphUnion(CyclicGraph(n*k),NullGraph(SymmetricGroup(n)));
    for i in [1..n] do
        for j in [1..k] do
            AddEdgeOrbit(g,[n*k+i,i+n*(j-1)]);
        od;
    od;
    return UnderlyingGraph(g);
end);
