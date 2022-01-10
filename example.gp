F=nfinit(y^2-11);\\The field Q(sqrt(11))
A=alginit(F, [y-3, -3]);\\Quaternion algebra (y-3,-3/F), which is split at one real place and has a>0 at that split place
area=algfdomarea(A);\\Computation of the area of the fundamental domain
dom=algfdom(A);\\Computation of the fundamental domain
elt=algmulvec(A, dom[1], [1,2,7,1,2,2,1,5,7]);\\Generating an element of the group (happens to be hyperbolic)
geod=algfdomrootgeodesic(elt, dom);\\Finding the geodesic in the fundamental domain
python_printfdom(dom, "fdom1");\\When printing fundamental domains, you must start the file name with "fd"
python_printarcs(geod[2], "geod1", 1, "fdom1");\\View the domain and the geodesic, ON WINDOWS SUBSYSTEM FOR LINUX ONLY (I think). On linux, you need to call "py geod1 fdom1" from the command line