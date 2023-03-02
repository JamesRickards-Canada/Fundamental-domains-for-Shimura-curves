print("\n\nType '?fdom' for help.\n\n");
addhelp(fdom, "This package can be used to compute fundamental domains for Shimura curves (the PARI code can be easily adapted to compute fundamental domains for any discrete subgroup of PSL(2, R)).\n For each subtopic ``P (p)'', call ?p to access a basic description and list of methods. Subtopics:\n Geometry (geo)\n Visualizing fundamental domains with Python (vfd)\n Quaternion methods (quat)");
parigp_version=version();
fdom_library=strprintf("./libfdom-%d-%d.so", parigp_version[1], parigp_version[2]);



default(parisize, "4096M");\\Must come at the end