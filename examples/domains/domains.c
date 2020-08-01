/*
 Copyright (c) 2020, Lawrence Livermore National Security, LLC. Produced at
 the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
 reserved. See files LICENSE and NOTICE for details.

 This file is part of CEED, a collection of benchmarks, miniapps, software
 libraries and APIs for efficient high-order finite element and spectral
 element discretizations for exascale applications. For more information and
 source code availability see http://github.com/ceed.

 The CEED research is supported by the Exascale Computing Project (17-SC-20-SC)
 a collaborative effort of two U.S. Department of Energy organizations (Office
 of Science and the National Nuclear Security Administration) responsible for
 the planning and preparation of a capable exascale ecosystem, including
 software, applications, hardware, advanced system engineering and early
 testbed platforms, in support of the nation's exascale computing imperative.
*/
#include <fms.h>
#include <fmsio.h>
#include <stdlib.h>
#include <stdio.h>

/* We'll make multiple domains and components/parts. */

/* All of the coordinates defined in a global ordering. */
double global_coords[] = {
0., 0.,
0.4, 0.,
0.65, 0.,
1., 0.,

0., 0.25,
0.25, 0.28,
0.65, 0.25,
1., 0.2,

0., 0.55,
0.16, 0.75,
0.5, 0.6,
1., 0.6,

0., 1.,
0.36, 1.,
0.65, 0.9,
1., 1.
};

/* Domain local index to global index */
int d0_global_index[] = {4,5,10,8,9,12};
int d1_global_index[] = {0,1,2,4,5,6,10};
int d2_global_index[] = {2,3,6,7,10,11,14,15};
int d3_global_index[] = {10,9,14,12,13,15};

/* Domain edges expressed in their local vertex ordering */
int d0_edge_vert[] = {
    0,1, /* 0 */
    0,3, /* 1 */
    1,3, /* 2 */
    1,4, /* 3 */
    1,2, /* 4 */
    2,4, /* 5 */
    3,4, /* 6 */
    4,5, /* 7 */
    3,5, /* 8 */
};
int d1_edge_vert[] = {
    0,3, /* 0 */
    0,1, /* 1 */
    1,3, /* 2 */
    1,4, /* 3 */
    1,5, /* 4 */
    1,2, /* 5 */
    2,5, /* 6 */
    3,4, /* 7 */
    4,5, /* 8 */
    5,6, /* 9 */
    6,4  /* 10 */
};
int d2_edge_vert[] = {
    0,2, /* 0 */
    0,3, /* 1 */
    0,1, /* 2 */
    1,3, /* 3 */
    2,3, /* 4 */
    2,4, /* 5 */
    3,4, /* 6 */
    3,5, /* 7 */
    4,5, /* 8 */
    4,6, /* 9 */
    5,6, /* 10 */
    6,7, /* 11 */
    5,7  /* 12 */
};
int d3_edge_vert[] = {
    0,1, /* 0 */
    0,2, /* 1 */
    1,2, /* 2 */
    1,3, /* 3 */
    1,4, /* 4 */
    2,4, /* 5 */
    2,5, /* 6 */
    3,4, /* 7 */
    4,5  /* 8 */
};

/* Domain triangles expressed in their local edge ordering */
int d0_tri_edge[] = {
    0,2,1, /* 0 */
    2,3,6, /* 1 */
    4,5,3, /* 2 */
    6,7,8  /* 3 */
};
int d1_tri_edge[] = {
    0,1,2, /* 0 */
    2,3,7, /* 1 */
    3,4,8, /* 2 */
    4,5,6, /* 3 */
    8,9,10 /* 4 */
};
int d2_tri_edge[] = {
    0,1,4, /* 0 */
    1,2,3, /* 1 */
    5,4,6, /* 2 */
    6,7,8, /* 3 */
    8,10,9,  /* 4 */
    10,12,11 /* 5 */
};
int d3_tri_edge[] = {
    0,1,2, /* 0 */
    3,4,7, /* 1 */
    4,2,5, /* 2 */
    5,6,8  /* 3 */
};

double *
append_domain_coords(int dom, double *ptr, const int *global_index, int nverts)
{
    int i;
    double pt[2];
    for(i = 0; i < nverts; i++)
    {
        int gidx = global_index[i];
        pt[0] = global_coords[2*gidx];
        pt[1] = global_coords[2*gidx+1];

        *ptr++ = pt[0];
        *ptr++ = pt[1];

        /* Print value to Point3D output. */
        printf("%lg %lg 0. %d\n", pt[0],pt[1], dom);
    }
    return ptr;
}

double *
append_domain_coords_edge(int dom, double *ptr, const int *global_index, int verts,
    const int *edges, int nedges)
{
    int i, v0, v1;
    double xoffset, yoffset, pt[2], travel = 0.02;

    for(i = 0; i < nedges; ++i)
    {
        /* Get the i'th edge and get the global indices of its vertices. */
        v0 = global_index[edges[2*i]];
        v1 = global_index[edges[2*i+1]];

        /* Figure out how much to perturb the point. */
        switch((v0+v1) % 8)
        {
        case 0: xoffset = -travel; yoffset = 0;       break;
        case 1: xoffset = -travel; yoffset = -travel; break;
        case 2: xoffset = 0.;      yoffset = -travel; break;
        case 3: xoffset = travel;  yoffset = -travel; break;
        case 4: xoffset = travel;  yoffset = 0;       break;
        case 5: xoffset = travel;  yoffset = travel;  break;
        case 6: xoffset = 0;       yoffset = travel;  break;
        case 7: xoffset = -travel; yoffset = travel;  break;
        }

        /* Compute the edge midpoint and perturb it. */
        pt[0] = 0.5 * (global_coords[v0*2+0] + global_coords[v1*2+0]) + xoffset;
        pt[1] = 0.5 * (global_coords[v0*2+1] + global_coords[v1*2+1]) + yoffset;

        *ptr++ = pt[0];
        *ptr++ = pt[1];

        /* Print value to Point3D output. */
        printf("%lg %lg 0. %d\n", pt[0],pt[1], dom);
    }

    return ptr;
}

int
main(int argc, char *argv[])
{
    const char *protocol = "ascii";
    FmsInt order = 1;
    if(argc > 1)
    {
        protocol = argv[1];
        if(argc > 2)
            order = atoi(argv[2]);
    }

    /* Create a mesh */
    FmsMesh mesh;
    FmsMeshConstruct(&mesh);

    /* Make 4 domains. */
    FmsDomain *domains;
    FmsMeshAddDomains(mesh, "domain", 4, &domains);

    /* Set the number of vertices for each domain */
    int nverts[4];
    nverts[0] = sizeof(d0_global_index)/sizeof(int);
    nverts[1] = sizeof(d1_global_index)/sizeof(int);
    nverts[2] = sizeof(d2_global_index)/sizeof(int);
    nverts[3] = sizeof(d3_global_index)/sizeof(int);
    FmsDomainSetNumVertices(domains[0], nverts[0]);
    FmsDomainSetNumVertices(domains[1], nverts[1]);
    FmsDomainSetNumVertices(domains[2], nverts[2]);
    FmsDomainSetNumVertices(domains[3], nverts[3]);

    /* Set the edges for each domain */
    int nedges[4];
    nedges[0] = sizeof(d0_edge_vert)/sizeof(int)/2;
    nedges[1] = sizeof(d1_edge_vert)/sizeof(int)/2;
    nedges[2] = sizeof(d2_edge_vert)/sizeof(int)/2;
    nedges[3] = sizeof(d3_edge_vert)/sizeof(int)/2;
    FmsDomainSetNumEntities(domains[0], FMS_EDGE, FMS_INT32, nedges[0]);
    FmsDomainSetNumEntities(domains[1], FMS_EDGE, FMS_INT32, nedges[1]);
    FmsDomainSetNumEntities(domains[2], FMS_EDGE, FMS_INT32, nedges[2]);
    FmsDomainSetNumEntities(domains[3], FMS_EDGE, FMS_INT32, nedges[3]);
    FmsDomainAddEntities(domains[0], FMS_EDGE, NULL, FMS_INT32, d0_edge_vert, nedges[0]);
    FmsDomainAddEntities(domains[1], FMS_EDGE, NULL, FMS_INT32, d1_edge_vert, nedges[1]);
    FmsDomainAddEntities(domains[2], FMS_EDGE, NULL, FMS_INT32, d2_edge_vert, nedges[2]);
    FmsDomainAddEntities(domains[3], FMS_EDGE, NULL, FMS_INT32, d3_edge_vert, nedges[3]);

    /* Set the triangles for each domain */
    int ntri[4];
    ntri[0] = sizeof(d0_tri_edge)/sizeof(int)/3;
    ntri[1] = sizeof(d1_tri_edge)/sizeof(int)/3;
    ntri[2] = sizeof(d2_tri_edge)/sizeof(int)/3;
    ntri[3] = sizeof(d3_tri_edge)/sizeof(int)/3;
    FmsDomainSetNumEntities(domains[0], FMS_TRIANGLE, FMS_INT32, ntri[0]);
    FmsDomainSetNumEntities(domains[1], FMS_TRIANGLE, FMS_INT32, ntri[1]);
    FmsDomainSetNumEntities(domains[2], FMS_TRIANGLE, FMS_INT32, ntri[2]);
    FmsDomainSetNumEntities(domains[3], FMS_TRIANGLE, FMS_INT32, ntri[3]);
    FmsDomainAddEntities(domains[0], FMS_TRIANGLE, NULL, FMS_INT32, d0_tri_edge, ntri[0]);
    FmsDomainAddEntities(domains[1], FMS_TRIANGLE, NULL, FMS_INT32, d1_tri_edge, ntri[1]);
    FmsDomainAddEntities(domains[2], FMS_TRIANGLE, NULL, FMS_INT32, d2_tri_edge, ntri[2]);
    FmsDomainAddEntities(domains[3], FMS_TRIANGLE, NULL, FMS_INT32, d3_tri_edge, ntri[3]);

    /* Make a component called surface that comprises all domains. */
    FmsComponent surface;
    FmsMeshAddComponent(mesh, "surface", &surface);
    FmsComponentAddDomain(surface, domains[0]);  /* adds all cells in domain to a part? */
    FmsComponentAddDomain(surface, domains[1]);
    FmsComponentAddDomain(surface, domains[2]);
    FmsComponentAddDomain(surface, domains[3]);

    FmsMeshFinalize(mesh);
    FmsMeshValidate(mesh);

    /* Make a data collection to contain the mesh. */
    FmsDataCollection dc;
    FmsDataCollectionCreate(mesh, "data collection", &dc);

    /* Create a field descriptor for the coordinates that live on "surface". */
    FmsFieldDescriptor coords_fd;
    FmsDataCollectionAddFieldDescriptor(dc, "coords descriptor", &coords_fd);
    FmsFieldDescriptorSetComponent(coords_fd, surface); /* define coords over surface */
    FmsFieldDescriptorSetFixedOrder(coords_fd, FMS_CONTINUOUS,
                                    FMS_NODAL_GAUSS_CLOSED, order);

    /* Since all domains were added to the "surface" component and the field 
       is defined over the "surface" component, we make a new field 
       that contains values from all domains, with their local vertex ordering.
     */
    int ndofs[4] = {0,0,0,0};
    if(order == 1)
    {
        ndofs[0] = nverts[0];
        ndofs[1] = nverts[1];
        ndofs[2] = nverts[2];
        ndofs[3] = nverts[3];
    }
    else if(order == 2)
    {
        ndofs[0] = nverts[0] + nedges[0];
        ndofs[1] = nverts[1] + nedges[1];
        ndofs[2] = nverts[2] + nedges[2];
        ndofs[3] = nverts[3] + nedges[3];
    }        
    int ntotal_dofs = ndofs[0] + ndofs[1] + ndofs[2] + ndofs[3];
    double *coord_data = (double *)malloc(ntotal_dofs * sizeof(double)*2);
    double *cptr = coord_data;

    /* Write a Point3D header since we will write coord values to stdout in a plottable way. */
    printf("X Y Z domain\n");

    /* Write the coordinate data. */
    cptr = append_domain_coords(0, cptr, d0_global_index, nverts[0]);
    if(order == 2)
        cptr = append_domain_coords_edge(0, cptr, d0_global_index, nverts[0], d0_edge_vert, nedges[0]);
    cptr = append_domain_coords(1, cptr, d1_global_index, nverts[1]);
    if(order == 2)
        cptr = append_domain_coords_edge(1, cptr, d1_global_index, nverts[1], d1_edge_vert, nedges[1]);
    cptr = append_domain_coords(2, cptr, d2_global_index, nverts[2]);
    if(order == 2)
        cptr = append_domain_coords_edge(2, cptr, d2_global_index, nverts[2], d2_edge_vert, nedges[2]);
    cptr = append_domain_coords(3, cptr, d3_global_index, nverts[3]);
    if(order == 2)
        cptr = append_domain_coords_edge(3, cptr, d3_global_index, nverts[3], d3_edge_vert, nedges[3]);

    /*Q: suppose that the coordinate order was 2+, would the extra dofs go here 
         or would they appear after their individual domain loops above?
     */
    FmsField coords;
    FmsDataCollectionAddField(dc, "coords", &coords);
    FmsFieldSet(coords, coords_fd, 2, FMS_BY_VDIM, FMS_DOUBLE, coord_data);

    /* Set the coordinates for the component. */
    FmsComponentSetCoordinates(surface, coords);

    /* Write the data out. */
    FmsIOWrite("domains.fms", protocol, dc);

    return 0;
}
