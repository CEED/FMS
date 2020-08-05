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
#include <math.h>

/* All of the coordinates defined in a global ordering. */
double global_coords[] = {
0., 0.,     /* 0 */
0.4, 0.,    /* 1 */
0.65, 0.,   /* 2 */
1., 0.,     /* 3 */
0., 0.25,   /* 4 */
0.25, 0.28, /* 5 */
0.65, 0.25, /* 6 */
1., 0.2,    /* 7 */
0., 0.55,   /* 8 */
0.16, 0.75, /* 9 */
0.5, 0.6,   /* 10 */
1., 0.6,    /* 11 */
0., 1.,     /* 12 */
0.36, 1.,   /* 13 */
0.65, 1.,   /* 14 */
1., 1.      /* 15 */
};

/* Edges */
int edge_vert[] = {
   /*vertical edges */
   0,4,   /* 0 */
   1,5,   /* 1 */
   2,6,   /* 2 */
   3,7,   /* 3 */
   4,8,   /* 4 */
   5,9,   /* 5 */
   6,10,  /* 6 */
   7,11,  /* 7 */
   8,12,  /* 8 */
   9,13,  /* 9 */
   10,14, /* 10 */
   11,15, /* 11 */
   /* horizontal edges */
   0,1,   /* 12 */
   1,2,   /* 13 */
   2,3,   /* 14 */
   4,5,   /* 15 */
   5,6,   /* 16 */
   6,7,   /* 17 */
   8,9,   /* 18 */
   9,10,  /* 19 */
   10,11, /* 20 */
   12,13, /* 21 */
   13,14, /* 22 */
   14,15, /* 23 */
};

/* Quads */
int quad_edge[] = {
    0,12,1,15,  /* 0 */
    1,13,2,16,  /* 1 */
    2,14,3,17,  /* 2 */
    4,15,5,18,  /* 3 */
    5,16,6,19,  /* 4 */
    6,17,7,20,  /* 5 */
    8,18,9,21,  /* 6 */
    9,19,10,22, /* 7 */
    10,20,11,23 /* 8 */
};

/* Adds coordinate dofs at vertices */
double *
append_coords_vert(double *ptr, const double *coords, int nverts)
{
    int i;
    for(i = 0; i < nverts; i++)
    {
        *ptr++ = coords[2*i];
        *ptr++ = coords[2*i+1];
    }
    return ptr;
}

/* Adds coordinate dofs along edges */
double *
append_coords_edge(int order, double *ptr, const double *coords, 
    const int *edges, int nedges, double travel)
{
    int i, j, v0, v1;
    double t, xoffset, yoffset, pt[2];

    if(order > 1)
    {
        for(i = 0; i < nedges; ++i)
        {
            /* Get the i'th edge so we know its vertices. */
            v0 = edges[2*i];
            v1 = edges[2*i+1];

            for(j = 0; j < order-1; ++j)
            {
                /* Figure out how much to perturb the point. */
                switch((i+j) % 8)
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

                t = ((double)(j+1)) / ((double)order);

                /* Compute a point along the edge and perturb it. */
                pt[0] = ((t) * coords[v0*2+0]) + ((1.-t) * coords[v1*2+0]) + xoffset;
                pt[1] = ((t) * coords[v0*2+1]) + ((1.-t) * coords[v1*2+1]) + yoffset;

                *ptr++ = pt[0];
                *ptr++ = pt[1];
            }
        }
    }

    return ptr;
}

/* Adds coordinate dofs on quad faces */
double *
append_coords_face(int order, double *ptr, const double *coords,
    const int *edges, const int *quads, int nquads, double travel)
{
    int i,j,edge,v0,v1;
    double pt[2], xoffset, yoffset;
    if(order >= 2)
    {
        for(i = 0; i < nquads; ++i)
        {
            /* NOTE: we assume edges point the same way. */
            int e0 = quads[i*4+0];
            int e2 = quads[i*4+2];

            /* 3012 seems to be the ordering to get interior dofs 
               in the right order. */
            int v0 = edges[e0*2+1];
            int v1 = edges[e0*2+0];
            int v2 = edges[e2*2+0];
            int v3 = edges[e2*2+1];

            /* Make interior points by blending the quad vertices */
            for(int jj = 1; jj < order; ++jj)
            for(int ii = 1; ii < order; ++ii)
            {
                float r = ((float)ii)/((float)order);
                float s = ((float)jj)/((float)order);

                pt[0] = (1.-r)*(1.-s)*coords[v0*2+0] + 
                        r     *(1.-s)*coords[v1*2+0] + 
                        r     *s     *coords[v2*2+0] + 
                        (1.-r)*s     *coords[v3*2+0];
                pt[1] = (1.-r)*(1.-s)*coords[v0*2+1] + 
                        r     *(1.-s)*coords[v1*2+1] + 
                        r     *s     *coords[v2*2+1] + 
                        (1.-r)*s     *coords[v3*2+1];

                switch((i+ii+jj)%4)
                {
                case 0: xoffset = -travel; yoffset = 0.;      break;
                case 1: xoffset = 0.;      yoffset = -travel; break;
                case 2: xoffset = travel;  yoffset = 0.;      break;
                case 3: xoffset = 0.;      yoffset = travel;  break;
                }

                /* perturb the point. */
                pt[0] += xoffset;
                pt[1] += yoffset;

                *ptr++ = pt[0];
                *ptr++ = pt[1];
            }
        }
    }
    return ptr;
}

/**
Determines the number of dofs on quad cells.
*/
int dofs_for_order(int order, int nverts, int nedges, int nquads)
{
    return nverts + (order-1)*nedges + (order-1)*(order-1)*nquads;
}

/**
Makes a radial scalar raised to a power.
@note coords are the order 3 coordinates.
*/
void
add_radial_scalar(FmsDataCollection dc, FmsComponent surface, 
    const char *name, const double *gc, int order, int nverts, int nedges, int nquads)
{
    int i, ndofs;
    float *data;
    double *coords, *cptr, x, y;
    char fdname[100];

    /* Create a field descriptor for the coordinates that live on "surface". */
    sprintf(fdname, "%s descriptor", name);
    FmsFieldDescriptor fd;
    FmsDataCollectionAddFieldDescriptor(dc, fdname, &fd);
    FmsFieldDescriptorSetComponent(fd, surface);
    FmsFieldDescriptorSetFixedOrder(fd, FMS_CONTINUOUS,
                                    FMS_NODAL_GAUSS_CLOSED, order);

    /* First, make coordinates of the target order. */
    ndofs = dofs_for_order(order, nverts, nedges, nquads);
    coords = (double *)malloc(ndofs * sizeof(double)*2);
    cptr = coords;
    cptr = append_coords_vert(cptr, gc, nverts);
    cptr = append_coords_edge(order, cptr, gc, edge_vert, nedges, 0.);
    cptr = append_coords_face(order, cptr, gc, edge_vert, quad_edge, nquads, 0.);

    /* Make the radial data. */
    data = (float *)malloc(ndofs * sizeof(float));
    for(i = 0; i < ndofs; ++i)
    {
        x = coords[i*2];
        y = coords[i*2+1];
        data[i] = (float)sqrt(x*x + y*y);
    }

#if 0
    /* Write coordinate dofs as Point3D so we can plot them. */
    char tmp[100];
    sprintf(tmp, "%s.3D", name);
    FILE *fp = fopen(tmp, "wt");
    fprintf(fp, "X Y Z dof\n");
    for(i = 0; i < ndofs; ++i)
        fprintf(fp, "%lg %lg 0. %f\n", coords[i*2],coords[i*2+1], data[i]);
    fclose(fp);
#endif

    /* Add the coordinates to the data collection. */
    FmsField f;
    FmsDataCollectionAddField(dc, name, &f);
    FmsFieldSet(f, fd, 1, FMS_BY_VDIM, FMS_FLOAT, data);

    /* Add some field metadata. (optional) */
    FmsMetaData mdata;
    FmsFieldAttachMetaData(f, &mdata);
    FmsMetaDataSetString(mdata, "units", "cm");

    free(data);
    free(coords);
}

void
add_doforder_scalar(FmsDataCollection dc, FmsComponent surface, 
    const char *name, int nverts, int nedges, int nquads)
{
    int i, order = 3, ndofs, nedge_data, nface_data;
    float *data, *ptr;
    double x, y;
    char fdname[100];

    /* Create a field descriptor for the coordinates that live on "surface". */
    sprintf(fdname, "%s descriptor", name);
    FmsFieldDescriptor fd;
    FmsDataCollectionAddFieldDescriptor(dc, fdname, &fd);
    FmsFieldDescriptorSetComponent(fd, surface);
    FmsFieldDescriptorSetFixedOrder(fd, FMS_CONTINUOUS,
                                    FMS_NODAL_GAUSS_CLOSED, order);

    /* Make the data. */
    nedge_data = (order-1)*nedges;
    nface_data = (order-1)*(order-1)*nquads;
    ndofs = dofs_for_order(order, nverts, nedges, nquads);
    data = (float *)malloc(ndofs * sizeof(float));
    ptr = data;
    for(i = 0; i < nverts; i++)
        *ptr++ = 1.f;
    for(i = 0; i < nedge_data; i++)
        *ptr++ = 2.f;
    for(i = 0; i < nface_data; i++)
        *ptr++ = 3.f;

    /* Add the coordinates to the data collection. */
    FmsField f;
    FmsDataCollectionAddField(dc, name, &f);
    FmsFieldSet(f, fd, 1, FMS_BY_VDIM, FMS_FLOAT, data);

    /* Add some field metadata. (optional) */
    FmsMetaData mdata;
    FmsFieldAttachMetaData(f, &mdata);
    FmsMetaDataSetString(mdata, "units", "dof order");

    free(data);
}

int
main(int argc, char *argv[])
{
    int i, nverts, nedges, nquads;
    const char *protocol = "ascii";
    FmsInt order = 1;
    if(argc > 1)
    {
        protocol = argv[1];
        if(argc > 2)
        {
            order = atoi(argv[2]);
            if(order < 1 || order > 4)
                return -1;
        }
    }

    /* Create a mesh */
    FmsMesh mesh;
    FmsMeshConstruct(&mesh);

    /* Make a domain. */
    FmsDomain *domains;
    FmsMeshAddDomains(mesh, "domain", 1, &domains);

    /* Set the number of vertices */
    nverts = (sizeof(global_coords)/sizeof(double))/2;
    FmsDomainSetNumVertices(domains[0], nverts);

    /* Set the edges */
    nedges = sizeof(edge_vert)/sizeof(int)/2;
    FmsDomainSetNumEntities(domains[0], FMS_EDGE, FMS_INT32, nedges);
    FmsDomainAddEntities(domains[0], FMS_EDGE, NULL, FMS_INT32, edge_vert, nedges);

    /* Set the quads */
    nquads = sizeof(quad_edge)/sizeof(int)/4;
    FmsDomainSetNumEntities(domains[0], FMS_QUADRILATERAL, FMS_INT32, nquads);
    FmsDomainAddEntities(domains[0], FMS_QUADRILATERAL, NULL, FMS_INT32, quad_edge, nquads);

    /* Make a component called surface. */
    FmsComponent surface;
    FmsMeshAddComponent(mesh, "surface", &surface);
    FmsComponentAddDomain(surface, domains[0]);

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

    /* Make the coordinate data. */
    int ndofs = dofs_for_order(order, nverts, nedges, nquads);
    double *coord_data = (double *)malloc(ndofs * sizeof(double)*2);
    double *cptr = coord_data;
    cptr = append_coords_vert(cptr, global_coords, nverts);
    cptr = append_coords_edge(order, cptr, global_coords, edge_vert, nedges, 0.02);
    cptr = append_coords_face(order, cptr, global_coords, edge_vert, quad_edge, nquads, 0.01);

    /* Add the coordinates to the data collection. */
    FmsField coords;
    FmsDataCollectionAddField(dc, "coords", &coords);
    FmsFieldSet(coords, coords_fd, 2, FMS_BY_VDIM, FMS_DOUBLE, coord_data);

    /* Add some field metadata. (optional) */
    FmsMetaData coords_mdata;
    FmsFieldAttachMetaData(coords, &coords_mdata);
    FmsMetaDataSetString(coords_mdata, "units", "cm");

    /* Set the coordinates for the component. */
    FmsComponentSetCoordinates(surface, coords);

    /* Add some scalars based on different orders of coordinates. */
    add_radial_scalar(dc, surface, "r1", global_coords, 1, nverts, nedges, nquads);
    add_radial_scalar(dc, surface, "r2", global_coords, 2, nverts, nedges, nquads);
    add_radial_scalar(dc, surface, "r3", global_coords, 3, nverts, nedges, nquads);
    add_radial_scalar(dc, surface, "r4", global_coords, 4, nverts, nedges, nquads);

    /* Add an order 3 scalar that includes the dof order. */
    add_doforder_scalar(dc, surface, "doforder", nverts, nedges, nquads);

    /* Add some metadata. (optional) */
    FmsMetaData mdata;
    int cycle = 1234;
    double dtime = 1.2345678;
    FmsDataCollectionAttachMetaData(dc, &mdata);
    FmsMetaData *mdata_objects;
    FmsMetaDataSetMetaData(mdata, "Info", 3, &mdata_objects);
    FmsMetaDataSetIntegers(mdata_objects[0], "cycle", FMS_INT32, 1, (void**)&cycle);
    FmsMetaDataSetScalars(mdata_objects[1], "time", FMS_DOUBLE, 1, (void**)&dtime);
    FmsMetaDataSetString(mdata_objects[2], "Description", "Quad example file");

    /* Write the data out. */
    FmsIOWrite("quads.fms", protocol, dc);

    /* Write coordinate dofs as Point3D so we can plot them. */
    printf("X Y Z dof\n");
    for(i = 0; i < ndofs; ++i)
        printf("%lg %lg 0. %d\n", coord_data[i*2],coord_data[i*2+1], i);

    free(coord_data);

    return 0;
}
