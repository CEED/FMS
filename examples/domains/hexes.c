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

/*#define DEBUG_PRINT 1*/

#define ADD_CELLS 1

/**
Determines the number of dofs on hex cells.
*/
int dofs_for_order(int order, int nverts, int nedges, int nfaces, int ncells)
{
    int O1 = order-1;
    return nverts + O1*nedges + O1*O1*nfaces + O1*O1*O1*ncells;
}

/**
Computes a 3D grid of vertices in a usual mesh vertex ordering
and all points get perturbed a little so we do not have a totally
regular grid.
*/
double *
compute_vertices(const int *dims, int *nverts)
{
    int i,j,k;
    double x,y,z,*coords, *ptr, px, py, pz;

    px = 0.1 / ((double)(dims[0]-1));
    py = 0.1 / ((double)(dims[1]-1));
    pz = 0.1 / ((double)(dims[2]-1));

    *nverts = dims[0]*dims[1]*dims[2];
    coords = ptr = (double *)malloc(*nverts * sizeof(double)*3);
    for(k = 0; k < dims[2]; k++)
    {
        double tz = ((double)k) / (double)(dims[2]-1);
        for(j = 0; j < dims[1]; j++)
        {
            double ty = ((double)j) / (double)(dims[1]-1);
            for(i = 0; i < dims[0]; i++)
            {
                x = ((double)i) / (double)(dims[0]-1);
                y = ty;
                z = tz;

                /* perturb points */
                x += px * cos(drand48() * 2 * M_PI);
                y += py * cos(drand48() * 2 * M_PI);
                z += pz * cos(drand48() * 2 * M_PI);

                *ptr++ = x;
                *ptr++ = y;
                *ptr++ = z;
            }
        }
    }
    return coords;
}

/**
Computes the edges that we will link into faces.
*/
int *
compute_edges(const int *dims, int *nedges)
{
    int i,j,k;
    int nxny = dims[0]*dims[1];
    int E1 = (dims[0]-1) * (dims[1])   * (dims[2]);
    int E2 = (dims[0])   * (dims[1]-1) * (dims[2]);
    int E3 = (dims[0])   * (dims[1])   * (dims[2]-1);
    *nedges = E1 + E2 + E3;
    int *e = (int *)malloc(*nedges * sizeof(int)*2);
    int *eptr = e;
 
    /* Make horizontal edges in XY planes. */
    for(k = 0; k < dims[2]; k++)
    for(j = 0; j < dims[1]; j++)
    {
        int knxy_plus_jnx = k*nxny + j*dims[0];
        for(i = 0; i < dims[0]-1; i++)
        {
            *eptr++ = knxy_plus_jnx + i;
            *eptr++ = knxy_plus_jnx + i + 1;
        }
    }

    /* Make vertical edges in XY planes. */
    for(k = 0; k < dims[2]; k++)
    for(j = 0; j < dims[1]-1; j++)
    {
        int knxy_plus_jnx = k*nxny + j*dims[0];
        int knxy_plus_j1nx = k*nxny + (j+1)*dims[0];
        for(i = 0; i < dims[0]; i++)
        {
            *eptr++ = knxy_plus_jnx + i;
            *eptr++ = knxy_plus_j1nx + i;
        }
    }

    /* Make Z edges. */
    for(k = 0; k < dims[2]-1; k++)
    for(j = 0; j < dims[1]; j++)
    {
        int knxy_plus_jnx = k*nxny + j*dims[0];
        int k1nxy_plus_jnx = (k+1)*nxny + j*dims[0];
        for(i = 0; i < dims[0]; i++)
        {
            *eptr++ = knxy_plus_jnx + i;
            *eptr++ = k1nxy_plus_jnx + i;
        }
    }
#ifdef DEBUG_PRINT
    for(i = 0; i < *nedges; i++)
        printf("edge %d: %d %d\n", i, e[2*i], e[2*i+1]);
#endif
    return e;
}

/**
Computes the faces that we'll use to make hexes. The ordering
of the edges within has been modified a bit to make dof ordering
for orders 3+ look right on the faces.
*/
int *
compute_faces(const int *dims, int *nfaces)
{
    int i,j,k,A,B,C,D, E1,E2,E3,F1,F2,F3, *faces, *f;

    /* number of edges. */
    E1 = (dims[0]-1) * (dims[1])   * (dims[2]);
    E2 = (dims[0])   * (dims[1]-1) * (dims[2]);
    E3 = (dims[0])   * (dims[1])   * (dims[2]-1);

    F1 = (dims[0]-1)*(dims[1]-1)*dims[2];     /* back/front */
    F2 = (dims[0]-1)*(dims[1])  *(dims[2]-1); /* bottom/top */
    F3 = (dims[0])  *(dims[1]-1)*(dims[2]-1); /* left/right */

    *nfaces = F1 + F2 + F3;

    faces = (int *)malloc(*nfaces * sizeof(int)*4);
    f = faces;

    /* back/front faces */
    for(k = 0; k < dims[2]; k++)
    for(j = 0; j < dims[1]-1; j++)
    for(i = 0; i < dims[0]-1; i++)
    {
        /*        C
              *------*
            D |      | B
              |      |
              *------*
                 A
         */

        /* left/right edges */
        A = k*(dims[0]-1)*(dims[1])  +j*(dims[0]-1)     + i;
        C = k*(dims[0]-1)*(dims[1])  +(j+1)*(dims[0]-1) + i;
        /* vertical edges */
        D = E1 + k*dims[0]*(dims[1]-1)    +j*dims[0]    + i;
        B = E1 + k*dims[0]*(dims[1]-1)    +j*dims[0]    + i + 1;       

        /* ADCB */
        *f++ = A;
        *f++ = D;
        *f++ = C;
        *f++ = B;
    }

    /* bottom/top faces */
    for(k = 0; k < dims[2]-1; k++)
    for(j = 0; j < dims[1]; j++)
    for(i = 0; i < dims[0]-1; i++)
    {
        /*          A
                *------*
            B  /      / D   (top view)
              /      /
             *------*
                 C
         */
        /* left/right edges */
        A = k*(dims[0]-1)*(dims[1])     +j*(dims[0]-1) + i;
        C = (k+1)*(dims[0]-1)*(dims[1]) +j*(dims[0]-1) + i;
        /* z edges */
        D = E1+E2 + k*dims[0]*dims[1]    +j*dims[0] + i + 1;
        B = E1+E2 + k*dims[0]*dims[1]    +j*dims[0] + i;       

        /* This edge order works well for order 3+ dof ordering */
        *f++ = D;
        *f++ = A;
        *f++ = B;
        *f++ = C;
    }

    /* left/right faces */
    for(k = 0; k < dims[2]-1; k++)
    for(j = 0; j < dims[1]-1; j++)
    for(i = 0; i < dims[0]; i++)
    {
        /*
                *
             C /|
              / | B
             *  | 
           D |  *
             | /
             |/ A
             *
         */
        /* z edges */
        A = E1+E2 + k*dims[0]*dims[1]    +j*dims[0] + i;       
        C = E1+E2 + k*dims[0]*dims[1]    +(j+1)*dims[0] + i;

        /* vertical edges */
        D = E1 + (k+1)*dims[0]*(dims[1]-1)  +j*dims[0] + i;
        B = E1 + k*dims[0]*(dims[1]-1)      +j*dims[0] + i;       

        /* ABCD */
        *f++ = A;
        *f++ = B;
        *f++ = C;
        *f++ = D;
    }

#ifdef DEBUG_PRINT
    printf("nfaces=%d\n", *nfaces);

    for(i = 0; i < *nfaces; i++)
        printf("face %d: %d %d %d %d\n", i, faces[4*i], faces[4*i+1], faces[4*i+2], faces[4*i+3]);
#endif

    return faces;
}

/**
Compute the hex cells from sets of faces.
*/
int *
compute_cells(const int *dims, int *ncells)
{
    int i,j,k;
    *ncells = (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    int *cells = (int *)malloc(*ncells * sizeof(int)*6);
    int *cptr = cells;

    int F1 = (dims[0]-1)*(dims[1]-1)*dims[2];     /* back/front */
    int F2 = (dims[0]-1)*(dims[1])  *(dims[2]-1); /* bottom/top */
    int F3 = (dims[0])  *(dims[1]-1)*(dims[2]-1); /* left/right */

    for(k = 0; k < dims[2]-1; k++)
    for(j = 0; j < dims[1]-1; j++)
    for(i = 0; i < dims[0]-1; i++)
    {
        int f0 = F1 + k*(dims[0]-1)*dims[1] + j*(dims[0]-1) + i;     /* bottom */
        int f1 = F1 + k*(dims[0]-1)*dims[1] + (j+1)*(dims[0]-1) + i; /* top */
        int f2 = (k+1)*(dims[0]-1)*(dims[1]-1) + j*(dims[0]-1) + i;  /* front */
        int f3 = k*(dims[0]-1)*(dims[1]-1) + j*(dims[0]-1) + i;      /* back  */
        int f4 = F1+F2 + k*dims[0]*(dims[1]-1) + j*dims[0] + i;      /* left */
        int f5 = F1+F2 + k*dims[0]*(dims[1]-1) + j*dims[0] + i+1;    /* right */

        /* This ordering works (see fms.h) */
        *cptr++ = f0;
        *cptr++ = f1;
        *cptr++ = f5;
        *cptr++ = f4;
        *cptr++ = f2;
        *cptr++ = f3;
    }
#ifdef DEBUG_PRINT
    for(i = 0; i < *ncells; i++)
        printf("cell %d: %d %d %d %d %d %d\n", i,
           cells[6*i+0],
           cells[6*i+1],
           cells[6*i+2],
           cells[6*i+3],
           cells[6*i+4],
           cells[6*i+5]
        );
#endif

    return cells;
}

/**
Adds the coordinate dofs for the vertices.
*/
double *
append_coord_dofs_verts(double *dest, const double *verts, int nverts, int component)
{
    int i;
    for(i = 0; i < nverts; ++i)
        dest[i] = verts[3*i + component];
    return dest + nverts;
}

/**
Adds coordinate dofs along edges.
*/
double *
append_coord_dofs_edge(int order, double *dest, const double *verts, 
    const int *edges, int nedges, double travel, int component)
{
    int i, j, v0, v1;
    double t, value, *ptr = dest;

    if(order > 1)
    {
        for(i = 0; i < nedges; ++i)
        {
            /* Get the i'th edge so we know its vertices. */
            v0 = edges[2*i];
            v1 = edges[2*i+1];

            for(j = 0; j < order-1; ++j)
            {
                t = ((double)(j+1)) / ((double)order);

                /* Compute a point along the edge and perturb it. */
                value = ((t) * verts[v0*3+component]) + ((1.-t) * verts[v1*3+component]);

                /* Perturb the point */
                value += travel * cos(drand48() * 2 * M_PI);

                *ptr++ = value;
            }
        }
    }

    return ptr;
}

/**
Adds coordinate dofs on quad faces.
*/
double *
append_coord_dofs_face(int order, double *dest, const double *verts,
    const int *edges, const int *faces, int nfaces, double travel,
    int component)
{
    int i,j,edge,v0,v1;
    double value, *ptr = dest;
    if(order >= 2)
    {
        for(i = 0; i < nfaces; ++i)
        {
            /* NOTE: we assume edges point the same way. */
            int e0 = faces[i*4+0];
            int e2 = faces[i*4+2];

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
                double r = ((double)ii)/((double)order);
                double s = ((double)jj)/((double)order);

                value = (1.-r)*(1.-s)*verts[v0*3+component] + 
                        r     *(1.-s)*verts[v1*3+component] + 
                        r     *s     *verts[v2*3+component] + 
                        (1.-r)*s     *verts[v3*3+component];

                /* Perturb the point */
                value += travel * cos(drand48() * 2 * M_PI);

                *ptr++ = value;
            }
        }
    }
    return ptr;
}

/* Adds coordinate dofs on cell interiors */
double *
append_coord_dofs_cell(int order, double *dest, const int *dims,
    const double *verts, double travel, int component)
{
    int i,j,k,ii,jj,kk,edge,v0,v1;
    double r,s,t,value, *ptr = dest;
#ifdef ADD_CELLS
    if(order >= 2)
    {
        /* NOTE: just do cells this way rather than tracing down through
                 faces/edges to get to vertices. */

        /* The permute array reorders the vertices into an order that
           will emit interior dofs in an order compatible with FMS/MFEM.
           This array is VERY IMPORTANT for getting dofs for orders 3+
           right.
         */
        int permute[] = {3,1,7,5,2,0,6,4};

        for(k = 0; k < dims[2]-1; ++k)
        for(j = 0; j < dims[1]-1; ++j)
        for(i = 0; i < dims[0]-1; ++i)
        {
            int v[8];
            v[permute[0]] = k*dims[0]*dims[1] + j*dims[0] + i;
            v[permute[1]] = k*dims[0]*dims[1] + j*dims[0] + i+1;
            v[permute[2]] = k*dims[0]*dims[1] + (j+1)*dims[0] + i;
            v[permute[3]] = k*dims[0]*dims[1] + (j+1)*dims[0] + i+1;
            v[permute[4]] = (k+1)*dims[0]*dims[1] + j*dims[0] + i;
            v[permute[5]] = (k+1)*dims[0]*dims[1] + j*dims[0] + i+1;
            v[permute[6]] = (k+1)*dims[0]*dims[1] + (j+1)*dims[0] + i;
            v[permute[7]] = (k+1)*dims[0]*dims[1] + (j+1)*dims[0] + i+1;

            /* Make interior points by blending the hex vertices */
            for(kk = 1; kk < order; ++kk)
            for(jj = 1; jj < order; ++jj)
            for(ii = 1; ii < order; ++ii)
            {
                r = ((double)ii)/((double)order);
                s = ((double)jj)/((double)order);
                t = ((double)kk)/((double)order);

                value = (1.-r)*(1.-s)*(1.-t)*verts[v[0]*3+component] + 
                        r     *(1.-s)*(1.-t)*verts[v[1]*3+component] + 
                        (1.-r)*s     *(1.-t)*verts[v[2]*3+component] + 
                        r     *s     *(1.-t)*verts[v[3]*3+component] +
                        (1.-r)*(1.-s)*t     *verts[v[4]*3+component] + 
                        r     *(1.-s)*t     *verts[v[5]*3+component] + 
                        (1.-r)*s     *t     *verts[v[6]*3+component] + 
                        r     *s     *t     *verts[v[7]*3+component];

                /* Perturb the point */
                value += travel * cos(drand48() * 2 * M_PI);

                *ptr++ = value;
            }
        }
    }
#endif
    return ptr;
}

/**
Makes a radial scalar.
*/
void
add_radial_scalar(FmsDataCollection dc, FmsComponent volume, 
    const char *name, int order, const int *dims,
    const double *verts, int nverts,
    const int *edges, int nedges, 
    const int *faces, int nfaces,
    int ncells)
{
    int i, ndofs;
    double *data, *coords, *cptr, *x, *y, *z;
    char fdname[100];

    /* Create a field descriptor for the coordinates that live on "surface". */
    sprintf(fdname, "%s descriptor", name);
    FmsFieldDescriptor fd;
    FmsDataCollectionAddFieldDescriptor(dc, fdname, &fd);
    FmsFieldDescriptorSetComponent(fd, volume);
    FmsFieldDescriptorSetFixedOrder(fd, FMS_CONTINUOUS,
                                    FMS_NODAL_GAUSS_CLOSED, order);

    /* First, make coordinates of the target order. */
    ndofs = dofs_for_order(order, nverts, nedges, nfaces, ncells);

    /* Make the coordinate data. Make it BYNODE xxxx...yyyy...zzzz... */
    coords = (double *)malloc(ndofs * sizeof(double)*3);
    x = coords;
    y = coords + ndofs;
    z = coords + 2*ndofs;

    /* Store coordinate dofs by node. */
    x = append_coord_dofs_verts(x, verts, nverts, 0);
    x = append_coord_dofs_edge(order, x, verts, edges, nedges, 0., 0);
    x = append_coord_dofs_face(order, x, verts, edges, faces, nfaces, 0., 0);
    x = append_coord_dofs_cell(order, x, dims, verts, 0., 0);

    y = append_coord_dofs_verts(y, verts, nverts, 1);
    y = append_coord_dofs_edge(order, y, verts, edges, nedges, 0., 1);
    y = append_coord_dofs_face(order, y, verts, edges, faces, nfaces, 0., 1);
    y = append_coord_dofs_cell(order, y, dims, verts, 0., 1);

    z = append_coord_dofs_verts(z, verts, nverts, 2);
    z = append_coord_dofs_edge(order, z, verts, edges, nedges, 0., 2);
    z = append_coord_dofs_face(order, z, verts, edges, faces, nfaces, 0., 2);
    z = append_coord_dofs_cell(order, z, dims, verts, 0., 2);

    /* Make the radial data. */
    x = coords;
    y = coords + ndofs;
    z = coords + 2*ndofs;
    data = (double *)malloc(ndofs * sizeof(double));
    for(i = 0; i < ndofs; ++i)
    {
        data[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
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
    FmsFieldSet(f, fd, 1, FMS_BY_VDIM, FMS_DOUBLE, data);

    /* Add some field metadata. (optional) */
    FmsMetaData mdata;
    FmsFieldAttachMetaData(f, &mdata);
    FmsMetaDataSetString(mdata, "units", "cm");

    free(data);
    free(coords);
}

/**
The main program.
*/
int
main(int argc, char *argv[])
{
    int i, ndofs, nverts, nedges, nfaces, ncells = 0;
    int *edges, *faces, *cells;
    double *verts, *coord_data, *x, *y, *z, px, py, pz;
    const char *protocol = "ascii";
    int dims[3] = {2,2,2};
    FmsInt order = 1;

    /* Handle some command line args. */
    if(argc > 1)
    {
        protocol = argv[1];
        if(argc > 2)
        {
            order = atoi(argv[2]);
            if(order < 1)
                return -1;
            if(argc >= 6)
            {
                dims[0] = atoi(argv[3]);
                dims[1] = atoi(argv[4]);
                dims[2] = atoi(argv[5]);
                if(dims[0] < 2 || dims[1] < 2 || dims[2] < 2)
                    return -2;
            }
        }
    }

    /* These vertices are just at the cell corners. We'll use them to make dofs later. */
    verts = compute_vertices(dims, &nverts);

    edges = compute_edges(dims, &nedges);
    faces = compute_faces(dims, &nfaces);
#ifdef ADD_CELLS
    cells = compute_cells(dims, &ncells);
#endif
    /* Create a mesh */
    FmsMesh mesh;
    FmsMeshConstruct(&mesh);

    /* Make a domain. */
    FmsDomain *domains;
    FmsMeshAddDomains(mesh, "domain", 1, &domains);

    /* Set the number of vertices */
    FmsDomainSetNumVertices(domains[0], nverts);

    /* Set the edges */
    FmsDomainSetNumEntities(domains[0], FMS_EDGE, FMS_INT32, nedges);
    FmsDomainAddEntities(domains[0], FMS_EDGE, NULL, FMS_INT32, edges, nedges);

    /* Set the faces */
    FmsDomainSetNumEntities(domains[0], FMS_QUADRILATERAL, FMS_INT32, nfaces);
    FmsDomainAddEntities(domains[0], FMS_QUADRILATERAL, NULL, FMS_INT32, faces, nfaces);
#ifdef ADD_CELLS
    /* Set the cells */
    FmsDomainSetNumEntities(domains[0], FMS_HEXAHEDRON, FMS_INT32, ncells);
    FmsDomainAddEntities(domains[0], FMS_HEXAHEDRON, NULL, FMS_INT32, cells, ncells);
#endif
    /* Make a component called surface. */
    FmsComponent volume;
    FmsMeshAddComponent(mesh, "volume", &volume);
    FmsComponentAddDomain(volume, domains[0]);

    FmsMeshFinalize(mesh);
    FmsMeshValidate(mesh);

    /* Make a data collection to contain the mesh. */
    FmsDataCollection dc;
    FmsDataCollectionCreate(mesh, "data collection", &dc);

    /* Create a field descriptor for the coordinates that live on "surface". */
    FmsFieldDescriptor coords_fd;
    FmsDataCollectionAddFieldDescriptor(dc, "coords descriptor", &coords_fd);
    FmsFieldDescriptorSetComponent(coords_fd, volume); /* define coords over surface */
    FmsFieldDescriptorSetFixedOrder(coords_fd, FMS_CONTINUOUS,
                                    FMS_NODAL_GAUSS_CLOSED, order);

    /* Make the coordinate data. Make it BYNODE xxxx...yyyy...zzzz... */
    ndofs = dofs_for_order(order, nverts, nedges, nfaces, ncells);
    coord_data = (double *)malloc(ndofs * sizeof(double)*3);
    x = coord_data;
    y = coord_data + ndofs;
    z = coord_data + 2*ndofs;

    /* These represent the amount we allow for perturbation of dofs to help 
       make things wavy. */
    px = 0.05 / ((double)(dims[0]-1));
    py = 0.05 / ((double)(dims[1]-1));
    pz = 0.05 / ((double)(dims[2]-1));

    /* Store coordinate dofs by node. */
    x = append_coord_dofs_verts(x, verts, nverts, 0);
    x = append_coord_dofs_edge(order, x, verts, edges, nedges, px, 0);
    x = append_coord_dofs_face(order, x, verts, edges, faces, nfaces, px, 0);
    x = append_coord_dofs_cell(order, x, dims, verts, px, 0);

    y = append_coord_dofs_verts(y, verts, nverts, 1);
    y = append_coord_dofs_edge(order, y, verts, edges, nedges, py, 1);
    y = append_coord_dofs_face(order, y, verts, edges, faces, nfaces, py, 1);
    y = append_coord_dofs_cell(order, y, dims, verts, py, 1);

    z = append_coord_dofs_verts(z, verts, nverts, 2);
    z = append_coord_dofs_edge(order, z, verts, edges, nedges, pz, 2);
    z = append_coord_dofs_face(order, z, verts, edges, faces, nfaces, pz, 2);
    z = append_coord_dofs_cell(order, z, dims, verts, pz, 2);

    /* Add the coordinates to the data collection. */
    FmsField coords;
    FmsDataCollectionAddField(dc, "coords", &coords);
    FmsFieldSet(coords, coords_fd, 3, FMS_BY_NODES, FMS_DOUBLE, coord_data);

    /* Add some field metadata. (optional) */
    FmsMetaData coords_mdata;
    FmsFieldAttachMetaData(coords, &coords_mdata);
    FmsMetaDataSetString(coords_mdata, "units", "cm");

    /* Set the coordinates for the component. */
    FmsComponentSetCoordinates(volume, coords);

    /* Add some scalars based on different orders of coordinates. */
    add_radial_scalar(dc, volume, "r1", 1, dims, verts, nverts, edges, nedges, faces, nfaces, ncells);
    add_radial_scalar(dc, volume, "r2", 2, dims, verts, nverts, edges, nedges, faces, nfaces, ncells);
    add_radial_scalar(dc, volume, "r3", 3, dims, verts, nverts, edges, nedges, faces, nfaces, ncells);

    /* How about a discontinuous scalar?...*/

    /* Add some metadata. (optional) */
    FmsMetaData mdata;
    int cycle = 1234;
    double dtime = 1.2345678;
    FmsDataCollectionAttachMetaData(dc, &mdata);
    FmsMetaData *mdata_objects;
    FmsMetaDataSetMetaData(mdata, "Info", 3, &mdata_objects);
    FmsMetaDataSetIntegers(mdata_objects[0], "cycle", FMS_INT32, 1, (void**)&cycle);
    FmsMetaDataSetScalars(mdata_objects[1], "time", FMS_DOUBLE, 1, (void**)&dtime);
    FmsMetaDataSetString(mdata_objects[2], "Description", "Hex example file");

    /* Write the data out. */
    FmsIOWrite("hex.fms", protocol, dc);

#if 0
    /* Write coordinate dofs as Point3D so we can plot them. */
    printf("X Y Z dof\n");
    for(i = 0; i < ndofs; ++i)
        printf("%lg %lg %lg %d\n", coord_data[i],coord_data[ndofs+i], coord_data[2*ndofs+i], i);
#endif

    free(coord_data);
    free(verts);
    free(edges);
    free(faces);
    free(cells);

    return 0;
}
