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
Determines the number of dofs on hex cells.
*/
int dofs_for_order(int order, int nverts, int nedges, int nfaces, int ncells)
{
    int O1 = order-1;
    return nverts + O1*nedges + O1*O1*nfaces + O1*O1*O1*ncells;
}

#if 0
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
#endif

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
#if 1
    for(i = 0; i < *nedges; i++)
        printf("edge %d: %d %d\n", i, e[2*i], e[2*i+1]);
#endif
    return e;
}

int *
compute_faces(const int *dims, int *nfaces)
{
    int i,j,k,A,B,C,D;
    /* number of edges. */
    int E1 = (dims[0]-1) * (dims[1])   * (dims[2]);
    int E2 = (dims[0])   * (dims[1]-1) * (dims[2]);
    int E3 = (dims[0])   * (dims[1])   * (dims[2]-1);

    int F1 = (dims[0]-1)*(dims[1]-1)*dims[2];     /* back/front */
    int F2 = (dims[0]-1)*(dims[1])  *(dims[2]-1); /* bottom/top */
    int F3 = (dims[0])  *(dims[1]-1)*(dims[2]-1); /* left/right */

    *nfaces = F1 + F2 + F3;
    printf("F1=%d\n", F1);
    printf("F2=%d\n", F2);
    printf("F3=%d\n", F3);
    printf("nfaces=%d\n", *nfaces);

    int *faces = (int *)malloc(*nfaces * sizeof(int)*4);
    int *f = faces;

/*
    For FMS_HEXAHEDRON, the faces (sides), "ABCDEF", the edges, "abcdefghijkl"
    and the vertices, "01234567", are ordered as follows:

              7--g--6
             /|    /|
            / l   / k   z=0      z=1      y=0      y=1      x=0      x=1
           h  |  f  |   bottom   top      front    back     left     right
          /   3-/c--2   2--c--3  7--g--6  4--e--5  6--g--7  7--h--4  5--f--6
         /   / /   /    |     |  |     |  |     |  |     |  |     |  |     |
        4--e--5   /     b  A  d  h  B  f  i  C  j  k  D  l  l  E  i  j  F  k
        |  d  |  b      |     |  |     |  |     |  |     |  |     |  |     |
        i /   j /       1--a--0  4--e--5  0--a--1  2--c--3  3--d--0  1--b--2
        |/    |/
        0--a--1

    For example, vertex "0" has coordinates (x,y,z)=(0,0,0), vertex "6" has
    coordinates (1,1,1), etc.
*/


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

        printf("A=%d, B=%d, C=%d, D=%d\n", A, B, C, D);
#if 1
        if(k == 0)
        {
            /* ADCB so normal faces out. */
            *f++ = A;
            *f++ = D;
            *f++ = C;
            *f++ = B;
        }
        else
#endif
        {
            /* ABCD */
            *f++ = A;
            *f++ = B;
            *f++ = C;
            *f++ = D;
        }
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

        printf("A=%d, B=%d, C=%d, D=%d\n", A, B, C, D);
#if 0
        if(j == 0)
        {
            /* ADCB so normal faces out. */
            *f++ = A;
            *f++ = D;
            *f++ = C;
            *f++ = B;
        }
        else
#endif
        {
            /* ABCD */
            *f++ = A;
            *f++ = B;
            *f++ = C;
            *f++ = D;
        }
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

        printf("A=%d, B=%d, C=%d, D=%d\n", A, B, C, D);
#if 1
        if(i == 0)
        {
            /* ADCB so normal faces out. */
            *f++ = A;
            *f++ = D;
            *f++ = C;
            *f++ = B;
        }
        else
#endif
        {
            /* ABCD */
            *f++ = A;
            *f++ = B;
            *f++ = C;
            *f++ = D;
        }
    }

#if 1
    printf("nfaces=%d\n", *nfaces);

    for(i = 0; i < *nfaces; i++)
        printf("face %d: %d %d %d %d\n", i, faces[4*i], faces[4*i+1], faces[4*i+2], faces[4*i+3]);
#endif

    return faces;
}

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

/*
    For FMS_HEXAHEDRON, the faces (sides), "ABCDEF", the edges, "abcdefghijkl"
    and the vertices, "01234567", are ordered as follows:

              7--g--6
             /|    /|
            / l   / k   z=0      z=1      y=0      y=1      x=0      x=1
           h  |  f  |   bottom   top      front    back     left     right
          /   3-/c--2   2--c--3  7--g--6  4--e--5  6--g--7  7--h--4  5--f--6
         /   / /   /    |     |  |     |  |     |  |     |  |     |  |     |
        4--e--5   /     b  A  d  h  B  f  i  C  j  k  D  l  l  E  i  j  F  k
        |  d  |  b      |     |  |     |  |     |  |     |  |     |  |     |
        i /   j /       1--a--0  4--e--5  0--a--1  2--c--3  3--d--0  1--b--2
        |/    |/
        0--a--1

    For example, vertex "0" has coordinates (x,y,z)=(0,0,0), vertex "6" has
    coordinates (1,1,1), etc.
*/

        *cptr++ = f0;
        *cptr++ = f1;
        *cptr++ = f3;
        *cptr++ = f2;
        *cptr++ = f4;
        *cptr++ = f5;

#if 0
/*This is the MFEM face ordering and it does not work. */
        *cptr++ = f0;
        *cptr++ = f2;
        *cptr++ = f5;
        *cptr++ = f3;
        *cptr++ = f4;
        *cptr++ = f1;
#endif
    }
#if 1
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

double *
append_coord_dofs_verts(double *dest, const double *src, int nverts, int component)
{
    int i;
    for(i = 0; i < nverts; ++i)
        dest[i] = src[3*i + component];
    return dest + nverts;
}

int
main(int argc, char *argv[])
{
    int i, nverts, nedges, nfaces, ncells;
    int *edges, *faces, *cells;
    double *verts;
    const char *protocol = "ascii";
    int dims[3] = {6,6,2};

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

    /* These vertices are just at the cell corners. We'll use them to make dofs later. */
    verts = compute_vertices(dims, &nverts);

    edges = compute_edges(dims, &nedges);
    faces = compute_faces(dims, &nfaces);
    cells = compute_cells(dims, &ncells);   

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
#if 1
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
    int ndofs = dofs_for_order(order, nverts, nedges, nfaces, ncells);

printf("nverts = %d\n", nverts);
printf("nedges = %d\n", nedges);
printf("nfaces = %d\n", nfaces);
printf("ncells = %d\n", ncells);

    double *coord_data = (double *)malloc(ndofs * sizeof(double)*3);
    double *x = coord_data;
    double *y = coord_data + ndofs;
    double *z = coord_data + 2*ndofs;

    /* Store order 1 coordinate dofs. */
    x = append_coord_dofs_verts(x, verts, nverts, 0);
    y = append_coord_dofs_verts(y, verts, nverts, 1);
    z = append_coord_dofs_verts(z, verts, nverts, 2);

#if 0
    double *cptr = coord_data;
    cptr = append_coords_vert(cptr, global_coords, nverts);
    cptr = append_coords_edge(order, cptr, global_coords, edge_vert, nedges, 0.02);
    cptr = append_coords_face(order, cptr, global_coords, edge_vert, quad_edge, nquads, 0.01);
    cptr = append_coords_cell(order, cptr, global_coords, edge_vert, quad_edge, nquads, 0.01);
#endif
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

#if 0
    /* Add some scalars based on different orders of coordinates. */
    add_radial_scalar(dc, surface, "r1", global_coords, 1, nverts, nedges, nquads);
    add_radial_scalar(dc, surface, "r2", global_coords, 2, nverts, nedges, nquads);
    add_radial_scalar(dc, surface, "r3", global_coords, 3, nverts, nedges, nquads);
    add_radial_scalar(dc, surface, "r4", global_coords, 4, nverts, nedges, nquads);

    /* Add an order 3 scalar that includes the dof order. */
    add_doforder_scalar(dc, surface, "doforder", nverts, nedges, nquads);
#endif

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

    /* Write coordinate dofs as Point3D so we can plot them. */
    printf("X Y Z dof\n");
    for(i = 0; i < ndofs; ++i)
        printf("%lg %lg %lg %d\n", coord_data[i],coord_data[ndofs+i], coord_data[2*ndofs+i], i);

    free(coord_data);

    return 0;
}
