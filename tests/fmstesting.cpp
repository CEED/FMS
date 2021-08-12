#include "fmstesting.hpp"
#include <cstdlib>

// #define PRINT_DATASETS
#ifdef PRINT_DATASETS
#include <fmsio.h>
#endif

int Construct2DData0(FmsDataCollection *out_dc) {
    if(!out_dc) return 1;
    // Taken from examples/meshing/simple_mesh.c
    FmsMesh mesh;
    FmsMeshConstruct(&mesh);
    FmsMeshSetPartitionId(mesh, 0, 1);

    FmsDomain *domain;
    FmsMeshAddDomains(mesh, "domains", 1, &domain);
    FmsDomainSetNumVertices(domain[0], 4);
    FmsDomainSetNumEntities(domain[0], FMS_EDGE, FMS_INT32, 4);
    FmsDomainSetNumEntities(domain[0], FMS_QUADRILATERAL, FMS_INT32, 1);

    const int edge_vert[] = {0,1, 1,2, 2,3, 3,0};
    FmsDomainAddEntities(domain[0], FMS_EDGE, NULL, FMS_INT32, edge_vert, 4);

    const int quad_edge[] = {0, 1, 2, 3};
    FmsDomainAddEntities(domain[0], FMS_QUADRILATERAL, NULL,
                        FMS_INT32, quad_edge, 1);

    FmsComponent volume;
    FmsMeshAddComponent(mesh, "volume", &volume);
    FmsComponentAddDomain(volume, domain[0]);

    FmsComponent boundary;
    FmsMeshAddComponent(mesh, "boundary", &boundary);
    FmsInt part_id;
    FmsComponentAddPart(boundary, domain[0], &part_id);
    const int bdr_edge[] = {2};
    FmsComponentAddPartEntities(boundary, part_id, FMS_EDGE, FMS_INT32,
                                FMS_INT32, FMS_INT32, NULL, bdr_edge, NULL, 1);

    FmsComponentAddRelation(volume, 1); // 1 = index of "boundary" component

    // Mesh Tags
    FmsTag material;
    FmsMeshAddTag(mesh, "material", &material);
    FmsTagSetComponent(material, volume);
    const int material_tags[] = {1};
    FmsTagSet(material, FMS_INT32, FMS_INT32, material_tags, 1);

    int tvs[1] = {1};
    const char * const descriptions[] = {"Test"};
    FmsTagAddDescriptions(material, FMS_INT32, tvs, descriptions, 1);

    FmsMeshFinalize(mesh);
    FmsMeshValidate(mesh);

    FmsDataCollection dc;
    FmsDataCollectionCreate(mesh, "data collection", &dc);

    FmsFieldDescriptor coords_fd;
    FmsDataCollectionAddFieldDescriptor(dc, "coords descriptor", &coords_fd);
    FmsFieldDescriptorSetComponent(coords_fd, volume);
    FmsInt coords_order = 1;
    FmsFieldDescriptorSetFixedOrder(coords_fd, FMS_CONTINUOUS,
                                    FMS_NODAL_GAUSS_CLOSED, coords_order);

    FmsField coords;
    FmsDataCollectionAddField(dc, "coords", &coords);
    const double coords_data[] = {
    0.,0.,
    1.,0.,
    1.,1.,
    0.,1.
    };
    FmsInt sdim = 2; // A 2D mesh embedded in "sdim"-dimensional space.
    FmsFieldSet(coords, coords_fd, sdim, FMS_BY_VDIM, FMS_DOUBLE, coords_data);

    FmsComponentSetCoordinates(volume, coords);
    *out_dc = dc;
#ifdef PRINT_DATASETS
    FmsIOWrite("2DData0.fms", "ascii", dc);
#endif
    return 0;
}

int Construct2DData1(FmsDataCollection *out_dc) {
    if(!out_dc) return 1;
    // Taken from examples/meshing/simple_mesh.c
    FmsMesh mesh;
    FmsMeshConstruct(&mesh);
    FmsMeshSetPartitionId(mesh, 0, 1);

    FmsDomain *domain = NULL;
    FmsMeshAddDomains(mesh, "domains", 1, &domain);
    FmsDomainSetNumVertices(domain[0], 3);
    FmsDomainSetNumEntities(domain[0], FMS_EDGE, FMS_INT32, 3);
    FmsDomainSetNumEntities(domain[0], FMS_TRIANGLE, FMS_INT32, 1);

    const int edge_vert[] = {0,1, 1,2, 2,0};
    FmsDomainAddEntities(domain[0], FMS_EDGE, NULL, FMS_INT32, edge_vert, 3);

    const int tri_edge[] = {0, 1, 2};
    FmsDomainAddEntities(domain[0], FMS_TRIANGLE, NULL,
                        FMS_INT32, tri_edge, 1);

    FmsComponent volume;
    FmsMeshAddComponent(mesh, "volume", &volume);
    FmsComponentAddDomain(volume, domain[0]);

    FmsComponent boundary;
    FmsMeshAddComponent(mesh, "boundary", &boundary);
    FmsInt part_id;
    FmsComponentAddPart(boundary, domain[0], &part_id);
    const int bdr_edge[] = {2};
    FmsComponentAddPartEntities(boundary, part_id, FMS_EDGE, FMS_INT32,
                                FMS_INT32, FMS_INT32, NULL, bdr_edge, NULL, 1);

    FmsComponentAddRelation(volume, 1); // 1 = index of "boundary" component

    // Mesh Tags
    FmsTag material;
    FmsMeshAddTag(mesh, "material", &material);
    FmsTagSetComponent(material, volume);
    const int material_tags[] = {1};
    FmsTagSet(material, FMS_INT32, FMS_INT32, material_tags, 1);

    int tvs[1] = {1};
    const char * const descriptions[] = {"Test"};
    FmsTagAddDescriptions(material, FMS_INT32, tvs, descriptions, 1);

    FmsMeshFinalize(mesh);
    FmsMeshValidate(mesh);

    FmsDataCollection dc;
    FmsDataCollectionCreate(mesh, "data collection", &dc);

    FmsFieldDescriptor coords_fd;
    FmsDataCollectionAddFieldDescriptor(dc, "coords descriptor", &coords_fd);
    FmsFieldDescriptorSetComponent(coords_fd, volume);
    FmsInt coords_order = 1;
    FmsFieldDescriptorSetFixedOrder(coords_fd, FMS_CONTINUOUS,
                                    FMS_NODAL_GAUSS_CLOSED, coords_order);

    FmsField coords;
    FmsDataCollectionAddField(dc, "coords", &coords);
    const double coords_data[] = {
        0.,0.,
        1.,0.,
        0.,1.
    };
    FmsInt sdim = 2; // A 2D mesh embedded in "sdim"-dimensional space.
    FmsFieldSet(coords, coords_fd, sdim, FMS_BY_VDIM, FMS_DOUBLE, coords_data);

    FmsComponentSetCoordinates(volume, coords);
    *out_dc = dc;
#ifdef PRINT_DATASETS
    FmsIOWrite("2DData1.fms", "ascii", dc);
#endif
    return 0;
}

// 1 element hex mesh
int Construct3DData0(FmsDataCollection *out_dc) {
    if(!out_dc) return 1;
    FmsMesh mesh;
    FmsMeshConstruct(&mesh);
    FmsMeshSetPartitionId(mesh, 0, 1);
    FmsDomain *domain;
    FmsMeshAddDomains(mesh, "domains", 1, &domain);
    FmsDomainSetNumVertices(domain[0], 8);
    const int edge_vert[] = {0,1, 2,3, 4,5, 6,7, 
                            0,2, 1,3, 4,6, 5,7,
                            0,4, 1,5, 2,6, 3,7};
    FmsDomainSetNumEntities(domain[0], FMS_EDGE, FMS_INT32, 12);
    FmsDomainAddEntities(domain[0], FMS_EDGE, NULL, FMS_INT32, edge_vert, 12);
    const int quad_edge[] = {0,5,1,4, 2,7,3,6, 0,9,2,8, 1,10,3,11, 4,8,6,10, 5,11,7,9};
    FmsDomainSetNumEntities(domain[0], FMS_QUADRILATERAL, FMS_INT32, 6);
    FmsDomainAddEntities(domain[0], FMS_QUADRILATERAL, NULL,
                        FMS_INT32, quad_edge, 6);
    const int hex_face[] = {0,1,2,3,4,5};
    FmsDomainSetNumEntities(domain[0], FMS_HEXAHEDRON, FMS_INT32, 1);
    FmsDomainAddEntities(domain[0], FMS_HEXAHEDRON, NULL, FMS_INT32, hex_face, 1);
    FmsComponent volume;
    FmsMeshAddComponent(mesh, "volume", &volume);
    FmsComponentAddDomain(volume, domain[0]);
    FmsComponent boundary;
    FmsMeshAddComponent(mesh, "boundary", &boundary);
    FmsInt part_id;
    FmsComponentAddPart(boundary, domain[0], &part_id);
    const int bdr_face[] = {0,2,4};
    FmsComponentAddPartEntities(boundary, part_id, FMS_QUADRILATERAL, FMS_INT32,
                                FMS_INT32, FMS_INT32, NULL, bdr_face, NULL, 3);
    FmsComponentAddRelation(volume, 1); // 1 = index of "boundary" component
    // Mesh Tags
    FmsTag material;
    FmsMeshAddTag(mesh, "material", &material);
    FmsTagSetComponent(material, volume);
    const int material_tags[] = {1};
    FmsTagSet(material, FMS_INT32, FMS_INT32, material_tags, 1);
    int tvs[1] = {1};
    const char * const descriptions[] = {"Test"};
    FmsTagAddDescriptions(material, FMS_INT32, tvs, descriptions, 1);
    FmsMeshFinalize(mesh);
    FmsMeshValidate(mesh);
    FmsDataCollection dc;
    FmsDataCollectionCreate(mesh, "data collection", &dc);
    FmsFieldDescriptor coords_fd;
    FmsDataCollectionAddFieldDescriptor(dc, "coords descriptor", &coords_fd);
    FmsFieldDescriptorSetComponent(coords_fd, volume);
    FmsInt coords_order = 1;
    FmsFieldDescriptorSetFixedOrder(coords_fd, FMS_CONTINUOUS,
                                    FMS_NODAL_GAUSS_CLOSED, coords_order);
    FmsField coords;
    FmsDataCollectionAddField(dc, "coords", &coords);
    const double coords_data[] = {
        -1,-1,-1, 1,-1,-1, -1,1,-1, 1,1,-1,
        -1,-1,1,  1,-1,1,  -1,1,1,  1,1,1};
    FmsInt sdim = 3;
    FmsFieldSet(coords, coords_fd, sdim, FMS_BY_VDIM, FMS_DOUBLE, coords_data);
    FmsComponentSetCoordinates(volume, coords);
    *out_dc = dc;
#ifdef PRINT_DATASETS
    FmsIOWrite("3DData0.fms", "ascii", dc);
#endif
    return 0;
}

// 1 element tet mesh
int Construct3DData1(FmsDataCollection *out_dc) {
    if(!out_dc) return 1;
    // Taken from examples/meshing/simple_mesh.c
    FmsMesh mesh;
    FmsMeshConstruct(&mesh);
    FmsMeshSetPartitionId(mesh, 0, 1);
    FmsDomain *domain;
    FmsMeshAddDomains(mesh, "domains", 1, &domain);
    FmsDomainSetNumVertices(domain[0], 4);
    const int edge_vert[] = {0,1, 0,2, 1,2, 0,3, 1,3, 2,3};
    FmsDomainSetNumEntities(domain[0], FMS_EDGE, FMS_INT32, 6);
    FmsDomainAddEntities(domain[0], FMS_EDGE, NULL, FMS_INT32, edge_vert, 6);
    const int tri_edge[] = {0,2,1, 1,3,5, 0,4,3, 2,5,4};
    FmsDomainSetNumEntities(domain[0], FMS_TRIANGLE, FMS_INT32, 4);
    FmsDomainAddEntities(domain[0], FMS_TRIANGLE, NULL,
                        FMS_INT32, tri_edge, 4);
    const int tet_faces[] = {0,1,2,3};
    FmsDomainSetNumEntities(domain[0], FMS_TETRAHEDRON, FMS_INT32, 1);
    FmsDomainAddEntities(domain[0], FMS_TETRAHEDRON, NULL, FMS_INT32, tet_faces, 1);
    FmsComponent volume;
    FmsMeshAddComponent(mesh, "volume", &volume);
    FmsComponentAddDomain(volume, domain[0]);
    FmsComponent boundary;
    FmsMeshAddComponent(mesh, "boundary", &boundary);
    FmsInt part_id;
    FmsComponentAddPart(boundary, domain[0], &part_id);
    const int bdr_face[] = {0,1,2};
    FmsComponentAddPartEntities(boundary, part_id, FMS_TRIANGLE, FMS_INT32,
                                FMS_INT32, FMS_INT32, NULL, bdr_face, NULL, 3);
    FmsComponentAddRelation(volume, 1); // 1 = index of "boundary" component
    // Mesh Tags
    FmsTag material;
    FmsMeshAddTag(mesh, "material", &material);
    FmsTagSetComponent(material, volume);
    const int material_tags[] = {1};
    FmsTagSet(material, FMS_INT32, FMS_INT32, material_tags, 1);
    int tvs[1] = {1};
    const char * const descriptions[] = {"Test"};
    FmsTagAddDescriptions(material, FMS_INT32, tvs, descriptions, 1);
    FmsMeshFinalize(mesh);
    FmsMeshValidate(mesh);
    FmsDataCollection dc;
    FmsDataCollectionCreate(mesh, "data collection", &dc);
    FmsFieldDescriptor coords_fd;
    FmsDataCollectionAddFieldDescriptor(dc, "coords descriptor", &coords_fd);
    FmsFieldDescriptorSetComponent(coords_fd, volume);
    FmsInt coords_order = 1;
    FmsFieldDescriptorSetFixedOrder(coords_fd, FMS_DISCONTINUOUS,
                                    FMS_NODAL_GAUSS_CLOSED, coords_order);
    FmsField coords;
    FmsDataCollectionAddField(dc, "coords", &coords);
    const double coords_data[] = {
        -1,-1,-1, 1,-1,-1, -1,1,-1, -1,-1,1};
    FmsInt sdim = 3;
    FmsFieldSet(coords, coords_fd, sdim, FMS_BY_VDIM, FMS_DOUBLE, coords_data);
    FmsComponentSetCoordinates(volume, coords);
    FmsMetaData mdata = NULL;
    FmsDataCollectionAttachMetaData(dc, &mdata);
    void *cycle = nullptr;
    FmsMetaDataSetIntegers(mdata, "cycle", FMS_UINT8, 1, &cycle);
    *static_cast<uint8_t*>(cycle) = 1u;
    *out_dc = dc;
#ifdef PRINT_DATASETS
    FmsIOWrite("3DData1.fms", "ascii", dc);
#endif
    return 0;
}

int ConstructAllDataCollections(FmsInt *size, FmsDataCollection **out_dcs) {
    if(!size) return 1;
    if(!out_dcs) return 2;
    FmsDataCollection *dcs = static_cast<FmsDataCollection*>(std::calloc(sizeof(FmsDataCollection), 4));
    if(Construct2DData0(&dcs[0])) {
        if(dcs[0]) FmsDataCollectionDestroy(&dcs[0]);
        std::free(dcs);
        return 3;
    }
    if(Construct2DData1(&dcs[1])) {
        FmsDataCollectionDestroy(&dcs[0]);
        if(dcs[1]) FmsDataCollectionDestroy(&dcs[1]);
        std::free(dcs);
        return 4;
    }
    if(Construct3DData0(&dcs[2])) {
        FmsDataCollectionDestroy(&dcs[0]);
        FmsDataCollectionDestroy(&dcs[1]);
        if(dcs[2]) FmsDataCollectionDestroy(&dcs[2]);
        std::free(dcs);
        return 5;
    }
    if(Construct3DData1(&dcs[3])) {
        FmsDataCollectionDestroy(&dcs[0]);
        FmsDataCollectionDestroy(&dcs[1]);
        FmsDataCollectionDestroy(&dcs[2]);
        if(dcs[3]) FmsDataCollectionDestroy(&dcs[3]);
        std::free(dcs);
        return 6;
    }
    *out_dcs = dcs;
    *size = 4;
    return 0;
}
