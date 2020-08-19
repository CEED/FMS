#include <fms.h>
#include <fmsio.h>
#include <gtest/gtest.h>

/***
 * Baseline data creation
 * ***/
int
Construct2DBaseline(FmsDataCollection *baseline) {
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

    *baseline = dc;
    return 0;
}

TEST(fmsio, WriteRead2DAscii) {
    FmsDataCollection baseline = NULL;
    ASSERT_EQ(Construct2DBaseline(&baseline), 0);
    ASSERT_TRUE(baseline);

    // No point in continuing if the write fails
    ASSERT_EQ(FmsIOWrite("TestWrite2DAscii.fms", "ascii", baseline), 0);

    FmsDataCollection dc = NULL;
    // No point in continuing if the read fails
    ASSERT_EQ(FmsIORead("TestWrite2DAscii.fms", "ascii", &dc), 0);
    ASSERT_TRUE(dc);
    EXPECT_EQ(FmsDataCollectionCompare(baseline, dc), 0);
}

TEST(fmsio, WriteRead2DYaml) {
    FmsDataCollection baseline = NULL;
    ASSERT_EQ(Construct2DBaseline(&baseline), 0);
    ASSERT_TRUE(baseline);

    // No point in continuing if the write fails
    ASSERT_EQ(FmsIOWrite("TestWrite2DAscii.fms", "yaml", baseline), 0);

    FmsDataCollection dc = NULL;
    // No point in continuing if the read fails
    ASSERT_EQ(FmsIORead("TestWrite2DAscii.fms", "yaml", &dc), 0);
    ASSERT_TRUE(dc);
    EXPECT_EQ(FmsDataCollectionCompare(baseline, dc), 0);
}