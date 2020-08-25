#include <fms.h>
#include <fmstesting.hpp>
#include <gtest/gtest.h>

TEST(Fms, Version) {
    FmsInt v = 0;
    ASSERT_EQ(FmsGetInterfaceVersion(&v), 0);
    EXPECT_EQ(FMS_INTERFACE_VERSION, v);
}

// Compare interface
TEST(Fms, DataCollectionCompare) {
    FmsDataCollection *dcs0 = NULL;
    FmsDataCollection *dcs1 = NULL;
    FmsInt size0 = 0, size1 = -1;
    ASSERT_EQ(ConstructAllDataCollections(&size0, &dcs0), 0);
    ASSERT_TRUE(dcs0);
    ASSERT_TRUE(size0);
    ASSERT_EQ(ConstructAllDataCollections(&size1, &dcs1), 0);
    ASSERT_TRUE(dcs1);
    ASSERT_TRUE(size1);
    ASSERT_EQ(size0, size1);
    for(FmsInt i = 0; i < size0; i++) {
        ASSERT_TRUE(dcs0[i]);
        EXPECT_EQ(FmsDataCollectionCompare(dcs0[i], dcs0[i]), 0);
        for(FmsInt j = 0; j < size1; j++) {
            ASSERT_TRUE(dcs1[j]);
            ASSERT_NE(dcs0[i], dcs1[j]);
            EXPECT_EQ(FmsDataCollectionCompare(dcs1[j], dcs1[j]), 0);
            if(i == j) {
                EXPECT_EQ(FmsDataCollectionCompare(dcs0[i], dcs1[j]), 0) 
                    << "DataCollection " << i << " compared not equal to " << j << ".";
            }
            else {
                EXPECT_NE(FmsDataCollectionCompare(dcs0[i], dcs1[j]), 0) 
                    << "DataCollection " << i << " compared equal to " << j << ".";
            }
        }
    }
    for(FmsInt i = 0; i < size0; i++) {
        EXPECT_EQ(FmsDataCollectionDestroy(&dcs0[i]), 0);
        EXPECT_EQ(FmsDataCollectionDestroy(&dcs1[i]), 0);
    }
    std::free(dcs0);
    std::free(dcs1);
}

TEST(Fms, MetaDataCompare) {
    // First we need to create a mesh and a datacollection to attach the metadata to
    FmsMesh mesh = NULL;
    ASSERT_EQ(FmsMeshConstruct(&mesh), 0);
    ASSERT_TRUE(mesh);
    ASSERT_EQ(FmsMeshFinalize(mesh), 0);
    FmsDataCollection dc = NULL;
    ASSERT_EQ(FmsDataCollectionCreate(mesh, "MetaDataCompareDC", &dc), 0);
    ASSERT_TRUE(dc);

    // Now start creating and comparing MetaData
    FmsMetaData toplevel_md = NULL;
    ASSERT_EQ(FmsDataCollectionAttachMetaData(dc, &toplevel_md), 0);
    ASSERT_TRUE(toplevel_md);
    FmsMetaData *mds = NULL;
    ASSERT_EQ(FmsMetaDataSetMetaData(toplevel_md, "MetaDataCompareDC", 8, &mds), 0);
    ASSERT_TRUE(mds);

    // Compare integer mds
    FmsInt *data0 = NULL;
    ASSERT_EQ(FmsMetaDataSetIntegers(mds[0], "Ints", FMS_INT_TYPE, 3, (void**)&data0), 0);
    ASSERT_TRUE(data0);
    data0[0] = 3; data0[1] = 5; data0[2] = 1;
    EXPECT_EQ(FmsMetaDataCompare(mds[0], mds[0]), 0);
    // First make sure that they compare equal when equal
    FmsInt *data1 = NULL;
    ASSERT_EQ(FmsMetaDataSetIntegers(mds[1], "Ints", FMS_INT_TYPE, 3, (void**)&data1), 0);
    data1[0] = 3; data1[1] = 5; data1[2] = 1;
    ASSERT_NE(mds[0], mds[1]);
    EXPECT_EQ(FmsMetaDataCompare(mds[0], mds[1]), 0);
    EXPECT_EQ(FmsMetaDataCompare(mds[1], mds[0]), 0);
    // Compare different names
    ASSERT_EQ(FmsMetaDataSetIntegers(mds[1], "Ints1", FMS_INT_TYPE, 3, (void**)&data1), 0);
    ASSERT_TRUE(data1);
    data1[0] = 3; data1[1] = 5; data1[2] = 1;
    EXPECT_NE(FmsMetaDataCompare(mds[0], mds[1]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[1], mds[0]), 0);
    // Compare different sizes
    ASSERT_EQ(FmsMetaDataSetIntegers(mds[1], "Ints", FMS_INT_TYPE, 4, (void**)&data1), 0);
    ASSERT_TRUE(data1);
    data1[0] = 3; data1[1] = 5; data1[2] = 1; data1[3] = 0;
    EXPECT_NE(FmsMetaDataCompare(mds[0], mds[1]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[1], mds[0]), 0);
    // Compare different contents
    ASSERT_EQ(FmsMetaDataSetIntegers(mds[1], "Ints", FMS_INT_TYPE, 3, (void**)&data1), 0);
    data1[0] = 2; data1[1] = 5; data1[2] = 1;
    EXPECT_NE(FmsMetaDataCompare(mds[0], mds[1]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[1], mds[0]), 0);
    data1[0] = 3; data1[1] = 0;
    EXPECT_NE(FmsMetaDataCompare(mds[0], mds[1]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[1], mds[0]), 0);
    data1[1] = 5; data1[2] = 0;
    EXPECT_NE(FmsMetaDataCompare(mds[0], mds[1]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[1], mds[0]), 0);
    data1[2] = 1;
    EXPECT_EQ(FmsMetaDataCompare(mds[0], mds[1]), 0);
    EXPECT_EQ(FmsMetaDataCompare(mds[1], mds[0]), 0);
    data1 = NULL;
    // Compare different int types
    uint16_t *data11 = NULL;
    ASSERT_EQ(FmsMetaDataSetIntegers(mds[1], "Ints", FMS_UINT16, 3, (void**)&data11), 0);
    ASSERT_TRUE(data11);
    data11[0] = 3; data11[1] = 5; data11[2] = 1;
    EXPECT_NE(FmsMetaDataCompare(mds[0], mds[1]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[1], mds[0]), 0);

    // Compare scalar mds
    double *data2 = NULL;
    ASSERT_EQ(FmsMetaDataSetScalars(mds[2], "Ints", FMS_DOUBLE, 4, (void**)&data2), 0);
    ASSERT_TRUE(data2);
    data2[0] = 3.; data2[1] = 5.; data2[2] = 1.; data2[3] = 10.;
    EXPECT_EQ(FmsMetaDataCompare(mds[2], mds[2]), 0);
    // Compare different metadata type
    EXPECT_NE(FmsMetaDataCompare(mds[0], mds[2]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[2], mds[0]), 0);
    ASSERT_EQ(FmsMetaDataSetScalars(mds[2], "Scalars", FMS_DOUBLE, 4, (void**)&data2), 0);
    ASSERT_TRUE(data2);
    data2[0] = 3.; data2[1] = 5.; data2[2] = 1.; data2[3] = 10.;
    // Compare different contents
    double *data3 = NULL;
    ASSERT_EQ(FmsMetaDataSetScalars(mds[3], "Scalars", FMS_DOUBLE, 4, (void**)&data3), 0);
    ASSERT_TRUE(data3);
    data3[0] = 13.; data3[1] = 5.; data3[2] = 1.; data3[3] = 10.;
    EXPECT_NE(FmsMetaDataCompare(mds[3], mds[2]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[2], mds[3]), 0);
    data3[0] = 3.; data3[1] = 13.;
    EXPECT_NE(FmsMetaDataCompare(mds[3], mds[2]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[2], mds[3]), 0);
    data3[1] = 5.; data3[2] = 13.;
    EXPECT_NE(FmsMetaDataCompare(mds[3], mds[2]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[2], mds[3]), 0);
    data3[2] = 1.; data3[3] = 13.;
    EXPECT_NE(FmsMetaDataCompare(mds[3], mds[2]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[2], mds[3]), 0);
    data3[3] = 10.;
    EXPECT_EQ(FmsMetaDataCompare(mds[3], mds[2]), 0);
    EXPECT_EQ(FmsMetaDataCompare(mds[2], mds[3]), 0);
    // Compare different sizes
    ASSERT_EQ(FmsMetaDataSetScalars(mds[3], "Scalars", FMS_DOUBLE, 3, (void**)&data3), 0);
    ASSERT_TRUE(data3);
    data3[0] = 3.; data3[1] = 5.; data3[2] = 1.;
    EXPECT_NE(FmsMetaDataCompare(mds[3], mds[2]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[2], mds[3]), 0);
    data3 = NULL;
    // Compare different scalar type
    float *data31 = NULL;
    ASSERT_EQ(FmsMetaDataSetScalars(mds[3], "Scalars", FMS_FLOAT, 4, (void**)&data31), 0);
    ASSERT_TRUE(data31);
    data31[0] = 3.f; data31[1] = 5.f; data31[2] = 1.f; data31[3] = 10.f;
    EXPECT_NE(FmsMetaDataCompare(mds[3], mds[2]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[2], mds[3]), 0);

    // Compare string mds
    const char *data4 = "Hello";
    ASSERT_EQ(FmsMetaDataSetString(mds[4], "String", data4), 0);
    EXPECT_EQ(FmsMetaDataCompare(mds[4], mds[4]), 0);
    // Compare equal strings
    const char *data5 = "Hello";
    ASSERT_EQ(FmsMetaDataSetString(mds[5], "String", data5), 0);
    EXPECT_EQ(FmsMetaDataCompare(mds[4], mds[5]), 0);
    EXPECT_EQ(FmsMetaDataCompare(mds[5], mds[4]), 0);
    // Compare different contents
    const char *data51 = "World";
    ASSERT_EQ(FmsMetaDataSetString(mds[5], "String", data51), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[4], mds[5]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[5], mds[4]), 0);

    FmsMetaData *data6 = NULL;
    ASSERT_EQ(FmsMetaDataSetMetaData(mds[6], "MetaDatas", 2, &data6), 0);
    ASSERT_TRUE(data6);
    int32_t *data61 = NULL;
    ASSERT_EQ(FmsMetaDataSetIntegers(data6[0], "Int", FMS_INT32, 1, (void**)&data61), 0);
    ASSERT_TRUE(data61);
    data61[0] = 7;
    double *data62 = NULL;
    ASSERT_EQ(FmsMetaDataSetScalars(data6[1], "Scalars", FMS_DOUBLE, 3, (void**)&data62), 0);
    ASSERT_TRUE(data62);
    data62[0] = 3.; data62[1] = 5.; data62[2] = 1.; 

    // Compare everything so far to unset md
    for(int i = 0; i < 7; i++) {
        EXPECT_NE(FmsMetaDataCompare(mds[i], mds[7]), 0);
        EXPECT_NE(FmsMetaDataCompare(mds[7], mds[i]), 0);
    }

    // Compare md of mds
    // Compare equal
    FmsMetaData *data7 = NULL;
    ASSERT_EQ(FmsMetaDataSetMetaData(mds[7], "MetaDatas", 2, &data7), 0);
    ASSERT_TRUE(data7);
    int32_t *data71 = NULL;
    ASSERT_EQ(FmsMetaDataSetIntegers(data7[0], "Int", FMS_INT32, 1, (void**)&data71), 0);
    ASSERT_TRUE(data71);
    data71[0] = 7;
    double *data72 = NULL;
    ASSERT_EQ(FmsMetaDataSetScalars(data7[1], "Scalars", FMS_DOUBLE, 3, (void**)&data72), 0);
    ASSERT_TRUE(data72);
    data72[0] = 3.; data72[1] = 5.; data72[2] = 1.; 
    EXPECT_EQ(FmsMetaDataCompare(mds[6], mds[7]), 0);
    EXPECT_EQ(FmsMetaDataCompare(mds[7], mds[6]), 0);
    data7 = NULL;
    data71 = NULL;
    data72 = NULL;
    // Compare different sizes
    ASSERT_EQ(FmsMetaDataSetMetaData(mds[7], "MetaDatas", 1, &data7), 0);
    ASSERT_TRUE(data7);
    ASSERT_EQ(FmsMetaDataSetIntegers(data7[0], "Int", FMS_INT32, 1, (void**)&data71), 0);
    ASSERT_TRUE(data71);
    data71[0] = 7;
    EXPECT_NE(FmsMetaDataCompare(mds[6], mds[7]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[7], mds[6]), 0);
    data7 = NULL;
    data71 = NULL;
    // Compare different content
    ASSERT_EQ(FmsMetaDataSetMetaData(mds[7], "MetaDatas", 2, &data7), 0);
    ASSERT_TRUE(data7);
    ASSERT_EQ(FmsMetaDataSetScalars(data7[0], "Scalars", FMS_DOUBLE, 3, (void**)&data72), 0);
    ASSERT_TRUE(data72);
    data72[0] = 3.; data72[1] = 5.; data72[2] = 1.;
    ASSERT_EQ(FmsMetaDataSetIntegers(data7[1], "Int", FMS_INT32, 1, (void**)&data71), 0);
    ASSERT_TRUE(data71);
    data71[0] = 7;
    EXPECT_NE(FmsMetaDataCompare(mds[6], mds[7]), 0);
    EXPECT_NE(FmsMetaDataCompare(mds[7], mds[6]), 0);
    // Make them equal again
    FmsMetaData temp = data7[0];
    data7[0] = data7[1]; data7[1] = temp;
    EXPECT_EQ(FmsMetaDataCompare(mds[6], mds[7]), 0);
    EXPECT_EQ(FmsMetaDataCompare(mds[7], mds[6]), 0);

    FmsDataCollectionDestroy(&dc);
}

TEST(Fms, FieldDescriptorCompare) {
    // We need a DataCollection to attach the FDs to
    FmsDataCollection dc = NULL;
    ASSERT_EQ(Construct2DData0(&dc), 0);
    ASSERT_TRUE(dc);
    FmsMesh mesh = NULL;
    ASSERT_EQ(FmsDataCollectionGetMesh(dc, &mesh), 0);
    ASSERT_TRUE(mesh);
    FmsComponent volume = NULL, boundary = NULL;
    ASSERT_EQ(FmsMeshGetComponent(mesh, 0, &volume), 0);
    ASSERT_EQ(FmsMeshGetComponent(mesh, 1, &boundary), 0);
    ASSERT_TRUE(volume); ASSERT_TRUE(boundary);
    ASSERT_NE(volume, boundary);

    // Compare equal FDs
    FmsDataCollection dc1 = NULL;
    ASSERT_EQ(Construct2DData0(&dc1), 0);
    ASSERT_TRUE(dc1);
    FmsInt dc_nfds = 0, dc1_nfds = -1;
    FmsFieldDescriptor *dc_fds = NULL, *dc1_fds = NULL;
    ASSERT_EQ(FmsDataCollectionGetFieldDescriptors(dc, &dc_fds, &dc_nfds), 0);
    ASSERT_EQ(FmsDataCollectionGetFieldDescriptors(dc1, &dc1_fds, &dc1_nfds), 0);
    ASSERT_TRUE(dc_fds); ASSERT_TRUE(dc1_fds);
    ASSERT_TRUE(dc_nfds); ASSERT_TRUE(dc1_nfds);
    ASSERT_EQ(dc_nfds, dc1_nfds);
    for(FmsInt i = 0; i < dc_nfds; i++) {
        ASSERT_TRUE(dc_fds[i]); ASSERT_TRUE(dc1_fds[i]);
        ASSERT_NE(dc_fds[i], dc1_fds[i]);
        EXPECT_EQ(FmsFieldDescriptorCompare(dc_fds[i], dc1_fds[i]), 0)
            << "The FD at index " << i << " did not compare equal.";
    }
    FmsDataCollectionDestroy(&dc1);

    // Compare different combinations of basis type
    // NOTE: Might need to use the two identical datacollection from above because it might 
    //  not always be okay to have two FDs with the same name in the same DataCollection
    for(int i = 0; i < 7; i++) {
        FmsFieldDescriptor fdi = NULL;
        ASSERT_EQ(FmsDataCollectionAddFieldDescriptor(dc, "FD", &fdi), 0);
        ASSERT_TRUE(fdi);
        ASSERT_EQ(FmsFieldDescriptorSetComponent(fdi, volume), 0);
        ASSERT_EQ(FmsFieldDescriptorSetFixedOrder(fdi, FMS_CONTINUOUS, static_cast<FmsBasisType>(i), 10), 0);
        for(int j = 0; j < 7; j++){
            FmsFieldDescriptor fdj = NULL;
            ASSERT_EQ(FmsDataCollectionAddFieldDescriptor(dc, "FD", &fdj), 0);
            ASSERT_EQ(FmsFieldDescriptorSetComponent(fdj, volume), 0);
            ASSERT_EQ(FmsFieldDescriptorSetFixedOrder(fdj, FMS_CONTINUOUS, static_cast<FmsBasisType>(j), 10), 0);
            if(i == j) {
                EXPECT_EQ(FmsFieldDescriptorCompare(fdi, fdj), 0)
                    << "FDs compared not equal with basistype " << i;
                EXPECT_EQ(FmsFieldDescriptorCompare(fdj, fdi), 0)
                    << "FDs compared not equal with basistype " << i;;
            }
            else {
                EXPECT_NE(FmsFieldDescriptorCompare(fdi, fdj), 0)
                    << "FDs compared equal with basistypes " << i << " " << j;
                EXPECT_NE(FmsFieldDescriptorCompare(fdj, fdi), 0)
                    << "FDs compared equal with basistypes " << i << " " << j;
            }
        }
    }

    // Compare different combinations of FieldType
    // TODO: Update this when HCURL is implemented
    for(int i = 0; i < 5; i++) {
        if(static_cast<FmsFieldType>(i) == FMS_HCURL) continue;
        FmsFieldDescriptor fdi = NULL;
        ASSERT_EQ(FmsDataCollectionAddFieldDescriptor(dc, "FD", &fdi), 0);
        ASSERT_TRUE(fdi);
        ASSERT_EQ(FmsFieldDescriptorSetComponent(fdi, volume), 0);
        ASSERT_EQ(FmsFieldDescriptorSetFixedOrder(fdi, static_cast<FmsFieldType>(i), FMS_NODAL_GAUSS_CLOSED, 10), 0);
        for(int j = 0; j < 5; j++){
            if(static_cast<FmsFieldType>(j) == FMS_HCURL) continue;
            FmsFieldDescriptor fdj = NULL;
            ASSERT_EQ(FmsDataCollectionAddFieldDescriptor(dc, "FD", &fdj), 0);
            ASSERT_EQ(FmsFieldDescriptorSetComponent(fdj, volume), 0);
            ASSERT_EQ(FmsFieldDescriptorSetFixedOrder(fdj, static_cast<FmsFieldType>(j), FMS_NODAL_GAUSS_CLOSED, 10), 0);
            ASSERT_NE(fdi, fdj);
            if(i == j) {
                EXPECT_EQ(FmsFieldDescriptorCompare(fdi, fdj), 0)
                    << "FDs compared not equal with fieldtype " << i;
                EXPECT_EQ(FmsFieldDescriptorCompare(fdj, fdi), 0)
                    << "FDs compared not equal with fieldtype " << i;;
            }
            else {
                EXPECT_NE(FmsFieldDescriptorCompare(fdi, fdj), 0)
                    << "FDs compared equal with fieldtypes " << i << " " << j;
                EXPECT_NE(FmsFieldDescriptorCompare(fdj, fdi), 0)
                    << "FDs compared equal with fieldtypes " << i << " " << j;
            }
        }
    }

    // Compare different orders
    for(int i = 1; i < 5; i++) {
        FmsFieldDescriptor fdi = NULL;
        ASSERT_EQ(FmsDataCollectionAddFieldDescriptor(dc, "FD", &fdi), 0);
        ASSERT_TRUE(fdi);
        ASSERT_EQ(FmsFieldDescriptorSetComponent(fdi, volume), 0);
        ASSERT_EQ(FmsFieldDescriptorSetFixedOrder(fdi, FMS_CONTINUOUS, FMS_NODAL_GAUSS_CLOSED, i), 0);
        for(int j = 1; j < 5; j++){
            if(static_cast<FmsFieldType>(j) == FMS_HCURL) continue;
            FmsFieldDescriptor fdj = NULL;
            ASSERT_EQ(FmsDataCollectionAddFieldDescriptor(dc, "FD", &fdj), 0);
            ASSERT_EQ(FmsFieldDescriptorSetComponent(fdj, volume), 0);
            ASSERT_EQ(FmsFieldDescriptorSetFixedOrder(fdj, FMS_CONTINUOUS, FMS_NODAL_GAUSS_CLOSED, j), 0);
            ASSERT_NE(fdi, fdj);
            if(i == j) {
                EXPECT_EQ(FmsFieldDescriptorCompare(fdi, fdj), 0)
                    << "FDs compared not equal with order " << i;
                EXPECT_EQ(FmsFieldDescriptorCompare(fdj, fdi), 0)
                    << "FDs compared not equal with order " << i;;
            }
            else {
                EXPECT_NE(FmsFieldDescriptorCompare(fdi, fdj), 0)
                    << "FDs compared equal with orders " << i << " " << j;
                EXPECT_NE(FmsFieldDescriptorCompare(fdj, fdi), 0)
                    << "FDs compared equal with orders " << i << " " << j;
            }
        }
    }

    // Compare different components
    FmsFieldDescriptor fdv = NULL, fdb = NULL;
    ASSERT_EQ(FmsDataCollectionAddFieldDescriptor(dc, "FD", &fdv), 0);
    ASSERT_EQ(FmsDataCollectionAddFieldDescriptor(dc, "FD", &fdb), 0);
    ASSERT_TRUE(fdv); ASSERT_TRUE(fdb);
    ASSERT_EQ(FmsFieldDescriptorSetComponent(fdv, volume), 0);
    ASSERT_EQ(FmsFieldDescriptorSetComponent(fdb, boundary), 0);
    ASSERT_EQ(FmsFieldDescriptorSetFixedOrder(fdv, FMS_CONTINUOUS, FMS_NODAL_GAUSS_CLOSED, 3), 0);
    ASSERT_EQ(FmsFieldDescriptorSetFixedOrder(fdb, FMS_CONTINUOUS, FMS_NODAL_GAUSS_CLOSED, 3), 0);
    ASSERT_NE(fdv, fdb);
    EXPECT_NE(FmsFieldDescriptorCompare(fdv, fdb), 0);
    EXPECT_NE(FmsFieldDescriptorCompare(fdb, fdv), 0);

    // Compare different names
    FmsFieldDescriptor diff_name = NULL;
    ASSERT_EQ(FmsDataCollectionAddFieldDescriptor(dc, "DiffName", &diff_name), 0);
    ASSERT_TRUE(diff_name);
    ASSERT_EQ(FmsFieldDescriptorSetComponent(diff_name, volume), 0);
    ASSERT_EQ(FmsFieldDescriptorSetFixedOrder(diff_name, FMS_CONTINUOUS, FMS_NODAL_GAUSS_CLOSED, 3), 0);
    EXPECT_NE(FmsFieldDescriptorCompare(fdv, diff_name), 0);
    EXPECT_NE(FmsFieldDescriptorCompare(diff_name, fdv), 0);

    FmsDataCollectionDestroy(&dc);
}
