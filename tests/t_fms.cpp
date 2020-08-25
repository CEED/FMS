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
    FmsMetaData md = NULL;
    FmsDataCollectionAttachMetaData()

    FmsDataCollectionDestroy(&dc);
}