#include <fms.h>
#include <fmsio.h>
#include <fmstesting.hpp>
#include <gtest/gtest.h>

TEST(FmsIO, WriteReadAscii) {
    FmsInt size = 0;
    FmsDataCollection *dcs = NULL;
    ASSERT_EQ(ConstructAllDataCollections(&size, &dcs), 0);
    ASSERT_TRUE(dcs);
    ASSERT_TRUE(size);
    for(FmsInt i = 0; i < size; i++) {
        ASSERT_TRUE(dcs[i]);
        ASSERT_EQ(FmsIOWrite("TestAscii.fms", "ascii", dcs[i]), 0);
        FmsDataCollection read_dc = NULL;
        ASSERT_EQ(FmsIORead("TestAscii.fms", "ascii", &read_dc), 0);
        ASSERT_TRUE(read_dc);
        EXPECT_EQ(FmsDataCollectionCompare(dcs[i], read_dc), 0);
        ASSERT_EQ(FmsDataCollectionDestroy(&dcs[i]), 0);
        ASSERT_EQ(FmsDataCollectionDestroy(&read_dc), 0);
    }
    std::free(dcs);
}

TEST(FmsIO, WriteReadYaml) {
    FmsInt size = 0;
    FmsDataCollection *dcs = NULL;
    ASSERT_EQ(ConstructAllDataCollections(&size, &dcs), 0);
    ASSERT_TRUE(dcs);
    ASSERT_TRUE(size);
    for(FmsInt i = 0; i < size; i++) {
        ASSERT_TRUE(dcs[i]);
        ASSERT_EQ(FmsIOWrite("TestYaml.fms", "yaml", dcs[i]), 0);
        FmsDataCollection read_dc = NULL;
        ASSERT_EQ(FmsIORead("TestYaml.fms", "yaml", &read_dc), 0);
        ASSERT_TRUE(read_dc);
        EXPECT_EQ(FmsDataCollectionCompare(dcs[i], read_dc), 0);
        ASSERT_EQ(FmsDataCollectionDestroy(&dcs[i]), 0);
        ASSERT_EQ(FmsDataCollectionDestroy(&read_dc), 0);
    }
    std::free(dcs);
}

TEST(FmsIO, WriteReadJson) {
    FmsInt size = 0;
    FmsDataCollection *dcs = NULL;
    ASSERT_EQ(ConstructAllDataCollections(&size, &dcs), 0);
    ASSERT_TRUE(dcs);
    ASSERT_TRUE(size);
    for(FmsInt i = 0; i < size; i++) {
        ASSERT_TRUE(dcs[i]);
        ASSERT_EQ(FmsIOWrite("TestJson.fms", "json", dcs[i]), 0);
        FmsDataCollection read_dc = NULL;
        ASSERT_EQ(FmsIORead("TestJson.fms", "json", &read_dc), 0);
        ASSERT_TRUE(read_dc);
        EXPECT_EQ(FmsDataCollectionCompare(dcs[i], read_dc), 0);
        ASSERT_EQ(FmsDataCollectionDestroy(&dcs[i]), 0);
        ASSERT_EQ(FmsDataCollectionDestroy(&read_dc), 0);
    }
    std::free(dcs);
}

TEST(FmsIO, WriteReadHdf5) {
    FmsInt size = 0;
    FmsDataCollection *dcs = NULL;
    ASSERT_EQ(ConstructAllDataCollections(&size, &dcs), 0);
    ASSERT_TRUE(dcs);
    ASSERT_TRUE(size);
    for(FmsInt i = 0; i < size; i++) {
        ASSERT_TRUE(dcs[i]);
        ASSERT_EQ(FmsIOWrite("TestHdf5.fms", "hdf5", dcs[i]), 0);
        FmsDataCollection read_dc = NULL;
        ASSERT_EQ(FmsIORead("TestHdf5.fms", "hdf5", &read_dc), 0);
        ASSERT_TRUE(read_dc);
        EXPECT_EQ(FmsDataCollectionCompare(dcs[i], read_dc), 0);
        ASSERT_EQ(FmsDataCollectionDestroy(&dcs[i]), 0);
        ASSERT_EQ(FmsDataCollectionDestroy(&read_dc), 0);
    }
    std::free(dcs);
}

#if 0
TEST(FmsIO, WriteReadSilo) {
    FmsInt size = 0;
    FmsDataCollection *dcs = NULL;
    ASSERT_EQ(ConstructAllDataCollections(&size, &dcs), 0);
    ASSERT_TRUE(dcs);
    ASSERT_TRUE(size);
    for(FmsInt i = 0; i < size; i++) {
        ASSERT_TRUE(dcs[i]);
        ASSERT_EQ(FmsIOWrite("TestSilo.fms", "silo", dcs[i]), 0);
        FmsDataCollection read_dc = NULL;
        ASSERT_EQ(FmsIORead("TestSilo.fms", "silo", &read_dc), 0);
        ASSERT_TRUE(read_dc);
        EXPECT_EQ(FmsDataCollectionCompare(dcs[i], read_dc), 0);
        ASSERT_EQ(FmsDataCollectionDestroy(&dcs[i]), 0);
        ASSERT_EQ(FmsDataCollectionDestroy(&read_dc), 0);
    }
    std::free(dcs);
}
#endif
