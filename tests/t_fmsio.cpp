#include <fms.h>
#include <fmsio.h>
#include "fmstesting.hpp"
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
    int isNotSupported = 0;
    for(FmsInt i = 0; i < size; i++) {
        ASSERT_TRUE(dcs[i]);
        // It's possible that Conduit is compiled without hdf5 support
        int retval = FmsIOWrite("TestHdf5.fms", "hdf5", dcs[i]);
        if(retval == 2) {
            isNotSupported = 1;
            break;
        }
        ASSERT_EQ(retval, 0);
        FmsDataCollection read_dc = NULL;
        ASSERT_EQ(FmsIORead("TestHdf5.fms", "hdf5", &read_dc), 0);
        ASSERT_TRUE(read_dc);
        EXPECT_EQ(FmsDataCollectionCompare(dcs[i], read_dc), 0);
        FmsDataCollectionDestroy(&read_dc);
    }
    for(FmsInt i = 0; i < size; i++) {
        FmsDataCollectionDestroy(&dcs[i]);
    }
    std::free(dcs);
    if(isNotSupported) GTEST_SKIP();
}

TEST(FmsIO, WriteReadConduitBin) {
    FmsInt size = 0;
    FmsDataCollection *dcs = NULL;
    ASSERT_EQ(ConstructAllDataCollections(&size, &dcs), 0);
    ASSERT_TRUE(dcs);
    ASSERT_TRUE(size);
    for(FmsInt i = 0; i < size; i++) {
        ASSERT_TRUE(dcs[i]);
        ASSERT_EQ(FmsIOWrite("TestConduitBin.fms", "conduit_bin", dcs[i]), 0);
        FmsDataCollection read_dc = NULL;
        ASSERT_EQ(FmsIORead("TestConduitBin.fms", "conduit_bin", &read_dc), 0);
        ASSERT_TRUE(read_dc);
        EXPECT_EQ(FmsDataCollectionCompare(dcs[i], read_dc), 0);
        ASSERT_EQ(FmsDataCollectionDestroy(&dcs[i]), 0);
        ASSERT_EQ(FmsDataCollectionDestroy(&read_dc), 0);
    }
    std::free(dcs);
}

TEST(FmsIO, WriteReadAutoDetect) {
    FmsDataCollection dc = NULL;
    ASSERT_EQ(Construct2DData0(&dc), 0);
    ASSERT_TRUE(dc);
    const char *protocols[5] = {
        "ascii", "json", "yaml", "hdf5", "conduit_bin"
    };
    for(int i = 0; i < 5; i++) {
        int retval = FmsIOWrite("TestAuto.fms", protocols[i], dc);
        // Sometimes hdf5 will not be supported by conduit
        if(retval == 2 && strcmp("hdf5", protocols[i]) == 0) {
            continue;
        }
        ASSERT_EQ(retval, 0);
        FmsDataCollection read_dc = NULL;
        ASSERT_EQ(FmsIORead("TestAuto.fms", NULL, &read_dc), 0);
        ASSERT_TRUE(read_dc);
        EXPECT_EQ(FmsDataCollectionCompare(dc, read_dc), 0);
        EXPECT_EQ(FmsDataCollectionCompare(read_dc, dc), 0);
        FmsDataCollectionDestroy(&read_dc);
    }
    FmsDataCollectionDestroy(&dc);
}
