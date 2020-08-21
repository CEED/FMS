#include <fms.h>
#include <fmsio.h>
#include <fmstesting.hpp>
#include <gtest/gtest.h>

TEST(FmsIO, WriteRead2DAscii) {
    FmsDataCollection dc0 = NULL;
    ASSERT_EQ(Construct2DData0(&dc0), 0);
    ASSERT_TRUE(dc0);

    // No point in continuing if the write fails
    ASSERT_EQ(FmsIOWrite("TestWrite2DAscii.fms", "ascii", dc0), 0);

    FmsDataCollection read_dc0 = NULL;
    // No point in continuing if the read fails
    ASSERT_EQ(FmsIORead("TestWrite2DAscii.fms", "ascii", &read_dc0), 0);
    ASSERT_TRUE(read_dc0);
    EXPECT_EQ(FmsDataCollectionCompare(dc0, read_dc0), 0);

    FmsDataCollectionDestroy(&dc0);
    FmsDataCollectionDestroy(&read_dc0);
}

TEST(FmsIO, WriteRead2DYaml) {
    FmsDataCollection dc0 = NULL;
    ASSERT_EQ(Construct2DData0(&dc0), 0);
    ASSERT_TRUE(dc0);

    // No point in continuing if the write fails
    ASSERT_EQ(FmsIOWrite("TestWrite2DYaml.fms", "yaml", dc0), 0);

    FmsDataCollection read_dc0 = NULL;
    // No point in continuing if the read fails
    ASSERT_EQ(FmsIORead("TestWrite2DYaml.fms", "yaml", &read_dc0), 0);
    ASSERT_TRUE(read_dc0);
    EXPECT_EQ(FmsDataCollectionCompare(dc0, read_dc0), 0);

    FmsDataCollectionDestroy(&dc0);
    FmsDataCollectionDestroy(&read_dc0);
}

TEST(FmsIO, WriteRead2DJson) {
    FmsDataCollection dc0 = NULL;
    ASSERT_EQ(Construct2DData0(&dc0), 0);
    ASSERT_TRUE(dc0);

    // No point in continuing if the write fails
    ASSERT_EQ(FmsIOWrite("TestWrite2DJson.fms", "json", dc0), 0);

    FmsDataCollection read_dc0 = NULL;
    // No point in continuing if the read fails
    ASSERT_EQ(FmsIORead("TestWrite2DJson.fms", "json", &read_dc0), 0);
    ASSERT_TRUE(read_dc0);
    EXPECT_EQ(FmsDataCollectionCompare(dc0, read_dc0), 0);

    FmsDataCollectionDestroy(&dc0);
    FmsDataCollectionDestroy(&read_dc0);
}

TEST(FmsIO, WriteRead2DHdf5) {
    FmsDataCollection dc0 = NULL;
    ASSERT_EQ(Construct2DData0(&dc0), 0);
    ASSERT_TRUE(dc0);

    // No point in continuing if the write fails
    ASSERT_EQ(FmsIOWrite("TestWrite2DHdf5.fms", "hdf5", dc0), 0);

    FmsDataCollection read_dc0 = NULL;
    // No point in continuing if the read fails
    ASSERT_EQ(FmsIORead("TestWrite2DHdf5.fms", "hdf5", &read_dc0), 0);
    ASSERT_TRUE(read_dc0);
    EXPECT_EQ(FmsDataCollectionCompare(dc0, read_dc0), 0);

    FmsDataCollectionDestroy(&dc0);
    FmsDataCollectionDestroy(&read_dc0);
}

#if 0
TEST(FmsIO, WriteRead2DSilo) {
    FmsDataCollection dc0 = NULL;
    ASSERT_EQ(Construct2DData0(&dc0), 0);
    ASSERT_TRUE(dc0);

    // No point in continuing if the write fails
    ASSERT_EQ(FmsIOWrite("TestWrite2DSilo.fms", "silo", dc0), 0);

    FmsDataCollection read_dc0 = NULL;
    // No point in continuing if the read fails
    ASSERT_EQ(FmsIORead("TestWrite2DSilo.fms", "silo", &read_dc0), 0);
    ASSERT_TRUE(read_dc0);
    EXPECT_EQ(FmsDataCollectionCompare(dc0, read_dc0), 0);

    FmsDataCollectionDestroy(&dc0);
    FmsDataCollectionDestroy(&read_dc0);
}
#endif
