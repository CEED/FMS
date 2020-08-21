#include <fms.h>
#include <fmstesting.hpp>
#include <gtest/gtest.h>

TEST(Fms, Version)
{
    FmsInt v = 0;
    ASSERT_EQ(FmsGetInterfaceVersion(&v), 0);
    EXPECT_EQ(FMS_INTERFACE_VERSION, v);
}

// Compare interface

TEST(Fms, DataCollectionCompare) {
    EXPECT_EQ(FmsDataCollectionCompare(NULL,NULL),0);

    FmsDataCollection dc0 = NULL;
    ASSERT_EQ(Construct2DData0(&dc0), 0);
    ASSERT_TRUE(dc0);
    EXPECT_EQ(FmsDataCollectionCompare(dc0, dc0), 0);
    EXPECT_NE(FmsDataCollectionCompare(dc0, NULL), 0);
    EXPECT_NE(FmsDataCollectionCompare(NULL, dc0), 0);

    FmsDataCollection dc00 = NULL;
    ASSERT_EQ(Construct2DData0(&dc00), 0);
    ASSERT_TRUE(dc00);
    EXPECT_EQ(dc00, dc00);

    ASSERT_NE(dc0, dc00);
    EXPECT_EQ(FmsDataCollectionCompare(dc0, dc00), 0);
    FmsDataCollectionDestroy(&dc00); 

    FmsDataCollection dc1 = NULL;
    ASSERT_EQ(Construct2DData1(&dc1), 0);
    ASSERT_TRUE(dc1);
    EXPECT_NE(FmsDataCollectionCompare(dc0, dc1), 0);

    FmsDataCollection dc11 = NULL;
    ASSERT_EQ(Construct2DData1(&dc11), 0);
    ASSERT_TRUE(dc11);
    EXPECT_EQ(FmsDataCollectionCompare(dc1, dc11), 0);
    FmsDataCollectionDestroy(&dc11);

    FmsDataCollection dc2 = NULL;
    ASSERT_EQ(Construct3DData0(&dc2), 0);
    ASSERT_TRUE(dc2);
    EXPECT_NE(FmsDataCollectionCompare(dc0, dc2), 0);
    EXPECT_NE(FmsDataCollectionCompare(dc1, dc2), 0);

    FmsDataCollection dc3 = NULL;
    ASSERT_EQ(Construct3DData1(&dc3), 0);
    ASSERT_TRUE(dc3);
    EXPECT_NE(FmsDataCollectionCompare(dc0, dc3), 0);
    EXPECT_NE(FmsDataCollectionCompare(dc1, dc3), 0);
    EXPECT_NE(FmsDataCollectionCompare(dc2, dc3), 0);

    FmsDataCollectionDestroy(&dc0);
    FmsDataCollectionDestroy(&dc1);
    FmsDataCollectionDestroy(&dc2);
    FmsDataCollectionDestroy(&dc3);
}