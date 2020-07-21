#include <gtest/gtest.h>
#include <fms.h>

TEST(fms, version)
{
    FmsInt v = 0;
    ASSERT_EQ(FmsGetInterfaceVersion(&v), 0);
    EXPECT_EQ(FMS_INTERFACE_VERSION, v);
}
