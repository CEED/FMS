#include <fms.h>
#include <stdio.h>

int main(int, char**)
{
    FmsMesh m;
    auto a = FmsMeshConstruct(&m);
    printf("Hello %d\n", a);
    return a;
}