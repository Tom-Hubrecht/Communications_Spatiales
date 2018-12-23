#include <stdio.h>

typedef unsigned int uint;

typedef struct primaryHeader {
    uint versionNumber : 3;
    uint type : 1;
    uint secondaryFlag : 1;
    uint apid : 11;
    uint sequenceFlags : 2;
    uint sequenceCount : 14;
} PrimaryHeader;

int main()
{

    PrimaryHeader header;
    header.versionNumber = 0;
    printf("%d\n", header.sequenceCoun);
    return 0;
}