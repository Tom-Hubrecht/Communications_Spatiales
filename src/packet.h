/*
Defines the standard space packet as in CCSDS 113.0-B-1
*/

typedef struct SPACE_PACKET SPACE_PACKET;
typedef struct SPACE_PACKET_PRIMARY_HEADER SPACE_PACKET_PRIMARY_HEADER;
typedef struct PACKET_DATA_FIELD PACKET_DATA_FIELD;
typedef unsigned int uint;


struct SPACE_PACKET_PRIMARY_HEADER
{
    uint PACKET_VERSION_NUMBER : 3;
    // PACKET IDENTIFICATION
    uint PACKET_TYPE : 1;
    uint SECONDARY_HEADER_FLAG : 1;
    uint APID : 11;
    // PACKET SEQUENCE CONTROL
    uint SEQUENCE_FLAGS : 2;
    uint PACKET_SEQUENCE_COUNT : 14;
    uint PACKET_DATA_LENGTH : 16;
}


struct PACKET_DATA_FIELD
{

}


struct SPACE_PACKET
{
    SPACE_PACKET_PRIMARY_HEADER PRIMARY_HEADER;

}
