/*
Defines the standard space packet as in CCSDS 113.0-B-1
*/

#include "bit_array.h"
#include "bar.h"

typedef struct SPACE_PACKET SPACE_PACKET;
typedef struct SPACE_PACKET_PRIMARY_HEADER SPACE_PACKET_PRIMARY_HEADER;
typedef struct PACKET_DATA_FIELD PACKET_DATA_FIELD;
typedef bar DATA_FIELD;
typedef unsigned int uint;


struct SPACE_PACKET_PRIMARY_HEADER
{
    uint PACKET_VERSION_NUMBER : 3; // Should always be 000
    // PACKET IDENTIFICATION
    uint PACKET_TYPE : 1; // 0 for telemetry, 1 for telecommand
    uint SECONDARY_HEADER_FLAG : 1; // Always 0 for now
    uint APID : 11;
    // PACKET SEQUENCE CONTROL
    uint SEQUENCE_FLAGS : 2;
    uint PACKET_SEQUENCE_COUNT : 14;
    uint PACKET_DATA_LENGTH : 16;
};

struct SPACE_PACKET
{
    SPACE_PACKET_PRIMARY_HEADER PRIMARY_HEADER;
    DATA_FIELD DATA_FIELD;

};

void packet_request(SPACE_PACKET packet, uint APID);
SPACE_PACKET *create_packet(uint PACKET_TYPE, uint APID, DATA_FIELD DATA_FIELD);
