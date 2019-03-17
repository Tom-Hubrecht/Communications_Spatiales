#include <wiringPi.h>
#include "bit_array.h"
#include "bar.h"

void send(bar *message)
{
    wiringPiSetup();
    pinMode(0, OUTPUT);

    size_t n = barlen(message);
    for(size_t i = 0; i < n; i++)
    {
        digitalWrite(0, barget(message, i));
        delay(100);
    }

}