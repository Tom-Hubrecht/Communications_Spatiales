#include <stdio.h>
#include <string.h>


#include <wiringPi.h>
#include "bit_array.h"
#include "bar.h"

#include "rpi.h"

void send(bar *message)
{
    char m;
    char r;

    wiringPiSetup();
    pinMode(0, OUTPUT);
    pinMode(29, INPUT);

    size_t n = barlen(message);
    for(size_t i = 0; i < n; i++)
    {
        m = barget(message, i);
        digitalWrite(0, m);
        r = digitalRead(29);
//        printf("%d -- %d\n", m, r);
//        delay(100);
    }

}
