#include <sys/time.h>


struct tictac {
        struct timeval base;

        void (*start)(struct tictac *); 
        double (*stop)(struct tictac *);
};


void tictac_initialise(struct tictac * tictac);
