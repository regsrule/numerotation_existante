#ifndef __TIME_MANAGER_H
#define __TIME_MANAGER_H

#include <sys/time.h>
typedef struct timezone timezone_t;
typedef struct timeval timeval_t;

static struct timeval _t1, _t2;
static struct timezone _tz;
static unsigned long _temps_residuel = 0;
#define top1() gettimeofday(&_t1, &_tz)
#define top2() gettimeofday(&_t2, &_tz)
unsigned long cpu_time(void);
void init_cpu_time(void);


void top_(timeval_t* tv, timezone_t* tz);
unsigned long get_temp_residuel(void);
unsigned long cpu_time_(timeval_t _t1,timeval_t _t2);
#endif
