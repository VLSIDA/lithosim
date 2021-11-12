#include <time.h>
/*########################################################
  Easy to use timer.
  ########################################################*/

class my_timer_t {
 public:
  my_timer_t() { stime = clock(); }
  void print_time(char * str) { 
    etime=clock();
    printf(str);
    printf(" TIME: %7.3f sec\n",((double) (etime - stime)) / CLOCKS_PER_SEC);
    stime = clock();  /* as a convenience. */
  }
 private:
 clock_t stime,etime;
};

