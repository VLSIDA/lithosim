#ifndef FLTIMAGE_H_
#define FLTIMAGE_H_

#include "pbmimage.h"

class pbm_image_t;
class flt_image_t;

class flt_image_t {
public:
  flt_image_t() {width=height=0;};
  flt_image_t(char *filename);
  flt_image_t(int w, int h);
  ~flt_image_t();

  void save(char *filename,float threshold);
  void save_pgm(char *filename);
  void save_ascii(char *filename,float threshold);
  void make_jinc(int P, float scale);
  void convolve(pbm_image_t &src,flt_image_t &k);
  void pack(float *barray,int n,unsigned char *str,float threshold);
  float *unpack(unsigned char *str,int n);

  float get(int x, int y) const;
  int get_height() const { return(height); };
  int  get_width() const { return(width); };

private:
  int width,height;
  float ** image;
};
#endif
