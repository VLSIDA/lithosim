
#ifndef PBMIMAGE_H_
#define PBMIMAGE_H_

#include <set>

#include "fltimage.h"
class flt_image_t;
class pbm_image_t;

class pbm_image_t {
 public:
  pbm_image_t(int w, int h);
  pbm_image_t(char *filename);
  ~pbm_image_t();

  void flip(int x, int y);
  void flip(int x, int y, int wx, int wy);
  void set(int x, int y, int val);

  int diff(const pbm_image_t &p);

  std::set<std::pair<int,int > > contour_bits();

  void pick_adjacent_bit(int *x, int *y) const;
  bool is_iso_region(int x, int y, int wx, int wy) const;
  void pick_iso_region(int *x, int *y, int xw, int yw) const;
  void pick_iso_region(int *x, int *y, int wx, int wy,const pbm_image_t &diff, int range) const;
  void pick_on_bit(int *x, int *y) const;

  void invert();
  void single_convolve(const pbm_image_t &src,const flt_image_t &k, int x, int y);
  void incremental_convolve(const pbm_image_t &src,const flt_image_t &k, int x, int y, int wx, int wy);
  void convolve(const pbm_image_t &src,const flt_image_t &k);
  void save(char *filename);
  void diff(const pbm_image_t &p1,const pbm_image_t &p2);

  int count() const;
  int complexity() const;
  int  singleton(int x, int y) const;
  int neighbors(int x, int y) const;
  int get(int x, int y) const;
  int get_extended(int x, int y) const;
  int get_height() const { return(height); };
  int get_width() const { return(width); };

  pbm_image_t &operator=(const pbm_image_t& f);

 private:
  int width, height;
  int widthbytes;
  unsigned char ** image;
};
#endif
