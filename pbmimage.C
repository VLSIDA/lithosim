#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <set>

#include "pbmimage.h"
#include "helper.h"

pbm_image_t::pbm_image_t(int w, int h) { // malloc "image" as 2D array of float.
  width=w;
  widthbytes=(int)(ceil(width/8.0));
  height=h;
  image=(unsigned char**)malloc(sizeof(unsigned char*)*h);
    for(int i=0;i<h;i++) {
      /* 8 bit entries per char. 
	 if width is 8 we have 1 
	 if width is 9 we have 2.*/
     image[i] = (unsigned char*)malloc(sizeof(unsigned char)*(int)(ceil(w/8)));
    }  
}

pbm_image_t::pbm_image_t(char *filename) {  
  FILE * pbmfile = fopen(filename,"r");
  if(!pbmfile) { printf("Error: can't open %s\n",filename); exit(-1); }
  /* get the magic number*/
  char type[10];   /* read the pnm file of the type specified in the file.*/
  fscanf(pbmfile,"%s",(char*)&type);
  width=0;
  height=0;
  if(!strcmp(type,"P4")) { /* PBM binary*/
    fscanf(pbmfile, "%d %d", &width, &height); /* read the width and the height.*/
    widthbytes=(int)(ceil(width/8.0));
    eatspace(pbmfile);
    /* create a pointer to the array of image width for this height.*/
    image = (unsigned char**)malloc(sizeof(unsigned char*)*height);   
    /* make sure any extra bits will end up as 0's since they are don't care inputs */
    /*    int remainder=(widthbytes*8 - width);
	  unsigned char remainder_mask = (1<<remainder)-1;*/
    for(int i=0 ; i<height; i++) {
      image[i]=(unsigned char*)malloc(sizeof(unsigned char)*widthbytes);
      fread(image[i],1,widthbytes,pbmfile);  
      /*image[i][widthbytes-1]|=remainder_mask;
      for(int j=0;j<widthbytes;j++) 
      image[i][j]=~(image[i][j]);*/
      /* for(int j=0;j<widthbytes;j++) 
	printf("%02x;",image[i][j]);
	printf("\n");*/

    }
  } else {
    printf("Not a PBM binary file.\n");
    exit(-1);
  }
  fclose(pbmfile);
  printf("Loaded %d x %d pbm: %s\n",width,height,filename);
}


/* returns the number of neighbors that are a different color */

pbm_image_t::~pbm_image_t() {
  for(int i=0;i<height;i++) 
    free(image[i]);
  free(image);
}


pbm_image_t&
pbm_image_t::operator=(const pbm_image_t& f) {
  assert(width==f.width);
  assert(height==f.height);
  for(int i=0;i<height;i++)
    for(int j=0;j<width/8;j++)
      image[i][j]=f.image[i][j];
}

// returns the number of neighbors that are a different color 
int 
pbm_image_t::singleton(int x, int y) const {   
  int cnt=0;
  int val = get_extended(x,y);
  if (val == get_extended(x+1,y)) cnt++;
  if (val == get_extended(x-1,y)) cnt++;
  if (val == get_extended(x,y+1)) cnt++;
  if (val == get_extended(x,y-1)) cnt++;
  return(cnt);
}

int 
pbm_image_t::neighbors(int x, int y) const {   
  int cnt=0;
  int val = get_extended(x,y);
  if (val == get_extended(x+1,y)) cnt++;
  if (val == get_extended(x-1,y)) cnt++;
  if (val == get_extended(x+1,y+1)) cnt++;
  if (val == get_extended(x-1,y+1)) cnt++;
  if (val == get_extended(x+1,y-1)) cnt++;
  if (val == get_extended(x-1,y-1)) cnt++;
  if (val == get_extended(x,y+1)) cnt++;
  if (val == get_extended(x,y-1)) cnt++;

  return(cnt);
}

std::set<std::pair<int,int> > 
pbm_image_t::contour_bits() {
  std::set<std::pair<int,int> > bits;
  // a contour bit has at least one white and one black neighbor
  int cnt=0;
  for(int i=0 ; i<height; i++) {
    for(int j=0;j<width;j++) {
      cnt=singleton(j,i);
      if (cnt>0 && cnt<4)
	bits.insert(std::pair<int,int>(j,i));
    }
  }

  return(bits);
}

// pick a bit with random probability depending on neighbors
// out of 8
// 8 = 1/10
// 7 = 2/10
// 6 = 3/10
// 5 = 4/10
// 4 = 5/10
// 3 = 6/10
// 2 = 7/10
// 1 = 8/10
// 0 = 9/10
void
pbm_image_t::pick_adjacent_bit(int *x, int *y) const {
  do {
    *x=(int)rand()%width;
    *y=(int)rand()%height;

  } while (rand()<(9.0-neighbors(*x,*y)*0.1));
}

// pick an x,y such that the whole region is 1 or 0 for inverting
bool
pbm_image_t::is_iso_region(int x, int y, int xw, int yw) const {
  int val=get(x,y);
  for(int xx=x;xx<x+xw;xx++) {
    for(int yy=y;yy<y+yw;yy++) {
      if (get_extended(xx,yy)!=val)
	return(false);
    }
  }
  return(true);
}

// pick a random on bit
void
pbm_image_t::pick_on_bit(int *x, int *y) const {
  do {
    *x=(int)rand()%width;
    *y=(int)rand()%height;
  }    while (!get(*x,*y));
}

// pick a random region that is
// 1. of size xw * yw
// 2. all set to the same value (on or off)
// 3. within range of a bit that is in error
void
pbm_image_t::pick_iso_region(int *x, int *y, int xw, int yw, const pbm_image_t &diff, int range) const {
  int deltax, deltay;
  do {
    deltax=(int)(rand()%range + range/2);
    deltay=(int)(rand()%range + range/2);
    diff.pick_on_bit(x,y);

  } while ((*x+deltax >= width) || (*y+deltay >= height) || (*x+deltax<0) || (*y+deltay<0) || !is_iso_region(*x,*y,xw,yw));
}

// pick a random region that is
// 1. of size xw * yw
// 2. all set to the same value (on or off)
void
pbm_image_t::pick_iso_region(int *x, int *y, int xw, int yw) const {
  do {
    *x=(int)(rand()%(width-xw));
    *y=(int)(rand()%(height-yw));

  } while (!is_iso_region(*x,*y,xw,yw));
}


void
pbm_image_t::invert() {
  for(int i=0 ; i<height; i++) {
    for(int j=0;j<widthbytes;j++) 
      image[i][j]=~(image[i][j]);
  }
}

// this updates the convolution of a single bit 
void
pbm_image_t::single_convolve(const pbm_image_t &src,const flt_image_t &k, int x, int y) {

  int xradius=k.get_width()/2;
  int yradius=k.get_height()/2;

      // initialize result for the intensity
      float temp=0.0;
      // for all the bits in the kernel 
      for(int i2=-yradius;i2<=yradius;i2++) {
	for(int j2=-xradius;j2<=xradius;j2++) {
          if (src.get_extended(x+j2,y+i2))
	    temp += k.get(j2+xradius-1,i2+yradius-1);
	}
      }

      //temp=get(src,j,i);
      // must do square to get intensity 
      if (temp*temp > 0.2)
	set(x,y,1);
      else
	set(x,y,0);


}

// this updates the entire region affected by this one bit 
// also the bits +w in x and y direction
// w should be >= 1
void
pbm_image_t::incremental_convolve(const pbm_image_t &src,const flt_image_t &k, int x, int y, int wx, int wy) {
  int xradius=k.get_width()/2;
  int yradius=k.get_height()/2;

  int xmin=(x-xradius>0)?x-xradius:0;
  int xmax=(x+xradius+wx<src.get_width())?x+xradius+wx:src.get_width();
  int ymin=(y-yradius>0)?y-yradius:0;
  int ymax=(y+yradius+wy<src.get_height())?y+yradius+wy:src.get_height();

  for(int i=ymin;i<ymax;i++) 
    for (int j=xmin;j<xmax;j++) 
      single_convolve(src,k,j,i);
}


void
pbm_image_t::convolve(const pbm_image_t &src,const flt_image_t &k) {
  assert(width==src.get_width());
  assert(height=src.get_height());
//   // We should consider the size of the mask, in calculating the divider. 
//   int y_divider = 4;   // split the input array into this many.  e.g. 256 squares. 
//   int x_divider = 2;   // split the input array into this many.  e.g. 256 squares. 
//   int y_part_size, y_remain_size, x_part_size, x_remain_size; 
//   y_part_size = src.get_height() / y_divider;  y_remain_size = src.get_height() % y_divider; 
//   x_part_size = src.get_width() /  x_divider;  x_remain_size = src.get_width()  % x_divider; 

//   // To partition into 256 problems - we split it into 16 x 16 squares. 
//   // To partition into 8 problems we do:  4 x 2 rectangles. 
//   // Create the partition #'s: 
//   int x_start = 0; 
//   int x_end = -1; 
//   int y_start = 0;  
//   int y_end = -1; 
  
//   for (int xcnt=1; xcnt<=x_divider; xcnt+=1) { 
//      x_start = x_end + 1; 
//      x_end += x_part_size;  
//      for (int ycnt=1; ycnt<=x_divider; ycnt+=1) {  
//         y_start = y_end + 1; 
//         y_end += y_part_size;  
//         printf("launch: x:(%d,%d), y:(%d,%d)\n", x_start, x_end, y_start, y_end); 
//      } 
//      if (y_remain_size) { 
//         printf("Deal with the y remainders\n"); 
//      } 
//   }; 

//   if (x_remain_size) { 
//         printf("Deal with the x remainders\n"); 
//   } 

  // for all the bits in the array 
  for(int i=0;i<src.get_height();i++) {
    for (int j=0;j<src.get_width();j++) {
      single_convolve(src,k,j,i);
    }
  }
}



void
pbm_image_t::diff(const pbm_image_t &p1, const pbm_image_t &p2) {
  unsigned char val;
  for(int i=0 ; i<height ; i++) 
  for(int j=0 ; j<widthbytes ; j++) 
    image[i][j]=p1.image[i][j]^p2.image[i][j];
}

void
pbm_image_t::save(char *filename) {
  FILE * out = fopen(filename,"w");
  if(!out) {
    printf("Error: can't open %s\n",filename);
    exit(-1);
  }
  
  fprintf(out,"P4\n");
  fprintf(out, "%d %d\n", width, height);
  for(int i=0 ; i<height ; i++) {
    fwrite(image[i],1,widthbytes,out);
  }
  fclose(out);
  printf("Saved %d x %d raw pbm: %s\n",width,height,filename);
}


int
pbm_image_t::count() const {
  int total=0;
  for(int i=0 ; i<height ; i++) 
    for(int j=0 ; j<width ; j++) {
      if(get(j,i))
	total++;
    }
  return(total);
}

// this computes the complexity heuristic of the image
//to reduce the number of isolated or protruding shapes.
// it computes the number of x and y transitions  
int
pbm_image_t::complexity() const {
  int val=0;
  for(int i=1 ; i<height ; i++) 
    for(int j=1 ; j<width ; j++) {
      val+=(4-singleton(j,i));
    }
  //  printf("xtrans %d ytrans %d\n",xtrans,ytrans);
  return(val/2);
}



// returns a 1 or a 0.
int 
pbm_image_t::get(int x, int y) const {   
  unsigned char abyte = image[y][x/8];	// get the corresponding byte.
  int val=(abyte >> (7 - (x % 8))) & 1; 	// get the bit out of the byte. 
  return(val);
}

int 
pbm_image_t::get_extended(int x, int y) const {   
  // extends the edge values before smallest index, and beyond largest index.
  int realx=x;
  if (x<0) 		realx=0;
  else if (x>=width)	realx=width-1;
  int realy=y;
  if (y<0) 		realy=0;
  else if (y>=height)	realy=height-1;
  unsigned char abyte = image[realy][realx/8];	// get the corresponding byte.
  int val=(abyte >> (7 - (realx % 8))) & 1; 	// get the bit out of the byte. 
  return(val);
}


// flips polarity of a given bit region 
// to the opposite of the base bit (x,y)
// note that w should be >= 1
// 1 means a single pixel
void
pbm_image_t::flip(int x, int y, int wx, int wy) {   
    for (int xx=x;xx<x+wx;xx++) 
        for (int yy=y;yy<y+wy;yy++) 
            if( ((width > xx) && (height > yy)) )
            //if( ((this.width > xx) && (this.height > yy)) )
                flip(xx,yy);
}

// flips the polarity of a single bit
void
pbm_image_t::flip(int x, int y) {   
      int bitnum = 7 - (x % 8);
      unsigned char bitmask = 1<<bitnum;
      image[y][x/8] ^= bitmask;
}

// flips polarity of a given bit 
void
pbm_image_t::set(int x, int y, int val) { 
  int bitnum = 7 - (x % 8);
  unsigned char bitmask = 1<<bitnum;
  if (val)
    image[y][x/8] |= bitmask;
  else
    image[y][x/8] &= ~bitmask;
}

int
pbm_image_t::diff(const pbm_image_t &b) {
  int diff=0;
  for(int i=0 ; i<height ; i++) 
    for(int j=0 ; j<width ; j++) 
      if (get(j,i)!=b.get(j,i)) {
	diff++;
      }
  return(diff);
}
