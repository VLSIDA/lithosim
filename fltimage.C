#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "helper.h"
#include "fltimage.h"

flt_image_t::flt_image_t(int w, int h) {
  height=h;
  width=w;
  image=(float**)malloc(sizeof(float*)*h);
  for(int i=0;i<h;i++) {
    image[i] = (float*)malloc(sizeof(float)*w);
  }
}

flt_image_t::flt_image_t(char *filename) {
  FILE * pbmfile = fopen(filename,"r");
  
  if(!pbmfile) {
    printf("Error: can't open %s\n",filename);
    exit(-1);
  }

  // get the magic number
  char type[10];
  fscanf(pbmfile,"%s",(char*)&type);
  if(!strcmp(type,"P4")) { // PBM binary
    fscanf(pbmfile, "%d %d", &width, &height);
    eatspace(pbmfile);
    image = (float**)malloc(sizeof(float*)*height);
    int numbytes=(int)(ceil(width/8.0));
    unsigned char buff[numbytes];

    for(int i=0 ; i<height; i++) {
      fread(buff,1,numbytes,pbmfile);
      //      for(int j=0;j<numbytes;j++)
      //	printf("%02x;",buff[j]);
      //	printf("\n");
      image[i]=unpack((unsigned char*)buff,numbytes);
    }
  } else {
    printf("Not a PBM binary file.\n");
    exit(-1);
  }
  fclose(pbmfile);
  printf("Loaded %d x %d pbm: %s\n",width,height,filename);
}

flt_image_t::~flt_image_t() {
  for(int i=0;i<height;i++) 
    free(image[i]);
  free(image);
}

void
flt_image_t::save_pgm(char *filename) {
  FILE * out = fopen(filename,"w");
  
  if(!out) {
    printf("Error: can't open %s\n",filename);
    exit(-1);
  }

  float max=0;
  float min=0; // only care about negatives
  for(int i=0 ; i<height ; i++) 
    for(int j=0;j<width;j++) {
      max=image[i][j]>max?image[i][j]:max;
      min=image[i][j]<min?image[i][j]:min;
    }
  
  fprintf(out,"P2\n");
  fprintf(out, "%d %d\n", width, height);
  fprintf(out,"255\n");
  for(int i=0 ; i<height ; i++,fprintf(out,"\n")) {
    for(int j=0;j<width;j++) {
      fprintf(out,"%d ",(int)((( (image[i][j]+(fabs(min))) / (max-min)))*255));
    }
  }

  fclose(out);
  printf("Saved %d x %d raw pbm: %s\n",width,height,filename);
}

void
flt_image_t::save_ascii(char *filename,float threshold) {
  FILE * out = fopen(filename,"w");
  
  if(!out) {
    printf("Error: can't open %s\n",filename);
    exit(-1);
  }
  
  fprintf(out,"P1\n");
  fprintf(out, "%d %d\n", width, height);
  for(int i=0 ; i<height ; i++, fprintf(out,"\n")) {
    for(int j=0 ; j<width  ; j++) {
      if (j!=0&&j%70==0) fprintf(out,"\n");
      fprintf(out, "%d",(image[i][j]>threshold));
    }
  }
  fclose(out);
  printf("Saved %d x %d ascii pbm: %s\n",width,height,filename);
}


void
flt_image_t::make_jinc(int P, float scale) {
  if (P%2==0) P++; // make sure it is odd so there is a center 
  width=height=P;
  image=(float**)malloc(sizeof(float*)*height);
  for(int i=0;i<height;i++)
    image[i]=(float*)malloc(sizeof(float)*width);
  //     x = [-floor(P/2) : 1 : floor(P/2)];
  //  y = x;
  //  [xx,yy]=meshgrid(x,y);
  //  r = sqrt(xx.^2+yy.^2);
  int offset=P/2;
  for(int i=0;i<height;i++) 
    for(int j=0;j<width;j++) 
      image[i][j]=sqrt((i-offset)*(i-offset)+(j-offset)*(j-offset));


  //   z = ones(size(r));
  //  k = find(r);
  //  z(k) = besselj(1, 2*pi*r(k)*scale) ./ ( 2*pi*r(k)*scale ); 
  float twopi=2*3.14159;
  for(int i=0;i<height;i++)
    for(int j=0;j<width;j++)
      image[i][j]=bessj1(twopi*image[i][j]*scale) / (twopi*image[i][j]*scale);

  // z(ceil(P/2)+1, ceil(P/2)+1) = 0.5; 
  image[offset][offset]=0.5;

  //   H = z / sum(sum(z));
  //    normalize it
  float sum=0;
  for(int i=0;i<height;i++)
    for(int j=0;j<width;j++)
      sum+=image[i][j];

  for(int i=0;i<height;i++)
    for(int j=0;j<width;j++)
      image[i][j]=(image[i][j]/sum);

  //   for(int i=0;i<height;i++,printf("\n"))
  //   for(int j=0;j<width;j++)
  //    printf("%f ",image[i][j]);
  //    printf("\n");

}


void
flt_image_t::convolve(pbm_image_t &src,flt_image_t &k) {
  assert(width==src.get_width());
  assert(height=src.get_height());

  for(int i=0;i<src.get_height();i++) {
    for (int j=0;j<src.get_width();j++) {

      image[i][j]=0.0;

      // for every point in the filter range
      for(int i2=-k.get_height()/2;i2<k.get_height()/2+1;i2++) {
	for(int j2=-k.get_width()/2;j2<k.get_width()/2+1;j2++) {
	  // the corresponding source image is offset by half the size of the filter 
          if (src.get_extended(j+j2,i+i2)) {
	    image[i][j]+= k.get(j2+k.get_width()/2-1,i2+k.get_height()/2-1);
          };
	}
      }

      // power itensity is the square 
      image[i][j]=image[i][j]*image[i][j];

    }
  }
}

void
flt_image_t::pack(float *barray,int n,unsigned char *str,float threshold) {
  int cnt=0;
  int b=0;
  unsigned char cur;
  while(b<n) {
    cur=0;
    if (barray[b++]>threshold) cur+=128;
    if (b<n && barray[b++]>threshold) cur+=64;
    if (b<n && barray[b++]>threshold) cur+=32;
    if (b<n && barray[b++]>threshold) cur+=16;
    if (b<n && barray[b++]>threshold) cur+=8;
    if (b<n && barray[b++]>threshold) cur+=4;
    if (b<n && barray[b++]>threshold) cur+=2;
    if (b<n && barray[b++]>threshold) cur+=1;
    str[cnt]=cur;
    cnt++;
  }
}

float *
flt_image_t::unpack(unsigned char *str,int n) {
  float *barray=(float*)malloc(sizeof(float)*8*n);
  for(int b=0;b<8*n;b++)
	barray[b]=1.0; 
  int cnt=0;
  unsigned char cur;
  for(int b=0;b<n;b++) {
    cur=str[b];
     if (0x80 & cur) barray[cnt]=0.0;
     cnt++;
     if (0x40 & cur) barray[cnt]=0.0;
     cnt++;
     if (0x20 & cur) barray[cnt]=0.0;
     cnt++;
     if (0x10 & cur) barray[cnt]=0.0;
     cnt++;
     if (0x08 & cur) barray[cnt]=0.0;
     cnt++;
     if (0x04 & cur) barray[cnt]=0.0;
     cnt++;
     if (0x02 & cur) barray[cnt]=0.0;
     cnt++;
     if (0x01 & cur) barray[cnt]=0.0;
     cnt++;
  }
  return(barray);
}





float 
flt_image_t::get(int x, int y) const { 
  int realx=x;
  if (x<0) 
    realx=0;
  else if (x>=width)
    realx=width-1;

  int realy=y;
  if (y<0) 
    realy=0;
  else if (y>=height)
    realy=height-1;

  return(image[realy][realx]);
}

void
flt_image_t::save(char *filename,float threshold) {
  FILE * out = fopen(filename,"w");
  if(!out) {
    printf("Error: can't open %s\n",filename);
    exit(-1);
  }
  
  fprintf(out,"P4\n");
  fprintf(out, "%d %d\n", width, height);
  for(int i=0 ; i<height ; i++) {
    int numbytes=(int)(ceil(width/8.0));
    unsigned char buff[numbytes];
    pack(image[i],width,(unsigned char*)&buff,threshold);
    fwrite(buff,1,numbytes,out);
  }

  fclose(out);
  printf("Saved %d x %d raw pbm: %s\n",width,height,filename);
}
