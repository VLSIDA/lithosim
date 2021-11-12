#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include <set>
#include <algorithm>

#include "timer.h"
#include "pbmimage.h"
#include "fltimage.h"

#define INITIAL_TEMP 50000
#define BEGIN_SCALE 0.95
#define MID_SCALE 0.99
#define END_SCALE 0.98
#define FINAL_TEMP 1.0
// how much more fixing an error is vs complexity to mask
#define ERROR_WEIGHT 1 //2.0
#define COMPLEXITY_WEIGHT 3 //1.0
#define CLIMB_WEIGHT 50.0

// Ian Lee
#define DEBUG 0

double 
schedule(double temp){
    if (temp > 1000)
        return(BEGIN_SCALE*temp);
    else if (temp > 100)
        return(MID_SCALE*temp);
    else
        return(END_SCALE*temp);
}


void
greedy(pbm_image_t &opced_mask, pbm_image_t &target, flt_image_t &kernel) {
    pbm_image_t contour(target.get_width(),target.get_height());
    //    pbm_image *contour2=create_pbm_image(original->width,original->height);
    contour.convolve(opced_mask,kernel);  
    int new_cost,old_cost=contour.diff(target);
    int accept=0,reject=0,counter=0;
    // find a list of "edge" contour bits
    std::set<std::pair<int,int> > bits = opced_mask.contour_bits();
    std::set<std::pair<int,int> >::iterator bitit;
    int bitnum,x,y;
    while (!bits.empty()) {
        // pick a random bit, remove it
        bitnum=(int)rand()%bits.size();
        int i=0;
        for(bitit=bits.begin();i<bitnum;bitit++,i++);
        x=bitit->first;
        y=bitit->second;
        bits.erase(bitit);
    
        // flip the bit
        opced_mask.flip(x,y);  
        contour.incremental_convolve(opced_mask,kernel,x,y,1,1);  
        new_cost=contour.diff(target);
        // if it improves things, keep it 
        if (new_cost>old_cost) {
            reject++;
            opced_mask.flip(x,y);  
            contour.incremental_convolve(opced_mask,kernel,x,y,1,1);  
        }
        else {
        //	printf("accept %d %d\n",x,y);

        accept++;
        old_cost=new_cost;
        // subtract the old contour bits if they became contained
        int cnt1=opced_mask.singleton(x-1,y);
        int cnt2=opced_mask.singleton(x,y-1);
        int cnt3=opced_mask.singleton(x+1,y);
        int cnt4=opced_mask.singleton(x,y+1);
        if (cnt1==0 || cnt1==4) {
	        bits.erase(std::pair<int,int>(x-1,y));
	        //	printf("remove %d %d\n",x-1,y);
        }
        if (cnt2==0 || cnt2==4) {
	        bits.erase(std::pair<int,int>(x,y-1));
	        //	printf("remove %d %d\n",x,y-1);
        }
        if (cnt3==0 || cnt3==4) {
	        bits.erase(std::pair<int,int>(x+1,y));
	        //	printf("remove %d %d\n",x+1,y);
        }
        if (cnt4==0 || cnt4==4) {
	        bits.erase(std::pair<int,int>(x,y+1));
	        //	printf("remove %d %d\n",x,y+1);
        }


        // add its new adjacent bits if they aren't contained
        if (x-1>=0 && cnt1>0 && cnt1<4) {
	        bits.insert(std::pair<int,int>(x-1,y));
	        //	printf("add %d %d\n",x-1,y);
        }
        if (y-1>=0 && cnt2>0 && cnt2<4) {
	        bits.insert(std::pair<int,int>(x,y-1));
	    // printf("add %d %d\n",x,y-1);
        }
        if (x+1<opced_mask.get_width() && cnt3>0 && cnt3<4) {
	        bits.insert(std::pair<int,int>(x+1,y));	
	        //	printf("add %d %d\n",x+1,y);
        }

        if (y+1<opced_mask.get_height() && cnt4>0 && cnt4<4) {
	        bits.insert(std::pair<int,int>(x,y+1));	
	        //	printf("add %d %d\n",x,y+1);
        }


        }

        if (counter++==20) {
            printf("size: %d\taccept: %d\treject: %d",bits.size(),accept,reject);
            printf("\tError: %d\tComplx: %d \tCOST: %d \n",
	        contour.diff(target),opced_mask.complexity(),old_cost);
            //      char name[20];
            //      sprintf(name,"status-%d.pnm",accept);
            //      opced_mask.save(name);

            counter=0;
        }
    }
}

void
anneal(pbm_image_t &opced_mask, pbm_image_t &target, flt_image_t &kernel,float min_feature) {

    pbm_image_t contour(target.get_width(),target.get_height());
    contour.convolve(opced_mask,kernel);  
    // instead of unconvolving, just keep a copy of the old contour
    pbm_image_t old_contour(target.get_width(),target.get_height());
    old_contour=contour;
    pbm_image_t diffmap(target.get_width(),target.get_height());
    diffmap.diff(old_contour,target);

    printf("Performing OPC...\n");
    // initialize temp & placement 
    float OPTIMAL_COMPLEXITY=opced_mask.complexity();
    float INITIAL_ERROR=contour.diff(target);
    printf("Initial complexity %f\n",OPTIMAL_COMPLEXITY);
    float old_cost=ERROR_WEIGHT*contour.diff(target)/INITIAL_ERROR + COMPLEXITY_WEIGHT;
    float new_cost;
    double temperature=INITIAL_TEMP;
    int inner_loop=0;
    float delta_c=0;
    int x,y,xw,yw;
    int counter=20;
    int accept=0,climb=0,reject=0;
    int statcnt=0;
    int courseness;
    int course_cnt=0;
    while (temperature > FINAL_TEMP) {

        if (temperature>5000)
            courseness=4;
        else if (temperature>1000)
            courseness=3;
        else if (temperature>100)
            courseness=2;
        else 
            courseness=1;


        if (counter==20) {
            printf("\nTEMP: %e\tsize: %2d  accept: %5d  climb: %5d  rej: %5d  ", temperature, courseness, accept, climb, reject);
            printf("Err: %4d  Cmpl: %4d  COST: %5.3f \n", old_contour.diff(target), opced_mask.complexity(), old_cost);
            counter=0;
            accept=climb=reject=0;
            //char name[20];
            //sprintf(name,"status-%d.pnm",statcnt++);
            //opced_mask.save(name);
        }
        else
            counter++;

        printf(".");
        fflush(stdout);

        while (inner_loop < 1000) {
            // pick the size of a region we want to flip
            // this picks a constant area, but changes the aspect ratio
            //      xw=(int)(rand()%courseness)+1;
            //      yw=(int)(courseness / xw);
            // find such a region
            //      opced_mask.pick_iso_region(&x,&y,xw,yw);
            //      opced_mask.pick_iso_region(&x,&y,xw,yw,diffmap,2*min_feature);
            
            #if DEBUG
                printf("\n\nanneal/inner_loop--%4d--", inner_loop);
                printf("x=%d y=%d    ", x, y);
            #endif
            xw=3;
            yw=3;
            opced_mask.pick_adjacent_bit(&x,&y);

            #if DEBUG
                printf("Pick %d %d",x,y);
            #endif

            // flip the bits
            opced_mask.flip(x,y,xw,yw);
            #if DEBUG
                printf("Flip %d %d +%d +%d\n",x,y,xw,yw);
            #endif

            // only incremental convolve if we can benefit from it
            contour.incremental_convolve(opced_mask,kernel,x,y,xw,yw);

            new_cost=ERROR_WEIGHT*contour.diff(target)/INITIAL_ERROR + COMPLEXITY_WEIGHT*opced_mask.complexity()/OPTIMAL_COMPLEXITY;
            delta_c = new_cost-old_cost;
	
            float r=(rand()/(float)RAND_MAX);
            if (delta_c < 0) {
                accept++;
                #if DEBUG
                    printf("Accept\told %f new %f delta %f\n",old_cost,new_cost,delta_c);
                #endif
	            old_contour=contour;
	            old_cost=new_cost;
	            diffmap.diff(old_contour,target);

            }
            else if (r < exp(-CLIMB_WEIGHT*INITIAL_TEMP*delta_c/temperature)) {
	            climb++;
                #if DEBUG
	                printf("Climbed (%f < %f)\n",r,exp(-CLIMB_WEIGHT*INITIAL_TEMP*delta_c/temperature));
	                printf("old: %f new %f delta_c %f\n",old_cost,new_cost,delta_c);
                #endif
	            old_contour=contour;
	            old_cost=new_cost;
	            diffmap.diff(old_contour,target);
                #if DEBUG
	                printf("GOTCHA - %d x=%d y=%d", inner_loop,x,y);
                #endif
            }
            else {
	            reject++;
                #if DEBUG
	                opced_mask.save("pre.pnm");
	                printf("Reject (%f < %f)\n",r,exp(-CLIMB_WEIGHT*INITIAL_TEMP*delta_c/temperature));
	                printf("old: %f new %f delta_c %f\n",old_cost,new_cost,delta_c);
                #endif
	            // flip the bit back if it got worse
	            opced_mask.flip(x,y,xw,yw);
	            //	opced_mask.save("post.pnm");
	            //	exit(1);
	            contour=old_contour; // undo the convolution

            }
            inner_loop++;
        }
        
        inner_loop=0;    
        temperature = schedule(temperature);
    }

    printf("\n");
    printf("Final:\n");
    printf("TEMP: %e\tsize: %2d  accept: %5d  climb: %5d  rej: %5d  ", temperature, courseness, accept, climb, reject);
    printf("Err: %4d  Cmpl: %4d  COST: %5.3f \n", old_contour.diff(target), opced_mask.complexity(), old_cost);
    printf("\n");
}


int 
main(int argc, char **argv) {

    srand(2);

    if (argc<3) 
    {
        printf("USAGE: %s <nm per pixel> <input pbm> [<aerial image>] <output contours> [<opced mask>]\n",argv[0]);
        printf("Using 5 arguments will generate aerial image.\n");
        printf("Using 6 arguments will run OPC.\n");
        printf("Suggested nm per pixel (given the tests in test directory):\n");
        printf("16 is ~180nm\n12 is ~130nm\n8 is ~90nm\n6 is ~65nm\n4 is ~45nm\n");
        exit(1);
    }

    int Resolution = atoi(argv[1]);    // nm per pixel in the input image 
    float NA = 0.95;          // numerical aperture 
    float lambda = 193;       // wavelength 
    float min_feature=(lambda/NA)/Resolution;
    int Filter_Size = (int)(min_feature);   // Jinc filter size to see first order assist bars
    my_timer_t t;


    pbm_image_t original(argv[2]);
    t.print_time("IMAGE READ");

    // make the jinc function 
    flt_image_t kernel;
    kernel.make_jinc(Filter_Size,(NA/lambda)*Resolution);
  
    printf("Convolution kernel of size: %d x %d\n", kernel.get_width(), kernel.get_height());
    kernel.save_pgm("jinc.pnm"); 
    t.print_time("KERNEL MADE");

    printf("Creating lithography contour...\n");
    if (argc<=4) {
        pbm_image_t contour(original.get_width(),original.get_height());
        contour.convolve(original,kernel);  // bit vector based.
        t.print_time("CONVOLUTION COMPLETE");
        contour.save(argv[3]); 
        t.print_time("CONTOUR IMAGE SAVED");
    }
    else if (argc<6) {
        // only generate the aerial image if we want it saved 
        printf("Generating aerial image...\n");
        flt_image_t aerial(original.get_width(),original.get_height());
        aerial.convolve(original,kernel);  
        t.print_time("CONVOLUTION COMPLETE");
        aerial.save_pgm(argv[3]); 
        t.print_time("AERIAL IMAGE SAVED");
        aerial.save(argv[4],0.2);
        t.print_time("CONTOUR IMAGE SAVED");
    }  
    else if (argc==6) {

        pbm_image_t opced_mask(argv[2]);    // start with the target as the mask 
        anneal(opced_mask,original,kernel,min_feature);
        //greedy(opced_mask,original,kernel);


        opced_mask.save(argv[5]); 
        t.print_time("OPC'ED MASK GENERATED");

        printf("Generating aerial image...\n");
        flt_image_t aerial(original.get_width(),original.get_height());
        aerial.convolve(opced_mask,kernel);  
        aerial.save_pgm(argv[3]); 
        t.print_time("AERIAL IMAGE SAVED");

        pbm_image_t contour(original.get_width(),original.get_height());
        contour.convolve(opced_mask,kernel);  
        printf("Final error: %d\n",contour.diff(original));
        printf("Final count: %d\n",opced_mask.count());
        contour.save(argv[4]);
        t.print_time("CONTOUR IMAGE SAVED");

    }
}    

