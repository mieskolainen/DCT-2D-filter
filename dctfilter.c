// -----------------------------------------------------------------
// Overlapping 2D-DCT (Discrete Cosine Transform) image filter
//
// Image model: z = y + n, where z is the measured image mxn
//                               y is the noise free image
//                               n is Gaussian noise
// 
// Compile with: g++ -O3 dctfilter.c dct.c -o dctfilter
// 
// Mikael Mieskolainen, {2011, 2021}
// ----------------------------------------------------------------- 

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <cstring>

#include "dct.h" // 2D-DCT Library

// This abstracts the image matrix as a linear array
// access image matrix pixel at row i and column j, with N columns
// The first for "->" and the second for "." access
#define pixel(data,i,j)   ( (data)->image[(i)*((data)->N) + (j)] )
#define pixel_h(data,i,j) ( (data).image[(i)*((data).N) + (j)] )

// This is the image object
struct IMAGE {

  IMAGE(int m, int n) {
    M = m;
    N = n;
  };

  double *image; // Data
  int M;         // Number of Rows
  int N;         // Number of Columns
};

// Subfunction prototypes
void ones(struct IMAGE*);
void zeros(struct IMAGE*);
void printValues(struct IMAGE*, int jump=1);
void pointDiv(struct IMAGE*, struct IMAGE*, struct IMAGE*);
void pointMul(struct IMAGE*, struct IMAGE*, struct IMAGE*);
bool dctfilter(struct IMAGE*, struct IMAGE*, int step, double thresh, int window_size);
double thrfunc(double* outblock, int window_size, double thresh);


// Main program
int main(int argc, char* argv[]) {

  char* inputfile;
  char* outputfile;

  // Image dimensions fixed here (no header in the binary files)
  int M = 512;              // Rows
  int N = 512;              // Cols
  
  // Filter parameters
  int step         = 1;     // Local window hopping step
  double threshold = 40;    // Frequency domain threshold
  int bytes        = 8;     // 64 bit (double) input expected (more formats TBD!)
  bool quiet       = false; // Print output

  int opt;
  while ((opt = getopt(argc, argv, "i:o:b:m:n:s:t:q")) != -1) {  
    switch(opt) {
      case 'i':
        inputfile = strdup(optarg);
        break;
      case 'o':
        outputfile = strdup(optarg);
        break;
      case 'b':
        bytes = atoi(optarg);
        if (bytes < 1) {
          printf("-b needs to be int >= 1 \n");
          return EXIT_FAILURE;
        }
        break;
      case 'm':
        M = atoi(optarg);
        if (M < 1) {
          printf("-m needs to be integer >= 1 \n");
          return EXIT_FAILURE;
        }
        break;
      case 'n':
        N = atoi(optarg);
        if (N < 1) {
          printf("-n needs to be integer >= 1 \n");
          return EXIT_FAILURE;
        }
        break;
      case 's':
        step = atoi(optarg);
        if (step < 1) {
          printf("-s needs to be integer >= 1 \n");
          return EXIT_FAILURE;
        }
        break;
      case 't':
        threshold = atof(optarg);
        if (threshold < 0) {
          printf("-t needs to be float >= 0 \n");
          return EXIT_FAILURE;
        }
        break;
      case 'q':
        quiet = true;
        break;
      case ':':
        printf("option needs a value \n");
        break;
      case '?':
        fprintf(stderr, "Usage: %s [-i inputfile] [-o outputfile] [-m nrows] [-n ncols] [-s step] [-t threshold] [-q quiet] \n", argv[0]);
        return EXIT_FAILURE;
    }
  }

  if (!quiet) {
    printf("-i inputfile:  %s \n", inputfile);  
    printf("-o outputfile: %s \n", outputfile);
    printf("-b bytes:      %d \n", bytes);
    printf("-m nrows:      %i \n", M);
    printf("-n ncols:      %i \n", N);
    printf("-s step:       %i \n", step); 
    printf("-t threshold:  %0.1f \n", threshold);
  }

  const int window_size = 8;  // Keep it at 8 (2D-DCT function is fixed 8x8)

  // ------------------------------------------------
  // Create image matrix for the original (noisy) image
  // and denoised image

  struct IMAGE z(M,N);
  struct IMAGE y_est(M,N);

  // ------------------------------------------------
  // Allocate memory

  z.image = (double*) malloc(z.M * z.N * sizeof(double));
  y_est.image = (double*) malloc(y_est.M * y_est.N * sizeof(double));
  
  if (z.image == NULL) {
    printf("dctfilter: Out of memory with z \n");
    return EXIT_FAILURE;
  }
  if (y_est.image == NULL) {
    printf("dctfilter: Out of memory with y_est \n");
    return EXIT_FAILURE;
  }

  // Init matrices with zeros
  zeros(&z);
  zeros(&y_est);
  
  // ------------------------------------------------
  // Read the image
  
  FILE *fid;
  fid = fopen(inputfile, "rb");

  if (fid) {
    int ret = fread(z.image, bytes, z.M * z.N, fid);
    fclose(fid);
  } else {
    printf("dctfilter: Error opening inputfile <%s> ! \n", inputfile);
    return EXIT_FAILURE;
  }
  
  // ------------------------------------------------
  // Filter the image

  if (dctfilter(&z, &y_est, step, threshold, window_size)) {
    // fine
  } else {
    printf("dctfilter: Problem in filtering \n");
    return EXIT_FAILURE;
  }
  
  // ------------------------------------------------
  // Write the image
  fid = fopen(outputfile, "wb");
  if (fid) { 
    int ret = fwrite(y_est.image, bytes, y_est.M * y_est.N, fid);
    fclose(fid);
  } else {
    printf("dctfilter: Error opening outputfile <%s> ! \n", outputfile);
    return EXIT_FAILURE;
  }

  // Free the image matrices
  free(z.image);
  free(y_est.image);
  
  if (!quiet) { printf("dctfilter: Filtering done \n"); }

  return EXIT_SUCCESS;
}


// 2D-filter function
bool dctfilter(struct IMAGE* z, struct IMAGE* y_est, int step, double thresh, int window_size) {
  
  // Boundary around image matrix
  const int b = window_size / 2;

  // 2D-DCT 8x8 buffers (matrix row by row as linear array)
  double inblock[(window_size*window_size)]  = {0};
  double outblock[(window_size*window_size)] = {0};

  // Allocate the aggregation matrix and weight matrix
  struct IMAGE A(z->M, z->N);
  struct IMAGE W(z->M, z->N);

  A.image = (double*) malloc(A.M * A.N * sizeof(double));
  W.image = (double*) malloc(W.M * W.N * sizeof(double));
  
  if (A.image == NULL) { printf("dctfilter: Out of memory with A \n"); return false; }
  if (W.image == NULL) { printf("dctfilter: Out of memory with W \n"); return false; }
  
  // Init matrices with zeros
  zeros(&A);
  zeros(&W);
  
  // Over all rows and columns of the image
  for (int i = b; i < z->M - b + 2; i += step) {
    for (int j = b; j < z->N - b + 2; j += step) {

      double* VP = &inblock[0]; // Current pointer

      // Accumulate the local window
      for (int m = i-b, k = 0; m < i+b; ++m) {
        for (int n = j-b; n < j+b; ++n, ++k) {

          // Check boundary
          *(VP + k) = (m < z->M && n < z->N) ? pixel(z,m,n) : 0.0;
        }
      }

      // 2D-DCT for the block
      fdct(VP, outblock);
      
      // Filtering in frequency domain & weight calculation
      const double weight = thrfunc(outblock, window_size, thresh);
      
      // Inverse 2D-DCT for the block
      idct(outblock, VP);

      // Aggregate the block and weight
      for (int m = i-b, k = 0; m < i+b; ++m) {
        for (int n = j-b; n < j+b; ++n, ++k) {

          if (m < z->M && n < z->N) { // Check boundary
            pixel_h(A,m,n) += weight * (*(VP + k));
            pixel_h(W,m,n) += weight;
          }
        }
      }
    }
  }
  
  // Weighted mean of the aggregation buffer
  for (int i = 0; i < y_est->M; ++i) {
    for (int j = 0; j < y_est->N; ++j) {
      pixel(y_est,i,j) = pixel_h(A,i,j) / pixel_h(W,i,j);
    }  
  }
  
  free(A.image);
  free(W.image);

  return true;
}

// Hard-thresholding and block weight function
double thrfunc(double* outblock, int window_size, double thresh) {

  double sum  = 0.0;
  double sum2 = 0.0;
  for (int i = 0; i < (window_size*window_size); ++i, ++outblock) {

    const double c = fabs(*outblock);

    if (c < thresh) {
      *outblock = 0;
    } else {
      sum  += c;   // Accumulate coefficient values
      sum2 += c*c;
    }
  }
  const double n   = window_size*window_size;
  const double var = (sum2 - sum*sum/n) / n;
  
  // Weight for this block
  const double W = var > 0 ? 1.0/var : 1.0;

  return W;
}

// Point wise division: y = a ./ b
void pointDiv(struct IMAGE *y, struct IMAGE *a, struct IMAGE *b) {
  
  for (int i = 0; i < y->M; ++i) {
    for (int j = 0; j < y->N; ++j) {
      pixel(y,i,j) = pixel(a,i,j) / pixel(b,i,j);
    }  
  }   
}

// Point wise multiplication: y = a .* b
void pointMul(struct IMAGE *y, struct IMAGE *a, struct IMAGE *b) {
  
  for (int i = 0; i < y->M; ++i) {
    for (int j = 0; j < y->N; ++j) {
      pixel(y,i,j) = pixel(a,i,j) * pixel(b,i,j);
    }
  }
}

// Init with ones
void ones(struct IMAGE *x) {
  
  for (int i = 0; i < x->M; ++i) {
    for (int j = 0; j < x->N; ++j) {
      pixel(x,i,j) = 1.0;
    }
  }
}

// Init with zeros
void zeros(struct IMAGE *x) {
  
  for (int i = 0; i < x->M; ++i) {
    for (int j = 0; j < x->N; ++j) {
      pixel(x,i,j) = 0.0;
    }
  }
}

void printValues(struct IMAGE *x, int jump) {
  
  for (int i = 0; i < x->M; i += jump) {
    for (int j = 0; j < x->N; j += jump) {
      printf("Value at [%d,%d] is %0.2f \n", i, j, pixel(x,i,j) );
    }
  }
}
