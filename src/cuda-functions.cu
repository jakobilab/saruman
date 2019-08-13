/*
 * cuda.cu
 *
 * CUDA Module
 *
 *  Created on: Mar 10, 2009
 *      Author: Tobias Jakobi
 */

// system includes

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>
#include <argtable2.h>
#include <pthread.h>
#include <math.h>

// cuda includes
#include <cuda.h>
//#include <cutil.h>

#include <helper_cuda.h>
#include <helper_gl.h>
//#include <helper_cuda_gl.h>
#include <helper_cuda_drvapi.h>
#include <helper_functions.h>
#include <helper_image.h>
#include <helper_math.h>
#include <helper_string.h>
#include <helper_timer.h>


//#include <cuda_runtime_api.h>
//#include <cuda_runtime.h>

// project includes


// not be be included - definition glitches
#include "util.h"
#include "saruman.h"
#include "hash.h"

extern "C"
{
  //uint3 const threadIdx;
  //uint3 const blockIdx;
  dim3 const blockDim;
  dim3 const gridDim;
  //int const warpSize;
}

/*!Swaps the chars  \a x and \a y. */
#define SWAP_CHAR( x, y ) {char c; c = x; x = y; y = c;}

#define NR_OF_BLOCKS 8

float alignment_time[2];

// internal defintion (differs from definition in header file)




/*
 * Process alignment in CUDA
 */



// ---------- CUDA KERNEL FUNCTIONS - EXECUTED ON GPU ONLY ---------------------

__device__ void
reverse_sequence (char* sequence, int sequence_length)
{
  int i, j;
  sequence_length = sequence_length - 1;

  for (i = 0, j = sequence_length; i < j; i++, j--)
    {
      SWAP_CHAR (sequence[i], sequence[j]);
    }

}

__device__ short int
minimize (short int a, short int b, short int c)
{

  short int min = a;

  if (min > b)
    {
      min = b;
        }
  if (min > c)
        {
      min = c;
        }
  return min;

    }


__device__ void
get_alignment (int match, int mismatch,
                  char *read_sequences, char* genome_sequences,
                  char *read_alignments, char *genome_alignments,
                  short int* edit_matrices_on_device_gpu,
                  int max_position, int edit_matrix_start,
                  int read_sequence_length, int read_start,
                  int genome_start, int genome_sequence_length,
                  int error, int end_offset, int start_offset, int threadID
                  )
{

  // set needed values
  int position = max_position + edit_matrix_start;
  int alignment_position = genome_start+error*threadID;
  int mismatch_correction = 0;

   // fprintf (stderr,"max length %d\n", genome_sequence_length+error);


  // generate alignment only for genome_sequence length
  while (position >= read_sequence_length+1)
    {

      short int diagonal;
      short int left;
      short int value = edit_matrices_on_device_gpu[position];
      short int top;

      // make sure to not exceed the borders of the matrix (i.e. set to high scores)
      if (position % (read_sequence_length + 1) != 0)
        {
          diagonal = edit_matrices_on_device_gpu[position - read_sequence_length - 2];
          left = edit_matrices_on_device_gpu[position - 1];
        }
      else
        {
          diagonal = 999;
          left = 999;
        }


      if (position != 0)
        {
          top = edit_matrices_on_device_gpu[position - read_sequence_length - 1];
        }
      else
        {
          top = 999;
        }

      // converting from 1d to 2d matrix and interpolating heigth and width

      short int width = 0;
      short int width_tmp = ((position % (read_sequence_length + 1)) - 1);
      if (width_tmp > 0) // make sure we don't exceed the matrix borders
        {
          width = width_tmp;
          width_tmp++;
        }
      else
        {
          width_tmp = 0;
        }

      short int heigth_tmp = ((position - edit_matrix_start - (width)) / (read_sequence_length + 1)) - 1;
      short int heigth = 0;
      if (heigth_tmp > 0) // make sure we don't exceed the matrix borders
        {
          heigth = heigth_tmp;
        }

      // if we have end gaps we don't want to score them
            if ((width == (read_sequence_length -1) && (heigth >= (genome_sequence_length-2*end_offset))) // means we have cost free gaps at the end of the matrix
                 || ((start_offset)&&(width == 0 && (heigth <= ((start_offset+1)*2)))) // cost free gaps at upper left matrix area
                    )
        {
          mismatch_correction = mismatch;
          //fprintf (stderr,"mismatch correction on\n");
        }else{
          mismatch_correction = 0;
          //fprintf (stderr,"mismatch correction off w%d  h%d \n",width,heigth);
          }
/*

      fprintf (stderr," === %c (%d) | %c (%d) [%d] <%d>\n", read_sequences[width + read_start],
             width,
              genome_sequences[heigth + genome_start],
              heigth,
              alignment_position-genome_start-error*threadID,end_offset);
     fprintf (stderr,"AP: %d\n",alignment_position-genome_start-error*(threadID));
*/


      // match and mismatch are prefered over indels here


/*
            fprintf (stderr,"my enviroment:\n", value);
            fprintf (stderr,"%d %d\n", diagonal,top);
            fprintf (stderr,"%d %d\n", left,value);
            fprintf (stderr,"costs: %d\n\n", mismatch-mismatch_correction);
*/

            //fprintf (stderr,"pos: %d\n", alignment_position);

      // mismatch
if (diagonal + mismatch == value)
        {
          read_alignments[alignment_position] = read_sequences[width + read_start];
          genome_alignments[alignment_position] = genome_sequences[heigth + genome_start];
          // fprintf (stderr,"=> mismatch\n");
          position = position - read_sequence_length - 2;
          //fprintf (stderr," === %c | %c \n", read_alignments[alignment_position],genome_alignments[alignment_position]);
          alignment_position++;
          continue;
        }
      //match
      else if (diagonal + match == value && (genome_sequences[genome_start + heigth] == read_sequences[read_start + width]))
        {
          read_alignments[alignment_position] = read_sequences[width + read_start];
          genome_alignments[alignment_position] = genome_sequences[heigth + genome_start];
          position = position - read_sequence_length - 2;
          //fprintf (stderr," === %c | %c \n", read_alignments[alignment_position],genome_alignments[alignment_position]);
          alignment_position++;
          //fprintf (stderr,"=> match %d\n",alignment_position);

          continue;

        }
      // deletion
      else if (top + mismatch - mismatch_correction == value)
        {
          read_alignments[alignment_position] = '_';
          genome_alignments[alignment_position] = genome_sequences[heigth + genome_start];
          position = position - read_sequence_length-1;
          //                   fprintf (stderr," === %c | %c \n", read_alignments[alignment_position],genome_alignments[alignment_position]);
          alignment_position++;
          //fprintf (stderr,"=> deletion\n");
          continue;

        }
      // insertion
      else if (left + mismatch == value)
        {
          read_alignments[alignment_position] = read_sequences[width + read_start];
          genome_alignments[alignment_position] = '_';
          //fprintf (stderr,"=> insertion\n");
          position = position - 1;
          //         fprintf (stderr," === %c | %c \n", read_alignments[alignment_position],genome_alignments[alignment_position]);
          alignment_position++;

          continue;
        }
      else
        // stop here
        {
/*
          fprintf (stderr,"=> end (%d)\n",mismatch_correction);
          //break;
          read_alignments[alignment_position] = '_';
          genome_alignments[alignment_position] = genome_sequences[heigth + genome_start];
          position = position - read_sequence_length-1;
          //                   fprintf (stderr," === %c | %c \n", read_alignments[alignment_position],genome_alignments[alignment_position]);
          alignment_position++;
          //fprintf (stderr,"=> deletion\n");
*/
          break;
        }



  //assert(alignment_position < max_size);

    }

/*
  fprintf (stderr,"alignment position %d\n", alignment_position);
  fprintf (stderr,"sequence start %d\n", genome_start);
  fprintf (stderr,"alignment start %d\n", genome_start+error);
  fprintf (stderr,"length %d\n", alignment_position-genome_start-error*(threadID));
  fprintf (stderr,"printing %d chars\n",alignment_position-genome_start+error);

  fprintf (stderr,"genome\t%.*s\n", alignment_position-genome_start-error*(threadID), &genome_alignments[genome_start+error*threadID]);
  fprintf (stderr,"read\t%.*s\n\n", alignment_position-genome_start-error*(threadID), &read_alignments[genome_start+error*threadID]);
*/

/*
    fprintf (stderr,"would print  %d\n", genome_start+error*threadID);

  fprintf (stderr,"genome\t%c\n", genome_alignments[genome_start+error*threadID]);
  fprintf (stderr,"read\t%c\n\n", read_alignments[genome_start+error*threadID]);
*/


  // reverse strings
  //fprintf (stderr,"AP\t%d\n",alignment_position);

  reverse_sequence (&read_alignments[genome_start+error*threadID], alignment_position-genome_start-error*(threadID));
  reverse_sequence (&genome_alignments[genome_start+error*threadID], alignment_position-genome_start-error*(threadID));

/*
  fprintf (stderr,"genome\t%.*s\n", alignment_position-genome_start-error*(threadID)+1, &genome_alignments[genome_start+error*threadID]);
  fprintf (stderr,"read\t%.*s\n\n", alignment_position-genome_start-error*(threadID)+1, &read_alignments[genome_start+error*threadID]);

  fprintf (stderr,"====== (%d/%d) ========\n",start_offset,end_offset);
*/
}

__global__ void
compute_score (int match, int mismatch,
               char *read_sequences, char* genome_sequences,
               int *read_offsets, int* genome_offsets,
               char* read_alignments, char* genome_alignments,
               int sample_size,
               short int* edit_matrices_on_device_gpu,
               int* edit_offsets,
               int* device_score_array,
               char* direction_array,
               int error
               )
{


  // get data for this thread

  int threadID = blockIdx.x * blockDim.x + threadIdx.x;

  // check if this thread has to be executed
  if (threadID < sample_size)
    {
          // get needed data
      short unsigned int genome_sequence_length = 0;
      short unsigned int read_sequence_length = 0;
      int read_start = -1;
      int genome_start = -1;
      int edit_matrix_start = 0;

      if (threadID != 0)
        {
          // get length

          read_sequence_length = read_offsets[threadID] - read_offsets[threadID - 1];
          genome_sequence_length = genome_offsets[threadID] - genome_offsets[threadID - 1];

          // gather start positions
          read_start = read_offsets[threadID - 1];
          genome_start = genome_offsets[threadID - 1];
          edit_matrix_start = edit_offsets[threadID];
        }
      else
        {
          // gather length, start position is 0
          read_sequence_length = read_offsets[threadID];
          genome_sequence_length = genome_offsets[threadID];
          edit_matrix_start = edit_offsets[threadID];
        }



      int max_score_position = 0;

      // initalize the memory

      if (threadID == 0)
        {
          read_start++;
          genome_start++;
        }



      short int start_offset = 0; // todo: is this processing needed?
      short int end_offset = 0;

      switch (direction_array[threadID])
        {
        case 1:
          end_offset = error;
          break;

        case 2:
          start_offset = error;
          break;

        case 3:
          start_offset = error;
          end_offset = error;
          break;

        case 5:
          end_offset = error;
          break;

        case 6:
          start_offset = error;
          break;

        case 7:
          start_offset = error;
          end_offset = error;
          break;
    }
                  //fprintf (stderr,"%d %d %d\n", end_offset,read_sequence_length,genome_sequence_length);

      for (int init_counter_1 = 0; init_counter_1 <= read_sequence_length; ++init_counter_1)
        {
          edit_matrices_on_device_gpu [edit_matrix_start + init_counter_1] = init_counter_1;
    }


      for (int init_counter_2 = 1; init_counter_2 <= genome_sequence_length; ++init_counter_2)
        {
          if (init_counter_2 < 2 * start_offset + 1)
            {
              edit_matrices_on_device_gpu [edit_matrix_start + (init_counter_2) * (read_sequence_length + 1)] =
                      edit_matrices_on_device_gpu [edit_matrix_start + (init_counter_2 - 1) * (read_sequence_length + 1)];
            }
          else
            {
              edit_matrices_on_device_gpu [edit_matrix_start + (init_counter_2) * (read_sequence_length + 1)] =
              edit_matrices_on_device_gpu [edit_matrix_start + (init_counter_2 - 1) * (read_sequence_length + 1)] + 1;
            }
/*
          else
            {
              //printf("init: %d\n",edit_matrix_start + (counter) * (read_sequence_length + 1));
              edit_matrices_on_device_gpu [edit_matrix_start + (counter) * (read_sequence_length + 1)] =
                      edit_matrices_on_device_gpu [edit_matrix_start + (counter - 1) * (read_sequence_length + 1)];
        }
*/
        }


      short int mismatch_correction = 0;

      for (short int heigth = 0; heigth < genome_sequence_length; ++heigth)
        {
          for (short int width = 0; width < read_sequence_length; ++width)
            {

            if ((width == (read_sequence_length -1) && (heigth >= (genome_sequence_length-2*end_offset))))
              {
                //fprintf (stderr,"mismatch costs off for %d|%d\n",heigth,width);
                mismatch_correction = mismatch;
            }else{
                 mismatch_correction = 0;
}


              short int from_top = 0;
              short int from_diagonal = 0;
              short int from_left = 0;

              // mismatch

             //fprintf (stderr,"%d:%d => %c:%c => ",heigth,width,genome_sequences[genome_start + heigth],read_sequences[read_start + width]);

              if (genome_sequences[genome_start + heigth] != read_sequences[read_start + width])
                {
                  from_diagonal = edit_matrices_on_device_gpu[edit_matrix_start + (heigth) *
                          (read_sequence_length + 1) + width] + mismatch;

              // match
                }
              else
                {
                  from_diagonal = edit_matrices_on_device_gpu[edit_matrix_start + (heigth) *
                          (read_sequence_length + 1) + width] + match;
                }

              // indels

              from_left = edit_matrices_on_device_gpu[edit_matrix_start + (heigth + 1) *
                      (read_sequence_length + 1) + width]
                      + mismatch;

              from_top = edit_matrices_on_device_gpu[edit_matrix_start + (heigth) *
                      (read_sequence_length + 1) + width + 1]
                      + mismatch - mismatch_correction;

              // calculate best score and save it in the matrix
              edit_matrices_on_device_gpu [edit_matrix_start + (heigth + 1) * (read_sequence_length + 1) + width + 1] =
              minimize (from_diagonal, from_left, from_top);
              // fprintf (stderr,"-> %d %d %d (%d)\n",from_diagonal,from_left,from_top,minimize (from_diagonal, from_left, from_top, this_algorithm_become_skynet_cost));
            }
        }
      max_score_position = (read_sequence_length + 1)*(genome_sequence_length + 1) - 1;

      device_score_array[threadID] = edit_matrices_on_device_gpu[edit_matrix_start + max_score_position];






      if (device_score_array[threadID]<= error){
          if (threadID == 0){
              error =0;
            }

/*
      for (short int heigth = 0; heigth <= genome_sequence_length; ++heigth)
        {

          for (short int width = 0; width <= read_sequence_length; ++width)
            {
              fprintf (stderr,"%*d ",3, (edit_matrices_on_device_gpu[edit_matrix_start + (heigth) * (read_sequence_length + 1) + width] < 900 ? edit_matrices_on_device_gpu[edit_matrix_start + (heigth) * (read_sequence_length + 1) + width] : 0));
            }
          fprintf (stderr,"\n");
        }

       fprintf (stderr,"\n\n");
*/
      get_alignment (match, mismatch,
                        read_sequences, genome_sequences,
                        read_alignments, genome_alignments,
                        edit_matrices_on_device_gpu,
                        max_score_position,
                        edit_matrix_start,
                        read_sequence_length,
                        read_start,
                        genome_start,
                        genome_sequence_length,
                        error,end_offset, start_offset, threadID
                        );
        }
    }
}


// ---------- HOST (CPU) CODE --------------------------------------------------

extern "C"
void
cuda_allocate_memory(int do_reverse) {

    //cudaError_t error_code; // holds the current error code

    cudaSetDevice(use_gpu_value);

    if (do_reverse==4) {

        checkCudaErrors(cudaMalloc((void**) & main_alloc_ptr_rev, cuda_memory_per_thread));
        //fprintf(stderr, "-> allocating %lu kb of GPU memory for reverse thread: %s\n", cuda_memory_per_thread / 1024 , cudaGetErrorString(error_code));

    } else {

        checkCudaErrors(cudaMalloc((void**) & main_alloc_ptr_for, cuda_memory_per_thread));
        //fprintf(stderr, "-> allocating %lu kb of GPU memory for forward thread: %s\n", cuda_memory_per_thread / 1024 , cudaGetErrorString(error_code));

    }


}

extern "C"
void
cuda_free_memory(int do_reverse) {

    if (do_reverse) {
    checkCudaErrors(cudaFree(main_alloc_ptr_rev));
    } else {
    checkCudaErrors(cudaFree(main_alloc_ptr_for));
    }


}
unsigned int
pad_memory(unsigned int unpadded) {


    return (unsigned int) (ceil((float)unpadded/512)*512);

}

extern "C"
void
cuda_needleman_wunsch (int match, int mismatch,
                     char* genome_sequences_on_host, char* read_sequences_on_host,
                     int * genome_offsets, int * read_offsets, int * edit_offsets,
                     int sample_size, int* read_id_array, int* genome_pos_array,
                     char* direction_array, int error, size_t total_genome_sequence_array_size,
                     size_t total_read_sequence_array_size, size_t total_edit_matrix_size,
                     int do_reverse
                     )
{


  // create a timer
  unsigned int timer = 0;
  //CUT_SAFE_CALL (cutCreateTimer (&timer));
  //CUT_SAFE_CALL (cutStartTimer (timer));

  fprintf (stderr,"\t\t-> starting CUDA kernel with %d samples from %s thread\n",sample_size,(do_reverse ? "reverse" : "forward"));

    if (do_reverse) {
        reverse_alignment_candidates += sample_size;
    } else {
        forward_alignment_candidates += sample_size;
    }

  // declare and init device variables

  char* main_alloc_ptr = NULL;

  // aligments on host
  char* read_alignments_on_host = NULL;
  char* genome_alignments_on_host = NULL;


//  char* edit_matrix_on_host = NULL; // FIXME  needed? only 0-filled matrix?

  // arrays holding the alignment scores
  //int* alignment_scores_on_device = NULL;
  int* alignment_scores_on_host = (int*) malloc (sample_size * sizeof (int));

  // allocate space needed for alignments on host
  read_alignments_on_host = (char *) calloc ((total_genome_sequence_array_size), sizeof (char));
  if (read_alignments_on_host == NULL)
    {
      fprintf (stderr, "ERROR: could not allocate memory for read alignments on host!\n");
      return;
    }

  genome_alignments_on_host = (char *) calloc ((total_genome_sequence_array_size), sizeof (char));
  if (genome_alignments_on_host == NULL)
    {
      fprintf (stderr, "ERROR: could not allocate memory for genome alignments on host\n");
      return;
    }



  //cudaError_t error_code; // holds the current error code

  // allocating device memory

  // summing up all needed memory and allocating it as one large chunk in memory

  size_t device_memory_needed = sample_size * sizeof (int) + // edit matrix offsets
                                sample_size * sizeof (int) + // genome offsets
                                sample_size * sizeof (int) + // read offsets
                                sample_size * sizeof (char) + // direction array;
                                total_edit_matrix_size * sizeof (short int) + // edit matrices
                                (total_genome_sequence_array_size-sample_size*error) * sizeof (char) + // genome sequences
                                total_read_sequence_array_size * sizeof (char) + // read sequences
                                sample_size * sizeof (int) + // score array
                                (total_genome_sequence_array_size)* sizeof (char) + // genome alignment
                                (total_genome_sequence_array_size)* sizeof (char);// read alignment


  if (be_verbose_value)
    {
  fprintf(stderr,"operation requires %d kb GPU memory\n",(int)device_memory_needed/1024);


    }


  // getting correct memory address

    if (do_reverse)
    {
        main_alloc_ptr = main_alloc_ptr_rev;
    }
  else
    {
        main_alloc_ptr = main_alloc_ptr_for;
    }


    checkCudaErrors(cudaMemset(main_alloc_ptr, 0, cuda_memory_per_thread));
    if (be_verbose_value) {
        //fprintf(stderr, "Memsetting device memory with 0: %s\n", cudaGetErrorString(error_code));
    }

    // some pointer arithmetic fun


  unsigned int genome_offset_ptr    = 0; // allocation starts here
  unsigned int read_offset_ptr      = sample_size * sizeof (int);
  unsigned int edit_offset_ptr      = sample_size * sizeof (int) * 2;
  unsigned int score_ptr            = sample_size * sizeof (int) * 3;
  unsigned int edit_matrix_ptr      = sample_size * sizeof (int) * 4;



  unsigned int direction_ptr        = edit_matrix_ptr + total_edit_matrix_size * sizeof(short int);
  //printf("pointer of edit matrix:  = %d\n",direction_ptr - edit_matrix_ptr);

  unsigned int read_seq_ptr         = direction_ptr + sample_size * sizeof(char);

  unsigned int gen_seq_ptr          = read_seq_ptr + total_read_sequence_array_size * sizeof (char);
  unsigned int read_aln_ptr         = gen_seq_ptr + (total_genome_sequence_array_size-sample_size*error) * sizeof (char);
  unsigned int gen_aln_ptr          = read_aln_ptr + (total_genome_sequence_array_size)* sizeof (char);


  int*  genome_offsets_on_device = (int*) &main_alloc_ptr[genome_offset_ptr];
  int*  read_offsets_on_device = (int*)  &main_alloc_ptr[read_offset_ptr];
  int*  edit_offsets_on_device =  (int*) &main_alloc_ptr[edit_offset_ptr];
  char* direction_array_on_device = &main_alloc_ptr[direction_ptr];

  short int* edit_matrices_on_device =  (short int*) &main_alloc_ptr[edit_matrix_ptr];

  char* read_sequences_on_device = &main_alloc_ptr[read_seq_ptr];
  char* genome_sequences_on_device = &main_alloc_ptr[gen_seq_ptr];

  char* read_alignments_on_device = &main_alloc_ptr[read_aln_ptr];
  char* genome_alignments_on_device = &main_alloc_ptr[gen_aln_ptr];

   int* alignment_scores_on_device =(int*) &main_alloc_ptr[score_ptr];



  // setup kernel parameters

  int block_size = NR_OF_BLOCKS; // change depending on used card and available multiprocessors

  // number of blocks actual used to get all alignments in parallel
  int blocks_used = sample_size / block_size + (sample_size % block_size == 0 ? 0 : 1);


  // now copying data from host to device

    //genome sequences from host to device
  checkCudaErrors(cudaMemcpy (genome_sequences_on_device, genome_sequences_on_host, (total_genome_sequence_array_size-sample_size*error) * sizeof (char), cudaMemcpyHostToDevice));
  if (be_verbose_value)
    {
      //fprintf (stderr,"Copying genome sequences (host->device): %s\n", cudaGetErrorString (error_code));
    }


  //read sequences from host to device
  checkCudaErrors(cudaMemcpy (read_sequences_on_device, read_sequences_on_host, total_read_sequence_array_size * sizeof (char), cudaMemcpyHostToDevice));
  if (be_verbose_value)
    {
      //fprintf (stderr,"Copying read sequences (host->device): %s\n", cudaGetErrorString (error_code));
    }


  //read sequence offsets from host to device
  checkCudaErrors(cudaMemcpy (read_offsets_on_device, read_offsets, sample_size * sizeof (int), cudaMemcpyHostToDevice));
  if (be_verbose_value)
    {
      //fprintf (stderr,"Copying read sequence offsets (host->device): %s\n", cudaGetErrorString (error_code));
    }

   //genome sequence offsets from host to device
  checkCudaErrors(cudaMemcpy (genome_offsets_on_device, genome_offsets, sample_size * sizeof (int), cudaMemcpyHostToDevice));
  if (be_verbose_value)
    {
      //fprintf (stderr,"Copying genome sequence offsets (host->device): %s\n", cudaGetErrorString (error_code));
    }

  //edit matrix offsets from host to device
  checkCudaErrors(cudaMemcpy (edit_offsets_on_device, edit_offsets, sample_size * sizeof (int), cudaMemcpyHostToDevice));
  if (be_verbose_value)
    {
      //fprintf (stderr,"Copying edit matrix offsets (host->device): %s\n", cudaGetErrorString (error_code));
    }

  //edit matrix offsets from host to device
  checkCudaErrors(cudaMemcpy (direction_array_on_device, direction_array, sample_size * sizeof (char), cudaMemcpyHostToDevice));
  if (be_verbose_value)
    {
      //fprintf (stderr,"Copying edit matrix offsets (host->device): %s\n", cudaGetErrorString (error_code));
    }

  if (be_verbose_value)
    {
      fprintf (stderr,"Launching CUDA kernel\n");
    }

  // to boldy go where no kernel has gone before...
  compute_score <<< blocks_used, block_size >>> (match, mismatch,
                                                 read_sequences_on_device,
                                                 genome_sequences_on_device,
                                                 read_offsets_on_device,
                                                 genome_offsets_on_device,
                                                 read_alignments_on_device,
                                                 genome_alignments_on_device,
                                                 sample_size,
                                                 edit_matrices_on_device,
                                                 edit_offsets_on_device,
                                                 alignment_scores_on_device,
                                                 direction_array_on_device,
                                                 error
                                                 );

  // stop timer, print time
  //CUT_SAFE_CALL (cutStopTimer (timer));

  // check if kernel execution generated an error
  //CUT_CHECK_ERROR ("Kernel execution failed: ");

  // copy result from device to host

  //genome sequence alignments from host to device
  cudaMemcpy (genome_alignments_on_host, genome_alignments_on_device, total_genome_sequence_array_size * sizeof (char), cudaMemcpyDeviceToHost);
  if (be_verbose_value)
    {
      //fprintf (stderr,"Copying genome sequence alignments (device->host): %s\n", cudaGetErrorString (error_code));
    }

  //read sequence alignments from host to device
  cudaMemcpy (read_alignments_on_host, read_alignments_on_device,total_genome_sequence_array_size * sizeof (char), cudaMemcpyDeviceToHost);
  if (be_verbose_value)
    {
      //fprintf (stderr,"Copying read sequence alignments (device->host): %s\n", cudaGetErrorString (error_code));
    }

  cudaMemcpy (alignment_scores_on_host, alignment_scores_on_device, sample_size * sizeof (int), cudaMemcpyDeviceToHost);
  if (be_verbose_value)
    {
      //fprintf (stderr,"Copying array of alignment scores (device->host): %s\n", cudaGetErrorString (error_code));
    }


  // printing some garbage on the screen
  for (int var = 0; var < sample_size; ++var)
    {

          short unsigned int tmp_direction = 0;

          if (direction_array[var]> 3){
              tmp_direction=1;
            }

      if (alignment_scores_on_host[var] <= error) // hits with an error under our threshold are always printed out
        {

          short unsigned int genome_sequence_length = 0;
          int genome_start = -1;
          int error_buffer = 0;

          if (var > 0)
            {
              // get length
              genome_sequence_length = genome_offsets[var] - genome_offsets[var - 1];
              // gather start positions
              genome_start = genome_offsets[var - 1];
              //fprintf (stderr,"genomestart: %d\n",genome_start);
              error_buffer=error*(var);


            }
          else
            {
              // gather length, start position is 0
              genome_sequence_length = genome_offsets[var];
              genome_start = 0;
              //error_tmp = error;
          }


          // compute start / stop positions

          unsigned int start = genome_pos_array[var];
          unsigned int stop  = genome_pos_array[var] + genome_sequence_length - 1;

          unsigned short int deletion = 0;
          unsigned short int insertion = 0;

          unsigned short int start_correction = 0;
          unsigned short int stop_correction = 0;

          // correct positions due to opening and closing gaps

          // opening gaps
          int i = 0;
          while (read_alignments_on_host[genome_start+error_buffer + i] == '_')
            {
              start++;
              start_correction++;

             // fprintf (stderr,"start++ (%c)\n",read_alignments_on_host[genome_start+error_buffer + i]);
              i++;
            }

          // closing gaps
          i = 0;





        // printf("last char is (%c)\n",read_alignments_on_host[genome_start + error_buffer + genome_sequence_length + i + error-1]);

           while (read_alignments_on_host[genome_start + error_buffer + genome_sequence_length + i + error-1] == '_'
                 || read_alignments_on_host[genome_start + error_buffer + genome_sequence_length + i + error-1] == 0)
            {
              if (read_alignments_on_host[genome_start + error_buffer + genome_sequence_length + i + error-1] == '_')
                {
                  stop--;
                  stop_correction++;
                  //fprintf (stderr,"stop-- (%c)\n", read_alignments_on_host[genome_start + error_buffer + genome_sequence_length + i + error-1]);
                }
              i--;
            };


          i = 0;
          int z = 0;
          short int end_of_read=0;

           //fprintf (stderr,"would scan %d, doing only %d\n",  genome_sequence_length,genome_sequence_length - stop_correction - start_correction);

          int alignment_error_detection = 0; // looks for alignments which are inconsistent to their edit distance
                                             // happens e.g. for repeat regions where alignments are shifted by X bases and therefore to many indels are
                                             // added
                                             // fixed with version 1.0.3

            while (i < genome_sequence_length - stop_correction - start_correction + z) {
               // fprintf (stderr,"scanning %d (%c)\n", i,read_alignments_on_host[genome_start + error_buffer + start_correction + i]);
                if (read_alignments_on_host[genome_start + error_buffer + start_correction + i] == 0) {

                    // stop here, end reached, kill this loop
                    end_of_read = i;
                    break;

                }

                if (read_alignments_on_host[genome_start + error_buffer + start_correction + i] == '_') {
                    deletion++;
                    //fprintf(stderr, "deletion @ %d\n", i);
                    alignment_error_detection++;
                    z++;
                    i++;
                    continue;

                }
                if (genome_alignments_on_host[genome_start + error_buffer + start_correction + i] == '_') {
                    insertion++;
                  // fprintf(stderr, "insertion @ %d\n", i);
                    alignment_error_detection++;
                    z++;
                    i++;
                    continue;

                }
                if (read_alignments_on_host[genome_start + error_buffer + start_correction + i] != genome_alignments_on_host[genome_start + error_buffer + start_correction + i]) {
                   // fprintf(stderr, "mismatch @ %d\n", i);
                    alignment_error_detection++;
                     i++;
                    continue;


                }
                i++;

            }


/*
          short int indel_correction = 0;

          if (insertion>deletion){
              indel_correction=insertion-deletion;
            } else {
              indel_correction=insertion;
              }

*/
          short int chars_to_print = 0;

          if (end_of_read > 0){

              chars_to_print = end_of_read-stop_correction;

          } else {
             chars_to_print = read_length_value+deletion;
          }

        // printf("start correction: %d | stop correction: %d | genome: %d | read %d | indel correction %d | flag %d | length %d\n",start_correction, stop_correction, insertion,deletion, chars_to_print, direction_array[var],genome_sequence_length );

          struct read_name_hash *hash_entry = NULL;
          HASH_FIND_INT (read_names, &read_id_array[var], hash_entry);

          printf ("%s\t%d\t%d\t%s\t%.*s\t%.*s\t%d\n",
                  hash_entry->name,
                  start,
                  chars_to_print-insertion+start-1,
                  (tmp_direction ? "<<" : ">>"),
                  chars_to_print, &read_alignments_on_host[genome_start +error_buffer+ start_correction],
                  chars_to_print, &genome_alignments_on_host[genome_start+error_buffer + start_correction],
                  alignment_scores_on_host[var]);


          if (tmp_direction == 0){
              forward_output_control[read_id_array[var]] = 1;
              forward_alignments_done++;
          }else{
              reverse_output_control[read_id_array[var]] = 1;
              reverse_alignments_done++;
          }

        }
    }
 //fprintf (stderr,"\t\t-> %s CUDA kernel finished after %4.2f ms\n", (do_reverse ? "reverse" : "forward"), cutGetTimerValue (timer));

  if (do_reverse)
    {
      //alignment_time[1] += cutGetTimerValue (timer);
    }
  else
    {
      //alignment_time[0] += cutGetTimerValue (timer);
    }


  //CUT_SAFE_CALL (cutDeleteTimer (timer));


  // cleaning up

    // GPU

/*
    checkCudaErrors (cudaFree (main_alloc_ptr));
    checkCudaErrors (cudaFree (read_alignments_on_device));
    checkCudaErrors (cudaFree (genome_alignments_on_device));
    checkCudaErrors (cudaFree (read_sequences_on_device));
    checkCudaErrors (cudaFree (genome_sequences_on_device));
    checkCudaErrors (cudaFree (alignment_scores_on_device));
    checkCudaErrors (cudaFree (direction_array_on_device));
    checkCudaErrors (cudaFree (edit_offsets_on_device));
    checkCudaErrors (cudaFree (genome_offsets_on_device));
    checkCudaErrors (cudaFree (read_offsets_on_device));
*/

    // host

    free (alignment_scores_on_host);
    free (read_alignments_on_host);
    free (genome_alignments_on_host);


  //cudaThreadExit (); /seems to kill cuda all cuda jobs (including active ones)
  //  cudaThreadSynchronize();
  //  sleep(10);
}

extern "C"
unsigned long int
get_cuda_infos(void){

   unsigned long int free, total;
   int gpuCount, i, best_gpu;
   long int most_vram_free;
   cudaDeviceProp deviceProp;
   CUresult res;
   CUdevice dev;
   CUcontext ctx;
   best_gpu = 0;
   most_vram_free = 0;



   cuInit(0);
   cuDeviceGetCount(&gpuCount);
   fprintf (stderr,"\n=== CUDA hardware ===\n");
   fprintf (stderr,"Detected %d GPU(s):\n\n",gpuCount);

       for (i=0; i<gpuCount; i++)
       {
       cuDeviceGet(&dev,i);
       cuCtxCreate(&ctx, 0, dev);
       cudaGetDeviceProperties(&deviceProp,dev);
       res = cuMemGetInfo (&free, &total);
       if (free > most_vram_free){

	most_vram_free=free;
                //fprintf (stderr, "free: %d vram: %d\n", free,most_vram_free);

	best_gpu=i;
       }
       fprintf (stderr, "GPU %d: %s\t", i,deviceProp.name);


  if (res == CUDA_SUCCESS)
    {

      fprintf (stderr, "[%*ld/%*ld MB available]\n",4,(unsigned long int)free/(1024*1024),4,(unsigned long int)total/(1024*1024));

      cuCtxDetach (ctx);

    }
  else
    {

      fprintf (stderr, "Could not get memory information from GPU (CUDA status code = %x)\n", res);
      fprintf (stderr, "Probably not enough free memory on the GPU for CUDA startup.\n");
      fprintf (stderr, "Try to disable secondary monitors or suspend 3d desktop effects.\n");
      fprintf (stderr, "Shutting down.\n");
      exit(-1);
    }

       cuCtxDetach(ctx);
       }
       
       if (most_vram_free < 128000000){
            fprintf (stderr,"\n-> SARUMAN was not able to find a GPU with enough VRAM available! (available: %d MB)\n",(unsigned long int)most_vram_free/(1024*1024));
            exit(1);
       } else if (use_gpu_value <= i ){
       fprintf (stderr,"\n-> SARUMAN was forced to run on GPU %d\n",use_gpu_value);
       } else {
       fprintf (stderr,"\n-> SARUMAN will run on GPU %d\n",best_gpu);
       use_gpu_value=best_gpu;
       }

  return most_vram_free;
}
