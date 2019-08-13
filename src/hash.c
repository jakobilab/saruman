/*
 * 
 * hash.c
 *
 *  Created on: May 4, 2009
 *      Author: Tobias Jakobi
 */

#include "util.h"   // imports from util class
#include "hash.h"
#include "saruman.h"
#include "cuda-functions.h"   // hash prototypes

#include <stdio.h>
#include <string.h>  // printf etc
#include <stdlib.h> // std. stuff
#include <math.h>   // ceil()
#include <unistd.h>

#include <string.h>
#include <argtable2.h>
#include <pthread.h>

// global and shared variables

struct qgram_hash *qgrams = NULL; /*!< holds the qgram index */

int start_pos_buffer = 100;


struct qgram_hash *
qgram_lookup (char* qgram)
{
  struct qgram_hash *s = NULL;
  //HASH_FIND_STR (qgrams, qgram, s);
  HASH_FIND (hh,qgrams, qgram, qgram_length ,s);
  return s;
}

unsigned short int
add_start_position (unsigned int read_id, int start_position, char direction_and_qgrams, struct result_hash **results, unsigned int * alignment_queue, int possible)
{

  //fprintf (stderr,"add start position: %d (%d)\n",start_position,exact_tresh);
  
  if (possible == exact_tresh){

  struct result_hash *result = NULL;
  int* start_positions = NULL;
  char* direction_and_qgram_array = NULL;
  char direction;


  if (direction_and_qgrams > 3)
    {
      direction = 1; //reverse
    }
  else
    {
      direction = 0; //forward
    }

  // lookup key == start
  HASH_FIND_INT (*results, &read_id, result);

  // no entry, create new one
  if (result == NULL)
    {
      // allocate memory for new entry in qgram hash
      result = (struct result_hash*) malloc (sizeof (struct result_hash));

      result->read_id = read_id;

      start_positions = malloc (start_pos_buffer * sizeof (int));
      direction_and_qgram_array = malloc (start_pos_buffer * sizeof (char));

      start_positions[0] = start_position;
      direction_and_qgram_array[0] = direction_and_qgrams;

      result->next_pos = 1;
      result->positions = start_positions;
      result->direction_and_matching_qgrams = direction_and_qgram_array;

      //printf ("new entry: %d (%d)\n",start_position,direction_and_qgrams);

      HASH_ADD_INT (*results, read_id, result);
      if(start_position > 0){
      *alignment_queue=*alignment_queue+1;
        }

      return 0;

      // read exists - start position for this read exists, too
    }
  else
    {

      int pos = result->next_pos;

      start_positions = result->positions; // get old pointer
      direction_and_qgram_array = result->direction_and_matching_qgrams;

      for (int current_position = 0; current_position < pos; current_position++)
        {
          char tmp_direction;
          if (direction_and_qgram_array[current_position] > 3)
            {
              tmp_direction = 1; //reverse
            }
          else
            {
              tmp_direction = 0; //forward
            }
          // start position and direction match

          //printf ("would update matching qgrams from %d to %d\n", direction_and_qgram_array[current_position], direction_and_qgrams);

          for (int i = maximal_errors_allowed * (-1); i <= maximal_errors_allowed; i++)
            {
              if ((start_position + i == start_positions[current_position] || start_position + i == start_positions[current_position]*(-1)) &&
                  direction_and_qgram_array[current_position] <= direction_and_qgrams)
                {
                  //printf ("rejecting %d because %d lies in error (%d -> %d) bounds, not adding\n",read_id, start_position,
                  //start_position-maximal_errors_allowed,start_position+maximal_errors_allowed);
                  return 1;
                }
            }

            // got a better hit in terms of flanking qgrams, update hit to save some alignment space later
          if ((start_positions[current_position] == start_position && tmp_direction == direction && direction_and_qgram_array[current_position] > direction_and_qgrams))
            {

              result->direction_and_matching_qgrams[current_position] = direction_and_qgrams;
              //printf ("updating matching qgrams from %d to %d\n", direction_and_qgram_array[current_position], direction_and_qgrams);
              return 2;
            }

          else if ((start_positions[current_position] == start_position && tmp_direction == direction) ||
                   // start position and direction match AND this read already has been aligned
                   (start_positions[current_position] == start_position * (-1) && tmp_direction == direction))
            {
              //printf("position %d exists\n",start_position);
              return 1;
            }


        }
      //realloc memory
      //printf("pos is %d|%d\n",pos,read_id);

      if (pos % start_pos_buffer == 0){
          //printf("old size: %d, new size: %d (%d)\n",((pos/start_pos_buffer)),((pos/start_pos_buffer)+1),pos);

          start_positions = (int *) realloc (start_positions, sizeof (int) * ((pos/start_pos_buffer)+1)*start_pos_buffer);
          direction_and_qgram_array = (char *) realloc (direction_and_qgram_array, sizeof (char) * ((pos/start_pos_buffer)+1)*start_pos_buffer);
      }


      // set new value
      start_positions[pos] = start_position;
      direction_and_qgram_array[pos] = direction_and_qgrams;

      result->next_pos++;

      result->positions = start_positions;
      result->direction_and_matching_qgrams = direction_and_qgram_array;

      if(start_position > 0){
      *alignment_queue=*alignment_queue+1;
        }
      return 2;
    }

    }else {
      return 1;
      }

}

unsigned short int
add_qgram_to_index (char *qgram, int start_position_paramater)
{
  struct qgram_hash *qgram_pointer = NULL;
  unsigned int* start_positions = NULL;

  // lookup key
  qgram_pointer = qgram_lookup (qgram);

  // no entry, create new one
  if (qgram_pointer == NULL)
    {

      // allocate memory for new entry in qgram hash
      qgram_pointer = (struct qgram_hash*) malloc (sizeof (struct qgram_hash));

      strcpy (qgram_pointer->qgram, qgram);

      start_positions = calloc (1, sizeof (int));

      start_positions[0] = start_position_paramater;
      pos_counter++;

      qgram_pointer->next_pos = 1;
      qgram_pointer->positions = start_positions;

      // add pointer to a char array (aka qgram) as key
      // handle name, head, key_ptr, key_len, item_ptt

      //HASH_ADD_STR (qgrams, qgram, qgram_pointer);
      HASH_ADD (hh,qgrams, qgram, qgram_length ,qgram_pointer);


      //printf ("%s %d\n",qgram,start_position_paramater);

      return 0;

      //entry exists
    }
  else
    {
      //printf ("%s %d\n",qgram,start_position_paramater);

      //fetch secondary qgram

      start_positions = qgram_pointer->positions; // get old pointer

      //realloc memory
      start_positions = realloc (start_positions, sizeof (int) *((qgram_pointer->next_pos) + 2));
      pos_counter++;

      // set new value
      start_positions[qgram_pointer->next_pos] = start_position_paramater;

      qgram_pointer->positions = start_positions;

      qgram_pointer->next_pos++;
      return 1;
    }

}

void
print_qgram_index_statistics (int qgram_size)
{


  hash_members = HASH_COUNT (qgrams);
  float memory = pos_counter * sizeof (int) + hash_members * sizeof (struct qgram_hash) +
          hash_members * sizeof (char*) + hash_members * sizeof (int*) + hash_members*qgram_size;


  fprintf (stderr, "\nqgram index statistics:\n");
  fprintf (stderr, "======================================\n");
  fprintf (stderr, "%d different qgrams.\n", hash_members);
  fprintf (stderr, "%d start positions total.\n", pos_counter);
  fprintf (stderr, "estimated memory usage: %.0f MB\n", memory / 1024 / 1024);
  fprintf (stderr, "======================================\n\n");

}

short int
perform_hashtable_identity_check (void)
{
  struct qgram_hash *s;

  int count = 0;
  int total = 0;

  // iterate through the complete hash and perform a lookup for each key
  for (s = qgrams; s != NULL; s = s->hh.next)
    {
      s = qgram_lookup (s->qgram);
      if (s != NULL)
        {
          count++; // hit
        }
      total++;
    }

  // 100% -> test passed else  -> bullshit

  printf ("qgram index indentity self check: %s.\n", ((float) 100 * count / total == 100 ? "passed" : "FAIL"));

  if ((float) 100 * count / total == 100)
    {
      return 0;
    }
  else
    {
      return -1;
    }

}

void *
cuda_daemon(void *threadarg) {


   int *id_ptr, do_reverse;
   id_ptr = (int *) threadarg;
   do_reverse = *id_ptr;

   cuda_allocate_memory(do_reverse);

   

    while (1) {



        if (do_reverse == 4 && reverse_thread_check == 1){

            reverse_thread_check = 2;
                   fprintf(stderr, "\t-> processing reverse alignments...\n");
            qgram_to_read_mapping(reads, read_length_value, qgram_length,
                    qgram_positions_in_read, number_of_position_in_read,
                    rest, do_reverse, chunk_value, number_of_reads_read,
                    maximal_errors_allowed);
                        fprintf(stderr, "\t-> finished reverse alignments...\n");

            reverse_thread_check = 3;


        } else if (do_reverse == 0 && forward_thread_check == 1) {


            forward_thread_check = 2;
                  fprintf(stderr, "\t-> processing forward alignments...\n");
            qgram_to_read_mapping(reads, read_length_value, qgram_length,
                    qgram_positions_in_read, number_of_position_in_read,
                    rest, do_reverse, chunk_value, number_of_reads_read,
                    maximal_errors_allowed);
            fprintf(stderr, "\t-> finished forward alignments...\n");

            forward_thread_check = 3;

        } if (reverse_thread_check == -1 || forward_thread_check == -1){
            cuda_free_memory(do_reverse);
            break;
        }
            sleep(1);

    }
    return NULL;
}


void *
qgram_to_read_mapping (char * read_array, unsigned short int readlength,   unsigned short int qgram_length,
    unsigned short int* qgram_start_positions, unsigned short int number_of_start_positions ,
    unsigned short int rest ,unsigned short int do_reverse, int chunksize, int read_number,
    int error)
{

  unsigned int alignment_queue = 0; /*!< holds the number of alignments in queue to be aligned */

  // run selfcheck before mapping
  //perform_hashtable_identity_check ();

  
  // offset = 2d matrix for each read holding all start positions per qgram
  // unsigned int offset = number_of_start_positions * max_counter;

  // compute number of runs needed to process all reads for a given chunksize
  unsigned short int max_runs = (unsigned short int) ceil (read_number / (float) chunksize);

  //printf("\n\nchunksize is %d, with %d reads total -> splitting into %d runs\n", chunksize, read_number, max_runs);

  unsigned int * data_array;

  struct result_hash *results = NULL; // initialize result hash

  unsigned int * data_array_offsets;

  for (unsigned int run = 0; run < max_runs; run++)
    {
      unsigned int data_array_index = 1;
      unsigned int allocation_buffer_default = chunksize*10; // initial size
      unsigned int allocation_buffer = allocation_buffer_default;

      //fprintf(stderr,"allocating space for  %d items\n",chunksize * number_of_start_positions);


      data_array_offsets = calloc (chunksize * number_of_start_positions, sizeof (unsigned int));
      data_array = calloc (allocation_buffer, sizeof (unsigned int));

      // no memory available
      if (data_array == NULL)
        {
          fprintf (stderr, "Could not allocate enough memory for data array\n");
          exit(-1);
        }

      //printf("chunk %d of %d %s\n\n", run + 1, max_runs, (do_reverse ? "(reverse: <<<)" : "(forward: >>>)"));

      for (unsigned int i = run * chunksize, p = 0; i < (run + 1) * chunksize; i++, p+=number_of_start_positions)
        {

      //printf("offset of read %d is %d \n",p,data_array_offsets[p]);



          if (i < (unsigned int)read_number)
            {
              char actual_read[readlength + 1];

              strncpy (actual_read, &read_array[i * readlength], readlength);
              actual_read[readlength] = '\0';

              if (do_reverse == 4)
                {
                  reverse_complement (actual_read, readlength);
                }

              //printf("read %d: %s\n", i + 1, actual_read);
 
              for (unsigned int x = 0; x < number_of_start_positions; x++)
                {
                  data_array_offsets[p+x]=data_array_index-1;

                  char qgram[qgram_length + 1];
                  strncpy (qgram, &actual_read[qgram_start_positions[x]], qgram_length);
                  qgram[qgram_length] = '\0';

                  struct qgram_hash *s;
                  s = qgram_lookup (qgram);

                  //printf("qgram #%d: %s - ", x + 1, qgram);
                  //printf("%s ", (s ? "found @ " : "not found"));

                    if (s) {
                        unsigned int * positions = s->positions;
                        int maxpos = s->next_pos;

                        for (int z = 0; z < maxpos; z++) {
                            //printf("%d (%d) /", positions[z],items_to_allocate-1);
                            data_array[data_array_index - 1] = positions[z];
                            data_array_index++;
                            allocation_buffer--;

                            if (allocation_buffer == 0) {
                                data_array = (unsigned int*) realloc(data_array, (data_array_index + allocation_buffer_default) * sizeof (unsigned int));
                                memset (&data_array[data_array_index],0,allocation_buffer_default * sizeof (unsigned int));
                                allocation_buffer = allocation_buffer_default;
                            }
                        }
                    }
                    data_array[data_array_index - 1] = 0;
                    //printf(" %d (%d)", 0,items_to_allocate-1);
                    data_array_index++;
                    allocation_buffer--;
                    
                    if (allocation_buffer == 0) {
                        data_array = (unsigned int*) realloc(data_array, (data_array_index + allocation_buffer_default) * sizeof (unsigned int));
                        memset (&data_array[data_array_index],0,allocation_buffer_default * sizeof (unsigned int));
                        allocation_buffer = allocation_buffer_default;
                    }
                    //printf("\n");
                }
            }
        }

      results = run_seashore (chunksize, run, number_of_start_positions, qgram_length, error, rest, do_reverse, data_array, results, &alignment_queue,data_array_offsets);
      free (data_array);
      free (data_array_offsets);
    }

  prepare_cuda_alignments (sw_mismatch_value, sw_match_value, error, results,&alignment_queue,do_reverse);
  clear_result_hash (results);
  return NULL;
}

struct result_hash *
run_seashore (int chunksize, int run, int number_of_start_positions, int qgram_length, int error, int rest, int do_reverse, unsigned int * data_array, struct result_hash *results, unsigned int * alignment_queue, unsigned int * data_array_offsets)
{

  int number_of_qgram_positions = number_of_start_positions;

  unsigned char first_qgram_matches = 0;

  short int mapped = 0;

  int read_number = number_of_reads_read;

 // unsigned short int max_runs = (unsigned short int) ceil (read_number / (float) chunksize);

  int limit = read_number;

  if (((run + 1) * chunksize) > limit){

      limit = read_number % chunksize;
            //fprintf (stderr,"limiting to %d! %d %d %d %d\n",limit,read_number,chunksize,run, max_runs);
  }

  for (int read_id = run * chunksize, counter = 0, counter2=0; (counter2<limit) && (read_id < (run + 1) * chunksize); read_id++, counter2++, counter+=number_of_qgram_positions)
    { // for each read


         //fprintf (stderr,"# reads: %d\n",strlen(reads)/read_length_value);


          char tmp_read[read_length_value + 1];
          strncpy (tmp_read, &reads[read_length_value * read_id], read_length_value);
          tmp_read[read_length_value] = '\0';
          //fprintf (stderr,"strcpy from %d to %d (%s | %d)\n",read_length_value * read_id,read_length_value * read_id+read_length_value,tmp_read,counter2);

          unsigned int actual_qgram_start_position = 0;
          unsigned int start_position_offset = 0;
          unsigned int in_read_start_position = 0;
          int genome_start = 0;

      // 2 qgrams of the first set match
      for (int actual_qgram = 0; actual_qgram < (number_of_qgram_positions / 2); actual_qgram++)
        { // for each of the X qgrams

          if (actual_qgram == 0)
            {
              first_qgram_matches = 1;
            }

          // computing start position in the read
          in_read_start_position = actual_qgram * qgram_length + 1; // compute startposition

          start_position_offset = 0;

          // for each genome start position of this qgram

          // find actual qgram start position
         actual_qgram_start_position = data_array[data_array_offsets[counter + actual_qgram] + start_position_offset];
                        //fprintf(stderr,"qgram is: %d|%d\n", actual_qgram,do_reverse);


          // position MUST be > 0; 0 means there are now other positions following in this set
          while (actual_qgram_start_position > 0) // 0 if there are no other positions left
            {
              //fprintf (stderr,"new cycle, pos is %d/%d/%d\n",actual_qgram_start_position, in_read_start_position,actual_qgram);

              // compute start position of read in genome
              genome_start = actual_qgram_start_position - in_read_start_position;
              // printf ("actual qgram pos is %d, actual in read start is %d\n", actual_qgram_start_position, in_read_start_position);

              //fprintf (stderr,"pos will be: %d \n",genome_start + 1 - genome_chunk_run*genome_chunk_value);

              // check if this will be a valid start point (has to be> 0)
              if (((genome_start + 1 - genome_chunk_run*(int)genome_chunk_value) > 0))
                {
                  // filter out all perfect hits

                  // gather read and genome sequence
                  char tmp_genome[read_length_value + 1];
                  //fprintf (stderr,"pos is: %d (%d,%d)\n",genome_start + 1 - genome_chunk_run*genome_chunk_value, genome_start +1,genome_chunk_run*genome_chunk_value);
                  strncpy (tmp_genome, &genome_sequence[genome_start + 1 - genome_chunk_run*genome_chunk_value], read_length_value);
                  tmp_genome[read_length_value] = '\0';
                  char tmp_read[read_length_value + 1];
                  strncpy (tmp_read, &reads[read_length_value * read_id], read_length_value);
                  tmp_read[read_length_value] = '\0';

                  // reverse complement of the sequence
                  if (do_reverse == 4)
                    {
                      reverse_complement (tmp_read, read_length_value);
                    }

                  // check if this is a perfect match
                  int cmp = strcmp (tmp_genome, tmp_read);

                  // find the read name
                  struct read_name_hash *hash_entry = NULL;
                  HASH_FIND_INT (read_names, &read_id, hash_entry);

                  // perfect if cmp == 0
                  if (cmp == 0)
                    {
                      if (add_start_position (read_id, (genome_start + 1)*(-1), (char) 0, &results, alignment_queue,exact_tresh) != 1)
                        {

                          unsigned int stop = genome_start + read_length_value;

                          // print out perfect reads
                          printf ("%s\t%d\t%d\t%s\t%s\t%s\t%d\n", hash_entry->name, genome_start + 1, stop, (do_reverse ? "<<" : ">>"), tmp_genome, tmp_read, 0);
                          mapped++;

                        }



                      // the next position is 0 aka empty -> break this complete loop and go to the next qgram
                      if (actual_qgram_start_position == 0)
                        {
                          // get the next start position
                          start_position_offset++; // one step forward
                          actual_qgram_start_position = data_array[data_array_offsets[counter + actual_qgram] + start_position_offset];
                          break;
                        }
                  }
                    //not perfect hit, mapping it manually
                  
                    

                      int possible = 1;
                      // perfect hits filtered, below only normal hits should be found

                     // fprintf(stderr,"from %d to %d\n", actual_qgram + 1,(number_of_qgram_positions / 2));


                      for (int next_qgram_number = actual_qgram + 1; next_qgram_number < (number_of_qgram_positions / 2); next_qgram_number++)
                        { // for each qgram following the actual
                           
                          //fprintf(stderr,"testing qgram pair  %d <> %d (first pass)\n", actual_qgram+1,next_qgram_number+1);


                          unsigned int read_start = (next_qgram_number) * qgram_length + 1; // start position of the next qgram
                          // shake left <> right to allow errors

                          // fprintf(stderr,"error shake from %d to %d\n",  abs ((error - (actual_qgram + 1) - 1))*(-1),abs ((error - (actual_qgram + 1) - 1)));
                          for (int error_zone = (error *(-1)); error_zone <= error ; error_zone++)

                            {
/*
                             fprintf(stderr,"counter: %d\n",counter);
                             fprintf(stderr,"next qgram number: %d\n",next_qgram_number);
                             fprintf(stderr,"start_position_offset: %d\n",start_position_offset);
                             fprintf(stderr,"data array offsets: %d\n",data_array_offsets[counter + next_qgram_number] + start_position_offset);
*/


                              // check if that value matches a start position of the next gram
                              unsigned int next_qgram_start_position = data_array[data_array_offsets[counter + next_qgram_number] + start_position_offset];
                              unsigned int next_qgram_start_position_offset = 0;
                              //fprintf(stderr,"lala: %d\n",data_array_offsets[counter + next_qgram_number] + start_position_offset);
                              while (next_qgram_start_position > 0)
                                {
                                  char flag = -1;
                                 //fprintf(stderr,"start position is now: %d, should be %d ()\n",genome_start + read_start + error_zone,next_qgram_start_position);

                                  if (next_qgram_start_position == genome_start + read_start + error_zone)
                                    {
                                      // check which qgrams matched and flag appopriate
                                      //printf("actual qgram %d, lal: %d\n",next_qgram_number,(number_of_qgram_positions / 2) - 1);
                                      if ((next_qgram_number == (number_of_qgram_positions / 2) - 1) && first_qgram_matches)
                                        {
                                          flag = 1 + do_reverse;
                                          //printf ("(first pass) first and last matching [%d]\n", genome_start + 1 - error);
                                          possible++;
                                          mapped++;
                                          add_start_position (read_id, genome_start + 1, flag, &results, alignment_queue,possible);
                                        }
                                      else if (actual_qgram == (number_of_qgram_positions / 2) - 1)
                                        {

                                          if ((int) (actual_qgram_start_position + qgram_length) > read_length_value)
                                            {
                                              flag = 3 + do_reverse;
                                              //printf ("(first pass) last matching [%d]\n", genome_start + 1 - error);
                                              possible++;
                                              mapped++;
                                              add_start_position (read_id, genome_start - error + 1, flag, &results, alignment_queue,possible);

                                            }
                                          else
                                            {
                                              flag = 2 + do_reverse;
                                              //printf ("(first pass) last matching [%d]\n", genome_start + 1 - error);
                                              possible++;
                                              mapped++;
                                              add_start_position (read_id, genome_start - error + 1, flag, &results, alignment_queue,possible);
                                            }

                                        }
                                      else if (first_qgram_matches)
                                        {
                                          flag = 1 + do_reverse;
                                          //printf ("(first pass) first matching [%d]\n", genome_start + 1);
                                          possible++;
                                          mapped++;
                                          add_start_position (read_id, genome_start + 1, flag, &results, alignment_queue,possible);
                                        }
                                      else
                                        {
                                          flag = 3 + do_reverse;
                                          //printf ("(first pass) nothing matching [%d]\n", genome_start + 1);
                                          possible++;
                                          mapped++;
                                          add_start_position (read_id, genome_start - error + 1, flag, &results, alignment_queue,possible);
                                        }

                                    }
                                  next_qgram_start_position_offset++;
                                  next_qgram_start_position = data_array[data_array_offsets[counter + next_qgram_number] + next_qgram_start_position_offset];
                                }
                            }
                        }

                      // 1 qgram from first set and 1 qgram from second set match
                      for (int next_qgram_number = 0; next_qgram_number < (number_of_qgram_positions / 2); next_qgram_number++)
                        { // for each qgram following the actual

                         // fprintf(stderr,"testing qgram pair  %d <> %d (second pass)\n", actual_qgram+1,next_qgram_number+3);


                          // shake left <> right to allow errors

                              for (int err = abs ((number_of_qgram_positions / 2) -1 - next_qgram_number - actual_qgram)*(-1); err <= abs (((number_of_qgram_positions / 2) -1 -next_qgram_number - actual_qgram)); err++)
                            {

                              int actual_error = err;
                              if (next_qgram_number - actual_qgram < 0)
                                {
                                  actual_error = err * (-1);
                                }

                              // start position of the next qgram
                              unsigned int next_qgram_in_read_start_position = ((number_of_qgram_positions / 2) - next_qgram_number - 1) * qgram_length + 1 + rest;

                              // check if that value matches a start position of the next gram

                              unsigned int start_position_offset = 0;
                              unsigned int next_qgram_start_position = data_array[data_array_offsets[counter  + ((number_of_qgram_positions / 2) + next_qgram_number)] + start_position_offset];

                              while (next_qgram_start_position > 0)
                                {
                                  //fprintf(stderr,"start position is now: %d, should be %d ()\n",genome_start + v + actual_error,next_qgram_start_position);

                                  //fprintf(stderr,"start position is now: %d, should be (%d) (%d/%d) \n",genome_start + next_qgram_in_read_start_position + actual_error,next_qgram_start_position,err,actual_error);


                                  //printf ("phase two, pos is %d\n",next_qgram_start_position);
                                  //printf ("genome pos is %d\n",genome_start);
                                  if (next_qgram_start_position == genome_start + next_qgram_in_read_start_position + actual_error)
                                    {
                                              // printf ("phase two, pos is %d\n",next_qgram_start_position);
                                  //printf ("genome pos is %d\n",genome_start);

                                      char flag = -1;
                                      // check which qgrams matched and flag appopriate
                                      if (next_qgram_number == 0 && first_qgram_matches)
                                        {
                                          flag = 1 + do_reverse;
                                          //printf ("(second pass) first and last matching [%d]\n", genome_start + 1 - error);
                                          possible++;mapped++;
                                          add_start_position (read_id, genome_start + 1, flag, &results, alignment_queue,possible);
                                        }
                                      else if (next_qgram_number == 0)
                                        {

                                         // printf ("qgram: %.*s\n", qgram_length, &tmp_read[v]);
                                         // printf ("length:%d\n", (int) (v + qgram_length));
                                          if ((int) (next_qgram_in_read_start_position + qgram_length) > read_length_value)
                                            {
                                              flag = 3 + do_reverse;
                                              //printf ("(second pass) last matching [%d]\n", genome_start + 1 - error);
                                              possible++;
                                              mapped++;
                                              add_start_position (read_id, genome_start - error + 1, flag, &results, alignment_queue,possible);

                                            }
                                          else
                                            {
                                              flag = 2 + do_reverse;
                                             // printf ("(second pass) last matching [%d]\n", genome_start + 1 - error);
                                              possible++;
                                              mapped++;
                                              add_start_position (read_id, genome_start - error + 1, flag, &results, alignment_queue,possible);
                                            }

                                        }

                                      else if (first_qgram_matches)
                                        {
                                          flag = 1 + do_reverse;
                                          //printf ("(second pass) first matching [%d]\n", genome_start + 1);
                                          possible++;
                                          mapped++;
                                          add_start_position (read_id, genome_start + 1, flag, &results, alignment_queue,possible);
                                        }
                                      else
                                        {
                                          flag = 3 + do_reverse;
                                          //printf ("(second pass) nothing matching [%d]\n", genome_start + 1 - error);
                                          possible++;
                                          mapped++;
                                          add_start_position (read_id, genome_start - error + 1, flag, &results, alignment_queue,possible);
                                        }
                                      // add canidate to result list

                                    }
                                  start_position_offset++;
                                  next_qgram_start_position = data_array[data_array_offsets[counter + ((number_of_qgram_positions / 2) + next_qgram_number)] + start_position_offset];

                                }

                            }

                        }

                      start_position_offset++; // one step forward
                      actual_qgram_start_position = data_array[data_array_offsets[counter + actual_qgram] + start_position_offset];



                  
                }
              else
                {

                  start_position_offset++; // one step forward
                  actual_qgram_start_position = data_array[data_array_offsets[counter + actual_qgram] + start_position_offset];
                }
            }
            first_qgram_matches = 0;
        }
       // }
        // output not aligned reads if not yet printed in any other thread
        if (mapped == 0 && print_non_matching_reads == 1) {

            struct read_name_hash *hash_entry = NULL;
            HASH_FIND_INT(read_names, &read_id, hash_entry);

            // test for null to avoid segfault
            if (hash_entry != NULL) {
                printf("%s\t*\t*\t%s\t%s\t*\t%d\n", hash_entry->name, (do_reverse ? "<<" : ">>"), tmp_read, read_length_value);

                 // forward thread
                if (do_reverse == 0) {
                    forward_output_control[read_id] = 1; // this read has not been mapped by the forward thread (as forward)

                 // reverse thread
                } else if (do_reverse == 4) {
                    reverse_output_control[read_id] = 1; // this read has not been mapped by the reverse thread (as reverse)
                }
            }
        } else {
            mapped = 0;
        }
        cuda_queue_check(sw_mismatch_value, sw_match_value, error, results, alignment_queue, do_reverse);
    }
  return results;

}

void
cuda_queue_check (int mismatch, int match, int error, struct result_hash *results, unsigned int * alignment_queue, int do_reverse)
{

  if (*alignment_queue > (unsigned)chunk_value)
    {
      prepare_cuda_alignments (mismatch, match, error, results,alignment_queue, do_reverse);
      *alignment_queue = 0;
    }

}

void
clear_result_hash (struct result_hash *results) {

  struct result_hash *current_entry;

  while(results) {
    current_entry = results;          /* copy pointer to first item     */
    HASH_DEL(results,current_entry);  /* delete; users advances to next */
    free(current_entry->positions);
    free(current_entry->direction_and_matching_qgrams);
    free(current_entry);              /* optional- if you want to free  */
  }

 }

void
clear_name_hash (void)
{

  struct read_name_hash *current_name;

  while (read_names)
    {
      current_name = read_names; /* copy pointer to first item     */
      HASH_DEL (read_names, current_name); /* delete; users advances to next */
      free (current_name); /* optional- if you want to free  */
    }

}

void
clear_qgram_index (void)
{

  struct qgram_hash *current_name;
  while (qgrams)
    {
      current_name = qgrams; /* copy pointer to first item     */
      HASH_DEL (qgrams, current_name); /* delete; users advances to next */
      free(current_name->positions);
      free (current_name); /* optional- if you want to free  */
    }

}

void
prepare_cuda_alignments (int mismatch, int match, int error, struct result_hash *results, unsigned int * alignment_queue, int do_reverse)

{


  if (*alignment_queue > 0)

    {
      struct result_hash *s; // holds the result of the hash lookup
      int counter = 0; //... counts...
      unsigned int start_offset = 0;
      unsigned int end_offset = 0;
      char tmp_direction = 0;
      size_t genome_array_size = 0;
      size_t edit_matrix_size = 0;
      int edit_size = 0;

      char* genome_array = malloc (chunk_value*10*(read_length_value+2*error) * sizeof (char)); // holds the genome sequences to be aligned
      char* read_array = malloc (chunk_value*10*(read_length_value) * sizeof (char)); // holds the read sequences to be aligned
      int* id_array = malloc (chunk_value*10*sizeof (int)); // holds the unique read ids
      int* pos_array = malloc (chunk_value*10*sizeof (int)); // holds the start positions in the genome
      char* direction_array = malloc (chunk_value*10*sizeof (char)); // forward or backwards
      int* genome_offsets = malloc (chunk_value*10*sizeof (int)); // stop positions of the genome sequences
      int* read_offsets = malloc (chunk_value*10*sizeof (int)); // stop positions of the read sequences
      int* edit_offsets = malloc (chunk_value*10*sizeof (int)); // start positions of the edit matrices

      // for each read in the result hash
      for (s = results; s != NULL; s = s->hh.next)
        {
          int id = s->read_id;
          //for each start position of this read
          //int* positions = s->positions; // load positions into this array
          int nextpos = s->next_pos;
          for (int i = 0; i < nextpos; i++)
            {
              // check if this position should be aligned
              if (s->positions[i] > 0 ) // remove -> here to improve performance
                {
                  //printf("%d\n",s->direction_and_matching_qgrams[i]);
                  switch (s->direction_and_matching_qgrams[i])
                    {
                    case 0:
                      tmp_direction = 0;
                      break;

                    case 1:
                      tmp_direction = 0;
                      end_offset = error;
                      break;

                    case 2:
                      tmp_direction = 0;
                      start_offset = error;
                      break;

                    case 3:
                      tmp_direction = 0;
                      start_offset = error;
                      end_offset = error;
                      break;

                    case 4:
                      tmp_direction = 1;
                      break;

                    case 5:
                      tmp_direction = 1;
                      end_offset = error;
                      break;

                    case 6:
                      tmp_direction = 1;
                      start_offset = error;
                      break;

                    case 7:
                      tmp_direction = 1;
                      start_offset = error;
                      end_offset = error;
                      break;
                    }

                  // do not try to align "out of genome" sequences
                  if(s->positions[i]+read_length_value + start_offset + end_offset > (unsigned)total_bases_read){
                      continue;
                  }

                  genome_array_size += start_offset + read_length_value + end_offset + error;
                  //genome_array = realloc (genome_array, genome_array_size);

                  // printf ("allocated %d for genome array\n",genome_array_size);

                  if (tmp_direction == 0)
                    { //forward direction


                      //printf ("position is: %d\n",s->positions[i]);
                      // extract genome sequence
                      strncpy (&genome_array[genome_array_size - read_length_value - start_offset - end_offset - (counter + 1) * error],
                               &genome_sequence[s->positions[i]-genome_chunk_run*genome_chunk_value], read_length_value + start_offset + end_offset);

                      //char tmp[read_length_value + 1];
                      //char tmp2[read_length_value + start_offset + end_offset + 1];


                      strncpy (&read_array[counter * (read_length_value)], &reads[(id) * read_length_value], read_length_value);
                      //strncpy (tmp2, &genome_sequence[s->positions[i]], read_length_value + start_offset + end_offset);

                      //strncpy (tmp, &reads[(id) * read_length_value], read_length_value);
                      //tmp[read_length_value] = '\0';
                      //tmp2[read_length_value + start_offset + end_offset] = '\0';
                      //printf ("genome\t%s (%d,%d,%d,%d)\n", tmp2, start_offset, end_offset,s->positions[i],s->direction_and_matching_qgrams[i]);
                      //printf ("read\t%s\n", tmp);
                    }
                  else
                    { // reverse direction

                      // hold the original read
                      char tmp[read_length_value + 1];
                      //char tmp2[read_length_value + start_offset + end_offset + 1];



                      // extract genome sequence
                      strncpy (&genome_array[genome_array_size - read_length_value - start_offset - end_offset - (counter + 1) * error],
                               &genome_sequence[s->positions[i]-genome_chunk_run*genome_chunk_value], read_length_value + start_offset + end_offset);

                      // extract read sequence
                      strncpy (tmp, &reads[(id) * read_length_value], read_length_value);

                      // build reverse complement
                      reverse_complement (tmp, read_length_value);
                      //strncpy (tmp2, &genome_sequence[s->positions[i]], read_length_value + start_offset + end_offset);

                      // extract read sequence
                      strncpy (&read_array[counter * (read_length_value)], tmp, read_length_value);


                      //tmp[read_length_value] = '\0';
                      //tmp2[read_length_value + start_offset + end_offset] = '\0';
                      //printf ("genome\t%s (%d,%d,%d,%d)\n", tmp2, start_offset, end_offset,s->positions[i],s->direction_and_matching_qgrams[i]);
                      //printf ("read\t%s\n", tmp);

                    }

                  // unique read id
                  id_array[counter] = id;
                  pos_array[counter] = s->positions[i];
/*
                  if (genome_chunk_run == 0){
                      pos_array[counter]++;
                    }
*/

                  direction_array[counter] = s->direction_and_matching_qgrams[i];

                  edit_matrix_size += (read_length_value + start_offset + end_offset + 1)*(read_length_value + 1);
                  genome_offsets[counter] = genome_array_size - (counter + 1) * error;
                  read_offsets[counter] = (counter + 1)*(read_length_value);
                  edit_offsets[counter] = edit_size;
                  edit_size += (read_length_value + 1)*((start_offset + read_length_value + end_offset) + 1);


                  counter++;
                  // allocate enough memory for the next cycle
/*
                  read_array = realloc (read_array, (counter + 1)*(read_length_value));
                  id_array = realloc (id_array, (counter + 1) * sizeof (int));
                  pos_array = realloc (pos_array, (counter + 1) * sizeof (int));



                  genome_offsets = realloc (genome_offsets, (counter + 1) * sizeof (int));
                  read_offsets = realloc (read_offsets, (counter + 1) * sizeof (int));
                  edit_offsets = realloc (edit_offsets, (counter + 1) * sizeof (int));

                  direction_array = realloc (direction_array, (counter + 1) * sizeof (char));
*/
                  s->positions[i] = s->positions[i]*(-1); // negative numbers are not aligned in the next run
                }
              end_offset = 0;
              start_offset = 0;
            }
        }
      // submit to cuda module an start the alignment


/*
      for (int i=0; i<counter;i++){

          printf("%d\t%d\n",i,read_offsets[i]);

      }
*/

      if (counter > 0){
      cuda_needleman_wunsch (match, mismatch, genome_array, read_array, genome_offsets, read_offsets,
                           edit_offsets, counter, id_array, pos_array, direction_array, error,
                           genome_array_size, counter*read_length_value, edit_matrix_size, do_reverse);
    
    }

          free (read_array);
    free (genome_array);
    free (direction_array);
    free (id_array);
    free (pos_array);
    free (read_offsets);
    free (genome_offsets);
    free (edit_offsets);
    *alignment_queue = 0;
    }
}
