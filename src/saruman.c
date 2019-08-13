/*
 * saruman.c
 *
 *  Created on: Apr 7, 2009
 *      Author: Tobias Jakobi
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>


#include "hash.h"       // hash stuff
#include "argtable2.h"	// command line parsing
#include "util.h"	// utility functions
#include "saruman.h"	// main functions
#include "cuda-functions.h" // cuda stuff

#define NUM_THREADS     3

struct thread_data thread_data_array[NUM_THREADS];

int
main(int argc, char * argv[]) {
    clock_t start_clock = clock();

    pthread_t threads[NUM_THREADS];

    // check input
    int check = check_user_input(argc, argv);

    // everthing ok, start programm
    if (check == 2) {
        line_offset = 0;

        // set maximal allowed errors
        maximal_errors_allowed = floor(read_length_value * error_rate_value);


        // set qgram length for this run
        qgram_length = ceil((float) read_length_value / (maximal_errors_allowed + 1)) - 1;


        // rest of the modulo operartion
        rest = read_length_value % qgram_length;

        // create index of lookup position


        // compute number of positions in array
        number_of_position_in_read = (int) (read_length_value / qgram_length);

        if (rest != 0) {
            number_of_position_in_read = number_of_position_in_read * 2;
        }

        //allocate array:

        qgram_positions_in_read = calloc(number_of_position_in_read, sizeof (int));


        unsigned int buffer_character_position;

        for (buffer_character_position = 0; buffer_character_position < (int) (read_length_value / qgram_length); buffer_character_position++) {
            qgram_positions_in_read[buffer_character_position] = buffer_character_position*qgram_length;
        }

        for (unsigned short int z = 1; z <= (unsigned short int) (read_length_value / qgram_length); z++, buffer_character_position++) {
            qgram_positions_in_read[buffer_character_position] = read_length_value - (z * qgram_length);
        }

        alignment_time_forward = 0;
        alignment_time_reverse = 0;

        unsigned long int cuda_memory = get_cuda_infos();

        cuda_memory_per_thread = cuda_memory*0.95 / 2;

        //cuda_memory = 0.9*cuda_memory;

        unsigned int system_set_cuda_chunksize = 0;

        unsigned int bytes_needed_per_read_on_GPU = (read_length_value + 2 * maximal_errors_allowed + 1)*(read_length_value + 1) * sizeof (short int) + // edit matrix for one alignment
                sizeof (int) * 4 + // auxilary data
                sizeof (char) + // structures
                (read_length_value + 3 * maximal_errors_allowed)*3 * sizeof (char) +// maximal storage needed for the alignment (read & genome) + storage needed for the genome sequence
                read_length_value * sizeof (char); // one read


        system_set_cuda_chunksize = cuda_memory / (bytes_needed_per_read_on_GPU + 100) / 2;

        if (chunk_value == 0) {

            chunk_value = system_set_cuda_chunksize;

        }


        fprintf(stderr, "\n");
        fprintf(stderr, "running saruman with the following parameters:\n");
        fprintf(stderr, "===================================================\n");
        fprintf(stderr, "alignment match cost: %d \n", sw_match_value);
        fprintf(stderr, "alignment mismatch cost: %d \n", sw_mismatch_value);
        fprintf(stderr, "qgram positions: %d\n", number_of_position_in_read);
        fprintf(stderr, "qgram start positions @ ");


        for (buffer_character_position = 0; buffer_character_position < number_of_position_in_read; buffer_character_position++) {
            fprintf(stderr, "%d ", qgram_positions_in_read[buffer_character_position]);
        }

        fprintf(stderr, "\n");

        exact_tresh = floor((float) read_length_value / (float) qgram_length) - maximal_errors_allowed + 1;

        fprintf(stderr, "read length: %d bp\n", read_length_value);
        fprintf(stderr, "read error rate: %f\n", error_rate_value);
        fprintf(stderr, "parallel alignments per direction on GPU: %d (%d bytes / alignment)\n", chunk_value, bytes_needed_per_read_on_GPU);
        fprintf(stderr, "filter threshold value: %d qgrams needed\n", exact_tresh);
        fprintf(stderr, "maximal errors allowed: %d\n", maximal_errors_allowed);
        fprintf(stderr, "qgram length: %d\n", qgram_length);
        fprintf(stderr, "===================================================\n");

        //get_cuda_infos ();

        fprintf(stderr, "\n-> acquiring data...\n");

        // read genome into string

        genome_sequence = NULL;

        //initialize variables

        FILE * file_pointer = NULL;
        unsigned long chars_to_read = genome_chunk_value + read_length_value;
        char * buffer;
        size_t buffer_chars_read = 0;
        total_bases_read = 0;

        file_pointer = fopen(genome_file_value, "r");
        if (file_pointer == NULL) {
            fprintf(stderr, "Could not load genome file '%s'. Check if the given file exists.\n", genome_file_value);
            exit(1);
        }

        char first_run = 1; // 1 = first run, filter fasta header, 0 = no header to filter
        bases_read = 0; // genome bases read

        // initiating
        line_offset = 0; // line offset for read input
        hash_members = 0; // number of hash members - number of unique qgrams
        pos_counter = 0; // number of (redundant) start positions in genome
        genome_chunk_run = 0;

        buffer_character_position = 0;

        buffer = calloc(chars_to_read, sizeof (char));
        if (buffer == NULL) {
            fputs("Could not allocate memory for buffer", stderr);
            exit(2);
        }

        genome_sequence = calloc(chars_to_read + 1, sizeof (char));
        if (genome_sequence == NULL) {
            fputs("Could not allocate memory for genome sequence", stderr);
            exit(2);
        }

        long file_position = 0;
        int *taskids[NUM_THREADS];

        int thread_return_value;
        void * thread_end_status;
        unsigned int genome_sequence_position = 0;

        do // read genome until end of file
        {

            //fprintf(stderr,"recovering file position: %d.\n",file_position);

            // setting file position

           genome_sequence_position = 0; // reset base pair counter

            if (!first_run) {
                fseek(file_pointer, file_position, SEEK_SET); // recover file position from previous run

                if (buffer_character_position < chars_to_read) {
                    fprintf(stderr, "reusing last buffer from pos %d\n", buffer_character_position);

                    for (unsigned int recycle_buffer_position = buffer_character_position; genome_sequence_position < genome_chunk_value && recycle_buffer_position < buffer_chars_read; recycle_buffer_position++) {

                        if (buffer[recycle_buffer_position] != '\n') { // skip new lines
                            genome_sequence[genome_sequence_position] = toupper(buffer[recycle_buffer_position]); // convert to uppercase
                            genome_sequence_position++;
                        }

                    }

                }

            }

            // copy the file into the buffer:
            buffer_chars_read = fread(buffer, 1, chars_to_read, file_pointer); // read maximal bytes_to_read

            //buffer[chars_to_read] = '\0'; // properly end buffer string for strlen

            unsigned int sequence_start = 0;

            // if first run, filter fasta header
            if (first_run) {
                // filter out fasta header
                char expression[] = "\n"; // sequence starts after first newline
                sequence_start = strcspn(buffer, expression); // start correction for fasta header

            }
            // processing fasta contents

            buffer_character_position = sequence_start; // set start to first character after header

            for (buffer_character_position = sequence_start; genome_sequence_position < genome_chunk_value && buffer_character_position < buffer_chars_read ; buffer_character_position++) {

                if (buffer[buffer_character_position] != '\n') { // skip new lines
                    genome_sequence[genome_sequence_position] = toupper(buffer[buffer_character_position]); // convert to uppercase
                    genome_sequence_position++;
                }

            } // ok, we have read X chars minus the chars of the header

            first_run = 0; // first run done - header was removed

            size_t genome_sequence_result = 1;

            while (genome_sequence_position < genome_chunk_value && genome_sequence_result == 1) {

                genome_sequence_result = fread(buffer, 1, 1, file_pointer); // read 1 element of 1 byte from file pointer into the buffer
                
                if (buffer[0] != '\n') {
                    genome_sequence[genome_sequence_position] = toupper(buffer[0]); // convert to upper case
                    genome_sequence_position++;
                }
            }

            bases_read += genome_sequence_position;
            total_bases_read += genome_sequence_position;
            // this is the position the next round will start

            file_position = ftell(file_pointer);

            // getting overhang

            int overhang_result = 0;

            while (overhang_result < read_length_value && genome_sequence_result == 1) {
                genome_sequence_result = fread(buffer, 1, 1, file_pointer);
                if (buffer[0] != '\n') {

                    genome_sequence[genome_sequence_position] = toupper(buffer[0]);
                    genome_sequence_position++;
                    overhang_result++;
                }

            }

            genome_sequence[genome_sequence_position] = '\0';

            fprintf(stderr, "-> fetchted %d (%d total) bp from %s.\n", bases_read, total_bases_read ,genome_file_value);

            // take time for qgram index creation
            clock_t qgram_start = clock();
            construct_qgram_index(genome_sequence, qgram_length, genome_chunk_run);
            clock_t qgram_stop = clock();

            forward_thread_check = 0;
            reverse_thread_check = 0;

            print_qgram_index_statistics(qgram_length);


            fprintf(stderr, "-> init forward working thread ...\n");



            taskids[0] = (int *) malloc(sizeof (int));
            *taskids[0] = 0;

            taskids[1] = (int *) malloc(sizeof (int));
            *taskids[1] = 4;


            thread_return_value = pthread_create(&threads[1], NULL, cuda_daemon, (void *) taskids[0]);

            if (thread_return_value) {
                fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", thread_return_value);
                exit(-1);
            }

            fprintf(stderr, "-> init reverse working thread ...\n");
            thread_return_value = pthread_create(&threads[2], NULL, cuda_daemon, (void *) taskids[1]);

            if (thread_return_value) {
                fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", thread_return_value);
                exit(-1);
            }

            // read reads into array

            while (((line_offset == 0) || (line_offset % read_number_value) == 0) &&
                    (reads = read_multiple_fasta_file(read_file_value, read_number_value, read_length_value))) {
                if (reads != NULL) {
                    
                    // setting up forward thread

                    forward_output_control = (char *) calloc(number_of_reads_read, sizeof (char));
                    alignment_time[0] = 0;
                    forward_thread_check = 1;

                    // setting up reverse thread

                    reverse_output_control = (char *) calloc(number_of_reads_read, sizeof (char));
                    alignment_time[1] = 0;
                    reverse_thread_check = 1;


                    sleep(2);


                    if (print_non_matching_reads == 1) {

                        for (unsigned int var = 0; var < number_of_reads_read; ++var) {

                            if (forward_output_control[var] == 0) {

                                struct read_name_hash *hash_entry = NULL;
                                HASH_FIND_INT(read_names, &var, hash_entry);

                                printf("%s\t%s\t%s\t%s\t%.*s\t%s\t%s\n",
                                        hash_entry->name,
                                        "*",
                                        "*",
                                        ">>",
                                        read_length_value, &reads[var * read_length_value],
                                        "*",
                                        "*");

                            } else if (reverse_output_control[var] == 0) {

                                struct read_name_hash *hash_entry = NULL;
                                HASH_FIND_INT(read_names, &var, hash_entry);

                                printf("%s\t%s\t%s\t%s\t%.*s\t%s\t%s\n",
                                        hash_entry->name,
                                        "*",
                                        "*",
                                        "<<",
                                        read_length_value, &reads[var * read_length_value],
                                        "*",
                                        "*");
                            }

                        }
                    }





                    while (forward_thread_check != 3 || reverse_thread_check != 3){
                        //fprintf(stderr, "-> Waiting for all threads to finish...\n");
                        sleep(2);
                    }

                    free(reads);
                    clear_name_hash();
                    free(reverse_output_control);
                    free(forward_output_control);
                }
            }

            forward_thread_check = -1;
            reverse_thread_check = -1;

            pthread_join(threads[1], &thread_end_status);
            fprintf(stderr, "-> forward thread finished with status %ld - spent %.0f ms for alignments\n", (long) thread_end_status, alignment_time[0]);


            pthread_join(threads[2], &thread_end_status);
            fprintf(stderr, "-> reverse thread finished with status %ld - spent %.0f ms for alignments\n", (long) thread_end_status, alignment_time[1]);

                     alignment_time_forward += alignment_time[0];
                     alignment_time_reverse += alignment_time[1];

            free(taskids[0]);
            free(taskids[1]);


            fprintf(stderr, "Mapping process complete.\nTotal execution time was %4.2f seconds (q-gram index creation: %4.2f seconds).\n",
                    ((double) clock() - start_clock) / CLOCKS_PER_SEC, ((double) qgram_stop - qgram_start) / CLOCKS_PER_SEC);

            fprintf(stderr, "Found %d (%d forward, %d reverse) candidate hits.\nGot %d (%d forward, %d reverse) suitables alignments.\nShutting down.\n", forward_alignment_candidates + reverse_alignment_candidates,forward_alignment_candidates, reverse_alignment_candidates, reverse_alignments_done + forward_alignments_done,reverse_alignments_done,forward_alignments_done);


            clear_qgram_index(); // clear the hash index
            pos_counter = 0;
            bases_read = 0;
            line_offset = 0;
            hash_members = 0;
            genome_chunk_run++;
        } while (genome_sequence_position >= genome_chunk_value); // while ends





        // terminate
        fclose(file_pointer);
        free(buffer);
        free(genome_sequence);
        free(qgram_positions_in_read);
        //free (read_file_value);


        //free_cuda_memory();

    }

    return (0);
}


