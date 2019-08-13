/*
 * saruman.h
 *
 *  Created on: Apr 7, 2009
 *      Author: Tobias Jakobi
 */

#ifndef SARUMAN_H_
#define SARUMAN_H_

// external variables

extern float alignment_time[2];
extern float error_rate_value;
extern int chunk_value;
extern int sw_match_value;
extern int sw_mismatch_value;
extern int read_number_value;
extern unsigned short int read_length_value;
extern unsigned short int be_verbose_value;
extern unsigned short int print_non_matching_reads;
extern const char* genome_file_value;
extern const char* read_file_value;
extern unsigned long int line_offset;
extern struct read_name_hash *read_names;
extern unsigned int number_of_reads_read;


/**
 *\brief work around struct for passing multiple arguments with pthread
 *
 * A structure needed to pass multiple arguments to the #qgram_to_read_mapping function called through pthread_start
 */
struct thread_data {
    char * read_array; /*!< array holding the reads for this thread*/
    unsigned short int readlength; /*!< holds the read length*/
    unsigned short int* qgram_start_positions; /*!< holds the qgram start positions in the read*/
    unsigned short int number_of_start_positions; /*!< holds the number of qgrams per read */
    unsigned short int qgram_length; /*!< holds the qgram length*/
    unsigned short int rest; /*!< holds the lenth of the overhang of the qgrams and the reads*/
    unsigned short int do_reverse; /*!< holds the direction*/
    int read_number;  /*!< holds number of reads to process per run*/
    int error; /*!< holds the maximal allowed error*/
    int chunksize; /*!< holds the chunksize*/
};

char* genome_sequence;
char* reads;
unsigned short int maximal_errors_allowed;
unsigned short int qgram_length;
int exact_tresh;
int genome_chunk_run;
int total_bases_read;
unsigned int bases_read;
char * forward_output_control;
char * reverse_output_control;
unsigned short int* qgram_positions_in_read;
unsigned short int number_of_position_in_read;
unsigned short int rest;
unsigned long int cuda_memory_per_thread;

unsigned int forward_alignments_done;
unsigned int reverse_alignments_done;

unsigned int forward_alignment_candidates;
unsigned int reverse_alignment_candidates;

float alignment_time_forward;
float alignment_time_reverse;

extern char* main_alloc_ptr_for;
extern char* main_alloc_ptr_rev;


/**
 * Main function. Starts the program and does some housekeeping before calling the mapping functions
 * @param argc number of CLI arguments
 * @param argv char array of the CLI arguments
 */
int main(int argc, char * argv[]);

#endif /* SARUMAN_H_ */
