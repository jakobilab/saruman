/*
 * hash.h
 *
 *  Created on: May 4, 2009
 *      Author: Tobias Jakobi
 */

#ifndef _HASH_H
#define	_HASH_H



// external variables 

extern int chunk_value;
extern int sw_match_value;
extern int sw_mismatch_value;
extern int genome_chunk_run;
extern int exact_tresh;
extern int total_bases_read;
extern unsigned int bases_read;
extern char * forward_output_control;
extern char * reverse_output_control;

extern char* genome_sequence;
extern unsigned short int read_length_value;
extern unsigned short int print_non_matching_reads;
extern char* reads;
extern struct read_name_hash *read_names;
extern unsigned short int maximal_errors_allowed;
extern unsigned short int* qgram_positions_in_read;
extern unsigned short int number_of_position_in_read;

// own variables

unsigned int pos_counter; /*!< counts number of start qgram start postions */
unsigned int hash_members; /*!< number of reads in the result array aka mapped reads */

short int forward_thread_check;
short int reverse_thread_check;


#include "util.h"
#include <uthash.h> // hash map

/**
 * \brief a structure needed to built a hash index of qgram start positions
 *
 * A structure needed to built a hash index of qgram start positions
 */
struct qgram_hash {
    char qgram[25]; /*!< \b KEY: unique qgram, hardcoded maximal qgram size */
    unsigned short int next_pos; /*!< next position in the following postions array */
    unsigned int* positions; /*!< holds the start positions of this qgram in the genome - dynamicly resized */
    UT_hash_handle hh; /*!< makes this structure hashable */
};

/**
 * \brief a structure needed to built a hash index of mapped reads
 *
 * A structure needed to built a hash index of mapped reads
 */
struct result_hash {
    unsigned int read_id; /*!< \b KEY: (unique) number of the read*/
    unsigned short int next_pos; /*!< next position in the following postions array */
    int* positions; /*!< holds the unique start positions of this read in the genome - dynamicly resized */
    char* direction_and_matching_qgrams; /*!< holds the direction of this hit together with the matching qgrams:<BR>
                       0:forward, start and end qgram match<BR>
                       1:forward, start qgram match<BR>
                       2:forward, end qgram match<BR>
                       3:forward, neither start or end match<BR>
                       4:reverse, start and end qgram match<BR>
                       5:reverse, start qgram match<BR>
                       6:reverse, end qgram match<BR>
                       7:reverse, neither start or end match
     */
    UT_hash_handle hh; /*!< makes this structure hashable */
};

/**
 * Checks if the given start position for the given qgram is a new one and stores
 * it in the corresponding struct of the #qgram_hash.
 * @param qgram the qgram
 * @param start_position_paramater start position to be added to the qgram
 * @return  0 if the qgram is new,
 *          1 if a read is updated with a new position
 */
unsigned short int add_qgram_to_index(char *qgram, int start_position_paramater);

/**
 * Checks if the given start position for the given read is a new one and stores
 * it in the corresponding struct of the #result_hash.
 * @param read_id unique number of the read
 * @param start_position start position to be added to the read
 * @param direction_and_qgrams holds the direction and the matching qgram(s)
 * @param results reference to the result hash of this thread
 * @param alignment_queue reference to the alignment queque of this thread
 * @return  0 if the read is new,
 *          1 if the start position already exists and
 *          2 if a read is updated with a new position
 */
unsigned short int
add_start_position (unsigned int read_id, int start_position, char direction_and_qgrams, struct result_hash **results, unsigned int * alignment_queue, int possible);

/**
 * Performs an self check of the #qgram_hash i.e. iterate over the hash and check
 * each key for existence
 * @return 0 if check is passed, 1 otherwise
 */
short int perform_hashtable_identity_check(void);

/**
 * Performs a lookup in the qgram hash table.
 * @param qgram qgram to look up
 * @return a pointer to the struct if found, NULL pointer else.
 */
struct qgram_hash *qgram_lookup(char* qgram);

/**
 * Tries to assign all needed qgrams for mapping to a each read
 * @param threadarg a struct holding all necessary information for the mapping and alignment phase
 * @return nothing
 */

void * qgram_to_read_mapping (char * read_array, unsigned short int readlength,   unsigned short int qgram_length,
    unsigned short int* qgram_start_positions, unsigned short int number_of_start_positions ,
    unsigned short int rest ,unsigned short int do_reverse, int chunksize, int read_number,
    int error);

/**
 * Prints some statistics about the #qgram_hash
 * @param qgram_size size of the used qgram
 * @return nothing
 */
void print_qgram_index_statistics(int qgram_size);

/**
 * Starts the SEASHORE mapping algorithm
 * @param chunksize number of reads to process in one run
 * @param run number of the actual run
 * @param number_of_start_positions number of qgram start positions per read
 * @param qgram_length length of one qgram
 * @param error maximal allowed error
 * @param rest rest of the read after beeing cut into X qgrams
 * @param do_reverse determines the direction of this run
 * @param data_array array holding the qgram to read mapping
 * @param results the result hash for the correspondig thread
 * @param alignment_queue holds the number of alignments in queue
 * @param data_array_offsets start positions for each qgram's start positions
 * @return nothing
 */
struct result_hash * run_seashore(int chunksize, int run, int number_of_start_positions, int qgram_length,
                                  int error, int rest, int do_reverse, unsigned int * data_array, struct result_hash *results,
                                  unsigned int * alignment_queue, unsigned int * data_array_offsets);

/**
 * Prepares and collects all necessary data for the CUDA alignment phase.
 * @param mismatch the mismatch costs
 * @param match the match costs
 * @param error maximal allowed error
 * @param results the result hash for the correspondig thread
 * @param alignment_queue holds the number of alignments in queue
 * @param do_reverse determines the direction of this run
 * @return nothing
 */
void prepare_cuda_alignments (int mismatch, int match, int error, struct result_hash *results,
                              unsigned int * alignment_queue, int do_reverse);


/**
 * Checks the number of alignments in queue and launches #prepare_cuda_alignments if the chunksize is reached
 * @param mismatch the mismatch costs
 * @param match the match costs
 * @param error maximal allowed error
 * @param results the result hash for the correspondig thread
 * @param alignment_queue holds the number of alignments in queue
 * @param do_reverse determines the direction of this run
 * @return nothing
 */
void cuda_queue_check (int mismatch, int match, int error, struct result_hash *results, unsigned int * alignment_queue, int do_reverse);

/**
 * Clears and frees the #result_hash
 * @param results reference to the hash to be cleared
 * @return nothing
 */
void clear_result_hash(struct result_hash *results);

/**
 * Clears and frees the #read_name_hash
 * @return nothing
 */
void clear_name_hash(void);

/**
 * Clears and frees the #qgram_hash
 * @return nothing
 */
void clear_qgram_index (void);


void *cuda_daemon(void *threadarg);


#endif	/* _HASH_H */

