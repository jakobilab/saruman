/*
 * cuda_functions.h
 *
 *  Created on: Apr 7, 2009
 *      Author: Tobias Jakobi
 */

#ifndef CUDA_FUNCTIONS_H_
#define CUDA_FUNCTIONS_H_

#include "saruman.h"



//external variables
extern struct result_hash *results;
extern struct read_name_hash *read_names;
extern unsigned short int be_verbose_value;
extern unsigned short int qgram_length;
extern unsigned short int use_gpu_value;

extern unsigned int forward_alignments_done;
extern unsigned int reverse_alignments_done;

extern unsigned int forward_alignment_candidates;
extern unsigned int reverse_alignment_candidates;

extern unsigned long int cuda_memory_per_thread;


char* main_alloc_ptr_for;
char* main_alloc_ptr_rev;


/**
 * Wrapper functions which allows to call the CUDA implemented Smith-Waterman
 * algorithm from other C files.
 * @param match costs for a match
 * @param mismatch costs for a mismatch
 * @param genome_sequences_on_host array of genome sequences
 * @param read_sequences_on_host array of read sequences
 * @param genome_offsets stop positions of each genome sequence in the genome array
 * @param read_offsets stop positions of each read sequence in the read array
 * @param edit_offsets stop positions of each edit matrix in the edit matrix array
 * @param genome_pos_array start positions per read in the reference genome
 * @param direction_array holds the direction and matching qgram(s) of each hit
 * @param error allowed error
 * @param total_genome_array_size total size of the genome array
 * @param total_read_sequence_array_size total size of the read array *
 * @param total_edit_matrix_size total size of the edit matrix array
 * @param do_reverse forward or reverse direction
 * @param read_id_array hods the id of each read in the set
 * @param sample_size size of the sample arrays
 */
extern void
cuda_needleman_wunsch(int match, int mismatch,
        char* genome_sequences_on_host, char* read_sequences_on_host,
        int * genome_offsets, int * read_offsets, int * edit_offsets,
        int sample_size, int* read_id_array, int* genome_pos_array,
        char* direction_array, int error, size_t total_genome_array_size,
        size_t total_read_sequence_array_size, size_t total_edit_matrix_size,
        int do_reverse
        );

extern
void
cuda_allocate_memory(int do_reverse);

extern
void
cuda_free_memory(int do_reverse);

/**
 * Reverses a given char array
 * @return amount of available CUDA memory in bytes
 */
extern unsigned long int
get_cuda_infos(void);


/**
 * Reverses a given char array
 * @param sequence the array/sequence to be reversed
 * @param k length of the array
 */
//__device__ void
//reverse_sequence(char* sequence, int k);

/**
 * Minimizes over four integer variables
 * @param a 1st integer
 * @param b 2nd integer
 * @param c 3rd integer
 * @param d 4th integer
 * @return maximum of the 4 values
 */
//__device__ short int
//minimize(short int a, short int b, short int c, short int d);


/**
 * CUDA function which constructs the Needleman-Wunsch alignments for the given sequences.
 * @param match costs for a match
 * @param mismatch costs for a mismatch
 * @param genome_sequences array of genome sequences
 * @param read_sequences array of read sequences
 * @param genome_alignments array holding the genomic alignment parts
 * @param read_alignments array holding the read part of the alignment
 * @param edit_matrices_on_device_gpu of edit matrices located on the device
 * @param max_position position of the maximal score (deprecated because result is always in the left bottom corner of the matrix)
 * @param edit_matrix_start start offset for the edit matrix of this alignment
 * @param read_sequence_length length of this read
 * @param read_start start offset of the read sequences of this alignment
 * @param genome_start start offset of the genomic sequences for this alignment
 * @param genome_sequence_length length of the genomic sample cut out
 * @param error allowed error
 * @param end_offset set, if this alignment got an overhang at the tail part
 * @param start_offset set, if this alignment got an overhand at the head part
 * @param threadID thread id of this thread

 */
//__device__ void
//get_alignment(int match, int mismatch,
//        char *read_sequences, char* genome_sequences,
//        char *read_alignments, char *genome_alignments,
//        char* edit_matrices_on_device_gpu,
//        int max_position, int edit_matrix_start,
//        int read_sequence_length, int read_start,
//        int genome_start, int genome_sequence_length,
//        int error, int end_offset, int start_offset, int threadID
//        );

/**
 * Actual CUDA function which computes the Needleman-Wunsch scores for the given sequenes.
 * @param match costs for a match
 * @param mismatch costs for a mismatch
 * @param genome_sequences array of genome sequences
 * @param read_sequences array of read sequences
 * @param genome_offsets stop positions of each genome sequence in the genome array
 * @param read_offsets stop positions of each read sequence in the read array
 * @param edit_offsets stop positions of each edit matrix in the edit matrix array
 * @param genome_alignments array holding the genomic alignment parts
 * @param read_alignments array holding the read part of the alignment
 * @param direction_array holds the direction and matching qgram(s) of each hit
 * @param error allowed error
 * @param edit_matrices_on_device_gpu of edit matrices located on the device
 * @param device_score_array holds the score for each alignment
 * @param sample_size size of the sample arrays
 */
//__global__ void
//compute_score (int match, int mismatch,
 //              char *read_sequences, char* genome_sequences,
 //              int *read_offsets, int* genome_offsets,
 //              char* read_alignments, char* genome_alignments,
 //              int sample_size,
 //              char* edit_matrices_on_device_gpu,
 //              int* edit_offsets,
 //              int* device_score_array,
 //              char* direction_array,
 //              int error
 //              );




#endif /* CUDA_FUNCTIONS_H_ */
