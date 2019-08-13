/*
 * util.h
 *
 *  Created on: Apr 7, 2009
 *      Author: Tobias Jakobi
 */

#ifndef UTIL_H_
#define UTIL_H_

#define MAX_READ_DESCRIPTION_LENGTH 60

// parameter declarations:


const char* genome_file_value; /*!< holds the path to the reference genome in fasta format */
const char* read_file_value; /*!< holds the path to the file in fasta format containing the reads*/
float error_rate_value; /*!< holds the error rate of the reads, e.g. 0.083 -> 8.3% */
int chunk_value; /*!< holds the number of reads to be processed in one run */
int sw_match_value; /*!< holds the costs for a match used by the CUDA Smith-Waterman module*/
int sw_mismatch_value; /*!< holds the costs for a match used by the CUDA Smith-Waterman module*/
int read_number_value; /*!< holds the number of total reads beeing mapped by the program*/
unsigned int genome_chunk_value; /*!< shows the version number */
unsigned short int use_gpu_value; /*!< selects the gpu*/
unsigned short int read_length_value; /*!< holds the length of each read in bp*/
unsigned short int be_verbose_value; /*!< generates more detailed output */
unsigned short int use_fastq; /*!< use fastq format */
unsigned short int print_non_matching_reads; /*!< print non matching reads */
unsigned short int help_value; /*!< shows the help*/
unsigned short int version_value; /*!< shows the version number */
unsigned long genome_file_offset;
unsigned long int line_offset;
unsigned int number_of_reads_read;

struct arg_dbl *error_rate_parameter; /*!< holds the error rate of the reads, e.g. 0.083 -> 8.3% */
struct arg_int *chunk_parameter; /*!< holds the number of reads to be processed in one run */
struct arg_int *sw_match_parameter; /*!< holds the costs for a match used by the CUDA Smith-Waterman module*/
struct arg_int *sw_mismatch_parameter; /*!< holds the costs for a match used by the CUDA Smith-Waterman module*/
struct arg_int *read_number_parameter; /*!< holds the number of total reads beeing mapped by the program*/
struct arg_int *read_length_parameter; /*!< holds the length of each read in bp*/
struct arg_file *genome_file_parameter; /*!< holds the path to the reference genome in fasta format */
struct arg_file *read_file_parameter; /*!< holds the path to the file in fasta format containing the reads*/
struct arg_int *genome_chunk_parameter; /*!< enables multithreading */
struct arg_lit *be_verbose_parameter; /*!< generates more detailed output */
struct arg_lit *help_parameter; /*!< shows the help*/
struct arg_lit *version_parameter; /*!< shows the version number */
struct arg_lit *fastq_parameter;  /*!< use fastq format for reads*/
struct arg_lit *cuda_parameter;  /*!< use fastq format for reads*/
struct arg_lit *print_non_matching_reads_parameter;  /*!< print out non matching reads*/
struct arg_int *use_gpu_parameter; /*!< selects the gpu */
struct arg_end *end; /*!< needed by the argtable library*/

#include <uthash.h> // hash map

struct read_name_hash *read_names; /*!< holds the results aka start positions in genome */

/**
 * \brief a structure needed to built a hash index which maps a read id to a human readable name
 *
 * A structure needed to built a hash index which maps a read id to a human readable name
 */
struct read_name_hash {
    int read_number; /*!< \b KEY: (unique) number of the read*/
    char name[MAX_READ_DESCRIPTION_LENGTH]; /*!< holds the name of this read */
    UT_hash_handle hh; /*!< makes this structure hashable */
};


void reverse_complement(char* sequence, int sequence_length);

/**
 * Checks the user input. Tests if all mandatory values are given.
 * @param argc number of CLI arguments
 * @param argv char array of the CLI arguments
 * @return  0 program starrted but directly terminated (show help, version etc)<BR>
 *          -1 missing mandaroty value<BR>
 *         2 all mandatory values ok, start mapping programm
 */
signed short int check_user_input(int argc, char * argv []);

/**
 * Reads and parses the given read file.
 * @param filename path to the file to read
 * @param max_entries number of reads to parse
 * @param readlength length of the reads
 * @return a char array containing all reads in order
 */
char* read_multiple_fasta_file(const char * filename, int max_entries, int readlength);

/**
 * Reads and parses the given genome file.
 * @param filename path to the file to read
 * @return a char array containing the genome
 */
char* read_genome_fasta_file(const char * filename);

/**
 * Constructs the main qgram index (#qgram_hash) holding all start positions for each qgram in the genome
 * @param genome_sequence the genome sequence returned by #read_genome_fasta_file
 * @param qgram_size the desired qgram size (computed automaticly as
 * <CODE>short int qgram_length = ceil(read_length / (maximal_errors_allowed+1))-1</CODE>)
 * @return nothing
 */
void construct_qgram_index(char * genome_sequence, int qgram_size, int run);


/**
 * Reverses a given string
 * @param sequence char sequence
 * @param sequence_length length of the sequence to be reversed
 * @return a char array containing the genome
 */
void reverse (char* sequence, int sequence_length);


#endif /* UTIL_H_ */
