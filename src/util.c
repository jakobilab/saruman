/*
 * util.c
 *
 *  Created on: Apr 7, 2009
 *      Author: Tobias Jakobi
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>

#include "argtable2.h"		// command line parsing
#include "util.h"		// utility functions
#include "saruman.h"            // main file
#include "hash.h"               // hash table stuff

#define VERSION "1.0.7 - compilation date: 23.08.2011"
#define SWAP_CHAR( x, y ) {char c; c = x; x = y; y = c;}

void
reverse (char* sequence, int sequence_length)
{
  int i, j;
  sequence_length--;

  for (i = 0, j = sequence_length; i < j; i++, j--)
    {
      SWAP_CHAR (sequence[i], sequence[j]);
    }

}

void
reverse_complement (char* sequence, int sequence_length)
{
  for (int i = 0; i < sequence_length; i++)
    {
      switch (sequence[i])
        {
        case 'A':
          sequence[i] = 'T';
          break;
        case 'C':
          sequence[i] = 'G';
          break;
        case 'G':
          sequence[i] = 'C';
          break;
        case 'T':
          sequence[i] = 'A';
          break;
        }
    }
  reverse (sequence, sequence_length);
}

signed short int
check_user_input (int argc, char * argv[])
{

  /* Define the allowable command line options, collecting them in argtable[] */

  // parameter has to occure exactly one time
  sw_match_parameter = arg_int0 ("m", "match", "<n>", "cost for a match used by the CUDA Needleman-Wunsch module (e.g. 0)");

  // parameter has to occure exactly one time
  sw_mismatch_parameter = arg_int0 ("i", "mismatch", "<n>", "cost for a mismatch used by the CUDA Needleman-Wunsch module (e.g. 1)");

  // parameter has to occure exactly one time
  error_rate_parameter = arg_dbl0 ("e", "error", "<n>", "maximal error rate tolerated in the reads (e.g. 0.06 ~ 6%)");

  // parameter has to occure exactly one time
  read_length_parameter = arg_int0 ("l", "length", "<n>", "read length in bp (e.g. 36)");

  // parameter has to occure exactly one time
  read_number_parameter = arg_int0 ("c", "count", "<n>", "number of reads to gather from input file in one chunk (e.g. 2000000)");

  // parameter has to occure exactly one time
  chunk_parameter = arg_int0 ("s", "chunksize", "<n>", "number of reads to align in one chunk (e.g. 400000)");

  // parameter has to occure exactly one time
  genome_file_parameter = arg_file0 ("g", "genomefile", "<file>", "the reference genome in fasta format");

  // parameter has to occure exactly one time
  read_file_parameter = arg_file0 ("r", "readfile", "<file>", "multiple fasta file containing reads");

  // parameter has to occure exactly one time
  genome_chunk_parameter = arg_int0 ("G", "genomechunksize" , "<n>" , "size of the genome chunks gathered in one pass (e.g. 5000000)");

  //standard help & version stuff
  be_verbose_parameter = arg_lit0 ("v", "verbose", "be verbose");

  //standard help & version stuff
  help_parameter = arg_lit0 ("h", "help", "print this help and exit");

  // print out program version
  version_parameter = arg_lit0 ("V", "version", "print version information and exit");

  // use fastq format instead of fasta
  fastq_parameter = arg_lit0 ("q", "fastq", "use fastq input format for reads instead of fasta");

  // use specified GPU
  use_gpu_parameter = arg_int0 ("u", "gpu" , "<n>" , "manually select the GPU to use starting from 0");

  // print non matching reads
  print_non_matching_reads_parameter = arg_lit0 ("p", "print", "print out non matching reads");

  end = arg_end (21);


  void* argtable[] = {error_rate_parameter, read_file_parameter, be_verbose_parameter,
    genome_file_parameter, help_parameter, version_parameter,
    read_length_parameter, read_number_parameter, chunk_parameter,
    sw_match_parameter, sw_mismatch_parameter, genome_chunk_parameter,
    fastq_parameter, print_non_matching_reads_parameter,use_gpu_parameter,end
  };

  const char* progname = "saruman";
  int exitcode = 0;
  int nerrors = 0;

  /* verify the argtable[] entries were allocated sucessfully */
  if (arg_nullcheck (argtable) != 0)
    {
      /* NULL entries were detected, some allocations must have failed */
      fprintf (stderr,"%s: insufficient memory\n", progname);
      exitcode = 1;
    }

  /* Parse the command line as defined by argtable[] */
  nerrors = arg_parse (argc, argv, argtable);


    // saving variables

    error_rate_value    = error_rate_parameter->dval[0];
    chunk_value         = chunk_parameter->ival[0];
    sw_match_value      = sw_match_parameter->ival[0];
    sw_mismatch_value   = sw_mismatch_parameter->ival[0];
    read_number_value   = read_number_parameter->ival[0];
    read_length_value   = read_length_parameter->ival[0];
    be_verbose_value    = be_verbose_parameter->count;
    read_file_value     = read_file_parameter->filename[0];
    genome_file_value   = genome_file_parameter->filename[0];
    genome_chunk_value  = genome_chunk_parameter->ival[0];
    use_fastq           = fastq_parameter->count;
    use_gpu_value       = use_gpu_parameter->ival[0];
    print_non_matching_reads  = print_non_matching_reads_parameter->count;

  /* special case: '--help' takes precedence over error reporting */
  if (help_parameter->count > 0)
    {
      fprintf (stderr,"'%s'- a GPU-based short read mapper - version %s\n\n", progname, VERSION);
      fprintf (stderr,"Usage: %s", progname);
      arg_print_syntax (stdout, argtable, "\n");
      arg_print_glossary (stdout, argtable, " %-25s %s\n");
      fprintf (stderr,"\nFor more information contact tjakobi@CeBiTec.Uni-Bielefeld.DE\n");
      fprintf (stderr,"\nor try http://www.cebitec.uni-bielefeld.de/brf/saruman/saruman.html\n\n");
      exit (-1);
    }

  /* special case: '--version' takes precedence error reporting */
  if (version_parameter->count > 0)
    {
      fprintf (stderr,"'%s' - version %s\n", progname, VERSION);
      //fprintf (stderr,"if you use this tool, please cite XXX\n");
      fprintf (stderr,"developed by Tobias Jakobi and Jochen Blom, [tjakobi|jblom]@CeBiTec.Uni-Bielefeld.DE\n");
      fprintf (stderr,"2010,2011 © CeBiTec Bielefeld University\n");

      exitcode = 0;
      exit (0);
    }

  /* If the parser returned any errors then display them and exit */
  if (nerrors > 0)
    {
      /* Display the error details contained in the arg_end struct.*/
      arg_print_errors (stdout, end, progname);
      fprintf (stderr,"Try '%s --help' for more information.\n", progname);
      exit (-1);
    }


  // checking if 2 files are given

    if (chunk_parameter->count == 0){
        chunk_value = 0;
    }

    if (use_gpu_parameter->count == 0){
        use_gpu_value = 5000;
    }

  if (!(read_file_parameter->count == 1 && genome_file_parameter ->count == 1))
    {
      fprintf (stderr,"A value for the genome reference genome is missing!\n\n");
      fprintf (stderr,"Usage: %s", progname);
      arg_print_syntax (stdout, argtable, "\n");
      arg_print_glossary (stdout, argtable, " %-25s %s\n");
      fprintf (stderr,"\nFor more information contact tjakobi@CeBiTec.Uni-Bielefeld.DE\n");
      fprintf (stderr,"\nor try http://www.cebitec.uni-bielefeld.de/brf/saruman/saruman.html\n\n");
      exit (-1);
    }

  // checking if length and number are given


  if (!(read_number_parameter->count == 1 && read_length_parameter ->count == 1))
    {
      fprintf (stderr,"You have to pass the read length and read number!\n\n");
      fprintf (stderr,"Usage: %s", progname);
      arg_print_syntax (stdout, argtable, "\n");
      arg_print_glossary (stdout, argtable, " %-25s %s\n");
      fprintf (stderr,"\nFor more information contact tjakobi@CeBiTec.Uni-Bielefeld.DE\n");
      fprintf (stderr,"\nor try http://www.cebitec.uni-bielefeld.de/brf/saruman/saruman.html\n\n");
      exit (-1);
    }



  if (!(error_rate_parameter->count == 1))
    {
      fprintf (stderr,"You have to pass an error rate!\n\n");
      fprintf (stderr,"Usage: %s", progname);
      arg_print_syntax (stdout, argtable, "\n");
      arg_print_glossary (stdout, argtable, " %-25s %s\n");
      fprintf (stderr,"\nFor more information contact tjakobi@CeBiTec.Uni-Bielefeld.DE\n");
      fprintf (stderr,"\nor try http://www.cebitec.uni-bielefeld.de/brf/saruman/saruman.html\n\n");
      exit (-1);
    }

/*
  if (!(chunk_parameter->count == 1))
    {
      //fprintf (stderr,"You have to pass a chunk size!\n");
      //exit (-1);
      chunk_value = 0;
    }
*/

   if (!(genome_chunk_parameter->count == 1) || (genome_chunk_parameter->ival[0]<  0))
    {
      fprintf (stderr,"You have to pass a genome chunk size!\n\n");
      fprintf (stderr,"Usage: %s", progname);
      arg_print_syntax (stdout, argtable, "\n");
      arg_print_glossary (stdout, argtable, " %-25s %s\n");
      fprintf (stderr,"\nFor more information contact tjakobi@CeBiTec.Uni-Bielefeld.DE\n");
      fprintf (stderr,"\nor try http://www.cebitec.uni-bielefeld.de/brf/saruman/saruman.html\n\n");
      exit (-1);
    }

    if (!(sw_match_parameter->count == 1))
    {
      fprintf (stderr,"You have to pass a valid match cost!\n\n");
      fprintf (stderr,"Usage: %s", progname);
      arg_print_syntax (stdout, argtable, "\n");
      arg_print_glossary (stdout, argtable, " %-25s %s\n");
      fprintf (stderr,"\nFor more information contact tjakobi@CeBiTec.Uni-Bielefeld.DE\n");
      fprintf (stderr,"\nor try http://www.cebitec.uni-bielefeld.de/brf/saruman/saruman.html\n\n");
      exit (-1);
    }

      if (!(sw_mismatch_parameter->count == 1))
    {
      fprintf (stderr,"You have to pass a valid mismatch cost!\n\n");
      fprintf (stderr,"Usage: %s", progname);
      arg_print_syntax (stdout, argtable, "\n");
      arg_print_glossary (stdout, argtable, " %-25s %s\n");
      fprintf (stderr,"\nFor more information contact tjakobi@CeBiTec.Uni-Bielefeld.DE\n");
      fprintf (stderr,"\nor try http://www.cebitec.uni-bielefeld.de/brf/saruman/saruman.html\n\n");
      exit (-1);
    }


  //everything is fine, start over
  if (read_file_parameter->count == 1 && genome_file_parameter->count == 1)
    {
      exitcode = 2;
    }

    // deallocate each non-null entry in argtable[]
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

  return exitcode;



}


char*
read_multiple_fasta_file (const char * filename, int max_entries, int readlength)
{

  unsigned int skipped_reads = 0;

  FILE *file_pointer;
  int read_counter = 0;
  int internal_counter=0;
  read_names = NULL; /*!< holds the results aka start positions in genome */

  file_pointer = fopen (filename, "r");
  if (file_pointer == NULL)
    {
      fprintf (stderr,"Could not load read file '%s'. Check if the given file exists.\n", filename);
      exit (1);
    }

  fprintf (stderr,"-> fetching reads, please wait... ");

  char* container;

  //length of input * size of a char on this system
  container = calloc (max_entries * (readlength + 1), sizeof (char));

  // no memory available for 1st dimension
  if (container == NULL)
    {
      fprintf (stderr, "not enough memory for fasta file");
      return NULL;
    }

  long int line_offset_tmp = line_offset+max_entries;
  long int line_offset_tmp2 = line_offset;

  //fprintf (stderr,"%d -> %d\n",line_offset,line_offset_tmp);
  //fprintf (stderr,"%d -> %d\n",internal_counter,line_offset_tmp);

    if (use_fastq)
    {
      fprintf (stderr, "(using fastq format)\n");
    }
  else
    {
      fprintf (stderr, "(using fasta format)\n");
    }


    char seqname[MAX_READ_DESCRIPTION_LENGTH];
    char seq[readlength+2];
    char qualname[MAX_READ_DESCRIPTION_LENGTH];
    char qual[readlength+2];


    while (fgets(seqname, sizeof(seqname), file_pointer) != NULL  && (internal_counter < line_offset_tmp) ) {

           fgets(seq, sizeof(seq), file_pointer);

               if (use_fastq)
    {
           fgets(qualname, sizeof(qualname), file_pointer);
           fgets(qual, sizeof(qual), file_pointer);
           qualname[strlen(qualname)-1] = '\0';
           qual[readlength] = '\0';
    }


/*
             fwrite(seqname, 1, strlen(seqname), stderr);
             fwrite(seq, 1, strlen(seq), stderr);
             //fwrite(qualname, 1, strlen(qualname), stderr);
             //fwrite(qual, 1, strlen(qual), stderr);

*/

           int length = strlen (seqname);
           seq[readlength] = '\0';
           seqname[length-1] = '\0';
           length--;
/*
           fprintf (stderr, "desc: %s %d\n",seqname,strlen(seqname));
           fprintf (stderr, "seq: %s\n",seq);
           fprintf (stderr, "qualname: %s %d\n",qualname,strlen(qualname));
           fprintf (stderr, "quality: %s\n",qual);
*/
      if (internal_counter >= line_offset_tmp2)
        {


          unsigned short int n_counter = 0;


          for (int i = 0; i < readlength; i++)
            {
              seq[i] = toupper (seq[i]);

              if (seq[i] == 'N')
                {
                  n_counter++;
                }
            }

          if (strlen(seq) < (unsigned) readlength){

                  // fprintf (stderr,"-> skipped malformed read %s %d\n", seq,strlen(seq));
                  skipped_reads++;
                  line_offset++;
          }


          else if (n_counter <= maximal_errors_allowed)
            {



              struct read_name_hash *hash_entry = NULL;
              // allocate memory for new entry in read_names hash
              hash_entry = (struct read_name_hash*) malloc (sizeof (struct read_name_hash));

              // copy to read container
              strncpy (&container[read_counter * readlength], seq, readlength);
              //strcpy (&container[read_counter * readlength], seq);

              // add id as key to the hash
              hash_entry->read_number = read_counter;

              // add value to the hash

              if (use_fastq){
              strncpy (hash_entry->name, &seqname[1],length);
              } else {
              strncpy (hash_entry->name, &seqname[1],length);
              }
              read_counter++;
              line_offset++;
              // commit changes to the hash
              HASH_ADD_INT (read_names, read_number, hash_entry);

            } else {
                  skipped_reads++;
                  line_offset++;
              }
        }
      internal_counter++;
      //fprintf (stderr,"counter %d\n",internal_counter);


    }

  // close file handle
  fclose (file_pointer);

  if (read_counter==0){
     container = read_multiple_fasta_file (filename, max_entries, readlength);
    }
    number_of_reads_read = read_counter;
    fprintf (stderr,"-> skipped %d malformed or quality filtered reads( >= %d Ns)\n", skipped_reads,maximal_errors_allowed);
    fprintf (stderr,"-> fetchted %d reads with %d bp starting at %ld from file %s\n", read_counter, readlength, line_offset_tmp2, filename);

  return container;

}

void
construct_qgram_index (char * genome_sequence, int qgram_size, int run)
{

  //clean_qgram_hash();

  int genome_length = strlen (genome_sequence);
  int max = genome_length - qgram_size;
  int percent = -1;
  char display[105];

  fprintf (stderr,"constructing qgram index, please wait...\n");

  if (genome_length <= qgram_size){
      fprintf (stderr,"last genome chunk shorter than qgram - abborting\n");
      exit(-1);
  }

  for (int genome_position = 0; genome_position <= max; genome_position++)
    {
      char* qgram;
      qgram = calloc ((qgram_size + 1), sizeof (char));
      //fprintf(stderr,"qram pos in real gneome would be: %d\n",genome_position+run*genome_chunk_value);
      strncpy (qgram, &genome_sequence[genome_position], qgram_size);
      qgram[qgram_size] = '\0';
      add_qgram_to_index (qgram, genome_position+run*genome_chunk_value);
      free (qgram);

      if ((genome_position * 100 / max)/2 > percent){

      percent ++;
      display[percent]='=';
      display[percent+1]='>';
      display[percent+2]='\0';

      fprintf (stderr,"%*d %% [%-*s]\r",2,genome_position * 100 / max,51,display);

      }
      }

  fprintf (stderr,"\ndone.\n");

  //print_qgram_index_statistics (qgram_size);
}



