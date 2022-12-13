#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<glib.h>
#include<regex.h>
#include<limits.h>
#include<stddef.h>
#include<unistd.h>
#include<ctype.h>
#include<stdbool.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<time.h>
#include<assert.h>

#ifndef READ_DATA_H_
#define READ_DATA_H_

#define MAX_CHROMOSOMES 23
#define TOTAL_MARKERS 500
#define BUFFER_SIZE 60000
#define MAX_HAPSEQS 18
#define MAX_INDIVIDUALS 11000
#define MAX_FILENAME 100
#undef DEBUG

typedef struct chromosome_all_info
{
  unsigned int chrom_id;
  double chrom_length;
  double chrom_recom_rate;
  unsigned int n_loci;
  double *markers;
  GHashTable* hashA;
  GHashTable* hashB;
  struct chromosome_all_info *next;
} chrom_data;

chrom_data* make_data_node(unsigned int index, double length, double recom_rate, unsigned int no_loci, double *loci);

void string_to_markers(char *string, unsigned int no_loci, double **markers);

void read_chrom(char *filename, GIOChannel *infile, unsigned int **chromosome_index, unsigned int **chromosome_no_loci, double **chromosome_length, double **chromosome_recom_rate, unsigned int *no_chromosomes, double **all_markers, unsigned int *no_total_markers);

chrom_data* data_linked_list_creation(unsigned int *chromosome_index, unsigned int *chromosome_no_loci, double *chromosome_length, double *chromosome_recom_rate, unsigned int no_chromosomes, double *chromosome_markers);

void hap_info_extract(char* string, unsigned int **sequence, double **frequency, unsigned int *no_sequence);

void iterator(gpointer key, gpointer value, gpointer user_data);

void read_pop(char *filename, GIOChannel *infile, chrom_data *head, char pop_name);

void read_classes(char *filename, GIOChannel *infile, char **gc,unsigned int *n_indv);

void free_a_hash_table_entry(gpointer key, gpointer value, gpointer user_data);

void free_all_key_value_entries (GHashTable *table);

void delete_chrom_data(chrom_data *head);

/* void output_filename(char *file_popA, char *file_popB, char *file_chrom, char *output); */

#endif
