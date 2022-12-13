#include "read_data.h"

GError *error = NULL;

chrom_data* make_data_node(unsigned int index, double length, double recom_rate, unsigned int no_loci, double *loci)
{
  chrom_data *new_data_node;
  if((new_data_node = malloc(sizeof(chrom_data))) == NULL)
    { fprintf(stderr, "Oops, out of memory!"); exit(1);}

  new_data_node->chrom_id = index;
  new_data_node->chrom_length = length;
  new_data_node->chrom_recom_rate = recom_rate;
  new_data_node->n_loci = no_loci;
  new_data_node->markers = malloc(no_loci * sizeof(double));
  assert(new_data_node->markers != NULL);
  for(int i = 0; i < no_loci; i++)
    {
      new_data_node->markers[i] = loci[i];
    }
  new_data_node->hashA = NULL;
  new_data_node->hashB = NULL;
  new_data_node->next = NULL;
  return(new_data_node);
}

void string_to_markers(char *string, unsigned int no_loci, double **markers)
{
  char *second_regexString = "([0-9]+)";
  char *second_string;                                                
  second_string = string; 
  regex_t *second_regexCompiled = (regex_t*)malloc(sizeof(regex_t));
  assert(second_regexCompiled != NULL);
  regmatch_t *pmatch;                                                 
  if(regcomp(second_regexCompiled, second_regexString, REG_EXTENDED|REG_NEWLINE))             
    {                                                                 
      printf("Could not compile regular expression.\n");              
      exit(1);                                                        
    }                                                                 
  pmatch = (regmatch_t*)malloc(sizeof(regmatch_t)*second_regexCompiled->re_nsub);
  assert(pmatch != NULL);
  if(regexec(second_regexCompiled, second_string, second_regexCompiled->re_nsub, pmatch, 0))  
    {                                                                 
      printf("Nothing matched with ""%s""\n", second_string);         
      exit(1);                                                        
    }   
  unsigned int n_loci = 0;                                            
  double *loci_position;                                              
  loci_position = malloc(TOTAL_MARKERS * sizeof(double));
  assert(loci_position != NULL);
  do {                                                                
    if (pmatch[0].rm_so != -1) {        /* The regex is matching part\
					   of a string */                                                       
      char *submatch;                                                 
      double val;                                                     
      size_t matchlen = pmatch[0].rm_eo - pmatch[0].rm_so;            
      submatch = (char*)malloc(matchlen+1);
      assert(submatch != NULL);
      strncpy(submatch, second_string + pmatch[0].rm_so, matchlen+1); 
      submatch[matchlen]='\0';                                        
      val = atof(submatch);                                           
      loci_position[n_loci] = val;                                    
      free(submatch);                                                 
    };                                                                
    second_string += pmatch[0].rm_eo;   /* Restart from last match */ 
    n_loci++;                                                         
  } while(!regexec(second_regexCompiled, second_string, second_regexCompiled->re_nsub, pmatch,0));
  *markers = malloc(n_loci * sizeof(double));
  assert(markers != NULL);
  for(int i = 0; i < n_loci; i++ )                                    
    {
      *(*markers + i) = loci_position[i];
#ifdef DEBUG
      printf("Locus[%d]: %lf\t %lf\n", i+1, loci_position[i],*(*markers + i));
#endif
    }  
  regfree(second_regexCompiled);                                      
  free(second_regexCompiled);                                         
  free(pmatch);
  free(loci_position);

}




void read_chrom(char *filename, GIOChannel *infile, unsigned int **chromosome_index, unsigned int **chromosome_no_loci, double **chromosome_length, double **chromosome_recom_rate, unsigned int *no_chromosomes, double **all_markers, unsigned int *no_total_markers)
{
  char *complete_file;
  if(g_io_channel_read_to_end(infile,&complete_file,NULL,&error) != G_IO_STATUS_NORMAL)
    {
      fprintf(stderr,"Found the file: '%s' but could not read the rest of the line\n ", filename);
      exit(1);
    }
  char *source;
  source = complete_file;
  char * regexString = "([0-9]+)[[:blank:]]+([0-9]+\\.?[0-9]*)[[:blank:]]+([0-9\
]+\\.?[0-9]*)[[:blank:]]+([0-9]+)[[:blank:]]+(.+)";
  size_t maxGroups = 6;
                                                                                
  regex_t regexCompiled;
  regmatch_t groupArray[maxGroups];

  unsigned int *chrom_index, *chrom_no_loci;
  double *chrom_length, *chrom_recom_rate;
  double *chrom_markers;
  chrom_index = calloc(MAX_CHROMOSOMES, sizeof(unsigned int));
  chrom_no_loci = calloc(MAX_CHROMOSOMES, sizeof(unsigned int));
  chrom_length = calloc(MAX_CHROMOSOMES, sizeof(double));
  chrom_recom_rate = calloc(MAX_CHROMOSOMES, sizeof(double));
  chrom_markers = calloc(MAX_CHROMOSOMES * TOTAL_MARKERS, sizeof(double));
  assert(chrom_index != NULL);
  assert(chrom_no_loci != NULL);
  assert(chrom_length != NULL);
  assert(chrom_recom_rate != NULL);
  assert(chrom_markers != NULL);
  
  if (regcomp(&regexCompiled, regexString, REG_EXTENDED|REG_NEWLINE))
    {
      printf("Could not compile regular expression.\n");
      exit(1);
    };

  unsigned int c = 0;
  unsigned int chrom_no = 0;
  unsigned int val_1, val_4;
  double val_2, val_3;
  int marker_index = 0;
  for(c = 0; c < MAX_CHROMOSOMES; c++)
    {
      if (regexec(&regexCompiled, source, maxGroups, groupArray, 0) == 0)
	{
	  unsigned int g = 0;
	  unsigned int offset = 0;
	  char sourceCopy[strlen(source) + 1];
	  strcpy(sourceCopy, source);
	  sourceCopy[groupArray[g].rm_eo] = 0;
	  for (g = 0; g < maxGroups; g++)
	    {
	      if (groupArray[g].rm_so == (size_t)-1)
		break;  // No more groups

	      if(g == 0)
		{
		  offset = groupArray[g].rm_eo;
		}
	      if(g == 1)
		{
		  char *submatch_1;
		  size_t matchlen_1 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_1 = (char*)malloc(matchlen_1+1);
		  assert(submatch_1 != NULL);
		  strncpy(submatch_1, sourceCopy + groupArray[g].rm_so, matchlen_1+1);
		  submatch_1[matchlen_1]='\0';
		  val_1 = atoi(submatch_1);
		  chrom_index[chrom_no] = val_1;
#ifdef DEBUG
		  printf("Chromosome ID: %u\n\n", val_1);
#endif
		  free(submatch_1);
		}
	      if(g == 2)
		{
		  char *submatch_2;
		  size_t matchlen_2 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_2 = (char*)malloc(matchlen_2+1);
		  assert(submatch_2 != NULL);
		  strncpy(submatch_2, sourceCopy + groupArray[g].rm_so, matchlen_2+1);
		  submatch_2[matchlen_2]='\0';
		  val_2 = atof(submatch_2);
		  chrom_length[chrom_no] = val_2;
#ifdef DEBUG
		  printf("Chromosome length (in Mb): %lf\n\n", val_2);
#endif
		  free(submatch_2);
		}
	      if(g == 3)
		{
		  char *submatch_3;
		  size_t matchlen_3 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_3 = (char*)malloc(matchlen_3+1);
		  assert(submatch_3 != NULL);
		  strncpy(submatch_3, sourceCopy + groupArray[g].rm_so, matchlen_3+1);
		  submatch_3[matchlen_3]='\0';
		  val_3 = atof(submatch_3);
		  chrom_recom_rate[chrom_no] = val_3;
#ifdef DEBUG
		  printf("Recombination rate (in cM/Mb): %lf\n\n", val_3);
#endif
		  free(submatch_3);
		}
	      if(g == 4)
		{
		  char *submatch_4;
		  size_t matchlen_4 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_4 = (char*)malloc(matchlen_4+1);
		  assert(submatch_4 != NULL);
		  strncpy(submatch_4, sourceCopy + groupArray[g].rm_so, matchlen_4+1);
		  submatch_4[matchlen_4]='\0';
		  val_4 = atoi(submatch_4);
		  chrom_no_loci[chrom_no] = val_4;
#ifdef DEBUG
		  printf("Number of markers: %u\n\n", val_4);
#endif
		  free(submatch_4);
		}
	      if(g == 5)
	      	{
		  char *submatch_5;
		  size_t matchlen_5 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_5 = (char*)malloc(matchlen_5+1);
		  assert(submatch_5 != NULL);
		  strncpy(submatch_5, sourceCopy + groupArray[g].rm_so, matchlen_5+1);
		  submatch_5[matchlen_5]='\0';
#ifdef DEBUG
		  printf("%s,\n\n",submatch_5);
#endif
		  double *loci;
		  string_to_markers(submatch_5, chrom_no_loci[chrom_no], &loci);
		  for(int i = 0; i < val_4; i++)
		    {
		      chrom_markers[marker_index + i] = loci[i];
		    }
		  free(loci);
		  free(submatch_5);
	      	}
	    }
	  marker_index = marker_index + val_4;	  
	  source = source + offset;
	  chrom_no++;
	}
    }
#ifdef DEBUG
  printf("# chrom:%u\n\n\n", chrom_no);
#endif
  
  *chromosome_index = malloc(chrom_no * sizeof(unsigned int));
  *chromosome_no_loci = malloc(chrom_no * sizeof(unsigned int));
  *chromosome_length = malloc(chrom_no * sizeof(double));
  *chromosome_recom_rate = malloc(chrom_no * sizeof(double));
  *all_markers = malloc(marker_index * sizeof(double));
  assert(chromosome_index != NULL);
  assert(chromosome_no_loci != NULL);
  assert(chromosome_length != NULL);
  assert(chromosome_recom_rate != NULL);
  assert(all_markers != NULL);

  for(int i = 0; i < chrom_no; i++)
    {
      *(*chromosome_index + i) = chrom_index[i];
      *(*chromosome_no_loci + i) = chrom_no_loci[i];
      *(*chromosome_length + i) = chrom_length[i];
      *(*chromosome_recom_rate + i) = chrom_recom_rate[i];
    }
  for(int i = 0; i < marker_index; i++)
    {
      *(*all_markers + i) = chrom_markers[i];
    }

  *no_chromosomes = chrom_no;
  *no_total_markers = marker_index;
  regfree(&regexCompiled);
  g_free(complete_file);
  free(chrom_index);
  free(chrom_no_loci);
  free(chrom_length);
  free(chrom_recom_rate);
  free(chrom_markers);
}

chrom_data* data_linked_list_creation(unsigned int *chromosome_index, unsigned int *chromosome_no_loci, double *chromosome_length, double *chromosome_recom_rate, unsigned int no_chromosomes, double *chromosome_markers)
{
  chrom_data* head, *current;
  int first = 0;
  int marker_index = 0;
  for(int i = 0; i < no_chromosomes; i++)
    {
      double *loci;
      loci = chromosome_markers + marker_index;
      chrom_data *next_data_node = make_data_node(chromosome_index[i], chromosome_length[i], chromosome_recom_rate[i], chromosome_no_loci[i], loci);
      marker_index = marker_index + chromosome_no_loci[i];
      if(first == 0)
	{
	  head = next_data_node;
	  current = head;
	  first = 1;
	}
      else
	{
	  current->next = next_data_node;
	  current = current->next;
	}
    }
  return(head);

}


void hap_info_extract(char* string, unsigned int **sequence, double **frequency, unsigned int *no_sequence)
{
  char *second_string;
  second_string = string;
  char *second_regexString = "([0-1]+:[0]*\\.[0-9]+)";
  regex_t *second_regexCompiled = (regex_t*)malloc(sizeof(regex_t));
  assert(second_regexCompiled != NULL);
  regmatch_t *pmatch;
  if(regcomp(second_regexCompiled, second_regexString, REG_EXTENDED|REG_NEWLINE))
    {                                                                            
      printf("Could not compile regular expression.\n");
      exit(1);
    }                                                                            
  pmatch = (regmatch_t*)malloc(sizeof(regmatch_t)*second_regexCompiled->re_nsub);
  assert(pmatch != NULL);
  if(regexec(second_regexCompiled, second_string, second_regexCompiled->re_nsub, pmatch, 0))
    {                                                                                       
      printf("Nothing matched with ""%s""\n", second_string);
      exit(1);
    }
  unsigned int n_hapseq = 0;
  unsigned int *hap_seq;
  double *hap_freq;
  hap_seq = malloc(MAX_HAPSEQS * sizeof(unsigned int));
  hap_freq = malloc(MAX_HAPSEQS * sizeof(double));
  assert(hap_seq != NULL);
  assert(hap_freq != NULL);
  do {                                                                          
    if (pmatch[0].rm_so != -1) {        /* The regex is matching part of a string */
      char *submatch;
      size_t matchlen = pmatch[0].rm_eo - pmatch[0].rm_so;
      submatch = (char*)malloc(matchlen+1);
      assert(submatch != NULL);
      strncpy(submatch, second_string + pmatch[0].rm_so, matchlen+1);
      submatch[matchlen]='\0';
#ifdef DEBUG
      printf("%s\n",submatch);
#endif
      char *source;
      source = submatch;
      char * regexString = "([0-1]+):([0]*\\.[0-9]+)";
      size_t maxGroups = 3;
      regex_t regexCompiled;
      regmatch_t groupArray[maxGroups];
      if (regcomp(&regexCompiled, regexString, REG_EXTENDED))                                      
        {                                                                                          
          printf("Could not compile regular expression.\n");
          exit(1);
        };
      if (regexec(&regexCompiled, source, maxGroups, groupArray, 0) == 0)
        {
          unsigned int g = 0;
          for (g = 0; g < maxGroups; g++)                                                          
            {                                                                                      
              if (groupArray[g].rm_so == (size_t)-1)                                               
                break;  // No more groups
              char sourceCopy[strlen(source) + 1];
              strcpy(sourceCopy, source);
              sourceCopy[groupArray[g].rm_eo] = 0;
	      if(g == 1)
		{
		  char *submatch_1;
		  long val_1;
		  size_t matchlen_1 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_1 = (char*)malloc(matchlen_1+1);
		  assert(submatch_1 != NULL);
		  strncpy(submatch_1, sourceCopy + groupArray[g].rm_so, matchlen_1+1);
		  submatch_1[matchlen_1]='\0';
		  val_1 = strtol(submatch_1, NULL, 2);
		  hap_seq[n_hapseq] = (unsigned int) val_1;
		  free(submatch_1);
		}
	      if(g == 2)
		{
		  char *submatch_2;
		  double val_2;
		  size_t matchlen_2 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_2 = (char*)malloc(matchlen_2+1);
		  assert(submatch_2 != NULL);
		  strncpy(submatch_2, sourceCopy + groupArray[g].rm_so, matchlen_2+1);
		  submatch_2[matchlen_2]='\0';
		  val_2 = atof(submatch_2);
		  hap_freq[n_hapseq] = val_2;
		  free(submatch_2);
		}
            }
        }
      regfree(&regexCompiled);
      free(submatch);
    };
    second_string += pmatch[0].rm_eo;   /* Restart from last match */
    n_hapseq++;
  } while(!regexec(second_regexCompiled, second_string, second_regexCompiled->re_nsub, pmatch, 0)); 
  *sequence = malloc(n_hapseq * sizeof(unsigned int));
  *frequency = malloc(n_hapseq * sizeof(double));
  assert(sequence != NULL);
  assert(frequency != NULL);
  for(int i = 0; i < n_hapseq; i++)
    {
      *(*sequence + i)  = hap_seq[i];
      *(*frequency + i) = hap_freq[i];
    }
  *no_sequence = n_hapseq;
  regfree(second_regexCompiled);
  free(second_regexCompiled);
  free(pmatch);
  free(hap_seq);
  free(hap_freq);
}

void iterator(gpointer key, gpointer value, gpointer user_data)  {
  printf(user_data, *(guint32*) key, *(gdouble*)value);
}

void read_pop(char *filename, GIOChannel *infile, chrom_data *head, char pop_name)
{
  chrom_data* current;
  current = head;
  char *complete_file;
  if(g_io_channel_read_to_end(infile,&complete_file,NULL,&error) != G_IO_STATUS_NORMAL)
    {
      fprintf(stderr,"Found the file: '%s' but could not read the rest of the line\n ", filename);
      exit(1);
    }
  char *source;
  source = complete_file;
  char * regexString = "([0-9]+)[[:blank:]]+(.+)";
  /* char * regexString = "([0-9]+)[[:blank:]]+([0-1]+:[0]*\\.[0-9]+)+"; */
  size_t maxGroups = 3;
                                                                                
  regex_t regexCompiled;
  regmatch_t groupArray[maxGroups];

  if (regcomp(&regexCompiled, regexString, REG_EXTENDED|REG_NEWLINE))
    {
      printf("Could not compile regular expression.\n");
      exit(1);
    };

  /* unsigned int c = 0; */
  unsigned int chrom_no = 0;
  unsigned int *temporary_key;
  gdouble* temporary_freq;
  while(current != NULL)
  /* for(c = 0; c < MAX_CHROMOSOMES; c++) */
    {
      if (regexec(&regexCompiled, source, maxGroups, groupArray, 0) == 0)
	{
	  unsigned int g = 0;
	  unsigned int offset = 0;
	  char sourceCopy[strlen(source) + 1];
	  strcpy(sourceCopy, source);
	  sourceCopy[groupArray[g].rm_eo] = 0;
	  for (g = 0; g < maxGroups; g++)
	    {
	      if (groupArray[g].rm_so == (size_t)-1)
		break;  // No more groups

	      if(g == 0)
		{
		  offset = groupArray[g].rm_eo;
		}
	      if(g == 1)
		{
		  char *submatch_1;
		  size_t matchlen_1 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_1 = (char*)malloc(matchlen_1+1);
		  assert(submatch_1 != NULL);
		  strncpy(submatch_1, sourceCopy + groupArray[g].rm_so, matchlen_1+1);
		  submatch_1[matchlen_1]='\0';
#ifdef DEBUG
		  printf("%s\n\n",submatch_1);
#endif
		  free(submatch_1);
		}
	      if(g == 2)
		{
		  if(pop_name == 'A')
		    {		 
		      char *submatch_2;
		      size_t matchlen_2 = groupArray[g].rm_eo - groupArray[g].rm_so;
		      submatch_2 = (char*)malloc(matchlen_2+1);
		      assert(submatch_2 != NULL);
		      strncpy(submatch_2, sourceCopy + groupArray[g].rm_so, matchlen_2+1);
		      submatch_2[matchlen_2]='\0';
#ifdef DEBUG
		      printf("~~ %s ~~\n\n", submatch_2);
#endif
		      unsigned int *hap_sequence, no_hap_sequence;
		      double *hap_frequency;
		      hap_info_extract(submatch_2, &hap_sequence, &hap_frequency, &no_hap_sequence);
#ifdef DEBUG
		      for(int i = 0; i < no_hap_sequence; i++)
			{
			  printf("Sequence: %u\tFrequency: %lf\n",hap_sequence[i], hap_frequency[i]);
			}
#endif
		      current->hashA = g_hash_table_new/* _full */(g_direct_hash, g_direct_equal/* , free, free */);
		      for(int i = 0; i < no_hap_sequence; i++)
			{
			  temporary_key = g_new(guint32,1);
			  *temporary_key = hap_sequence[i];
			  temporary_freq = g_new(gdouble,1);
			  *temporary_freq = hap_frequency[i];
			  assert(temporary_key != NULL);
			  assert(temporary_freq != NULL);
			  g_hash_table_insert(current->hashA, temporary_key, temporary_freq);
			}
		      
		      free(hap_sequence);
		      free(hap_frequency);
		      free(submatch_2);
		      
		    }
		  else
		    {
		      char *submatch_2;
		      size_t matchlen_2 = groupArray[g].rm_eo - groupArray[g].rm_so;
		      submatch_2 = (char*)malloc(matchlen_2+1);
		      assert(submatch_2 != NULL);
		      strncpy(submatch_2, sourceCopy + groupArray[g].rm_so, matchlen_2+1);
		      submatch_2[matchlen_2]='\0';
#ifdef DEBUG
		      printf("~~ %s ~~\n\n", submatch_2);
#endif
		      unsigned int *hap_sequence, no_hap_sequence;
		      double *hap_frequency;
		      hap_info_extract(submatch_2, &hap_sequence, &hap_frequency, &no_hap_sequence);
#ifdef DEBUG
		      for(int i = 0; i < no_hap_sequence; i++)
			{
			  printf("Sequence: %u\tFrequency: %lf\n",hap_sequence[i], hap_frequency[i]);
			}
#endif
		      current->hashB = g_hash_table_new/* _full */(g_direct_hash, g_direct_equal/* , free, free */);
		      for(int i = 0; i < no_hap_sequence; i++)
			{
			  temporary_key = g_new(guint32,1);
			  *temporary_key = hap_sequence[i];
			  temporary_freq = g_new(gdouble,1);
			  *temporary_freq = hap_frequency[i];
			  assert(temporary_key != NULL);
			  assert(temporary_freq != NULL);
			  g_hash_table_insert(current->hashB, temporary_key, temporary_freq);
			}
		      		      
		      free(hap_sequence);
		      free(hap_frequency);
		      free(submatch_2);

		    }

		}
	    }
	  source = source + offset;
	  chrom_no++;
	  /* printf("\n\n"); */
	}
      current = current->next;
    }

  regfree(&regexCompiled);
  g_free(complete_file);
}


void read_classes(char *filename, GIOChannel *infile, char **gc,unsigned int *n_indv)
{
  char *complete_file;
  if(g_io_channel_read_to_end(infile,&complete_file,NULL,&error) != G_IO_STATUS_NORMAL)
    {
      fprintf(stderr,"Found the file: '%s' but could not read the rest of the line\n ", filename);
      exit(1);
    }
  char *source;
  source = complete_file;
  char *regexString = "([a-z])";
  regex_t *regexCompiled = (regex_t*) malloc(sizeof(regex_t));
  assert(regexCompiled != NULL);
  regmatch_t *pmatch;
  if(regcomp(regexCompiled, regexString, REG_EXTENDED|REG_NEWLINE))
    {
      printf("Could not compile regular expression.\n");
      exit(1);
    }
  pmatch = (regmatch_t*)malloc(sizeof(regmatch_t) * regexCompiled->re_nsub);
  assert(pmatch != NULL);
  if(regexec(regexCompiled, source, regexCompiled->re_nsub, pmatch, 0))
    {
      printf("Nothing matched with ""%s""\n", source);
      exit(1);
    }
  unsigned int no_individuals = 0;
  char *genealogical_class;
  genealogical_class = malloc((MAX_INDIVIDUALS + 1) * sizeof(char));
  assert(genealogical_class != NULL);
  int first = 0;
  do {
    if(pmatch[0].rm_so != -1) {
      char *submatch;
      size_t matchlen = pmatch[0].rm_eo - pmatch[0].rm_so;
      submatch = (char*)malloc(matchlen+1);
      assert(submatch != NULL);
      strncpy(submatch,source + pmatch[0].rm_so,matchlen+1);
      submatch[matchlen] = '\0';
#ifdef DEBUG
      printf("%s\n",submatch);
#endif
      if(first == 0)
	{
	  strcpy(genealogical_class,submatch);
	  first = 1;
	}
      else
	{
	  strcat(genealogical_class,submatch);
	}
      free(submatch);
    }
    source = source + pmatch[0].rm_eo;
    no_individuals++;
  } while(!regexec(regexCompiled,source,regexCompiled->re_nsub,pmatch,0));
#ifdef DEBUG
  printf("Total number of individuals:%u\n",no_individuals);
#endif
  *n_indv = no_individuals;
  *gc = (char*)malloc(no_individuals+1);
  assert(gc != NULL);
  /* for(int i = 0; i < no_individuals; i++) */
  /*   { */
  /*     strcpy(gc[i],genealogical_class[i]); */
  /*     printf("Class:%c\t%c\n",genealogical_class[i],gc[i]); */
  /*   } */

#ifdef DEBUG
  printf("%s has length %ld\n",genealogical_class, strlen(genealogical_class));
#endif
  strcpy(*gc,genealogical_class);
#ifdef DEBUG
  printf("%s\n",*gc);
#endif
  regfree(regexCompiled);
  free(regexCompiled);
  free(pmatch);
  g_free(complete_file);
  free(genealogical_class);
}

void free_a_hash_table_entry(gpointer key, gpointer value, gpointer user_data)
{
  g_free(key);
  g_free(value);
}

void free_all_key_value_entries (GHashTable *table)
{
    g_hash_table_foreach (table, free_a_hash_table_entry, NULL);
    g_hash_table_destroy (table);
}

void delete_chrom_data(chrom_data *head)
{
  chrom_data *temp;
  while(head != NULL)
    {
      temp = head;
      head = head->next;
      free(temp->markers);
      free_all_key_value_entries (temp->hashA);
      free_all_key_value_entries (temp->hashB);
      /* g_hash_table_destroy(temp->hashA); */
      /* g_hash_table_destroy(temp->hashB); */
      free(temp);
    }
}

/* void output_filename(char *file_popA, char *file_popB, char *file_chrom, char *output) */
/* { */
/*   /\* char *source_A = "c10_m10_h10_an20_hc0.1.popA"; *\/ */
/*   /\* char *source_B = "c10_m10_h10_an20_hc1.popB"; *\/ */
/*   char * regexString = "(c[0-9]+)_(m[0-9]+)_(h[0-9]+)_(a[nu][0-9]+)_(hc(0\\.1|1|EQUAL))\\.pop[AB]"; */
/*   size_t maxGroups = 7; */


/*   regex_t regexCompiled; */
/*   regmatch_t groupArray[maxGroups]; */
/*   unsigned int m = 0; */
/*   char * cursor_A, *cursor_B; */

/*   printf("%s\t%s\n",file_popA,file_popB); /\* just checking *\/ */
/*   cursor_A = file_popA; */
/*   cursor_B = file_popB; */
/*   /\* strcat(cursor_A,file_popA); *\/ */
/*   /\* strcat(cursor_B,file_popB); *\/ */


/*   if (regcomp(&regexCompiled, regexString, REG_EXTENDED)) */
/*     { */
/*       printf("Could not compile regular expression.\n"); */
/*       exit(1); */
/*     }; */




/*   char filename[100]; */
/*   filename[0] = '\0';   */

/*   while(regexec(&regexCompiled, cursor_A, maxGroups, groupArray, 0) == 0 && regexec(&regexCompiled, cursor_B, maxGroups, groupArray, 0) == 0) */
/*     { */

/*       unsigned int g = 0; */
/*       unsigned int offset_A = 0; */
/*       unsigned int offset_B = 0; */

/*       for (g = 0; g < (maxGroups-1); g++) */
/*         { */
/*           if (groupArray[g].rm_so == (size_t)-1) */
/*             break;  // No more groups */
	  
/*           else if (g == 0){ */
/*             offset_A = groupArray[g].rm_eo; */
/* 	    offset_B = groupArray[g].rm_eo; */
/* 	  } */

/* 	  else if(g != 0){ */
/* 	    char cursorCopy_A[strlen(cursor_A) + 1]; */
/* 	    char cursorCopy_B[strlen(cursor_B) + 1]; */
/* 	    strcpy(cursorCopy_A, cursor_A); */
/* 	    strcpy(cursorCopy_B, cursor_B); */
/* 	    cursorCopy_A[groupArray[g].rm_eo] = 0; */
/* 	    cursorCopy_B[groupArray[g].rm_eo] = 0; */
/* 	    printf("Group %u: %s\tGroup %u: %s\n", */
/* 		   g,cursorCopy_A + groupArray[g].rm_so,g,cursorCopy_B + groupArray[g].rm_so); */


/* 	    if(strcmp(cursorCopy_A + groupArray[g].rm_so,cursorCopy_B + groupArray[g].rm_so) == 0){ */

/* 	      strcat(filename,cursorCopy_A + groupArray[g].rm_so); */
/* 	      strcat(filename,"_"); */
/* 	    } */
/* 	    else { */
/* 	      fprintf(stderr,"\nOops, parameter values does not match for two populations!\n"); */
/* 	      exit(1); */
/* 	    } */
/* 	  } */
/*         } */

      
/*       cursor_A += offset_A; */
/*       cursor_B += offset_B; */
/*       m++; */

/*     } */



/*   /\* char *source_chrom = "c10_m10_r1.chrom"; *\/ */
/*   char *regexString_chrom = "(c[0-9]+)_(m[0-9]+)_(r[0-9]+)\\.chrom"; */
/*   size_t maxGroups_chrom = 4; */

/*   regex_t regexCompiled_chrom; */
/*   regmatch_t groupArray_chrom[maxGroups_chrom]; */
/*   unsigned int m_chrom; */
/*   char *cursor_chrom; */

/*   if (regcomp(&regexCompiled_chrom, regexString_chrom, REG_EXTENDED)) */
/*     { */
/*       printf("Could not compile regular expression.\n"); */
/*       exit(1); */
/*     }; */


/*   m_chrom = 0; */
/*   cursor_chrom = file_chrom; */


/*   while(regexec(&regexCompiled_chrom, cursor_chrom, maxGroups_chrom, groupArray_chrom, 0) == 0) */
/*     { */

/*       unsigned int g_chrom = 0; */
/*       unsigned int offset_chrom = 0; */
/*       for (g_chrom = 0; g_chrom < maxGroups_chrom; g_chrom++) */
/*         { */
/*           if (groupArray_chrom[g_chrom].rm_so == (size_t)-1) */
/*             break;  // No more groups */

/*           if (g_chrom == 0) */
/*             offset_chrom = groupArray_chrom[g_chrom].rm_eo; */

/*           char cursorCopy_chrom[strlen(cursor_chrom) + 1]; */
/*           strcpy(cursorCopy_chrom, cursor_chrom); */
/*           cursorCopy_chrom[groupArray_chrom[g_chrom].rm_eo] = 0; */


/*           printf("Match %u, Group %u: [%2u-%2u]: %s\n", */
/*                  m_chrom, g_chrom, groupArray_chrom[g_chrom].rm_so, groupArray_chrom[g_chrom].rm_eo, */
/*                  cursorCopy_chrom + groupArray_chrom[g_chrom].rm_so); */
/*   	  if(g_chrom == (maxGroups_chrom-1)){ */
/*   	      strcat(filename,cursorCopy_chrom + groupArray_chrom[g_chrom].rm_so); */
/* 	      strcat(filename,".sim"); */
/*   	  } */
/*         } */
/*       cursor_chrom += offset_chrom; */
/*       m_chrom++; */
/*     } */


  

/*   printf("%s\n", filename); */
/*   strcpy(output,filename); */

  
/*   regfree(&regexCompiled); */
/*   regfree(&regexCompiled_chrom); */
/* } */

