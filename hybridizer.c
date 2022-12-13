#include "read_data.h"
/* #include<assert.h> */
/* #include<unistd.h> */
/* #include<ctype.h> */
/* #include<stdbool.h> */
/* #include <gsl/gsl_rng.h> */
/* #include <gsl/gsl_randist.h> */
/* #include<time.h> */

GError *e = NULL;
GIOChannel *infile_chrom = NULL;
GIOChannel *infile_A = NULL;
GIOChannel *infile_B = NULL;
GIOChannel *infile_gc = NULL;
GIOChannel *outfile = NULL;
/* GIOChannel *outfile_simulation = NULL; */

gsl_rng *r;

int model_a = 0;
int model_b = 0;
int model_c = 0;
int model_d = 0;
int model_e = 0;
int model_f = 0;

typedef struct node
{
  unsigned int state;
  double point;
  struct node* next;
}chromosome_state;

chromosome_state *make_node(unsigned int s,double l)
{
  chromosome_state *new_node;
  if((new_node = malloc(sizeof(chromosome_state)))==NULL)
    { fprintf(stderr,"Oops, out of memory!"); exit(1); }

  new_node->state = s;
  new_node->point = l;
  new_node->next = NULL;
  return(new_node);

}

int cmpfunc (const void * a, const void * b) /* Comparing function for sorting an array in ascending order */
{
  if (*(double*)a > *(double*)b)
    return 1;
  else if (*(double*)a < *(double*)b)
    return -1;
  else
    return 0;
}


unsigned int no_of_recombinations(double mu)
{
  unsigned int k = gsl_ran_poisson (r,mu);
  return(k);
}

void recombination_position(double **array, unsigned int n_recombination)
{
  if(n_recombination == 0)
    {
      *array = malloc(sizeof(double));
      assert(array != NULL);
      **array = 1;
#ifdef DEBUG
      printf("There was no recombination\n");
#endif
    }
  else
    {
      *array = malloc((n_recombination + 1) * sizeof(double));
      assert(array != NULL);
      for (int i = 0 ; i < n_recombination ; i++)
	{
	  *(*array + i) = gsl_rng_uniform_pos(r);
#ifdef DEBUG
	  printf ("%lf\t", *(*array+i));
#endif
	}
      *(*array + n_recombination) = 1; /* Adding 1 (lenght of chromosome) at the end of the array for convenience */
      /* printf("\n"); */
      qsort(*array,n_recombination+1, sizeof(double), cmpfunc);
#ifdef DEBUG
      printf("\nAfter sorting the recombinations over the chromosome in ascending order is: \n");
      
      for(int j = 0 ; j < n_recombination+1 ; j++ ) {
	printf("%lf\t", *(*array + j));
      }
#endif
    }
  /* printf("\n"); */

}


void delete_linked_list(chromosome_state *head)
{
  chromosome_state *temp;
  while(head != NULL)
    {
      temp = head;
      head = head->next;
      free(temp);
    }
}

chromosome_state* Linked_List(double *array, unsigned int n_recombination)
{
  
  chromosome_state *head;
  unsigned int k = gsl_ran_bernoulli(r,0.5);

  if(n_recombination == 0)
    {
      head = make_node(k,*array);
      head->next = NULL;
#ifdef DEBUG
      printf("The only state : %u\n",head->state);
#endif
    }
  else
    {
      chromosome_state *current;
      unsigned int temp_state;
      head = make_node(k,*array); /* Creating the first node (head node) */
      temp_state = k;
#ifdef DEBUG
      printf("\n First state: %u Time to first crossover: %lf\n",head->state,head->point);
#endif
      int first = 1;
      for(int i = 1; i < n_recombination+1 ; i++)
	{
	  if (first == 1) /* Creating the second node */
	    {
	      if(temp_state == 0)
		{
		  chromosome_state *next_node = make_node(temp_state + 1,*(array+i));
		  current = next_node;
#ifdef DEBUG
		  printf("\n State:%u Waiting time:%lf",current->state, current->point);
#endif
		  temp_state = temp_state + 1;
		  head->next = current;
		  first = 0;
		}
	      else
		{
		  chromosome_state *next_node = make_node(temp_state - 1,*(array+i));
		  current = next_node;
#ifdef DEBUG
		  printf("\n State:%u Waiting time:%lf",current->state, current->point);
#endif
		  temp_state = temp_state - 1;
		  head->next = current;
		  first = 0;
		}
	    }
	  else /* Creating the other nodes after the second one */
	    {
	      if( temp_state == 0)
		{
		  chromosome_state *next_node = make_node(temp_state + 1, *(array+i));
		  current->next = next_node;
#ifdef DEBUG
		  printf("\n State:%u Waiting time:%lf",current->next->state, current->next->point);
#endif
		  temp_state = temp_state + 1;
		  current = current->next;
		}
	      else
		{
		  chromosome_state *next_node = make_node(temp_state - 1, *(array+i));
		  current->next = next_node;
#ifdef DEBUG
		  printf("\n State:%u Waiting time:%lf",current->next->state, current->next->point);
#endif
		  temp_state = temp_state - 1;
		  current = current->next;
		}
	    }
	}
    }
  return(head);
}

unsigned int setBit(unsigned int n, unsigned int k)
{
  unsigned int temp;
  temp = n;
  temp = (temp | (1 << k));

  return(temp);
}


void recombinant_ancestry(unsigned int *I_A, unsigned int *I_B, chromosome_state *r, unsigned int n_loci, double *locus_position)
{
  unsigned int A, B;
  A = 0;
  B = 0;
  chromosome_state *head;
  head = r;
  for(int i = 0; i < n_loci; i++)
    {
      if(locus_position[i] <= head->point)
	{
	  if(head->state == 0)
	    {
	      A = setBit(A, n_loci-i-1);
	    }
	  else
	    {
	      B = setBit(B, n_loci-i-1);
	    }
	}
      else
	{
	  head = head->next;
	  while(locus_position[i] > head->point)
	    {
	      head = head->next;
	    }
	  if(head->state == 0)
	    {
	      A = setBit(A, n_loci-i-1);
	    }
	  else
	    {
	      B = setBit(B, n_loci-i-1);
	    }
	}
    }
  *I_A = A;
  *I_B = B;
  
}


unsigned int recom_hap(unsigned int C_A, unsigned int C_B, unsigned int I_A, unsigned int I_B)
{
  unsigned int temp1, temp2, final;
  temp1 = C_A & I_A ;
  temp2 = C_B & I_B ;
  final = temp1 | temp2 ;
  return(final) ;
}


unsigned int checkSet(unsigned int n, unsigned int k)
{
  unsigned int a;
  a = (n & (1 << k));
  if(a == 0)
    {
      return(0);
    }
  else
    {
      return(1);
    }
}

void simulate_model_a(chrom_data *head, unsigned int no_chromosomes, char *simulated_output[] /* char simulated_output[][BUFFER_SIZE] */)
{
  chrom_data *current;
  unsigned int /* table_size_A, */ table_size_B;
  unsigned int homolog_chromosomes[no_chromosomes][2];

  current = head;
  unsigned int j = 0;
  unsigned int total_markers = 0;
  while(current != NULL)
    {
      table_size_B = g_hash_table_size(current->hashB);
      gpointer *keys;
      keys = g_hash_table_get_keys_as_array(current->hashB, &table_size_B);
      double *values;
      values = malloc(table_size_B * sizeof(double));
      assert(values != NULL);
      for(int i = 0; i < table_size_B; i++)
	{
	  values[i] = *(gdouble*) g_hash_table_lookup(current->hashB,keys[i]);
#ifdef DEBUG
	  printf("Sequence[%d] = %u\n",i,*(guint32*)keys[i]);
	  printf("%u has frequency %lf\n",*(guint32*)keys[i], values[i]);
#endif
	}

      unsigned int *hapB_indicator_1, *hapB_indicator_2;
      hapB_indicator_1 = malloc(table_size_B * sizeof(unsigned int));
      hapB_indicator_2 = malloc(table_size_B * sizeof(unsigned int));
      assert(hapB_indicator_1 != NULL);
      assert(hapB_indicator_2 != NULL);
      gsl_ran_multinomial(r,table_size_B,1,values,hapB_indicator_1);
      gsl_ran_multinomial(r,table_size_B,1,values,hapB_indicator_2);
#ifdef DEBUG
      for(int i = 0; i < table_size_B; i++)
	{
	  printf("[%d]:%u\t%u\n",i,hapB_indicator_1[i],hapB_indicator_2[i]);
	}
#endif
      for(unsigned int i = 0; i < table_size_B; i++)
	{
	  if(hapB_indicator_1[i] == 1)
	    {
	      homolog_chromosomes[j][0] = *(guint32*)keys[i];
	    }
	}
      for(unsigned int i = 0; i < table_size_B; i++)
	{
	  if(hapB_indicator_2[i] == 1)
	    {
	      homolog_chromosomes[j][1] = *(guint32*)keys[i];
	    }
	}
      j++;
      free(hapB_indicator_1);
      free(hapB_indicator_2);
      free(values);
      g_free(keys);
      /* printf("\n\n"); */
      total_markers = total_markers + current->n_loci;
      current = current->next;
    }
#ifdef DEBUG
  for(int k = 0; k < no_chromosomes; k++)
    {
      printf("%u\t%u\n",homolog_chromosomes[k][0],homolog_chromosomes[k][1]);
    }
#endif
  j = 0;
  int k = 0;
  char indv[total_markers][5];
  current = head;
  while(current != NULL && j < total_markers && k < no_chromosomes)
    {
      for(int i = 0; i < current->n_loci; i++)
  	{
  	  sprintf(indv[j],"\t%u|%u",checkSet(homolog_chromosomes[k][0],current->n_loci-i-1),checkSet(homolog_chromosomes[k][1],current->n_loci-i-1));

  	  j++;
  	}
      k++;
      current = current->next;
    }
  
  for(int i = 0; i < total_markers; i++)
    {
      strcat(simulated_output[i],indv[i]);
      /* printf("%s\n",simulated_output[i]); */
    }
}


void simulate_model_d(chrom_data *head, unsigned int no_chromosomes, char *simulated_output[] /* char simulated_output[][BUFFER_SIZE] */)
{
  chrom_data *current;
  unsigned int table_size_A;
  unsigned int homolog_chromosomes[no_chromosomes][2];
  
  current = head;
  unsigned int j = 0;
  unsigned int total_markers = 0;
  while(current != NULL)
    {
      table_size_A = g_hash_table_size(current->hashA);
      gpointer *keys;
      keys = g_hash_table_get_keys_as_array(current->hashA, &table_size_A);
      double *values;
      values = malloc(table_size_A * sizeof(double));
      assert(values != NULL);
      for(int i = 0; i < table_size_A; i++)
	{
	  values[i] = *(gdouble*) g_hash_table_lookup(current->hashA,keys[i]);
#ifdef DEBUG
	  printf("Sequence[%d] = %u\n",i,*(guint32*)keys[i]);
	  printf("%u has frequency %lf\n",*(guint32*)keys[i], values[i]);
#endif
	}

      unsigned int *hapA_indicator_1, *hapA_indicator_2;
      hapA_indicator_1 = malloc(table_size_A * sizeof(unsigned int));
      hapA_indicator_2 = malloc(table_size_A * sizeof(unsigned int));
      assert(hapA_indicator_1 != NULL);
      assert(hapA_indicator_2 != NULL);

      gsl_ran_multinomial(r,table_size_A,1,values,hapA_indicator_1);
      gsl_ran_multinomial(r,table_size_A,1,values,hapA_indicator_2);
#ifdef DEBUG
      for(int i = 0; i < table_size_A; i++)
	{
	  printf("[%d]:%u\t%u\n",i,hapA_indicator_1[i],hapA_indicator_2[i]);
	}
#endif
      for(unsigned int i = 0; i < table_size_A; i++)
	{
	  if(hapA_indicator_1[i] == 1)
	    {
	      homolog_chromosomes[j][0] = *(guint32*)keys[i];
	    }
	}
      for(unsigned int i = 0; i < table_size_A; i++)
	{
	  if(hapA_indicator_2[i] == 1)
	    {
	      homolog_chromosomes[j][1] = *(guint32*)keys[i];
	    }
	}
      j++;
      free(hapA_indicator_1);
      free(hapA_indicator_2);
      free(values);
      g_free(keys);
      /* printf("\n\n"); */
      total_markers = total_markers + current->n_loci;
      current = current->next;
    }
#ifdef DEBUG
  for(int k = 0; k < no_chromosomes; k++)
    {
      printf("%u\t%u\n",homolog_chromosomes[k][0],homolog_chromosomes[k][1]);
    }
#endif
  j = 0;
  int k = 0;
  char indv[total_markers][5];
  current = head;
  while(current != NULL && j < total_markers && k < no_chromosomes)
    {
      for(int i = 0; i < current->n_loci; i++)
  	{
  	  sprintf(indv[j],"\t%u|%u",checkSet(homolog_chromosomes[k][0],current->n_loci-i-1),checkSet(homolog_chromosomes[k][1],current->n_loci-i-1));

  	  j++;
  	}
      k++;
      current = current->next;
    }
  
  for(int i = 0; i < total_markers; i++)
    {
      strcat(simulated_output[i],indv[i]);
      /* printf("%s\n",simulated_output[i]); */
    }

}

void simulate_model_c(chrom_data *head, unsigned int no_chromosomes, char *simulated_output[] /* char simulated_output[][BUFFER_SIZE] */)
{
  chrom_data *current;
  unsigned int table_size_A, table_size_B;
  unsigned int homolog_chromosomes[no_chromosomes][2];

  current = head;
  unsigned int j = 0;
  unsigned int total_markers = 0;
  while(current != NULL)
    {
      table_size_A = g_hash_table_size(current->hashA);
      table_size_B = g_hash_table_size(current->hashB);
      gpointer *keys_A, *keys_B;
      keys_A = g_hash_table_get_keys_as_array(current->hashA, &table_size_A);
      keys_B = g_hash_table_get_keys_as_array(current->hashB, &table_size_B);
      double *values_A, *values_B;
      values_A = malloc(table_size_A * sizeof(double));
      values_B = malloc(table_size_B * sizeof(double));
      assert(values_A != NULL);
      assert(values_B != NULL);
      for(int i = 0; i < table_size_A; i++)
	{
	  values_A[i] = *(gdouble*) g_hash_table_lookup(current->hashA,keys_A[i]);
#ifdef DEBUG
	  printf("Sequence[%d] = %u\n",i,*(guint32*)keys_A[i]);
	  printf("%u has frequency %lf\n",*(guint32*)keys_A[i], values_A[i]);
#endif
	}

      for(int i = 0; i < table_size_B; i++)
	{
	  values_B[i] = *(gdouble*) g_hash_table_lookup(current->hashB,keys_B[i]);
#ifdef DEBUG
	  printf("Sequence[%d] = %u\n",i,*(guint32*)keys_B[i]);
	  printf("%u has frequency %lf\n",*(guint32*)keys_B[i], values_B[i]);
#endif
	}

	  
      unsigned int *hapA_indicator, *hapB_indicator;
      hapA_indicator = malloc(table_size_A * sizeof(unsigned int));
      hapB_indicator = malloc(table_size_B * sizeof(unsigned int));
      assert(hapA_indicator != NULL);
      assert(hapB_indicator != NULL);
      gsl_ran_multinomial(r,table_size_A,1,values_A,hapA_indicator);
      gsl_ran_multinomial(r,table_size_B,1,values_B,hapB_indicator);
#ifdef DEBUG
      for(int i = 0; i < table_size_A; i++)
	{
	  printf("[%d]:%u\n",i,hapA_indicator[i]);
	}
      for(int i = 0; i < table_size_B; i++)
	{
	  printf("[%d]:%u\n",i,hapB_indicator[i]);
	}
#endif

      for(unsigned int i = 0; i < table_size_A; i++)
	{
	  if(hapA_indicator[i] == 1)
	    {
	      homolog_chromosomes[j][0] = *(guint32*)keys_A[i];
	    }
	}
      for(unsigned int i = 0; i < table_size_B; i++)
	{
	  if(hapB_indicator[i] == 1)
	    {
	      homolog_chromosomes[j][1] = *(guint32*)keys_B[i];
	    }
	}
      j++;
      free(hapA_indicator);
      free(hapB_indicator);
      free(values_A);
      free(values_B);
      g_free(keys_A);
      g_free(keys_B);
      /* printf("\n\n"); */
      total_markers = total_markers + current->n_loci;
      current = current->next;
    }
#ifdef DEBUG
  for(int k = 0; k < no_chromosomes; k++)
    {
      printf("%u\t%u\n",homolog_chromosomes[k][0],homolog_chromosomes[k][1]);
    }
#endif
  j = 0;
  int k = 0;
  char indv[total_markers][5];
  current = head;
  while(current != NULL && j < total_markers && k < no_chromosomes)
    {
      for(int i = 0; i < current->n_loci; i++)
  	{
  	  sprintf(indv[j],"\t%u|%u",checkSet(homolog_chromosomes[k][0],current->n_loci-i-1),checkSet(homolog_chromosomes[k][1],current->n_loci-i-1));

  	  j++;
  	}
      k++;
      current = current->next;
    }
  
  for(int i = 0; i < total_markers; i++)
    {
      strcat(simulated_output[i],indv[i]);
      /* printf("%s\n",simulated_output[i]); */
    }

  
}


void simulate_model_b(chrom_data *head, unsigned int no_chromosomes, char *simulated_output[] /* char simulated_output[][BUFFER_SIZE] */)
{
  chrom_data *current;
  unsigned int table_size_A, table_size_B;
  unsigned int homolog_chromosomes[no_chromosomes][2];
  current = head;
  unsigned int j = 0;
  unsigned int total_markers = 0;
  while(current != NULL)
    {
      double mu; /* mu is the expected number of recombinations for the length of chromosome under consideration*/
      unsigned int n_r; /* number of recombinations */ 
      mu = ( current->chrom_recom_rate * current->chrom_length / 100.0);
      n_r = no_of_recombinations(mu);
#ifdef DEBUG
      printf("Number of recombinations: %u\n",n_r);
#endif
      double *location_r;
      recombination_position(&location_r,n_r);
      chromosome_state *head_r;
      head_r = Linked_List(location_r,n_r);

      chromosome_state *temp_head_r;
      temp_head_r = head_r;
      while(temp_head_r != NULL)
	{
#ifdef DEBUG
	  printf("\nSTATE: '%u'\t Point: %lf",temp_head_r->state,temp_head_r->point);
#endif
	  temp_head_r = temp_head_r->next;
	}
      /* printf("\n"); */
      double *loci_position;
      loci_position = malloc(current->n_loci * sizeof(double));
      assert(loci_position != NULL);
      for(int i = 0; i < current->n_loci; i++)
	{
	  loci_position[i] = (current->markers[i] / (current->chrom_length * pow(10,6)));
	}
      unsigned int A_indicator, B_indicator;
      recombinant_ancestry(&A_indicator, &B_indicator, head_r, current->n_loci, loci_position);
#ifdef DEBUG
      printf("I_A = %u\tI_B = %u\n",A_indicator,B_indicator);
#endif

      table_size_A = g_hash_table_size(current->hashA);
      table_size_B = g_hash_table_size(current->hashB);
      gpointer *keys_A, *keys_B;
      keys_A = g_hash_table_get_keys_as_array(current->hashA, &table_size_A);
      keys_B = g_hash_table_get_keys_as_array(current->hashB, &table_size_B);
      double *values_A, *values_B;
      values_A = malloc(table_size_A * sizeof(double));
      values_B = malloc(table_size_B * sizeof(double));
      assert(values_A != NULL);
      assert(values_B != NULL);
      for(int i = 0; i < table_size_A; i++)
	{
	  values_A[i] = *(gdouble*) g_hash_table_lookup(current->hashA,keys_A[i]);
#ifdef DEBUG
	  printf("Sequence[%d] = %u\n",i,*(guint32*)keys_A[i]);
	  printf("%u has frequency %lf\n",*(guint32*)keys_A[i], values_A[i]);
#endif
	}

      for(int i = 0; i < table_size_B; i++)
	{
	  values_B[i] = *(gdouble*) g_hash_table_lookup(current->hashB,keys_B[i]);
#ifdef DEBUG
	  printf("Sequence[%d] = %u\n",i,*(guint32*)keys_B[i]);
	  printf("%u has frequency %lf\n",*(guint32*)keys_B[i], values_B[i]);
#endif
	}

	  
      unsigned int *hapA_indicator_1, *hapA_indicator_2, *hapB_indicator_2;
      hapA_indicator_1 = malloc(table_size_A * sizeof(unsigned int));
      hapA_indicator_2 = malloc(table_size_A * sizeof(unsigned int));
      hapB_indicator_2 = malloc(table_size_B * sizeof(unsigned int));
      assert(hapA_indicator_1 != NULL);
      assert(hapA_indicator_2 != NULL);
      assert(hapB_indicator_2 != NULL);
      gsl_ran_multinomial(r,table_size_A,1,values_A,hapA_indicator_1);
      gsl_ran_multinomial(r,table_size_A,1,values_A,hapA_indicator_2);
      gsl_ran_multinomial(r,table_size_B,1,values_B,hapB_indicator_2);

      for(unsigned int i = 0; i < table_size_A; i++)
	{
	  if(hapA_indicator_1[i] == 1)
	    {
	      homolog_chromosomes[j][0] = *(guint32*)keys_A[i];
	    }
	}
      unsigned int chrom_A, chrom_B;
      for(int i = 0; i < table_size_A; i++)
	{
	  if(hapA_indicator_2[i] == 1)
	    {
	      chrom_A = *(guint32*)keys_A[i];
	    }
	}

      for(int i = 0; i < table_size_B; i++)
	{
	  if(hapB_indicator_2[i] == 1)
	    {
	      chrom_B = *(guint32*)keys_B[i];
	    }
	}
#ifdef DEBUG
      printf("\n ChromA: %u\tChromB:%u\n",chrom_A,chrom_B);
#endif
      unsigned int recombinant_chromosome;
      recombinant_chromosome = recom_hap(chrom_A,chrom_B,A_indicator,B_indicator);
      homolog_chromosomes[j][1] = recombinant_chromosome;
      j++;
      free(location_r);
      free(loci_position);
      delete_linked_list(head_r);
      free(hapA_indicator_1);
      free(hapA_indicator_2);
      free(hapB_indicator_2);
      free(values_A);
      free(values_B);
      g_free(keys_A);
      g_free(keys_B);
      /* printf("\n\n"); */
      total_markers = total_markers + current->n_loci;
      current = current->next;
    }
#ifdef DEBUG
  for(int k = 0; k < no_chromosomes; k++)
    {
      printf("%u\t%u\n",homolog_chromosomes[k][0],homolog_chromosomes[k][1]);
    }
#endif
  j = 0;
  int k = 0;
  char indv[total_markers][5];
  current = head;
  while(current != NULL && j < total_markers && k < no_chromosomes)
    {
      for(int i = 0; i < current->n_loci; i++)
  	{
  	  sprintf(indv[j],"\t%u|%u",checkSet(homolog_chromosomes[k][0],current->n_loci-i-1),checkSet(homolog_chromosomes[k][1],current->n_loci-i-1));

  	  j++;
  	}
      k++;
      current = current->next;
    }
  
  for(int i = 0; i < total_markers; i++)
    {
      strcat(simulated_output[i],indv[i]);
      /* printf("%s\n",simulated_output[i]); */
    }


}

void simulate_model_e(chrom_data *head, unsigned int no_chromosomes, char *simulated_output[] /* char simulated_output[][BUFFER_SIZE] */)
{
  chrom_data *current;
  unsigned int table_size_A, table_size_B;
  unsigned int homolog_chromosomes[no_chromosomes][2];
  current = head;
  unsigned int j = 0;
  unsigned int total_markers = 0;
  while(current != NULL)
    {
      double mu; /* mu is the expected number of recombinations for the length of chromosome under consideration*/
      unsigned int n_r; /* number of recombinations */ 
      mu = ( current->chrom_recom_rate * current->chrom_length / 100.0);
      n_r = no_of_recombinations(mu);
#ifdef DEBUG
      printf("Number of recombinations: %u\n",n_r);
#endif
      double *location_r;
      recombination_position(&location_r,n_r);
      chromosome_state *head_r;
      head_r = Linked_List(location_r,n_r);

      chromosome_state *temp_head_r;
      temp_head_r = head_r;
      while(temp_head_r != NULL)
	{
#ifdef DEBUG
	  printf("\nSTATE: '%u'\t Point: %lf",temp_head_r->state,temp_head_r->point);
#endif
	  temp_head_r = temp_head_r->next;
	}
      /* printf("\n"); */
      double *loci_position;
      loci_position = malloc(current->n_loci * sizeof(double));
      assert(loci_position != NULL);
      for(int i = 0; i < current->n_loci; i++)
	{
	  loci_position[i] = (current->markers[i] / (current->chrom_length * pow(10,6)));
	}
      unsigned int A_indicator, B_indicator;
      recombinant_ancestry(&A_indicator, &B_indicator, head_r, current->n_loci, loci_position);
#ifdef DEBUG
      printf("I_A = %u\tI_B = %u\n",A_indicator,B_indicator);
#endif

      table_size_A = g_hash_table_size(current->hashA);
      table_size_B = g_hash_table_size(current->hashB);
      gpointer *keys_A, *keys_B;
      keys_A = g_hash_table_get_keys_as_array(current->hashA, &table_size_A);
      keys_B = g_hash_table_get_keys_as_array(current->hashB, &table_size_B);
      double *values_A, *values_B;
      values_A = malloc(table_size_A * sizeof(double));
      values_B = malloc(table_size_B * sizeof(double));
      assert(values_A != NULL);
      assert(values_B != NULL);
      for(int i = 0; i < table_size_A; i++)
	{
	  values_A[i] = *(gdouble*) g_hash_table_lookup(current->hashA,keys_A[i]);
#ifdef DEBUG
	  printf("Sequence[%d] = %u\n",i,*(guint32*)keys_A[i]);
	  printf("%u has frequency %lf\n",*(guint32*)keys_A[i], values_A[i]);
#endif
	}

      for(int i = 0; i < table_size_B; i++)
	{
	  values_B[i] = *(gdouble*) g_hash_table_lookup(current->hashB,keys_B[i]);
#ifdef DEBUG
	  printf("Sequence[%d] = %u\n",i,*(guint32*)keys_B[i]);
	  printf("%u has frequency %lf\n",*(guint32*)keys_B[i], values_B[i]);
#endif
	}

	  
      unsigned int *hapB_indicator_1, *hapA_indicator_2, *hapB_indicator_2;
      hapB_indicator_1 = malloc(table_size_B * sizeof(unsigned int));
      hapA_indicator_2 = malloc(table_size_A * sizeof(unsigned int));
      hapB_indicator_2 = malloc(table_size_B * sizeof(unsigned int));
      assert(hapB_indicator_1 != NULL);
      assert(hapA_indicator_2 != NULL);
      assert(hapB_indicator_2 != NULL);
      gsl_ran_multinomial(r,table_size_B,1,values_B,hapB_indicator_1);
      gsl_ran_multinomial(r,table_size_A,1,values_A,hapA_indicator_2);
      gsl_ran_multinomial(r,table_size_B,1,values_B,hapB_indicator_2);

      for(unsigned int i = 0; i < table_size_B; i++)
	{
	  if(hapB_indicator_1[i] == 1)
	    {
	      homolog_chromosomes[j][0] = *(guint32*)keys_B[i];
	    }
	}
      unsigned int chrom_A, chrom_B;
      for(int i = 0; i < table_size_A; i++)
	{
	  if(hapA_indicator_2[i] == 1)
	    {
	      chrom_A = *(guint32*)keys_A[i];
	    }
	}

      for(int i = 0; i < table_size_B; i++)
	{
	  if(hapB_indicator_2[i] == 1)
	    {
	      chrom_B = *(guint32*)keys_B[i];
	    }
	}
#ifdef DEBUG
      printf("\n ChromA: %u\tChromB:%u\n",chrom_A,chrom_B);
#endif
      unsigned int recombinant_chromosome;
      recombinant_chromosome = recom_hap(chrom_A,chrom_B,A_indicator,B_indicator);
      homolog_chromosomes[j][1] = recombinant_chromosome;
      j++;
      free(location_r);
      free(loci_position);
      delete_linked_list(head_r);
      free(hapB_indicator_1);
      free(hapA_indicator_2);
      free(hapB_indicator_2);
      free(values_A);
      free(values_B);
      g_free(keys_A);
      g_free(keys_B);
      /* printf("\n\n"); */
      total_markers = total_markers + current->n_loci;
      current = current->next;
    }
#ifdef DEBUG
  for(int k = 0; k < no_chromosomes; k++)
    {
      printf("%u\t%u\n",homolog_chromosomes[k][0],homolog_chromosomes[k][1]);
    }
#endif
  j = 0;
  int k = 0;
  char indv[total_markers][5];
  current = head;
  while(current != NULL && j < total_markers && k < no_chromosomes)
    {
      for(int i = 0; i < current->n_loci; i++)
  	{
  	  sprintf(indv[j],"\t%u|%u",checkSet(homolog_chromosomes[k][0],current->n_loci-i-1),checkSet(homolog_chromosomes[k][1],current->n_loci-i-1));

  	  j++;
  	}
      k++;
      current = current->next;
    }
  
  for(int i = 0; i < total_markers; i++)
    {
      strcat(simulated_output[i],indv[i]);
      /* printf("%s\n",simulated_output[i]); */
    }

  
}

void simulate_model_f(chrom_data *head, unsigned int no_chromosomes, char *simulated_output[] /* char simulated_output[][BUFFER_SIZE] */)
{
  chrom_data *current;
  unsigned int table_size_A, table_size_B;
  unsigned int homolog_chromosomes[no_chromosomes][2];
  current = head;
  unsigned int j = 0;
  unsigned int total_markers = 0;
  while(current != NULL)
    {
      double mu; /* mu is the expected number of recombinations for the length of chromosome under consideration*/
      unsigned int n_r_1, n_r_2; /* number of recombinations for chromosome 1 & chromosome 2 */ 
      mu = ( current->chrom_recom_rate * current->chrom_length / 100.0);
      n_r_1 = no_of_recombinations(mu);
      n_r_2 = no_of_recombinations(mu);
#ifdef DEBUG
      printf("Number of recombinations for Chromsome 1: %u\n",n_r_1);
      printf("Number of recombinations for Chromsome 2: %u\n",n_r_2);
#endif
      double *location_r_1, *location_r_2;
      recombination_position(&location_r_1,n_r_1);
      recombination_position(&location_r_2,n_r_2);
      chromosome_state *head_r_1, *head_r_2;
      head_r_1 = Linked_List(location_r_1,n_r_1);
      head_r_2 = Linked_List(location_r_2,n_r_2);

      chromosome_state *temp_head_r_1, *temp_head_r_2;
      temp_head_r_1 = head_r_1;
      temp_head_r_2 = head_r_2;
      while(temp_head_r_1 != NULL)
	{
#ifdef DEBUG
	  printf("\nSTATE: '%u'\t Point: %lf",temp_head_r_1->state,temp_head_r_1->point);
#endif
	  temp_head_r_1 = temp_head_r_1->next;
	}
      /* printf("\n"); */
      while(temp_head_r_2 != NULL)
	{
#ifdef DEBUG
	  printf("\nSTATE: '%u'\t Point: %lf",temp_head_r_2->state,temp_head_r_2->point);
#endif
	  temp_head_r_2 = temp_head_r_2->next;
	}
      /* printf("\n"); */
	  
      double *loci_position;
      loci_position = malloc(current->n_loci * sizeof(double));
      assert(loci_position != NULL);
      for(int i = 0; i < current->n_loci; i++)
	{
	  loci_position[i] = (current->markers[i] / (current->chrom_length * pow(10,6)));
	}
      unsigned int A_indicator_1, B_indicator_1, A_indicator_2, B_indicator_2;
      recombinant_ancestry(&A_indicator_1, &B_indicator_1, head_r_1, current->n_loci, loci_position);
      recombinant_ancestry(&A_indicator_2, &B_indicator_2, head_r_2, current->n_loci, loci_position);
#ifdef DEBUG
      printf("Chromosome 1: I_A = %u\tI_B = %u\n",A_indicator_1,B_indicator_1);
      printf("Chromosome 2: I_A = %u\tI_B = %u\n",A_indicator_2,B_indicator_2);
#endif
      table_size_A = g_hash_table_size(current->hashA);
      table_size_B = g_hash_table_size(current->hashB);
      gpointer *keys_A, *keys_B;
      keys_A = g_hash_table_get_keys_as_array(current->hashA, &table_size_A);
      keys_B = g_hash_table_get_keys_as_array(current->hashB, &table_size_B);
      double *values_A, *values_B;
      values_A = malloc(table_size_A * sizeof(double));
      values_B = malloc(table_size_B * sizeof(double));
      assert(values_A != NULL);
      assert(values_B != NULL);
      for(int i = 0; i < table_size_A; i++)
	{
	  values_A[i] = *(gdouble*) g_hash_table_lookup(current->hashA,keys_A[i]);
#ifdef DEBUG
	  printf("Sequence[%d] = %u\n",i,*(guint32*)keys_A[i]);
	  printf("%u has frequency %lf\n",*(guint32*)keys_A[i], values_A[i]);
#endif
	}

      for(int i = 0; i < table_size_B; i++)
	{
	  values_B[i] = *(gdouble*) g_hash_table_lookup(current->hashB,keys_B[i]);
#ifdef DEBUG
	  printf("Sequence[%d] = %u\n",i,*(guint32*)keys_B[i]);
	  printf("%u has frequency %lf\n",*(guint32*)keys_B[i], values_B[i]);
#endif
	}

	  
      unsigned int *hapA_indicator_1, *hapB_indicator_1, *hapA_indicator_2, *hapB_indicator_2;
      hapA_indicator_1 = malloc(table_size_A * sizeof(unsigned int));
      hapB_indicator_1 = malloc(table_size_B * sizeof(unsigned int));
      hapA_indicator_2 = malloc(table_size_A * sizeof(unsigned int));
      hapB_indicator_2 = malloc(table_size_B * sizeof(unsigned int));
      assert(hapA_indicator_1 != NULL);
      assert(hapB_indicator_1 != NULL);
      assert(hapA_indicator_2 != NULL);
      assert(hapB_indicator_2 != NULL);
      gsl_ran_multinomial(r,table_size_A,1,values_A,hapA_indicator_1);
      gsl_ran_multinomial(r,table_size_B,1,values_B,hapB_indicator_1);
      gsl_ran_multinomial(r,table_size_A,1,values_A,hapA_indicator_2);
      gsl_ran_multinomial(r,table_size_B,1,values_B,hapB_indicator_2);

      unsigned int chrom_A_1, chrom_B_1;
      for(int i = 0; i < table_size_A; i++)
	{
	  if(hapA_indicator_1[i] == 1)
	    {
	      chrom_A_1 = *(guint32*)keys_A[i];
	    }
	}
      for(int i = 0; i < table_size_B; i++)
	{
	  if(hapB_indicator_1[i] == 1)
	    {
	      chrom_B_1 = *(guint32*)keys_B[i];
	    }
	}
      unsigned int recombinant_chromosome_1;
      recombinant_chromosome_1 = recom_hap(chrom_A_1,chrom_B_1,A_indicator_1,B_indicator_1);

	  
      unsigned int chrom_A_2, chrom_B_2;
      for(int i = 0; i < table_size_A; i++)
	{
	  if(hapA_indicator_2[i] == 1)
	    {
	      chrom_A_2 = *(guint32*)keys_A[i];
	    }
	}

      for(int i = 0; i < table_size_B; i++)
	{
	  if(hapB_indicator_2[i] == 1)
	    {
	      chrom_B_2 = *(guint32*)keys_B[i];
	    }
	}
#ifdef DEBUG
      printf("\n ChromA_1: %u\tChromB_1:%u\n",chrom_A_1,chrom_B_1);
      printf("\n ChromA_2: %u\tChromB_2:%u\n",chrom_A_2,chrom_B_2);
#endif
      unsigned int recombinant_chromosome_2;
      recombinant_chromosome_2 = recom_hap(chrom_A_2,chrom_B_2,A_indicator_2,B_indicator_2);
      homolog_chromosomes[j][0] = recombinant_chromosome_1;
      homolog_chromosomes[j][1] = recombinant_chromosome_2;
      j++;
      free(location_r_1);
      free(location_r_2);
      free(loci_position);
      delete_linked_list(head_r_1);
      delete_linked_list(head_r_2);
      free(hapA_indicator_1);
      free(hapB_indicator_1);
      free(hapA_indicator_2);
      free(hapB_indicator_2);
      free(values_A);
      free(values_B);
      g_free(keys_A);
      g_free(keys_B);
      /* printf("\n\n"); */
      total_markers = total_markers + current->n_loci;
      current = current->next;
    }
#ifdef DEBUG
  for(int k = 0; k < no_chromosomes; k++)
    {
      printf("%u\t%u\n",homolog_chromosomes[k][0],homolog_chromosomes[k][1]);
    }
#endif
  j = 0;
  int k = 0;
  char indv[total_markers][5];
  current = head;
  while(current != NULL && j < total_markers && k < no_chromosomes)
    {
      for(int i = 0; i < current->n_loci; i++)
  	{
  	  sprintf(indv[j],"\t%u|%u",checkSet(homolog_chromosomes[k][0],current->n_loci-i-1),checkSet(homolog_chromosomes[k][1],current->n_loci-i-1));

  	  j++;
  	}
      k++;
      current = current->next;
    }
  
  for(int i = 0; i < total_markers; i++)
    {
      strcat(simulated_output[i],indv[i]);
      /* printf("%s\n",simulated_output[i]); */
    }

  
}


/* void delete_chrom_data(chrom_data *head) */
/* { */
/*   chrom_data *temp; */
/*   while(head != NULL) */
/*     { */
/*       temp = head; */
/*       head = head->next; */
/*       free(temp->markers); */
/*       g_hash_table_destroy(temp->hashA); */
/*       g_hash_table_destroy(temp->hashB); */
/*       free(temp); */
/*     } */
/* } */




int main(int argc, char **argv)
{
  int c;
  int pflag = 0, Aflag = 0, Bflag = 0, gflag = 0;
  char *filenames[4];
  for(int i = 0; i < 4; i++)
    {
      filenames[i] = calloc(MAX_FILENAME, sizeof(char));
      assert(filenames[i] != NULL);
    }
  while((c = getopt(argc, argv, ":p:A:B:g:")) != -1)
    {
      switch(c)
	{
	case 'p':
	  strcpy(filenames[0],optarg);
	  pflag = 1;
	  break;
	case 'A':
	  strcpy(filenames[1],optarg);
	  Aflag = 1;
	  break;
	case 'B':
	  strcpy(filenames[2],optarg);
	  Bflag = 1;
	  break;
	case 'g':
	  strcpy(filenames[3],optarg);
	  gflag = 1;
	  break;
	case ':':
	  if(optopt == 'p'){
	    fprintf(stderr,"Missing argument for option -%c\n", optopt);}
	  else if(optopt == 'A'){
	    fprintf(stderr,"Missing argument for option -%c\n", optopt);}
	  else if(optopt == 'B'){
	    fprintf(stderr,"Missing argument for option -%c\n", optopt);}
	  else if(optopt == 'g'){
	    fprintf(stderr,"Missing argument for option -%c\n", optopt);}
	  for(int i = 0; i < 4; i++)
	    {
	      free(filenames[i]);
	    }
	  exit(1);
	case '?':
	  if (isprint (optopt))
	    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
	  else
	    fprintf(stderr,"Unknown option character `\\x%x'.\n",optopt);
	  for(int i = 0; i < 4; i++)
	    {
	      free(filenames[i]);
	    }
	  exit(1);
	}
    }
  
  if(pflag == 0 || Aflag == 0 || Bflag == 0 || gflag ==0 )
    {
      if(pflag == 0)
	{
	  fprintf(stderr,"%s: missing -p option\n",argv[0]);
	}
      if(Aflag == 0)
	{
	  fprintf(stderr,"%s: missing -A option\n",argv[0]);
	}
      if(Bflag == 0)
	{
	  fprintf(stderr,"%s: missing -B option\n",argv[0]);

	}
      if(gflag == 0)
	{
	  fprintf(stderr,"%s: missing -g option\n",argv[0]);
	}
      for(int i = 0; i < 4; i++)
	{
	  free(filenames[i]);
	}
      exit(1);
    }


  if(pflag)
    {
      infile_chrom = g_io_channel_new_file(filenames[0],"r",&e);
      if(infile_chrom == NULL)
	{
	  fprintf(stderr,"Could not open file %s\n", filenames[0]);
	  exit(1);
	}
    }
  if(Aflag)
    {
      infile_A = g_io_channel_new_file(filenames[1],"r",&e);
      if(infile_A == NULL)
	{
	  fprintf(stderr,"Could not open file %s\n", filenames[1]);
	  exit(1);
	}

    }
  if(Bflag)
    {
      infile_B = g_io_channel_new_file(filenames[2],"r",&e);
      if(infile_B == NULL)
	{
	  fprintf(stderr,"Could not open file %s\n", filenames[2]);
	  exit(1);
	}
    }
  if(gflag)
    {
      infile_gc = g_io_channel_new_file(filenames[3],"r",&e);
      if(infile_gc == NULL)
	{
	  fprintf(stderr,"Could not open file %s\n", filenames[3]);
	  exit(1);
	}

    }


  /* char simulated_data_filename[100]; */
  /* output_filename(filenames[1],filenames[2],filenames[0],simulated_data_filename); */
  /* printf("### %s\n",simulated_data_filename); */


  outfile = g_io_channel_new_file(argv[9],"a+",&e);
  /* outfile = g_io_channel_new_file(simulated_data_filename,"a+",&e); */
  
  const gsl_rng_type * T;
     
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r,clock());
  

  unsigned int *chrom_index, *chrom_no_loci;
  unsigned int no_chrom;
  double *chrom_length, *chrom_recom_rate;
  double *loci;
  unsigned int total_markers;
  read_chrom(filenames[0],infile_chrom,&chrom_index,&chrom_no_loci,&chrom_length,&chrom_recom_rate,&no_chrom,&loci,&total_markers);

#ifdef DEBUG
  for(int i = 0; i < total_markers; i++)
    {
      printf("%lf\t",loci[i]);
    }
#endif
  
  chrom_data *head, *current;
  head = data_linked_list_creation(chrom_index, chrom_no_loci, chrom_length, chrom_recom_rate, no_chrom, loci);
  current = head;
  unsigned int n_chromosomes = 0;

  while(current != NULL)
    {
      /* printf("chromosome number:%u\tchromosome length:%lf\trate_of_recom:%lf\tnumber_of_loci:%u\n",current->chrom_id,current->chrom_length, current->chrom_recom_rate, current->n_loci); */
      /* for(int i = 0; i < current->n_loci; i++) */
      /* 	{ */
      /* 	  printf("Marker[%d]:%lf\n",i+1,current->markers[i]); */
      /* 	} */
      n_chromosomes++;
      current = current->next;
    }

  /* printf("\nNo. of chromosomes:%u\tTotal no.of markers:%u\n",n_chromosomes,total_markers); */

  char population_name[] = {'A','B','\0'};
  read_pop(filenames[1],infile_A,head,population_name[0]);

  /* current = head; */
  /* while(current != NULL) */
  /*   { */
  /*     g_hash_table_foreach(current->hashA, (GHFunc)iterator, "Haplotype Sequence: %u has Frequency %lf \n"); */
  /*     printf("\n\n"); */
  /*     current = current->next; */
  /*   } */


  read_pop(filenames[2],infile_B,head,population_name[1]);

  /* current = head; */
  /* while(current != NULL) */
  /*   { */
  /*     g_hash_table_foreach(current->hashB, (GHFunc)iterator, "Haplotype Sequence: %u has Frequency %lf \n"); */
  /*     printf("\n\n"); */
  /*     current = current->next; */
  /*   } */

  
  /* char output[total_markers][BUFFER_SIZE]; */
  char *output[total_markers];
  for(int i = 0; i < total_markers; i++)
    {
      output[i] = malloc(BUFFER_SIZE);
      assert(output[i] != NULL);
    }

  current = head;
  int j = 0;
  while(current != NULL && j < total_markers)
    {
      for(int i = 0; i < current->n_loci; i++)
  	{
  	  sprintf(output[j],"%u:%-10u",current->chrom_id,(unsigned int)current->markers[i]);
  	  j++;

  	}
      current = current->next;
    }

  char *class;
  unsigned int no_individuals;
  read_classes(filenames[3],infile_gc,&class,&no_individuals);

#ifdef DEBUG
  printf("%s has strlen:%lu and sizeof:%lu\n",class,strlen(class),sizeof(class));
#endif
  
  /* printf("\n---------------- Simulating Chromosomes -------------------\n"); */
  for(int i = 0; i < strlen(class); i++)
    {
      if(class[i] == 'a'){simulate_model_a(head,n_chromosomes,output);}
      else if(class[i] == 'b'){simulate_model_b(head,n_chromosomes,output);}
      else if(class[i] == 'c'){simulate_model_c(head,n_chromosomes,output);}
      else if(class[i] == 'd'){simulate_model_d(head,n_chromosomes,output);}
      else if(class[i] == 'e'){simulate_model_e(head,n_chromosomes,output);}
      else if(class[i] == 'f'){simulate_model_f(head,n_chromosomes,output);}
    }

  /* printf("---------------- FINAL OUTPUT -----------------\n\n"); */

  char *header = NULL;
  header = malloc(BUFFER_SIZE);
  strcpy(header,"chrom:pos");
  char *indv_header[no_individuals];
  for(int i = 0; i < no_individuals; i++)
    {
      indv_header[i] = malloc(BUFFER_SIZE);
      
    }
  for(int i = 0; i < no_individuals; i++)
    { 
      sprintf(indv_header[i],"\ti%d",i+1);
      strcat(header,indv_header[i]);
    }
  strcat(header,"\n");



  
  /* char header[BUFFER_SIZE] = "chrom:pos"; */
  /* char indv_header[no_individuals][BUFFER_SIZE]; */
  /* /\* strcat(header,"chrom:pos"); *\/ */
  /* for(int i = 0; i < no_individuals; i++) */
  /*   { */
  /*     sprintf(indv_header[i],"\ti%d",i+1); */
  /*     strcat(header,indv_header[i]); */
  /*   } */
  /* strcat(header,"\n"); */
  /* /\* printf("%s",header); *\/ */

  for(int i = 0; i < total_markers; i++)
    {
      strcat(output[i],"\n");
      /* printf("%s",output[i]); */
    }



  
  gsize bytes_written;
  g_io_channel_write_chars(outfile,header,strlen(header),&bytes_written,&e);
  for(int i = 0; i < total_markers; i++)
    {
      g_io_channel_write_chars(outfile,output[i],strlen(output[i]),&bytes_written,&e);
    }

  /* outfile_simulation = g_io_channel_new_file(argv[10],"a+",&e); */
  /* for(int i = 0; i < strlen(class); i++) */
  /*   { */
  /*     g_io_channel_write_unichar(outfile_simulation,class[i],&e); */

  /*   } */

  for(int i = 0; i < no_individuals; i++)
    {
      free(indv_header[i]);
    }
  free(header);

  for(int i = 0; i < total_markers; i++)
    {
      free(output[i]);
    }

  
  free(class);
  delete_chrom_data(head);
  free(chrom_index);
  free(chrom_no_loci);
  free(chrom_length);
  free(chrom_recom_rate);
  free(loci);

  gsl_rng_free(r);
  g_io_channel_unref (infile_chrom);
  g_io_channel_unref (infile_A);
  g_io_channel_unref (infile_B);
  g_io_channel_unref (infile_gc);
  g_io_channel_unref (outfile);
  /* g_io_channel_unref (outfile_simulation); */
  for(int i = 0; i < 4; i++)
    {
      free(filenames[i]);
    }

  
  return 0;
}




