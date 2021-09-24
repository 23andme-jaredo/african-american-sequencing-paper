/*  expected-correct.c calculates the expected proportion of correct genotypes ie. the mean of unphredded GQ

    Copyright (C) 2020 23andme

    Author: Jared O'Connell <jaredo@23andme.com>
*/

#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <math.h>
#include "bcftools.h"

bcf_hdr_t *in_hdr, *out_hdr;
int num_sample;
int num_gq=0;
int *f_gq=NULL;

float unphred(int pl) { return((float) 1. - pow(10., -pl / 10.)); }

const char *about(void)
{
  return "calculates the expected proportion of correct genotypes ie. the mean of unphredded GQ";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
  in_hdr  = in;
  out_hdr = out;
  num_sample =  bcf_hdr_nsamples(in_hdr);
  bcf_hdr_append(out_hdr, "##INFO=<ID=EXPECTED_CORRECT,Number=1,Type=Float,Description=\"expected proportion of correct genotypes\">");
  return 0;
}

void destroy(void)
{
  free(f_gq);
}

void kahan(float x,float *c,float *sum) {
  float y,t;
  y = x - *c;
  t = *sum + y;
  *c = (t - *sum) - y;
  *sum = t;
}

bcf1_t *process(bcf1_t *rec)
{
  int i=0,status;
  float expected_correct;
  float c,N=num_sample;

  status=bcf_get_format_int32(in_hdr, rec, "GQ", &f_gq, &num_gq);
  //  fprintf(stderr,"%d %d\n",status,num_sample);
  if (status != num_sample) error("%d must equal to %d\n", status, num_sample);

  expected_correct=0.;
  c=0.;
  for(i=0;i<status;i++) kahan(unphred(f_gq[i]),&c,&expected_correct);
  expected_correct/=N;

  bcf_update_info_float(out_hdr, rec, "EXPECTED_CORRECT", &expected_correct, 1);

  return(rec);
}
