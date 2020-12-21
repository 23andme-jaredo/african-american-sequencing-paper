/*  inbreeding.c returns 1-Observed[nhets]/Expected[nhets]

    Copyright (C) 2019 23andMe

    Author: Jared O'Connell <jaredo@23andme.com>
*/

#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <math.h>

bcf_hdr_t *in_hdr, *out_hdr;
float *i_inbreeding=NULL,*i_af=NULL;
int *i_achet=NULL;
int n_achet=0,n_inbreeding=0,n_af=0;

const char *about(void)
{
  return "calculates the InbreedingCoef 1-Observed[nhets]/Expected[nhets]";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
  in_hdr  = in;
  out_hdr = out;
  bcf_hdr_append(out_hdr, "##INFO=<ID=InbreedingCoef,Number=A,Type=Float,Description=\"1-Observed[nhets]/Expected[nhets]\">");
  return 0;
}

void destroy(void)
{
  free(i_achet);
  free(i_inbreeding);
  free(i_af);
}

bcf1_t *process(bcf1_t *rec)
{
  int i=0,status,an,n_an=1;
  int *i_an=&an;
  i_inbreeding=(float *)realloc(i_inbreeding,(rec->n_allele-1)*sizeof(float));

  status =  bcf_get_info_float(in_hdr,rec,"AF",&i_af,&n_af);
  assert(status==(rec->n_allele-1));
  status =  bcf_get_info_int32(in_hdr,rec,"AC_Het",&i_achet,&n_achet);
  assert(status==(rec->n_allele-1));
  status =  bcf_get_info_int32(in_hdr,rec,"AN",&i_an,&n_an);
  assert(status==1);

  for(i=1;i<rec->n_allele;i++) {
    if(i_af[i-1]>0 && i_af[i-1]<1) {
      i_inbreeding[i-1] = 1. - (float)i_achet[i-1] / ( (float)an * i_af[i-1] * (1-i_af[i-1]) );
    } else {
      i_inbreeding[i-1]=0.;
    }
    //    fprintf(stderr,"%d %f %f = %f\n",i_achet[i-1],(float)an,i_af[i-1], i_inbreeding[i-1]);     //debug
  }
  bcf_update_info_float(out_hdr, rec, "InbreedingCoef", i_inbreeding, rec->n_allele-1);
  return(rec);
}
