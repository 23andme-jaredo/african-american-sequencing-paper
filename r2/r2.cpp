#include <getopt.h>
#include <math.h>

#include <utility>
#include <map>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <set>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

extern "C" {
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
}

int strsplit(const std::string &input, const char split, std::vector<std::string> &out)
{
  std::istringstream ss(input);
  out.clear();
  std::string tmp;
  int count = 0;
  while (std::getline(ss, tmp, split))    {
    out.push_back(tmp);
    count++;
  }
  return (count);
}

std::vector<float> parse_af_bins(const std::string & input) {
  std::vector<float> ret;
  std::vector<std::string> splits;
  strsplit(input,',',splits);
  for(auto it=splits.begin();it!=splits.end();it++) {
    ret.push_back(std::atof(it->c_str()));
  }
  std::sort(ret.begin(),ret.end());  
  return(ret);
}

std::string record2string(bcf_hdr_t const *header, bcf1_t *record)
{
  bcf_unpack(record, BCF_UN_ALL);
  std::stringstream ss;
  ss << bcf_hdr_id2name(header, record->rid) << ":" << record->pos + 1 << ":" << record->d.allele[0];
  for (int i = 1; i < record->n_allele; i++) {
    ss << ":" << record->d.allele[i];
  }
  return(ss.str());
}

void die(const std::string &s)
{
  std::cerr << "ERROR: " << s << "\nExiting..." << std::endl;
  exit(1);
}

double cor_from_table(std::vector< std::vector<int > > & K) {
  double n=0,x_ss=0,y_ss=0,y_mean=0,x_mean=0,xy_ss=0;
  for(size_t i=0;i<3;i++) {
    for(size_t j=0;j<3;j++) {
      double ct = K[i][j];
      x_mean += i*ct;
      y_mean += j*ct;
      n+=ct;
    }
  }
  x_mean/=n;
  y_mean/=n;
  for(size_t i=0;i<3;i++) {
    for(size_t j=0;j<3;j++) {
      double ct = K[i][j];
      x_ss +=  pow(i-x_mean,2.) * ct;
      y_ss +=  pow(j-y_mean,2.) * ct;
      xy_ss += (i-x_mean)*(j-y_mean)*ct;
    }
  }
  return(pow(xy_ss,2.)/(x_ss*y_ss));
}

void print_matrix(std::vector< std::vector<int > > & X) {
  for(auto it1=X.begin();it1!=X.end();it1++) {
    for(auto it2=it1->begin();it2!=it1->end();it2++) {
      std::cout<<*it2<<"\t";
    }
    std::cout<<std::endl;
  }
}

std::string get_intersect(bcf_hdr_t *hdr1,bcf_hdr_t *hdr2) {
  std::set<std::string> s1,s2;
  for(size_t i=0;i<bcf_hdr_nsamples(hdr1);i++) {
    s1.insert(hdr1->samples[i]);
  }
  for(size_t i=0;i<bcf_hdr_nsamples(hdr2);i++) {
    s2.insert(hdr2->samples[i]);
  }
  std::string ret="";
  for(auto it=s2.begin();it!=s2.end();it++){
    if(s1.count(*it)) {
      if(ret.size()>0) {
	ret = (*it)+","+ret;
      } else {
	ret=*it;
      }
    }
  }
  return(ret);
}

std::vector<int> match(bcf_hdr_t *hdr1,bcf_hdr_t *hdr2) {
  std::vector<int> ret(bcf_hdr_nsamples(hdr1));
  for(size_t i=0;i<bcf_hdr_nsamples(hdr1);i++) {
    ret[i]=0;
    char *query=hdr1->samples[i];
    while(ret[i]<bcf_hdr_nsamples(hdr2) && strcmp(hdr2->samples[ret[i]],query)!=0) {
      ret[i]++;
    }
    if(ret[i]==bcf_hdr_nsamples(hdr2)) ret[i]=-1;
  }
  return(ret);
}

int get_dosage(int *genotypes,int sample_index) {
  int g0=genotypes[sample_index*2];
  int g1=genotypes[sample_index*2+1];
  
  if(sample_index==-1||bcf_gt_is_missing(g0)||bcf_gt_is_missing(g1)){
    return(3);
  } else {
    return(bcf_gt_allele(g0)+bcf_gt_allele(g1));
  }
}

bcf_srs_t *build_reader(char *guess_vcf,char *truth_vcf,const std::string & regions) {
  bcf_srs_t *bcf_reader=bcf_sr_init();
  bcf_sr_set_opt(bcf_reader, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_EXACT);
  bcf_sr_set_opt(bcf_reader, BCF_SR_REQUIRE_IDX);
  bcf_sr_set_threads(bcf_reader, 4);
  if(!regions.empty()) {
    bcf_sr_set_regions(bcf_reader, regions.c_str(),0);
  }

  if(!(bcf_sr_add_reader(bcf_reader, guess_vcf))) die("problem opening "+(std::string)guess_vcf);
  if(!(bcf_sr_add_reader(bcf_reader, truth_vcf))) die("problem opening "+(std::string)truth_vcf);
  return(bcf_reader);
}

int get_af_bin(float a,std::vector<float> & x) {
  int ret=0;
  if(a==x.back()) return(x.size()-2);
  assert(x[0]==0.0);
  assert(x.back()==1.0);
  if(!(a>=0.0 && a<=1.0)) {
    return(-1);
  }
  while(a>=x[ret]) ret++;
  return(ret-1);
}

int r2_main(char *guess_vcf,char *truth_vcf,bool assume_homref,const std::string & regions,std::string af_bins_string,const std::string & af_tag)
{
  std::vector<float> af_bins= parse_af_bins(af_bins_string);
  //  for(auto it=af_bins.begin();it!=af_bins.end();it++) std::cerr<<*it<<std::endl;//debug
  auto lg = spdlog::get("stderr");
  lg->info("Guess: {}",guess_vcf);
  lg->info("Truth: {}",truth_vcf);
  bcf_srs_t *bcf_reader=build_reader(guess_vcf,truth_vcf,regions);
  bcf_hdr_t *hdr1=bcf_reader->readers[0].header;
  bcf_hdr_t *hdr2=bcf_reader->readers[1].header;
  std::string intersecting_ids = get_intersect(hdr1,hdr2);
  
  if(intersecting_ids.size()==0) {
    die("no intersecting sample ids found");
  }
  else {
    std::cerr<<std::count(intersecting_ids.begin(),intersecting_ids.end(),',')+1<<" intersecting ids found"<<std::endl;
  }
  
  bcf_sr_set_samples(bcf_reader,intersecting_ids.c_str(),0); //this doesnt seem to work.
  std::vector<int> idx=match(hdr2,hdr1);
  int *gt_guess=NULL,*gt_truth=NULL,num_gt_guess=0,num_gt_truth=0;
  
  int num_af_bins =  af_bins.size()-1;
  //  std::vector< std::vector<int> > empty(4,std::vector<int>(4,0));
  std::map< int, std::vector< std::vector< std::vector<int> > > > confusion;
  std::vector< std::vector< std::vector<int> > > empty(num_af_bins,std::vector< std::vector<int> >(4, std::vector<int>(4,0)));
  std::map< int, std::vector<std::pair<int,float> > > avg_af;

  int num_truth_variants=0,num_guess_variants=0;
  while ( bcf_sr_next_line(bcf_reader) ) {
    
    if(bcf_sr_has_line(bcf_reader,1)) { //line is in the truth data
      bcf1_t *line_truth = bcf_sr_get_line(bcf_reader,1);

      if(line_truth->n_allele!=2) {
	lg->warn("Truth at {} has >2 alleles. Ignored",record2string(hdr2,line_truth));
      }
      else {
	num_truth_variants++;
	if(bcf_sr_has_line(bcf_reader,0) || assume_homref) { //line is in the estimated data
	  if(bcf_sr_has_line(bcf_reader,0)) {
	    num_guess_variants++;
	    bcf1_t *line_guess = bcf_sr_get_line(bcf_reader,0);
	    bcf_get_genotypes(hdr1, line_guess, &gt_guess, &num_gt_guess);
	    if(line_guess->n_allele!=2) {
	      lg->error("Guess at {} has >2 alleles. This should not be possible for imputed data",record2string(hdr1,line_guess));
	    }
	  } else { //No variant in guess - punish by filling its GTs a 0/0
	    lg->warn("Guess at {} had no variant, filling as 0/0.",record2string(hdr2,line_truth));
	    num_gt_guess=2*bcf_hdr_nsamples(hdr1);
	    gt_guess=(int *)realloc(gt_guess,2*num_gt_guess*sizeof(int));
	    std::fill(gt_guess,gt_guess+2*num_gt_guess,bcf_gt_unphased(0));
	  }
	  bcf_get_genotypes(hdr2, line_truth, &gt_truth, &num_gt_truth);	  
	  int num_concordant=0;
	  int num_non_missing=0;
	  int num_af=1;
	  float af;
	  float *dummy=&af; //#htsliblife
	  if(bcf_get_info_float(hdr2,line_truth,af_tag.c_str(),&dummy,&num_af)!=1) {
	    lg->warn("Problem fetching INFO/"+af_tag+" at  {}. Skipping. ",record2string(hdr2,line_truth));
	  } else {
	    int bin=get_af_bin(af,af_bins);
	    if(bin==-1||bin>=num_af_bins) {
	      lg->error("Bad AF={} at {}",af,record2string(hdr2,line_truth));
	      exit(1);
	    }
	    int variant=bcf_get_variant_types(line_truth);
	    if(!avg_af.count(variant)) avg_af[variant].assign(num_af_bins,std::pair<int,float>(0,0.));
	    avg_af[variant][bin].first++;
	    avg_af[variant][bin].second+=af;

	    for(size_t i=0;i<bcf_hdr_nsamples(hdr2);i++) {
	      if(idx[i]!=-1) {
		int guess=get_dosage(gt_guess,idx[i]);
		int truth=get_dosage(gt_truth,i);
		if(guess != -3 && truth != -3) {
		  if(guess==truth) num_concordant++;
		  num_non_missing++;
		}
		assert(guess<4&&truth<4&&truth>=0&&guess>=0);
		if(!confusion.count(bcf_get_variant_types(line_truth))) {
		  confusion[bcf_get_variant_types(line_truth)] = empty;
		}
		assert(bin<confusion[variant].size());
		confusion[variant][bin][guess][truth]++;
	      }
	    }
	  }
	}
      }
    }
  }

  lg->info("{} out of {} true variants were in the guess set",num_guess_variants,num_truth_variants);
  //  print_matrix(confusion);

  if(confusion.count(VCF_SNP))  {
    for(size_t i=0;i<num_af_bins;i++) {
      int n = avg_af[VCF_SNP][i].first;
      std::cout<<"SNP"<<"\t"<<n<<"\t"<<avg_af[VCF_SNP][i].second/n<<"\t"<<cor_from_table(confusion[VCF_SNP][i])<<std::endl;
    }
  }
  if(confusion.count(VCF_SNP))  {
    for(size_t i=0;i<num_af_bins;i++) {
      int n = avg_af[VCF_INDEL][i].first;
      std::cout<<"INDEL"<<"\t"<<n<<"\t"<<avg_af[VCF_INDEL][i].second/n<<"\t"<<cor_from_table(confusion[VCF_INDEL][i])<<std::endl;
    }
  }

  //  if(confusion.count(VCF_INDEL))  std::cout<<"INDEL r2="<<cor_from_table(confusion[VCF_INDEL])<<std::endl;
  bcf_sr_destroy(bcf_reader);  
  free(gt_guess);
  free(gt_truth);
  return(0);
}

int usage() {
  std::cerr << "\nAbout:   Calculates Pearson's correlation between genotype dosages." << std::endl;
  std::cerr <<   "Usage:   r2 guess.bcf truth.bcf" << std::endl;
  std::cerr << "" << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "    --missing-as-homref   treat missing variants as homref call" << std::endl;
  std::cerr << "    -af-bins <list>      allele frequency bins eg. 0.0,0.1,0.5,1" << std::endl;
  std::cerr << "    -af-tag <string>     allele frequency tag to use, default is AF" << std::endl;
  return(1);
}
 
int main(int argc, char **argv) {
  if(argc<3) {
    return(usage());
  }
  static struct option loptions[] = {
    {"missing-as-homref",0,0,'m'},
    {"region",1,0,'r'},
    {"af-bins",1,0,'b'},
    {"af-tag",1,0,'a'},
    {0,0,0,0}
  };
  std::string af_tag="AF";
  std::string af_bins="0,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0";
  std::string regions="";
  int c;
  bool assume_homref=false;
  auto lg = spdlog::stderr_color_mt("stderr");
  spdlog::set_pattern(" [%c] [%l] %v");  
  while ((c = getopt_long(argc, argv, "mr:b:a:",loptions,NULL)) >= 0)  {
    switch (c)
      {
      case 'm': assume_homref=true; break;
      case 'r': regions = (optarg); break;
      case 'b': af_bins = (optarg); break;
      case 'a': af_tag = (optarg); break;
      default: 
	if(optarg!=NULL) {lg->error("Unknown argument: {}",optarg);}
	else {return(usage());};
      }
  }
  if(optind>=argc-1) return(usage());
  r2_main(argv[optind],argv[optind+1],assume_homref,regions,af_bins,af_tag.c_str());
  lg->info("Done");
  spdlog::drop_all();  
}
