#include "CorrelationCalculator.hh"

void die(const std::string &s)
{
  std::cerr << "ERROR: " << s << "\nExiting..." << std::endl;
  exit(1);
}

//TODO: to be super cautious we should switch to Kahan summation here. That said, there are only a few values.
double cor_from_table(std::vector<std::vector<int>> &K, size_t dosage_resolution, bool count_missing_as_homref)
{
  auto lg = spdlog::get("stderr");
  double n = 0, x_ss = 0, y_ss = 0, y_mean = 0, x_mean = 0, xy_ss = 0;
  size_t J=dosage_resolution;
  if(count_missing_as_homref) J++;    
  for (size_t i = 0; i < 3; i++)
  {
    for (size_t j = 0; j <= J; j++)
    {
      double dosage = 2. * (double)j / (double)dosage_resolution;
      if(j==(1+dosage_resolution)) {
        dosage = 0.;
        lg->warn("Will treat {} missing genotypes in estimates as homozygous reference",K[i][j]);
      }
      double ct = K[i][j];
      x_mean += i * ct;
      y_mean += dosage * ct;
      n += ct;
    }
  }
  x_mean /= n;
  y_mean /= n;
  for (size_t i = 0; i < 3; i++)
  {
    for (size_t j = 0; j <= J; j++)
    {
      double dosage = 2. * (double)j / (double)dosage_resolution;
      if(j==(1+dosage_resolution)) dosage = 0.; 
      double ct = K[i][j];
      x_ss += pow(i - x_mean, 2.) * ct;
      y_ss += pow(dosage - y_mean, 2.) * ct;
      xy_ss += (i - x_mean) * (dosage - y_mean) * ct;
    }
  }
  return (pow(xy_ss, 2.) / (x_ss * y_ss));
}

void print_matrix(const std::vector<std::vector<int>> &X,const std::string & row_name="")
{
  for (auto it1 = X.begin(); it1 != X.end(); it1++)
  {
    std::cout << row_name << "\t";
    for (auto it2 = it1->begin(); it2 != it1->end(); it2++)
    {
      std::cout << *it2 << "\t";
    }
    std::cout << std::endl;
  }
}

std::string get_intersect(bcf_hdr_t *hdr1, bcf_hdr_t *hdr2)
{
  std::set<std::string> s1, s2;
  for (size_t i = 0; i < bcf_hdr_nsamples(hdr1); i++)
  {
    s1.insert(hdr1->samples[i]);
  }
  for (size_t i = 0; i < bcf_hdr_nsamples(hdr2); i++)
  {
    s2.insert(hdr2->samples[i]);
  }
  std::string ret = "";
  for (auto it = s2.begin(); it != s2.end(); it++)
  {
    if (s1.count(*it))
    {
      if (ret.size() > 0)
      {
        ret = (*it) + "," + ret;
      }
      else
      {
        ret = *it;
      }
    }
  }
  return (ret);
}

std::vector<int> match(bcf_hdr_t *hdr1, bcf_hdr_t *hdr2)
{
  std::vector<int> ret(bcf_hdr_nsamples(hdr1));
  for (size_t i = 0; i < bcf_hdr_nsamples(hdr1); i++)
  {
    ret[i] = 0;
    char *query = hdr1->samples[i];
    while (ret[i] < bcf_hdr_nsamples(hdr2) && strcmp(hdr2->samples[ret[i]], query) != 0)
    {
      ret[i]++;
    }
    if (ret[i] == bcf_hdr_nsamples(hdr2))
      ret[i] = -1;
  }
  return (ret);
}

int get_dosage(int *genotypes, int sample_index)
{
  int g0 = genotypes[sample_index * 2];
  int g1 = genotypes[sample_index * 2 + 1];
  auto lg = spdlog::get("stderr");
  if (sample_index<0) 
  {
    lg->error("bad sample index");
    exit(1);
  }
  else if (bcf_gt_is_missing(g0) || bcf_gt_is_missing(g1))
  {
    return (3);
  }
  else
  {
    return (bcf_gt_allele(g0) + bcf_gt_allele(g1));
  }
}

float get_dosage(float *genotypes, int sample_index)
{
  if (sample_index == -1)
  {
    return (3);
  }
  else
  {
    return (genotypes[sample_index]);
  }
}

bcf_srs_t *build_reader(char *guess_vcf, char *truth_vcf, const std::string &regions)
{
  bcf_srs_t *bcf_reader = bcf_sr_init();
  bcf_sr_set_opt(bcf_reader, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_EXACT);
  bcf_sr_set_opt(bcf_reader, BCF_SR_REQUIRE_IDX);
  bcf_sr_set_threads(bcf_reader, 4);
  if (!regions.empty())
  {
    bcf_sr_set_regions(bcf_reader, regions.c_str(), 0);
  }

  if (!(bcf_sr_add_reader(bcf_reader, guess_vcf)))
    die("problem opening " + (std::string)guess_vcf);
  if (!(bcf_sr_add_reader(bcf_reader, truth_vcf)))
    die("problem opening " + (std::string)truth_vcf);
  return (bcf_reader);
}

int get_af_bin(float a, std::vector<float> &x)
{
  int ret = 0;
  if (a == x.back())
    return (x.size() - 2);
  assert(x[0] == 0.0);
  assert(x.back() == 1.0);
  if (!(a >= 0.0 && a <= 1.0))
  {
    return (-1);
  }
  while (a >= x[ret])
    ret++;
  return (ret - 1);
}

int strsplit(const std::string &input, const char split, std::vector<std::string> &out)
{
  std::istringstream ss(input);
  out.clear();
  std::string tmp;
  int count = 0;
  while (std::getline(ss, tmp, split))
  {
    out.push_back(tmp);
    count++;
  }
  return (count);
}

std::vector<float> parse_af_bins(const std::string &input)
{
  std::vector<float> ret;
  std::vector<std::string> splits;
  strsplit(input, ',', splits);
  for (auto it = splits.begin(); it != splits.end(); it++)
  {
    ret.push_back(std::atof(it->c_str()));
  }
  std::sort(ret.begin(), ret.end());
  return (ret);
}

void CorrelationCalculator::SetAfBins(const std::string &af_bins_string)
{
  _af_bins = parse_af_bins(af_bins_string);
  _num_af_bins = _af_bins.size() - 1;
  _empty.assign(_num_af_bins, std::vector<std::vector<int>>(4, std::vector<int>(_dosage_resolution + 2, 0)));
}

std::string record2string(bcf_hdr_t const *header, bcf1_t *record)
{
  bcf_unpack(record, BCF_UN_ALL);
  std::stringstream ss;
  ss << bcf_hdr_id2name(header, record->rid) << ":" << record->pos + 1 << ":" << record->d.allele[0];
  for (int i = 1; i < record->n_allele; i++)
  {
    ss << ":" << record->d.allele[i];
  }
  return (ss.str());
}

CorrelationCalculator::CorrelationCalculator(char *guess_vcf, char *truth_vcf, const std::string &regions, bool assume_homref)
{
  _dosage_resolution = 200;
  _assume_homref = assume_homref;
  _af_tag = "AF";
  _af_bins = {0, 1};
  //  for(auto it=af_bins.begin();it!=af_bins.end();it++) std::cerr<<*it<<std::endl;//debug
  auto lg = spdlog::get("stderr");
  lg->info("Guess: {}", guess_vcf);
  lg->info("Truth: {}", truth_vcf);
  _bcf_reader = build_reader(guess_vcf, truth_vcf, regions);
  _hdr1 = _bcf_reader->readers[0].header;
  _hdr2 = _bcf_reader->readers[1].header;
  std::string intersecting_ids = get_intersect(_hdr1, _hdr2);

  if (intersecting_ids.size() == 0)
  {
    lg->error("no intersecting sample ids found");
    exit(1);
  }
  else
  {
    lg->info("{} intersecting ids found", std::count(intersecting_ids.begin(), intersecting_ids.end(), ',') + 1);
  }

  bcf_sr_set_samples(_bcf_reader, intersecting_ids.c_str(), 0); //this doesnt seem to work.
  _idx = match(_hdr2, _hdr1);
  _gt_guess = NULL;
  _ds_guess = NULL;
  _gt_truth = NULL;
  _num_gt_guess = 0;
  _num_gt_truth = 0;
}

int CorrelationCalculator::GetGuessDosage()
{
  int status;
  auto lg = spdlog::get("stderr"); //sloppy
  if (bcf_sr_has_line(_bcf_reader, 0))
  {
    _num_guess_variants++;
    bcf1_t *_line_guess = bcf_sr_get_line(_bcf_reader, 0);
    if (_use_gt)
    {
      status = bcf_get_genotypes(_hdr1, _line_guess, &_gt_guess, &_num_gt_guess);
    }
    else
    {
      status = bcf_get_format_float(_hdr1, _line_guess, _ds_tag.c_str(), &_ds_guess, &_num_gt_guess);
      if (status != bcf_hdr_nsamples(_hdr1))
      {
        lg->error("Problem fetching FORMAT/{}", _ds_tag);
        exit(1);
      }
    }

    if (_line_guess->n_allele != 2)
    {
      lg->error("Guess at {} has >2 alleles. This should not be possible for imputed data", record2string(_hdr1, _line_guess));
    }
  }
  else
  { //No variant in guess - punish by filling its GTs a 0/0
    //lg->warn("Guess at {} had no variant, filling as 0/0.", record2string(_hdr2, _line_truth));
    if (_use_gt)
    {
      _num_gt_guess = 2 * bcf_hdr_nsamples(_hdr1);
      _gt_guess = (int *)realloc(_gt_guess, 2 * _num_gt_guess * sizeof(int));
      //std::fill(_gt_guess, _gt_guess + 2 * _num_gt_guess, bcf_gt_unphased(0));
      std::fill(_gt_guess, _gt_guess + 2 * _num_gt_guess, bcf_gt_missing);
    }
    else
    {
      _num_gt_guess = bcf_hdr_nsamples(_hdr1);
      _ds_guess = (float *)realloc(_ds_guess, _num_gt_guess * sizeof(float));
      //std::fill(_ds_guess, _ds_guess + _num_gt_guess, 0.);
      std::fill(_ds_guess, _ds_guess + _num_gt_guess, 3.);
    }
  }
  return (status);
}

void CorrelationCalculator::main()
{
  auto lg = spdlog::get("stderr");
  if (_ds_tag != "GT")
  {
    _use_gt = false;
    lg->info("Using {} for dosages", _ds_tag);
  }
  else
  {
    _use_gt = true;
    lg->info("Using GT for dosages");
  }

  int num_truth_variants = 0;
  _num_guess_variants = 0;
  while (bcf_sr_next_line(_bcf_reader))
  {
    if (bcf_sr_has_line(_bcf_reader, 1))
    { //line is in the truth data
      _line_truth = bcf_sr_get_line(_bcf_reader, 1);
      if (_line_truth->n_allele != 2)
      {
        lg->warn("Truth at {} has >2 alleles. Ignored", record2string(_hdr2, _line_truth));
      }
      else
      {
        num_truth_variants++;        
        if (bcf_sr_has_line(_bcf_reader, 0) || _assume_homref)
        { //line is in the estimated data
          GetGuessDosage();          
          bcf_get_genotypes(_hdr2, _line_truth, &_gt_truth, &_num_gt_truth);
          int num_concordant = 0;
          int num_non_missing = 0;
          int num_af = 1;
          float af;
          float *dummy = &af; //#htsliblife
          if (bcf_get_info_float(_hdr2, _line_truth, _af_tag.c_str(), &dummy, &num_af) != 1 || af==bcf_float_missing || isnan(af))
          {
            lg->warn("Problem fetching INFO/" + _af_tag + " at  {}. Skipping. ", record2string(_hdr2, _line_truth));
          }
          else
          {
            int bin = get_af_bin(af, _af_bins);
            if (bin == -1 || bin >= _num_af_bins)
            {
              lg->error("Bad AF={} at {}", af, record2string(_hdr2, _line_truth));
              exit(1);
            }
            int variant = bcf_get_variant_types(_line_truth);
            if (!_avg_af.count(variant))
              _avg_af[variant].assign(_num_af_bins, std::pair<int, float>(0, 0.));
            _avg_af[variant][bin].first++;
            _avg_af[variant][bin].second += af;

            for (size_t i = 0; i < bcf_hdr_nsamples(_hdr2); i++)
            {
              if (_idx[i] != -1)
              {
                double guessf;
                if (_use_gt)
                  guessf = get_dosage(_gt_guess, _idx[i]);
                else
                  guessf = get_dosage(_ds_guess, _idx[i]);
                int guess = (guessf == 3) ? _dosage_resolution + 1 : round(_dosage_resolution * guessf / 2.);
                int truth = get_dosage(_gt_truth, i);
                if (guess < _dosage_resolution && truth < 3)
                {
                  if (guess == truth)
                    num_concordant++;
                  num_non_missing++;
                }
                assert(guess < (_dosage_resolution + 2) && truth < 4 && truth >= 0 && guess >= 0);
                if (!_confusion.count(bcf_get_variant_types(_line_truth)))
                {
                  _confusion[bcf_get_variant_types(_line_truth)] = _empty;
                }
                assert(bin < _confusion[variant].size());
                _confusion[variant][bin][truth][guess]++;
              }
            }
          }
        }
      }
    }
  }
  lg->info("{} out of {} true variants were in the guess set", _num_guess_variants, num_truth_variants);
}

void CorrelationCalculator::printResultSummary()
{
  auto lg = spdlog::get("stderr");
  if (_confusion.count(VCF_SNP))
  {
    for (size_t i = 0; i < _num_af_bins; i++)
    {
      int n = _avg_af[VCF_SNP][i].first;
      std::cout << "SNP"
                << "\t" << n << "\t" << _avg_af[VCF_SNP][i].second / n << "\t" << cor_from_table(_confusion[VCF_SNP][i], _dosage_resolution,_assume_homref) << std::endl;
    }
  }
  else
  {
    lg->info("No SNPs observed");
  }
  if (_confusion.count(VCF_INDEL))
  {
    for (size_t i = 0; i < _num_af_bins; i++)
    {
      int n = _avg_af[VCF_INDEL][i].first;
      std::cout << "INDEL"
                << "\t" << n << "\t" << _avg_af[VCF_INDEL][i].second / n << "\t" << cor_from_table(_confusion[VCF_INDEL][i], _dosage_resolution,_assume_homref) << std::endl;
    }
  }
  else
  {
    lg->info("No indels observed");
  }
}

std::vector<std::vector<int>> CorrelationCalculator::collapseConfusion(std::vector<std::vector<std::vector<int>>> & confusion) {
  std::vector<std::vector<int>> result(std::vector< std::vector<int> >(3, std::vector<int>(4, 0)));
  for(size_t i = 0; i < _num_af_bins; i++) {
    for(size_t j=0;j<3;j++) {
      for(size_t k=0;k<(_dosage_resolution+2);k++) {
        if(k==(_dosage_resolution+1))
          result[j][3]+=confusion[i][j][1+_dosage_resolution];
        else
          result[j][round(2. * (double)k/(double)_dosage_resolution)]+=confusion[i][j][k];
      }
    }
  }
  return(result);
}

void CorrelationCalculator::printConfusion()
{  
  auto lg = spdlog::get("stderr");
  if (_confusion.count(VCF_SNP)) {
    auto m=collapseConfusion(_confusion[VCF_SNP]);
    print_matrix(m,"SNP\t");
  }
  else  {
    lg->info("No SNPs observed");
  }
  if (_confusion.count(VCF_INDEL)) {
    auto m=collapseConfusion(_confusion[VCF_INDEL]);
    print_matrix(m,"INDEL\t");
  }
  else  {
    lg->info("No indels observed");
  }
}

CorrelationCalculator::~CorrelationCalculator()
{
  bcf_sr_destroy(_bcf_reader);
  free(_gt_guess);
  free(_gt_truth);
  free(_ds_guess);
}
