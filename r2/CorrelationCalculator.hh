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

extern "C"
{
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
}

class CorrelationCalculator
{

public:
  CorrelationCalculator(char *guess_vcf, char *truth_vcf, const std::string &regions,bool assume_homref);
  ~CorrelationCalculator();
  void SetAfBins(const std::string &af_bins_string);
  void printResultSummary();
  void printConfusion();
  void main();
  void SetAfTag(const std::string &af_tag) { _af_tag = af_tag; };
  void SetDsTag(const std::string &ds_tag) { _ds_tag = ds_tag; };

private:
  int _num_guess_variants;
  float *_ds_guess;
  bool _use_gt;
  bcf1_t *_line_truth,*_line_guess;
  bcf_srs_t *_bcf_reader;
  bcf_hdr_t *_hdr1, *_hdr2;
  std::vector<float> _af_bins;
  std::vector<int> _idx;
  int *_gt_guess, *_gt_truth, _num_gt_guess, _num_gt_truth;
  size_t _num_af_bins;
  std::map<int, std::vector<std::vector<std::vector<int>>>> _confusion;
  std::map<int, std::vector<std::pair<int, float>>> _avg_af;
  std::string _af_tag,_ds_tag;
  std::vector<std::vector<std::vector<int>>> _empty;
  int GetGuessDosage();
  std::vector<std::vector<int>> collapseConfusion(std::vector<std::vector<std::vector<int>>> & confusion);
  bool _assume_homref;
  size_t _dosage_resolution;
};
