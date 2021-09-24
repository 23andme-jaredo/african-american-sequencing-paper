#include <getopt.h>
#include "CorrelationCalculator.hh"

int usage()
{
  std::cerr << "\nAbout:   Calculates Pearson's correlation between genotype dosages." << std::endl;
  std::cerr << "Usage:   r2 guess.bcf truth.bcf" << std::endl;
  std::cerr << "" << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "    --missing-as-homref   treat missing variants as homref call" << std::endl;
  std::cerr << "    --confusion           output confusion matrices" << std::endl;
  std::cerr << "    --af-bins <list>      allele frequency bins eg. 0,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0" << std::endl;
  std::cerr << "    --af-tag <string>     allele frequency tag to use, default is AF" << std::endl;
  std::cerr << "    --dosage <string>     FORMAT value to use as genotype dosage (defaults to GT)" << std::endl;
  return (1);
}

int main(int argc, char **argv)
{
  if (argc < 3)
  {
    return (usage());
  }
  static struct option loptions[] = {
      {"missing-as-homref", 0, 0, 'm'},
      {"confusion", 0, 0, 'c'},
      {"region", 1, 0, 'r'},
      {"af-bins", 1, 0, 'b'},
      {"dosage", 1, 0, 'd'},
      {"af-tag", 1, 0, 'a'},
      {0, 0, 0, 0}};
  std::string ds_tag = "GT";
  std::string af_tag = "AF";
  std::string af_bins = "0,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0";
  std::string regions = "";
  int c;
  bool assume_homref = false;
  bool confusion = false;
  auto lg = spdlog::stderr_color_mt("stderr");
  spdlog::set_pattern(" [%c] [%l] %v");
  while ((c = getopt_long(argc, argv, "cmr:b:a:d:", loptions, NULL)) >= 0)
  {
    switch (c)
    {
    case 'm':
      assume_homref = true;
      break;
    case 'c':
      confusion = true;
      break;
    case 'r':
      regions = (optarg);
      break;
    case 'd':
      ds_tag = (optarg);
      break;
    case 'b':
      af_bins = (optarg);
      break;
    case 'a':
      af_tag = (optarg);
      break;
    default:
      if (optarg != NULL)
      {
        lg->error("Unknown argument: {}", optarg);
      }
      else
      {
        return (usage());
      };
    }
  }
  if (optind >= argc - 1)
    return (usage());
  CorrelationCalculator r2(argv[optind], argv[optind + 1], regions, assume_homref);
  r2.SetAfTag(af_tag);
  r2.SetDsTag(ds_tag);
  r2.SetAfBins(af_bins);
  r2.main();
  if(confusion)   r2.printConfusion();
  else r2.printResultSummary();
  lg->info("Done");
  spdlog::drop_all();
}
