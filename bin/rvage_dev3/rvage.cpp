//
//  main.cpp
//  mutanni
//
//  Created by pkalbers on 25/07/2016.
//  Copyright © 2016 pkalbers. All rights reserved.
//


#include <algorithm>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <exception>
#include <stdexcept>

#include "load_map.h"
#include "load_vcf.h"
#include "load_gen.h"
#include "load_hap.h"
#include "load_bin.h"
#include "join_bin.h"
#include "load_hmm.h"
#include "load_sim.h"

#include "detect_share.h"
#include "select_share.h"
#include "load_share.h"

#include "ibd_by_pair.h"
#include "ibd_by_site.h"
#include "print_ibd.h"
#include "count_share.h"

#include "infer_age.h"

#include "Command.hpp"
#include "Redirect.hpp"
#include "Clock.hpp"


//
// Main
//
int main(int argc, const char * argv[])
{
	Command::Line line(argc, argv);
	
	// common arguments
	Command::Value< std::string >  mode("Execution mode: 'IBD', 'AGE', or 'PSM' (capitalisation ignored)");
	Command::Value< std::string >  method('m', "method", "Inference method: 'DGT', 'FGT', or 'HMM' (capitalisation ignored)");
	Command::Value< std::string >  output('o', "out", "Prefix for generated output files");
	Command::Value< size_t >       thread('t', "treads", "Number of threads for parallel execution");
	Command::Value< size_t >       buffer('b', "buffer", "Memory buffer, upper limit (megabytes)");
	Command::Value< size_t >       seed("seed", "Set seed for random operations");
	
	// input arguments
	Command::Value< std::string >    input_bin_file('i', "input", "Pre-processed binary input file");
	Command::Vector< std::string >   input_bin_join("joinInputFiles", "Join  multiple pre-processed binary input files");
	Command::Value< std::string >    input_map_file("map", "Gen map input file (optionally GZIP compressed)");
	Command::Value< double >         input_rec_rate("rec", "Recombination rate, per site per generation");
	Command::Value< std::string >    input_vcf_file("vcf", "VCF input file (optionally GZIP compressed)");
	Command::Array< std::string, 2 > input_gen_file("gen", "GEN and SAMPLE input files (optionally GZIP compressed)");
	Command::Array< std::string, 2 > input_hap_file("hap", "HAP and SAMPLE input files (optionally GZIP compressed)");
	Command::Array< std::string, 2 > input_hmm_file("hmm", "Hidden Markov Model, input files containing empirical initial state and emission probabilties (optionally GZIP compressed)");
	Command::Value< std::string >    input_sim_file("sim", "Override IBD detection by including simulated results");
	Command::Bool                    input_compress("compress", "Enable genotype compression when generating binary file from text input files");
	Command::Value< size_t >         input_lines("window", "Number of lines to be buffered while processing input formats");
	Command::Array< size_t, 2 >      input_chunk("chunk", "Chunk of chromosome to be parsed, position interval in (from, last]");
	
	// share detection arguments
	Command::Array< double, 2 >   share_f_range('f', "Allele frequency interval for rare variant sharing");
	Command::Array< size_t, 2 >   share_k_range('k', "Allele count interval for rare variant sharing");
	Command::Vector< double >     share_f_list("fList", "List of allele frequency targets to detect rare variant sharing");
	Command::Vector< size_t >     share_k_list("kList", "List of allele count targets to detect rare variant sharing");
	Command::Value< size_t >      share_max_sites("randomSites", "Maximum number of shared sites (random subset)");
	Command::Value< size_t >      share_max_pairs("randomPairs", "Maximum number of sample pairs (random subset)");
	Command::Array< size_t, 2 >   share_pos_range("region", "Position interval (inclusive) in which to detect shared rare variants");
	Command::Value< std::string > share_select("positions", "File containing target positions (any white-space separation)");
	Command::Value< size_t >      share_position("position", "Target position");
	Command::Value< std::string > share_load('x', "index", "Sharing index file; see --saveShareIndex");
	
	
	// age estimation parameters
	
	Command::Value< double > max_missing("maxMissing", "Maximum missing rate per pairwise comparison");
	Command::Value< size_t > effective_size("Ne", "Effective population size");
	Command::Value< double > mutation_rate("mut", "Mutation rate, per site per generation");
	
	//Command::Value< double > theta_manual("setTheta", "Set theta value");
	//Command::Bool            theta_watterson("estTheta", "Theta is calculated using the Watterson estimator");
	
	//Command::Value< bool > age_use_mut_clock("useMutClock", "Include mutation clock in age estimation");
	//Command::Value< bool > age_use_rec_clock("useRecClock", "Include recombination clock in age estimation");
	//Command::Value< bool > age_use_hard_brks("useHardBreaks", "Assume certain recombination breakpoints in age estimation");
	//Command::Value< bool > age_use_post_prob("usePosteriors", "Include posterior probabilities in age estimation (HMM only)");
	
	//Command::Value< bool > age_filter_fixed("filterProp", "Filter pairs: fixed proportion");
	//Command::Value< bool > age_filter_detect("filterAuto", "Filter pairs: automatic threshold detection");
	Command::Value< bool > age_nearest_neighb("selectNN", "Select nearest neighbours");
	Command::Value< bool > age_treeconsistent("countConsDiff", "Count pairwise differences only if consistent with tree");
	
	Command::Value< size_t > age_limit_sharers("maxConcordant", "Number of concordant pairs to be randomly selected");
	Command::Value< size_t > age_outgroup_size("maxDiscordant", "Number of discordant pairs to be randomly selected");
	Command::Value< size_t > age_breakpt_range("considerBreakpoints", "Number of sites to be considered from breakpoint towards focal site");
	
	
	// other parameters
	
	Command::Bool print_cle_distr("cleDistr", "Write CLE distribution to output file");
	Command::Bool print_ccf_distr("ccfDistr", "Write CCF distribution to output file");
	Command::Bool print_ibd_info("writeIBD", "Write complete IBD information to output file");
	Command::Bool local_tmp_files("localTmpFiles", "Temporary files are stored in local directory when parsing input file");
	Command::Bool save_share_index("saveShareIndex", "Save sharing index to file");
	Command::Bool silent("silent", "Prevent output to console");
	Command::Bool hmmi("hmmi", "HMM, iterative approach");
	
	
	bool mode_is_ibd = false;
	bool mode_is_age = false;
	bool mode_is_psm = false;
	
	bool method_is_dgt = false;
	bool method_is_fgt = false;
	bool method_is_hmm = false;
	bool method_is_sim = false;
	
	bool share_is_good = false;
	bool share_is_detect = false;
	bool share_is_select = false;
	bool share_is_load   = false;
	
	
	// parse command line
	try
	{
		line.parse();
		
		line.get(silent, false);
		
		
		line.get(output, true);
		
		line.get(mode, false);
		
		line.get(thread, false, size_t(1));
		line.get(buffer, false, std::numeric_limits<size_t>::max()); // default: all individuals
		line.get(seed, false);
		
		
		if (!line.get(input_bin_file, false)) // BIN
		{
			if (!line.get(input_bin_join, false)) // join BINs
			{
				// genetic map
				line.get(input_map_file, false);
				
				// constant recombination rate
				line.get(input_rec_rate, false, double(1e-08)); // default: 1e-08
				
				
				int input_n = 0;
				
				if (line.get(input_vcf_file, false)) ++input_n; // VCF
				if (line.get(input_gen_file, false)) ++input_n; // GEN
				if (line.get(input_hap_file, false)) ++input_n; // HAP
				
				if (input_n == 0)
					throw std::runtime_error("No input format provided");
				if (input_n > 1)
					throw std::runtime_error("Different input formats provided");
				
				line.get(input_compress, false, false);
				line.get(input_lines, false, size_t(5e05)); // default: 0.5 million lines
				line.get(input_chunk, false);
				
				line.get(local_tmp_files, false, false);
			}
		}
		
		
		if (!method_is_sim)
		{
			const bool share_f = line.get(share_f_range, false);
			const bool share_k = line.get(share_k_range, false);
			
			if (share_f && share_k)
				throw std::invalid_argument("Different definitions provided for detection of sharing");
			
			if (share_f || share_k)
			{
				share_is_good   = true;
				share_is_detect = true;
			}
			else
			{
				const bool share_fl = line.get(share_f_list, false);
				const bool share_kl = line.get(share_k_list, false);
				
				if (share_fl && share_kl)
					throw std::invalid_argument("Different lists provided for detection of sharing");
				
				if (share_fl || share_kl)
				{
					share_is_good   = true;
					share_is_detect = true;
				}
				else
				{
					share_is_select = line.get(share_select, false);
					
					if (!share_is_select)
					{
						share_is_select = line.get(share_position, false);
					}
					
					share_is_load = line.get(share_load, false);
					
					if (share_is_select && share_is_load)
						throw std::invalid_argument("Cannot load sharing index when target positions are provided");
					
					if (share_is_select || share_is_load)
						share_is_good = true;
				}
			}
			
			if (share_is_good)
			{
				line.get(share_max_sites, false, size_t(0));
				line.get(share_max_pairs, false, size_t(0));
				line.get(share_pos_range, false);
				
				line.get(save_share_index, false, false);
			}
		}
		
		
		if (mode.good())
		{
			if (mode.value == "IBD" || mode.value == "ibd") mode_is_ibd = true; else
			if (mode.value == "AGE" || mode.value == "age") mode_is_age = true; else
			if (mode.value == "PSM" || mode.value == "psm") mode_is_psm = true; else
				throw std::invalid_argument("Unknown mode: " + mode.value);
			
			if (!mode_is_psm)
			{
				line.get(method, true);
				
				if (method.value == "DGT" || method.value == "dgt") method_is_dgt = true; else
				if (method.value == "FGT" || method.value == "fgt") method_is_fgt = true; else
				if (method.value == "HMM" || method.value == "hmm") method_is_hmm = true; else
				if (method.value == "SIM" || method.value == "sim") method_is_sim = true; else
					throw std::invalid_argument("Unknown method: " + method.value);
			}
			
			
			if (method_is_hmm)
			{
				line.get(input_hmm_file, false);
				line.get(hmmi, false);
			}
			
			if (method_is_sim)
			{
				line.get(input_sim_file, true);
				
				if (!mode_is_age)
					throw std::invalid_argument("Mode must be 'AGE' when method is 'SIM'");
			}
			
			if (mode_is_ibd)
			{
				line.get(print_ibd_info, false);
			}
			
			
			if (!share_is_good && !method_is_sim)
				throw std::invalid_argument("No sharing definition provided");
			

			line.get(effective_size, false, size_t(10000)); // default: 10000
			line.get(mutation_rate, false, double(1e-08)); // default: 1e-08
			line.get(max_missing, false, double(0.05)); // defaul: 5%
			
			//line.get(theta_watterson, false);
			//line.get(theta_manual, false);
			
			//line.get(age_use_hard_brks, false, false);
			//line.get(age_use_mut_clock, false, true);
			//line.get(age_use_rec_clock, false, true);
			
			//line.get(age_filter_fixed, false, false);
			//line.get(age_filter_detect, false, false);
			line.get(age_nearest_neighb, false, true);
			line.get(age_treeconsistent, false, true);
			
			//if (!method_is_sim)
			//	line.get(age_use_post_prob, false, true);
			
			line.get(age_limit_sharers, false);
			line.get(age_outgroup_size, false);
			line.get(age_breakpt_range, false);
			
			line.get(print_cle_distr, false, false);
			line.get(print_ccf_distr, false, false);
		}

		line.finish();
	}
	catch(const int & help)
	{
		return EXIT_SUCCESS;
	}
	catch(const std::exception & error)
	{
		std::cout << error.what() << std::endl;
		return EXIT_FAILURE;
	}
	
	
	// set seed
	if (seed.good())
	{
		set_random_seed(seed);
	}
	
	
	// redirect log/err output
	Redirect redirect_log(std::clog, output.value + ".log");
	Redirect redirect_err(std::cerr, output.value + ".err");
	
	if (silent)
	{
		Redirect redirect_out(std::cout, nullptr);
	}
	
	
	std::cout << std::endl;
	std::cout << "====================================================" << std::endl;
	std::cout << "             Rare variant age estimation            " << std::endl;
	std::cout << "====================================================" << std::endl;
	std::cout << ":: #######:: ##::: ##:::: ##::::: ######:: #######::" << std::endl;
	std::cout << ":: ##... ##: ##::: ##:: ##. ##:: ##....::: ##....:::" << std::endl;
	std::cout << ":: #######:: ##::: ##: ##::. ##: ##:: ###: #####::::" << std::endl;
	std::cout << ":: ##. ##:::. ##: ##:: ########: ##::. ##: ##..:::::" << std::endl;
	std::cout << ":: ##::. ##:::. ##:::: ##::: ##:. ######:: #######::" << std::endl;
	std::cout << "::..::::..:::::..:::::..::::..:::......:::.......:::" << std::endl;
	std::cout << "====================================================" << std::endl;
	std::cout << "                  Patrick K. Albers                 " << std::endl;
	std::cout << "====================================================" << std::endl;
	std::cout << std::endl;
	
	
	// print programme call to log
	line.print(std::clog);
	
	
	// begin timer
	Clock runtime;
	
	
	// main algorithm
	try
	{
		// process input data
		
		Gen::Grid::Data grid;
		
		if (input_bin_file.good())
		{
			grid = load_bin(input_bin_file);
		}
		else if (input_bin_join.good())
		{
			input_bin_file.value = join_bin(output, input_bin_join);
			
			grid = load_bin(input_bin_file);
			
			grid->print_sample(output.value + ".sample.txt");
			grid->print_marker(output.value + ".marker.txt");
		}
		else
		{
			Gen::Map gmap = (input_map_file.good()) ? load_map(input_map_file): load_map(input_rec_rate);
			
			if (!input_chunk.good())
			{
				input_chunk.value[0] = 0;
				input_chunk.value[1] = 0;
			}
			
			if (input_vcf_file.good())
				input_bin_file = load_vcf(input_vcf_file, gmap, output, input_lines, input_compress, local_tmp_files, input_chunk[0], input_chunk[1]);
			
			if (input_gen_file.good())
				input_bin_file = load_gen(input_gen_file[0], input_gen_file[1], gmap, output, input_lines, input_compress, local_tmp_files, input_chunk[0], input_chunk[1]);
			
			if (input_hap_file.good())
				input_bin_file = load_hap(input_hap_file[0], input_hap_file[1], gmap, output, input_lines, input_compress, local_tmp_files, input_chunk[0], input_chunk[1]);
			
			grid = load_bin(input_bin_file);
			
			grid->print_sample(output.value + ".sample.txt");
			grid->print_marker(output.value + ".marker.txt");
		}
		
		grid->cache(buffer);

		
		// execute pairwise sharing matrix
		if (mode_is_psm)
		{
			if (share_f_range.good())
			{
				const double size = static_cast<double>(grid->sample_size() * 2);
				share_k_range.value[0] = static_cast<size_t>(std::round(share_f_range[0] * size));
				share_k_range.value[1] = static_cast<size_t>(std::round(share_f_range[1] * size));
			}
			
			if (share_k_range.good() || share_f_range.good())
			{
				count_share(output, grid, share_k_range[0], share_k_range[1]);
			}
			
			if (share_f_list.good())
			{
				const double size = static_cast<double>(grid->sample_size() * 2);
				
				for (size_t i = 0, n = share_f_list.value.size(); i < n; ++i)
				{
					share_k_list.value.push_back( static_cast<size_t>(std::round(share_f_list[i] * size)) );
				}
			}
			
			if (share_k_list.good() || share_f_list.good())
			{
				count_share(output, grid, share_k_list);
			}
		}
		else
		if (share_is_good)
		{
			Gen::Share::Data share = std::make_shared<Gen::Share>();
			
			
			if (share_pos_range.good())
			{
				share->region.first  = share_pos_range.value[0];
				share->region.second = share_pos_range.value[1];
			}
			
			if (share_max_sites.good())
			{
				share->max_sites = share_max_sites.value;
			}
			
			if (share_max_pairs.good())
			{
				share->max_pairs = share_max_pairs.value;
			}
			
			
			if (share_is_detect)
			{
				if (share_f_range.good())
				{
					const double size = static_cast<double>(grid->sample_size() * 2);
					share_k_range.value[0] = static_cast<size_t>(std::round(share_f_range[0] * size));
					share_k_range.value[1] = static_cast<size_t>(std::round(share_f_range[1] * size));
				}
				
				if (share_k_range.good() || share_f_range.good())
				{
					share_is_good = detect_share(share, grid, share_k_range[0], share_k_range[1]);
				}
				
				if (share_f_list.good())
				{
					const double size = static_cast<double>(grid->sample_size() * 2);
					
					for (size_t i = 0, n = share_f_list.value.size(); i < n; ++i)
					{
						share_k_list.value.push_back( static_cast<size_t>(std::round(share_f_list[i] * size)) );
					}
				}
				
				if (share_k_list.good() || share_f_list.good())
				{
					share_is_good = detect_share(share, grid, share_k_list);
				}
			}
			
			if (share_is_select)
			{
				if (share_select.good())
					share_is_good = select_share(share, grid, share_select);
				
				if (share_position.good())
					share_is_good = select_share(share, grid, share_position);
			}
			
			if (share_is_load)
			{
				share_is_good = load_share(share_load, share, grid);
			}
			
			
			if (!share_is_good)
			{
				throw std::runtime_error("Sharing table is empty");
			}
			
			
			if (save_share_index.good())
			{
				save_share(share, output, grid);
			}
			
			
			if (mode.good())
			{
				IBD::DetectMethod method = IBD::DETECT_VOID;
				
				if (method_is_dgt) method = IBD::DETECT_DGT;
				if (method_is_fgt) method = IBD::DETECT_FGT;
				if (method_is_hmm) method = IBD::DETECT_HMM;
				if (method_is_sim) method = IBD::DETECT_SIM;
				
				
				// prepare hmm
				
				IBD::HMM::Model::Data hmm_model;
				
				if (method_is_hmm)
				{
					if (input_hmm_file.good())
						hmm_model = load_hmm(grid, input_hmm_file[0], input_hmm_file[1], effective_size, output);
					else
						hmm_model = load_hmm(grid, effective_size);
					
					if (hmmi.good())
						hmm_model->do_iterative = true;
				}
				
				
				// execute IBD detection only
				
				if (mode_is_ibd)
				{
					if (print_ibd_info.good())
					{
						print_ibd(method, output, share, grid, hmm_model);
					}
					else
					{
						if (share_select.good())
							ibd_by_site(method, max_missing, output, share, grid, runtime, hmm_model, thread);
						else
							ibd_by_pair(method, max_missing, output, share, grid, runtime, hmm_model, thread);
					}
				}
				
				
				// execute age estimation
				
				if (mode_is_age)
				{
					Age::Param::Data param = std::make_shared< Age::Param >(grid, effective_size, mutation_rate);
					
					//if (age_use_hard_brks.good())
					//	param->use_hard_brks = age_use_hard_brks.value;
					
					//if (age_use_mut_clock.good())
					//	param->use_mut_clock = age_use_mut_clock.value;
					
					//if (age_use_rec_clock.good())
					//	param->use_rec_clock = age_use_rec_clock.value;
					
					//if (age_use_post_prob.good())
					//	param->use_post_prob = age_use_post_prob.value;
					
					if (age_limit_sharers.good())
						param->limit_sharers = age_limit_sharers.value;
					
					if (age_outgroup_size.good())
						param->outgroup_size = age_outgroup_size.value;
					
					if (age_breakpt_range.good())
						param->breakpt_range = age_breakpt_range.value;
					
					//if (theta_manual.good())
					//	param->set_theta(theta_manual.value);
					
					//if (theta_watterson.good())
					//	param->est_theta(grid);
					
//					if (age_filter_fixed.good())
//						param->apply_filter_fixed = age_filter_fixed.value;
					
//					if (age_filter_detect.good())
//						param->apply_filter_detect = age_filter_detect.value;
					
					if (age_nearest_neighb.good())
						param->apply_nearest_neighb = age_nearest_neighb.value;
					
					if (age_treeconsistent.good())
						param->use_tree_consistency = age_treeconsistent.value;
					
					param->threads = thread;
					
					
					infer_age(param, method, max_missing, output, print_cle_distr, print_ccf_distr, share, grid, runtime, hmm_model, nullptr, thread);
				}
			}
		}
		
		
		if (method_is_sim && mode_is_age)
		{
			std::cout << "IBD detection override" << std::endl << std::endl;
			std::clog << "IBD detection override" << std::endl << std::endl;
			
			// load simulation results
			IBD::SIM::Result::Data sim_result = load_sim(input_sim_file.value);
			
			
			Age::Param::Data param = std::make_shared< Age::Param >(grid, effective_size, mutation_rate);
			IBD::DetectMethod method = IBD::DETECT_SIM;
			
			//if (age_use_hard_brks.good())
			//	param->use_hard_brks = age_use_hard_brks.value;
			
			//if (age_use_mut_clock.good())
			//	param->use_mut_clock = age_use_mut_clock.value;
			
			//if (age_use_rec_clock.good())
			//	param->use_rec_clock = age_use_rec_clock.value;
			
			if (age_breakpt_range.good())
				param->breakpt_range = age_breakpt_range.value;
			
			//if (theta_manual.good())
			//	param->set_theta(theta_manual.value);
			
			//if (theta_watterson.good())
			//	param->est_theta(grid);
			
			if (age_treeconsistent.good())
				param->use_tree_consistency = age_treeconsistent.value;
			
			
			infer_age(param, method, max_missing, output, print_cle_distr, print_ccf_distr, nullptr, grid, runtime, nullptr, sim_result, thread);
		}
	}
	catch(const std::exception & error)
	{
		std::cout << error.what() << std::endl;
		return EXIT_FAILURE;
	}
	
	
	std::cout << "This is the end!" << std::endl;
	std::clog << "This is the end!" << std::endl;
	runtime.print(std::cout);
	runtime.print(std::clog);
	
	
	return EXIT_SUCCESS;
}

