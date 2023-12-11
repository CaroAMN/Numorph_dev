function align = numorph_align(sample_name, output_directory)
NM_setup;
sample = sample_name;

config = NM_config('process',sample);
NM_process(config, "align",sample);