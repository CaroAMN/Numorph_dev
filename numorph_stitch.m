function stitch = stitching(sample_name)
NM_setup
sample = sample_name;
config = NM_config('process',sample);
NM_process(config, "stitch",sample);