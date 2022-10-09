from dedalus.tools import post
post.merge_process_files('analysis',cleanup=True)
post.merge_process_files('checkpoint',cleanup=True)
print('Done of merge process, go to download the data with original slurm number')