# Bond Lab phyluce workflow


### STEP 3: Match contigs to probes
After moving your contigs.fasta files to a single folder, run match contigs to probes. We use 80% minimum coverage and 80% minimum coverage. The UCEProbes_blended_Spider2Kv1+Arachnid11Kv1-WPM2020.fasta is the blended spider and arachnid probe list- this does not change.
```
phyluce_assembly_match_contigs_to_probes \
     --contigs /your/directory/to/contigs-only/ \
     --probes /your/directory/to/probe/file/UCEProbes_blended_Spider2Kv1+Arachnid11Kv1-WPM2020.fasta \
     --output /your/directory/probe_match_UCEs_spadescaffolds_min80_rapidtest \
     --min-coverage 80 \
     --min-identity 80 \
     --log-path /your/directory/log_files
 ```
