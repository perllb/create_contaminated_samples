# create_contaminated_samples
Simulate cross-sample-contamination by mixing fastq reads from two samples

### 1. Edit config.txt 
- Set full path to fastq files
- Sample A: Original sample
- Sample B: Contaminating sample
- If simulating 5% contamination, 95% of reads will be taken from sample A, and 5% from sample B.

### 2. run simulation
Simulate 5% contamination, with total read depth of 80M in output fastq.
```
bash simulate_contamination.sh -p 5 -r 80
```

