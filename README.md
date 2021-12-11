Computational DNA storage simulator serves as a tool to conduct trial experiments for deciding suitable encoding and decoding (with/without clustering) mechanisms and also assessing the requirement of logical and sequencing redundancy thereby decoding the data back with minimum error.

Please follow the procedure to integrate our codes with existing biological sequencer- Deepsimulator that mimics the process of nanopore sequencing.

1. Download deepsimulator from https://github.com/liyu95/DeepSimulator
2. Repace the file deep_simulator.sh with the file dna_storage_channel_simulator.sh present here.
3. Add the three files - case2_beta_non-parametric.py, case2_log-gamma_non-parametric.py, case2_normal_non_parametric.py inside util/
4. Execute the script deepsimulator.sh with two extra parameters d for specifying distribution  (1-normal,2-beta,3-loggamma) and Y for specifying coverage less than or equal to 60x
5. –n  –1 is to avoid the splitting of reads. 
6. e, f, g are the other parameters for specifying noise, frequency,and text signal file respectively. It can be varied based on our requirement of basecalling accuracy and need of the output signal file.
              
              E.g.: / dna_storage_channel_simulator.sh –i input_677reads.txt –n  –1 –d 1 –Y 60 

