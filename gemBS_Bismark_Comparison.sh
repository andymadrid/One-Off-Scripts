GemBS comparison to Bismark
# general set up for GemBS run
# working with one human Chiari sample to start and compare

# get into the working directory
cd /media/Data/WGBS/Chiari_Benny/testGEMBS

# prepare the reference
### edit the example files as seen fit
gemBS prepare -c example.conf -t example.csv 

# index the reference
gemBS index

# run the remaining pipeline
gemBS map
gemBS call
gemBS extract
gemBS map-report
gemBS call-report

# need to finish this at some point
