# Tapestation-Gel-Annotater
Tool to annotate gel images from Tapestation

#INSTALL THESE FIRST
pip install pandas
pip install Pillow
pip install numpy
pip install glob2

#REQUIRES IN FOLDER
Tapestation Gel PNG images
1 corresponding summary file from of current plate setup
1 corresponding Tapestation sampleTable.csv of current plate setup

#TO DO/QOL
Import to excel sheet in grid layout
Fix naming - extra PNG in filename

#ERROR TROUBLESHOOTING
Number of wells from sampleTable.csv have to match number of purple wells from images combined
Number of wells must be the same between images or else it throws error -> check line 229 for np.arange stop value
Spelling of names must be the same between sampletable.csv and dict created from matrix file
NaN error can occur when 1) there are no X,Y coordinates for mean group categories from dataframe (ie 1) when no purple pixels present in last image or 2) the csv sheet does not match the number of coordinates created/vice versa
Error can occur if image is missing, leading to less wells compared to those in the list
If error occurs, delete created annotated images first as they are a copy, subsequent run of script will use copies resulting in unordered gel annotations
