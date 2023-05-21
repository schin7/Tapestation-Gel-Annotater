import os
import shutil
from PIL import Image, ImageDraw, ImageFont
import numpy as np
import pandas as pd
import glob

def target_primer_db(target):
    '''Contains list of Targets + their Primers'''
    
    target1 = target
    tg1 = {'StartDB':['Start'],\
           'PrimerTarget':['SetName1','SetName2'],\
           'EndDB':['End']}
    
    return tg1[target1]


def createCopies():
    '''create copies of gel png files to work on, to avoid overwriting main files'''
    src = os.getcwd()
    for item in os.listdir(src):
        s = os.path.join(src, item)
        if s.endswith(".png"):
            shutil.copy(s, os.path.join(src, str(s) + "_annotated.png"))    
            
            #ISSUE - insert if exist statement so that rerunning script more than once doesnt create another set of pngs



def readSUM():    
    '''read matrix csv and delete columns which arent needed in nested dictionary'''
    SEQ_SUM = glob.glob("*Matrix.SEQ_SUM")
    seq = SEQ_SUM[0]
    dfSUM = pd.read_csv(seq,sep='\t', header=(0))
    #del dfSUM['Date (MM-DD-YY)']
    #del dfSUM['Time (HH:MM)']
    #del dfSUM['Version']
    #del dfSUM['Unnamed: 17'] #extra column created in dataframe
    dfSUM = dfSUM.drop(['Date (MM-DD-YY)','Time (HH:MM)','Version','Unnamed: 17'], axis=1, errors='ignore')

    #use regex to remove unneeded characters in headers - make clean Target IDS to parse later
    dfSUM.columns = dfSUM.columns.str.replace(r'[()"=]','').str.replace(r'\b(\d*CHAR10\d*)\b','').str.replace(r'&&','-')
    #create dynamic nested dictionary for sample ID -> targets -> call (key, key, value)
    #check if headers match

    d = dict(dfSUM.set_index('Sample ID').groupby(level = 0).apply(lambda x : x.to_dict(orient= 'list')))
    
    return d


def readTapestationList():
    '''read tapestation summary table for ordered list of wells:samples'''
    tapeTable = glob.glob("*_sampleTable.csv")
    plateOrder = tapeTable[0]
    #set columns for gel summary csv
    #col_list = ["Well", "Sample Description"]
    col_list = [0, 2]
    dfWells = pd.read_csv(plateOrder, usecols=col_list) 
    orderedList = dfWells.iloc[:,1]
    
    return orderedList

def f(x):
  '''drop NaN value pairs from dict created by using groupby function in imageCalc()'''
  y = x.dropna()
  return np.nan if y.empty else y

def imageCalc():
    '''create empty lists to store data from image calculations 
    data - big array of X,Y coordinates that are appended to in each image file iteration.  
    Array used to: 1) to determine number of purple wells on image and have a count for when to change image, files do not provide this information'''
    # 2) catch images that can be differently sized. Expected to see repeated coordinates in dataframe, unless image is different
    
    #initialize arrays
    data = pd.DataFrame([])
    xmax = []
    fnames = []
    
    images = glob.glob("*_annotated.png")
    for image in images:
        with open(image, 'rb') as file:
            fnames.append(image)
            #convert/ensure file is RGB
            imGel = Image.open(file).convert('RGB')
            #Look for all pixels that are purple
            im  = np.array(imGel)
            purple = [128,0,128]
            #numpy.where() method to retrieve a tuple indices of two arrays where the first array contains the x-coordinates of the pixels of color and the second array contains the y-coordinates
            indices  = np.where(np.all(im==purple,axis=2))
            #zip() method to get a list of tuples containing those points
            L = list(zip(indices[0], indices[1]))
            df=pd.DataFrame(L, columns = ['X','Y'])
            #rearrange columns in dataframe for correct coordinates (Y, X instead of X,Y for pixel drawing) - pillow uses different coordinate system
            cols = df.columns.tolist()
            cols = cols[-1:] + cols[:-1]
            df = df[cols]
            #Purple Well dimensions 62 pixels wide, 7 pixels high, 11 pixels in between each, 134 is center line of well - Add buffer pixels for placement of text - Pillow placement slightly off from actual pixel location
            #Groupby function used on Y coordinates as X are all relatively the same between 131-137. Used to categorize all dataframe points into ranges representing each well - take the mean of each well to get the center point along Y
            dfGroups = df.groupby(pd.cut(df["Y"],np.arange(81, 1350, 73))).mean().agg(f)
            #ADJUST LINE 229 ARANGE STOP VALUE IF WELL NUMBER IN IMAGE IS DIFFERENT
            #count the number of categories of max number of wells onp/n an image
            xcount, ycount = dfGroups.count()
            xmax.append(int(xcount))
            data = data.append(dfGroups, ignore_index=True)
        del imGel
        
    return data, xmax, fnames


def annotateImages(meanRange, maxWells, pngNames, sampList):
    '''read each line in the summary csv, if line = ladder skip, else split names and run through dictionary
    #obtain sample ID and target ID for nested dictioanry to get calls
    #draw texts for both calls and sample names both'''   
   
    #initialize counters
    dfcounter = 0
    LadderBuffer = 0
    wellCounter = 0
    maxcounter = 0  
    fcount= 0
    
    #read each line in sample table and split into target id and sample id
    for each_line in sampList:
        new_sample_id = ''
        Sample_ID = ''
        targetPrimer = ''
        sample_split = Sample_ID.split('-') 
        primer_finding_counter = 0 
        
        if each_line == 'Ladder':
            call_value = 'Ladder'
            new_sample_id = ''
            LadderBuffer = LadderBuffer + 75   #buffer x-axis for when a ladder is present 
            
        else:
            split_name =each_line.split('-')
            gel_target_name = split_name[0]
            if gel_target_name == 'Ebien' or gel_target_name == 'Enceph':
                gel_target_name = 'Micro'
            elif (gel_target_name == "tcdA") or (gel_target_name == "tcdB"):
                gel_target_name = "Cdiff"
            elif gel_target_name == "ELT":
                gel_target_name = "ETEC"
            elif gel_target_name == "NoroGII":
                gel_target_name = "Noro"
            elif gel_target_name == "NoroGI":
                gel_target_name = "Noro"
            elif gel_target_name == "SapoGI":
                gel_target_name = "Sapo"
            elif gel_target_name == "SapoGII":
                gel_target_name = "Sapo"
            elif gel_target_name == "SapoGIV":
                gel_target_name = "Sapo"
            elif gel_target_name == "SapoGV":
                gel_target_name = "Sapo" 
                
            Primer_list = target_primer_db(gel_target_name)    
            if gel_target_name == "Sapo":
                sample_split_rejoin = '-'.join(sample_split[0::]).upper()
            else:
                sample_split_rejoin = '-'.join(split_name[1::]).upper() 
          
            for i in range(len(Primer_list)):
                each_primer_item = Primer_list[i].upper()
                if (sample_split_rejoin.find(each_primer_item+"-")>=0) == True:
                    new_sample_id = sample_split_rejoin.replace(each_primer_item+'-','') 
                    new_sample_id1 = new_sample_id.split('-')
                    new_sample_id1.pop(-1)
                    Sample_ID = '-'.join(new_sample_id1)
                    Primer = Primer_list[i]
                    Target_ID = gel_target_name + "-" + Primer
                    targetPrimer = Target_ID #
                    primer_finding_counter +=1 
                    sample_result = seqcallDict[new_sample_id][Target_ID]
                    call_value = "".join(sample_result)  
                    #print "Call for "+ new_sample_id + " and " + Target_ID + " is: "+ call_value  

                #if primer_finding_counter == 0:
                    #print "Unknown primer detected"
                #if (new_sample_id == '') and (Sample_ID == ''):
                    #Primer = sample_split[1] 
                    #sample_split.pop(0) 
                    #sample_split.pop(-1) 
                    #Sample_ID = '-'.join(sample_split)
       
        wellCounter = wellCounter + 1   
        maxwellcount = maxWells[maxcounter]    
        pngFile = pngNames[fcount]
        
        #print "wellcounter: " + str(wellCounter)
        #print "maxcounter: " + str(maxcounter)
        #print "fcount: " + str(fcount)
        
        if wellCounter == maxwellcount:
            #counters to change to next file being analyzed in pillow
            maxcounter = maxcounter + 1        
            wellCounter = 0            
            fcount = fcount + 1            
            imGel = Image.open(pngFile).convert('RGB')  
        else:
            #print new_sample_id
            imGel = Image.open(pngFile).convert('RGB')            
        
        #create images in pillow
        font = ImageFont.truetype('C:\Windows\Fonts\calibri.ttf', 11)
        swidth, sheight = font.getsize(new_sample_id)
        tpwidth, tpheight = font.getsize(targetPrimer)
        image2 = Image.new('RGBA', (swidth, sheight), (0, 0, 0, 0)) #empty transparent canvas
        image3 = Image.new('RGBA', (tpwidth, tpheight), (0, 0, 0, 0)) # empty transparent canvas
        
        #d1 = sequencing call, d2 = sample id, d3 = target + primer id
        d1 = ImageDraw.Draw(imGel)
        d2 = ImageDraw.Draw(image2)
        d3 = ImageDraw.Draw(image3)
        
        #create a separate image and rotate text/entire canvas to be placed on original - cannot rotate text on original without rotating the canvas, will rotate gel
        d2.text((0, 0), new_sample_id, font=font, fill=(105,105,105))  
        d3.text((0, 0), targetPrimer, font=font, fill=(105,105,105))  
        image2 = image2.rotate(90, resample=Image.BICUBIC, expand=1)
        image3 = image3.rotate(90, resample=Image.BICUBIC, expand=1)
        
        #get X, Y coordinates and iterate through with a counter
        X, Y = meanRange.iloc[dfcounter] 
        
        
        
        print "Annotating: " + str(each_line)

        #if X == "NaN":
            #break

        #buffer pixels added to X, Y coordinates to account for Pillow coordinate system        
        if call_value == "POS":
            d1.text((int(X)-10,int(Y)-104), call_value, font=font, fill=(0,255,0)) #green
        elif call_value == "NEG":
            d1.text((int(X)-10,int(Y)-104), call_value, font=font, fill=(255,0,0)) #red
        elif call_value == "Poor POS":
            d1.text((int(X)-10,int(Y)-104), call_value, font=font, fill=(255,255,0)) #yellow
        else: #invalids, NT, checks, multiple calls
            d1.text((int(X)-10,int(Y)-104), call_value, font=font, fill=(0,0,0)) #black
            
        
        imGel.paste(image2, (int(X)+int(LadderBuffer),int(Y)+10), image2) 
        imGel.paste(image3, (int(X)-13+int(LadderBuffer),int(Y)+10), image3)   
        
        imGel.save(pngFile,'PNG',quality=95, subsampling=0, optimize=True)  
        
        #delete images after every loop for variable to clear
        imGel.close()    
        del d1, d2, d3
        LadderBuffer = 0
        dfcounter = dfcounter + 1

print "--------------------------------------------------------------------------------"
print "Analyzing Tapestation Images..."
print "--------------------------------------------------------------------------------"
print " "

#run functions
createCopies()
seqcallDict = readSUM()
meanRange, maxWells, pngNames = imageCalc()
sampList = readTapestationList()
annotateImages(meanRange, maxWells, pngNames, sampList)

print " "
print '----------------------------'
print 'The analysis is now complete'
print '----------------------------'
print " "

raw_input("Press enter to exit ")
