# -*- coding: utf-8 -*-
#%%
if __name__ == '__main__':
    import subprocess
    import sys
    import Bio
    import pandas as pd
    import os
    from configparser import ConfigParser
    
    os.chdir("C:/Users/James_Mann1/Desktop/TestBench")
    print('-------------------------------------------')
    print('--------------PrimedSherlock---------------')
    print('-------------Mann et al 2020---------------')
    print('-------------------------------------------')
    print('------------Major Dependencies-------------')
    print('------- Primed RPA : Higgins M et al ------')
    print('--------Cas-Offinder : Bae S et al---------')
    print('-------------------------------------------')
    print('----If First Time Use Check User Guide-----')
    print('-------------------------------------------')
    print('-------------------------------------------\n\n')
    
    master_dir = os.getcwd()
    ##############################################################################
    #This section goes over the config file and sets up the file for user editing
    ##############################################################################
    try:
        config = ConfigParser()
        config.read('config.ini')
        theta = config.get('main','Amplicon_Max_Size')
        print("[Primed Sherlock] Found valid Config File")
    except:
        print("Unable to find config.ini, please fill out newly created config.ini")
        config.read('config.ini')
        config.add_section('main')
        config.set('main', 'Amplicon_Max_Size', '375')
        config.set('main', 'Amplicon_Min_Size', '150')
        config.set('main', 'Max_Background_Identitiy', '80')
        config.set('main', 'Off_Target_Mismatch_Score', '7')
        config.set('main', 'Working_Directory', 'C:/Users/James_Mann1/Desktop/Testbench')
        config.set('main', 'Primed_RPA_Output', 'C:/Users/James_Mann1/Desktop/Testbench')
        config.set('main', 'Input.fna_Location','C:/Users/James_Mann1/Desktop/Testbench')
        config.set('main', 'Background_Blast_File', 'C:/Users/James_Mann1/Desktop/Testbench')
        config.add_section('Features')
        config.set('Features', 'Run_Cas-off', 'Yes')
        with open('config.ini', 'w') as f:
            config.write(f)
    working_dir = config.get('main', 'Working_Directory')
    print(os.getcwd())
    ###############################################################################
    #This section of the code contains housekeeping functions, such as deletion
    #of created file folders.
    ###############################################################################
    #%%
    import shutil as sh #imported for rmtree function
    delete_dir = [] #stores all the created dir for later
    working_dir = config.get('main', 'Working_Directory')
    catch = 0
    #Cas-Off Target Cas-Offinder Files
    #try:
       # os.chdir(working_dir)
        #sh.rmtree("Cas-OffInput")
       # print("[Primed Sherlock] Found existing Cas-OffInput Directory, deleting directory")
    #except:
      #  catch = 1
    
    #Cas-On
    try:
        os.chdir(working_dir)
        sh.rmtree("Cas-On")
        print("[Primed Sherlock] Found existing Cas-On Directory, deleting directory")
    except:
        catch =  1
    
    #Cas-On-Input
    try:
        os.chdir(working_dir)
        sh.rmtree("Cas-ONInput")
        print("[Primed Sherlock] Found existing Cas-OffInput Directory, deleting directory")
    except:
        catch = 1
    #Cas-On-Out
    try:
        os.chdir(working_dir)
        sh.rmtree("Cas-On-Out")
        print("[Primed Sherlock] Found existing Cas-On-Out Directory, deleting directory")
    except:
        catch = 1
    #Cas-On-Temp
    try:
        os.chdir(working_dir)
        sh.rmtree("Cas-On-Temp")
        print("[Primed Sherlock] Found existing Cas-On-Temp Directory, deleting directory")
    except:
        catch = 1
        
    #crRNA_Pairs
    try:
        os.chdir(working_dir)
        sh.rmtree("crRNA_Pairs")
        print("[Primed Sherlock] Found existing crRNA_Pairs Directory, deleting directory")
    except:
        catch = 1
        
    #Off
    try:
        os.chdir(working_dir)
        sh.rmtree("Off")
        print("[Primed Sherlock] Found existing Off Directory, deleting directory")
    except:
        catch = 1
    
    #On
    try:
        os.chdir(working_dir)
        sh.rmtree("On")
        print("[Primed Sherlock] Found existing On Directory, deleting directory")
    except:
        catch = 1
    
    if catch == 1:
        print("[Primed Sherlock] Several old directories were found and deleted, leaving directories or tampering with normal function can lead to program instability")
        print("[Primed Sherlock] If you know what you are doing, disable this section by deleting lines 51-115")
    catch = 0 #This value is returned to zero as it is reused elsewhere, several times.
    
    ###############################################################################
    #Extracting PrimedRPA File
    ###############################################################################
    working_dir = config.get('main', 'Working_Directory')
    os.chdir(working_dir)
    import shutil
    try:
        shutil.rmtree("crRNA_Pairs")
    except:
        print("[Primed Sherlock] Unable to remove crRNA Pairs, this is normal for first time operation")
    try:
        shutil.rmtree("PrimerComboFasta")
    except:
        catch = 0
    try:
        shutil.rmtree("TempAlignDir")
        shutil.rmtree("TempDIR")
    except:
        print("[Primed Sherlock] Manually delete all directories if issues persist")
    if str(config.get('Features', "Run_Cas-off")) == 'Yes':
        print("[Primed Sherlock] Running Cas-offinder using PrimedRPA results")
        PrimedOutputLoc = config.get('main','Primed_RPA_Output')
        os.chdir(PrimedOutputLoc)
        items = os.listdir(".")
        newlist = []
        for names in items:
            if names.endswith("Output_Sets.csv"):
                newlist.append(names)
        if len(newlist) >= 2:
            print("Ensure there is only one output file in the working directory, the first output file will be used otherwise")
        try:
            
            data = pd.read_csv(newlist[0])
            rpa_csv = open(newlist[0], "r")
        except:
            print("No Output file exists in directory: " + str(PrimedOutputLoc))
        
        theta = config.get('main','Amplicon_Max_Size')
        beta = config.get('main','Amplicon_Min_Size')
        maxback = config.get('main', 'Max_Background_Identitiy')
    
          
        columns = list(data.columns.values.tolist()) 
        Amp = data[data.columns[0]]
        FPrimer = data[data.columns[3]]
        RPrimer = data[data.columns[12]]
        MaxBackground = data[data.columns[6]]
        FPrimerSite = data[data.columns[1]]
        RPrimerSite = data[data.columns[10]]
        time_count = 0
        
        
        New_Amp = []
        New_FPrimer = []
        New_RPrimer = []
        New_MaxBackground = []
        NewFPrimerSite = []
        NewRPrimerSite = []
        for amplicons in Amp:
            if amplicons <= int(theta) and amplicons >= int(beta) :
               if MaxBackground[time_count] <=  int(maxback):
                New_Amp.append(Amp[time_count])
                New_FPrimer.append(FPrimer[time_count])
                New_RPrimer.append(RPrimer[time_count])
                New_MaxBackground.append(MaxBackground[time_count])
                NewFPrimerSite.append(FPrimerSite[time_count])
                NewRPrimerSite.append(RPrimerSite[time_count])
            time_count = time_count + 1
    ##############################################################################
    #Creates Multifastas of each respective primersite
    ##############################################################################
    try:
        Input_File = config.get('main','Input.fna_Location')
    except:
        print("[Primed Sherlock] Unable to locate Input.fna file")
    
    os.chdir(Input_File)
    items = os.listdir(".")
    newlist = []
    for names in items:
        if names.endswith("nput.fna"):
            newlist.append(names)
    
    input_fasta = open(newlist[0], "r")
    try:
        os.mkdir("PrimerComboFasta")
    except:
        print("[Primed Sherlock] Temporary Primer Combo Folder")
    from Bio import SeqIO
    n = 0
    sequence = []
    seq_files = []
    time_count = 0
    for primer in New_FPrimer:
        time_count = 0
        CurrentForward = primer
        CurrentReverse = New_RPrimer[n]
        CCF = NewFPrimerSite[n]
        CCR = NewRPrimerSite[n]
        Amplicon_Size = New_Amp[n]
        #print(CCF)
        #print(CCR)
        CutFix = CCR - CCF
        CutFix = int(CutFix) + len(CurrentReverse)
        CurrentBackgroundScore = New_MaxBackground[n]
        FileName = str(CurrentForward + "_" + CurrentReverse + "_" + str(Amplicon_Size))
        #print(FileName)
        os.chdir("PrimerComboFasta")
        output_file = open(FileName + '.fna', "w") 
        os.chdir(Input_File)
        inFile = open(newlist[0], "r")
        for seq_record in SeqIO.parse(inFile,'fasta'):
                seq_id = seq_record.id
                #sequence = str(seq_record.seq.ungap("-").upper())
                sequence = str(seq_record.seq.upper())
                #PrimerSite1 = sequence.find(CurrentForward)
                length = len(sequence) #length of code
                Cut1Done = sequence[(int(CCF)): (length)]
                finalsequence = Cut1Done[0:(int(CutFix))]
                output_file.write(">" + seq_record.id + "\n" + finalsequence + "\n")
        #print(time_count , " primer pair amplicons extracted for input fasta." + "\n")
        seq_files.append(FileName)
        n = n + 1
    ##############################################################################
    #Makes Consensus of primer specific extracted input multifasta
    ##############################################################################
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    try:
        os.chdir(Input_File)
        os.mkdir("TempAlignDir")
    except:
        catch = 1
    os.chdir("PrimerComboFasta")  
    com_dir = "PrimerComboFasta"
    tem_dir = "TempAlignDir"
    consensus_list = []
    def seq_align(Input_File,seq_files,com_dir,tem_dir):
        for amplicon_pairs in seq_files:
            os.chdir(Input_File)
            os.chdir(com_dir)
            seq_alignment = AlignIO.read(amplicon_pairs + ".fna", 'fasta')
            summary_align = AlignInfo.SummaryInfo(seq_alignment)
            os.chdir(Input_File)
            os.chdir(tem_dir)
            summary_align.dumb_consensus(0.5)
            finalconsequence = summary_align.dumb_consensus(0.5)
            consensus_list.append(str(amplicon_pairs) + "con" + '.fna')
            with open(str(amplicon_pairs) + "con" + '.fna', "w") as output_file:
                output_file.write(">" + str(amplicon_pairs) + "\n" + str(finalconsequence) + "\n")
    seq_align(Input_File,seq_files,com_dir,tem_dir)
    def getComplement(seq,
    				 reverse=False,
    				 rule='N2N'):
    	seqComp = ""
    	for base in seq.upper():
    		base = base.upper()
    		if base == "A":
    			seqComp += "T"
    		elif base == "T":
    			seqComp += "A"
    		elif base == "C":
    			seqComp += "G"
    		elif base == "G":
    			seqComp += "C"
    		elif base == 'U':
    			seqComp += "A"
    		elif base == "N":
    			if rule == 'N2-':
    				seqComp += '-'
    			elif rule == 'N2N':
    				seqComp += "N"
    		elif base == "-":
    			seqComp += "N"
    		# If Character not any of above assign it as N
    		else:
    			seqComp += "N"
    	if reverse:
    		return(seqComp[::-1])
    	else:
    		return(seqComp)
    
    ##############################################################################
    #Finds High Value crRNA targets
    ##############################################################################
    n = 0
    try:
        os.chdir(Input_File)
        os.mkdir("crRNA_Pairs")
    except:
        catch = 1
    os.chdir(Input_File)
    os.chdir(tem_dir)
    sequence_table = []
    crRNASiteLocation = []
    crRNAlist = [] 
    from Bio import SeqUtils
    for primers in consensus_list:
       #print("\n")
       #print("New Primer Pair " + str(primers))
       os.chdir(Input_File)
       os.chdir(tem_dir)
       FPrimer = New_FPrimer[n]
       RPrimer = New_RPrimer[n]
       n = n + 1
       inFile = open(primers,"r")
       crRNA = []
       crRNA2 = []
       leon = 0
       for seq_record in SeqIO.parse(inFile,'fasta'):
            os.chdir(Input_File)
            os.chdir("crRNA_Pairs")        
            output_file = open("TTT_" + primers + "_crRNA" + '.csv', "w")
            os.chdir(Input_File)
            os.chdir(tem_dir)
            output_file.write(str(primers + ","))
            output_file.write("\n")
            os.chdir(Input_File)
            os.chdir("crRNA_Pairs")        
            output_file = open("AAA_" + primers + "_crRNA" + '.csv', "w")
            os.chdir(Input_File)
            os.chdir(tem_dir)
            output_file.write(str(primers + ","))
            output_file.write("\n")
            sequence_table = []
            seq_id = seq_record.id
            amplicon_sequence = str(seq_record.seq.ungap("-").upper())
            amplicon_sequence_complement = getComplement(amplicon_sequence, False, 'N2N')
            pattern = "TTTNNNNNNNNNNNNNNNNNNNNN"
            crRNA = SeqUtils.nt_search(amplicon_sequence, pattern)
            pattern = "NNNNNNNNNNNNNNNNNNNNNAAA"        
            crRNA2 = SeqUtils.nt_search(amplicon_sequence, pattern)     
            try:
                crRNA.pop(0)
                UniquecrRNA = []
                crRNA_count = 0
                for targetsite in crRNA:
                    crRNAlocation = targetsite
                    crRNAlocation = crRNAlocation
                    length = len(amplicon_sequence) #length of code
                    Cut1Done = amplicon_sequence[(int(crRNAlocation)): (length)]
                    #crRNAactual.append(Cut1Done[0:24])
                    temp_crRNA = str(Cut1Done[0:24])
                    if temp_crRNA not in sequence_table:
                        crRNASiteLocation.append(targetsite)
                        sequence_table.append(temp_crRNA)
                        top_end = length
                        if crRNASiteLocation[crRNA_count] >= 30 and crRNASiteLocation[crRNA_count] <= top_end:
                            UniquecrRNA.append(temp_crRNA)
                        merged_list = [(sequence_table[i], crRNASiteLocation[i]) for i in range (0, len(sequence_table))]
                        #UniquecrRNA = list(dict.fromkeys(merged_list))
                        #UniquecrRNA = list(dict.fromkeys(crRNAactual))
                    os.chdir(Input_File)
                    os.chdir("crRNA_Pairs")
                    output_file = open("TTT_" + primers + "_crRNA" + '.csv', "a")
                    UniquecrRNA = list(dict.fromkeys(UniquecrRNA))
                    for values in UniquecrRNA:
                        output_file.write(str(values + ","))
                        output_file.write(str("\n"))
                    crRNA_count = crRNA_count + 1
            except: 
                print("\n")
                print('No crRNA site + ')  
            try:
                crRNA2.pop(0)
                UniquecrRNA = []
                crRNA_count = 0
                for targetsite in crRNA2:
                    crRNAlocation = targetsite
                    crRNAlocation = crRNAlocation
                    length = len(amplicon_sequence) #length of code
                    Cut1Done = amplicon_sequence[(int(crRNAlocation)): (length)]
                    #crRNAactual.append(Cut1Done[0:24])
                    temp_crRNA = str(Cut1Done[0:24])
                    
                    if temp_crRNA not in sequence_table:
                        crRNASiteLocation.append(targetsite)
                        sequence_table.append(temp_crRNA)
                        top_end = length
                        if crRNASiteLocation[crRNA_count] >= 30 and crRNASiteLocation[crRNA_count] <= top_end:
                            UniquecrRNA.append(temp_crRNA)
                        merged_list = [(sequence_table[i], crRNASiteLocation[i]) for i in range (0, len(sequence_table))]
                        #UniquecrRNA = list(dict.fromkeys(merged_list))
                        #UniquecrRNA = list(dict.fromkeys(crRNAactual))
                    os.chdir(Input_File)
                    os.chdir("crRNA_Pairs")
                    output_file = open("AAA_" + primers + "_crRNA" + '.csv', "a")
                    UniquecrRNA = list(dict.fromkeys(UniquecrRNA))
                    for values in UniquecrRNA:
                        output_file.write(str(values + ","))
                        output_file.write(str("\n"))
                    crRNA_count = crRNA_count + 1
            except: 
                print("\n")
                print('No crRNA site - ')  
    ###############################################################################
    #Grabs crRNA info
    ###############################################################################
    import os
    import pandas as pd
    import shutil as sh
    path = os.getcwd()
    crRNAfiles = os.listdir(path)
    try:
        os.chdir(Input_File)
        sh.rmtree("Cas-OffInput")
    except:
        catch = 1
    os.chdir(path)
    for value in crRNAfiles:
        value_csv = pd.read_csv(value)
        print(value_csv.head())
        if len(value_csv) >= 1:
            os.chdir(Input_File)
            try:
                os.mkdir("Cas-OffInput")
                os.chdir("Cas-OffInput")
            except:
                os.chdir("Cas-OffInput")
                catch = 1
            shortend_value = value.find("con.fna")
            short = value[0:int(shortend_value)]
            input_file = open(str(short) + "input" + ".txt", "w")
            input_file.write(str(working_dir) + "/off" + "\n")
            rows = value_csv[value_csv.columns[0]]
            if value.startswith("TTT_"):
                input_file.write(str("TTTNNNNNNNNNNNNNNNNNNNNNNNNN") + "\n")
            if value.startswith("AAA_"):
                input_file.write(str("NNNNNNNNNNNNNNNNNNNNNNNNNAAA") + "\n")       
            currentRNA = []
            for crRNA in rows:
                currentRNA.append(crRNA)
            currentRNA = list(dict.fromkeys(currentRNA)) 
            for crRNA in currentRNA:
                input_file.write(str(crRNA) + " 10" + "\n")
        os.chdir(path)
    ##############################################################################
    #Makes Cas-offinder off-target
    ##############################################################################
    #%%
    from Bio import SeqIO
    background_file_loc = config.get('main', 'Background_Blast_File')
    os.chdir(background_file_loc)
    items = os.listdir(".")
    newlist = []
    for names in items:
        try:
            if names.endswith("db.fasta"):
                newlist.append(names)
                Background_File = newlist[0]
        except:
            print("[Primed Sherlock] No blast_db.fasta file found in " + str(background_file_loc))
    if len(newlist) >= 2:
       print("Ensure there is only one output file in the working directory, the first output file will be used otherwise")
    print(newlist)
    
    print("[Primed Sherlock] Using " + str(Background_File) + " as background input fasta file.")
    sequences={}
        # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(Background_File, "fasta"):
        # Take the current sequence
        sequence = str(seq_record.seq).upper()
        # Check if the current sequence is according to the user parameters
        if (len(sequence) >= 500 and
            (float(sequence.count("N"))/float(len(sequence)))*100 <= 10):
        # If the sequence passed in the test "is it clean?" and it isn't in the
        # hash table, the sequence and its id are going to be in the hash
            if sequence not in sequences:
                sequences[sequence] = seq_record.id
          # If it is already in the hash table, we're just gonna concatenate the ID
          # of the current sequence to another one that is already in the hash table
            #else:
                #sequences[sequence] += "_" + seq_record.id
    
    with open("clear_" + Background_File, "w+") as output_file:
            # Just read the hash table and write on the file as a fasta format
            for sequence in sequences:
                output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")
    
    
    Background_File = "clear_" + Background_File
    for record in SeqIO.parse(Background_File,"fasta"):
        file = str(record.id)
        try:
            os.mkdir("Off")   
        except:
            catcha = []
        os.chdir("Off") 
        with open(str(file), "w") as handle: 
            SeqIO.write(record, handle, "fasta")
        os.chdir(background_file_loc)
    os.chdir(master_dir)
    
    os.remove("clear_blast_db.fasta")
    
    
    print("[Primed Sherlock] This next section may take several hours or days to complete. See Multithreading in Config for settings to speed this up.")
    
##############################################################################
#Runs Cas off-finder on background sequences
##############################################################################
# #%%
# import shutil as sh
# Input_File ='c:/Users/James_Mann1/Desktop/TestBench'

# import os
# try:
#     os.chdir(Input_File)
# except:
#     catch = 1
# try:
#     sh.rmtree("Cas-Out")
# except:
#     catch = 1
# try:
#     os.mkdir("Cas-Out")
# except:
#     catch = 1
# #print(os.getcwd())
# os.chdir("Cas-OffInput")
# path = os.getcwd()
# Input_Casoff = os.listdir(path)
# src = Input_File + "/cas-offinder.exe"
# dst = Input_File + "/Cas-OffInput/cas-offinder.exe"
# sh.copyfile(src, dst)
# for value in Input_Casoff:
#     subprocess.call("cas-offinder " + value + " G0 " + " " + value + "out.txt")
#     #os.system("cas-offinder " + value + " G0 " + " " + value + "out.txt")
#################################################################################
#Multiprocessing, Cas-Off
#Runs Cas-Off on background sequences
#################################################################################
#%%
import os
import shutil as sh
import subprocess
from multiprocessing import Pool

def cas_prep(Input_File, Output_Dir, Input_path):
    import subprocess
    import os
    import shutil as sh
    try:
        os.chdir(Input_File)
    except:
        catch = 1
    os.chdir(Output_Dir)
    Input_Casoff = os.listdir(Input_path)
    good_check = []
    for values in Input_Casoff:
        if values.endswith("input.txt"):
            good_check.append(values)
    Input_Casoff = good_check
    src = Input_File + "/cas-offinder.exe"
    dst = Output_Dir + "/" + "cas-offinder.exe"
    sh.copyfile(src, dst)
    return Input_Casoff

def Cas_off(x):
    subprocess.call("cas-offinder " + x + " G0 " + " " + x + "out.txt")
    #subprocess.call("cas-offinder " + x + " G0 " + " " + x + "out.txt")
    return x

if __name__ == '__main__':
    Input_File ='c:/Users/James_Mann1/Desktop/TestBench'
    Output_Dir = Input_File + "/" + "Cas-OffInput"
    Input_path = Output_Dir
    Input_Casoff = cas_prep(Input_File, Output_Dir, Input_path, )
    with Pool(32) as p:
        print(p.map(Cas_off, Input_Casoff))
        
print("[Primed Sherlock] Done processing off-targets with Cas-Off Finder. Please Cite Author")
##############################################################################
#After Primed RPA and Cas-OffFinder Is Ran, This checks for crRNA targets
#and finds those that have little to no matches to background sequences
##############################################################################
#%%
if __name__ == '__main__':
    import os
    import pandas as pd
    #working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' (REMOVE)
    working_dir = config.get('main', 'Working_Directory')
    os.chdir(working_dir)
    path = os.getcwd() #Checks for sectional code piece
    if path.endswith("Cas-OffInput") == True:
        print("[Primed Sherlock] Correct Working Directory Loaded")
    else:
        try:
            path = path + '\Cas-OffInput'
            os.chdir(path)
            files = os.listdir(path)
        except:
            print(path + "Fix working directory")
    files = os.listdir(path) #Searches all the files
    print("[Primed Sherlock] Coverting Aquired Cas Off-Target Data to CSV format")
    for file in files:
        if file.endswith('txtout.txt') == True:
            file_out = file + ".csv"
            with open(file) as infile, open(file_out, 'w') as outfile:
                outfile.write("Cas_off-Results" + "," + "SeqID" + "," + "Location" + "," + "Mismatch" + "," + "MismatchSeq" + "," + "Score" + "," + "\n")
                for line in infile:
                      outfile.write(" ".join(line.split()).replace(' ', ','))
                      outfile.write("," + "\n") # trailing comma shouldn't matter  
    print("[Primed Sherlock] Program has finished converting Cas-Off results to CSV for analysis") 
    
    
    
    ###############################################################################
    #This code segment sets a threshold for offtargets as well as discards crRNA
    #that Cas-Offinder has determined to potentially target background sequences
    ###############################################################################    
    import os
    import pandas as pd
    user_defined_value = config.get('main', 'Off_Target_Mismatch_Score') #This is the value for Cas-Offinder's mismatch, goes after crRNA so e.g. crRNA 7
    user_defined_value = int(user_defined_value)
    #user_defined_value = 7 # defines value for specificity, less then or equal to. (7 normally, used for Cas-off output) (Remove if works)
    #exception_count = 0
    pathtest = os.getcwd()
    path = os.getcwd()
    if pathtest.endswith("Cas-OffInput") == True:
        print("[line 503] Correct Working Directory")
    else:
        try:
            path = path + '\Cas-OffInput'
            os.chdir(path)
        except:
            print(pathtest + "Fix working directory")
    files = os.listdir(path)
    path_length = len(path)
    path_length = path_length - 13
    output_path = path[0:int(path_length)]
    output_check = 0
    primer_list = []
    primer_string = []
    for file in files:
        if file.endswith('.csv'):
           file1 = file   
           crRNA_List = []
           check_for_unique = 0
           Good_List = [] #List of good RPA Amplicons with crRNA and low off target
           if file.endswith('txt.csv'):
                file_name = file1
                print(file_name)
                file_data = pd.read_csv(file)
                #Grabs Primers + Amplicon Size
                Trim_loc = file1.find("_") + 1
                length = len(file1)
                file2 = file1[Trim_loc:length]
                #F-Primer
                Trim_loc = file2.find("_") + 1
                length2 = len(file2)
                FPrimer = file2[0:Trim_loc - 1]
                file2 = file2[Trim_loc:length2]    
                Trim_loc = file2.find("_") + 1
                length2 = len(file2)
                RPrimer = file2[0:Trim_loc - 1]
                file2 = file2[Trim_loc:length2]    
                length3 = len(file2)
                Trim_loc = file2.find("i")
                Amplicon = file2[0:Trim_loc]  
                file_details = []
                file_details.append(FPrimer)
                file_details.append(RPrimer)
                file_details.append(Amplicon)
                columns = list(file_data.columns.values.tolist()) 
                crRNA = file_data[file_data.columns[0]]
                SeqID = file_data[file_data.columns[1]]
                Loc = file_data[file_data.columns[2]]
                Score = file_data[file_data.columns[5]]
                #Grabs all Unique crRNA and puts into list
                for crRNA_Target in crRNA:
                    crRNA_List.append(crRNA_Target)
                Unique_crRNA_List = list(dict.fromkeys(crRNA_List))
                #Does the legwork
                #holder
                Good_crRNA_List = []
                Uniquedetails = []
                for Unique in Unique_crRNA_List:
                    count_sum = 0
                    bad_sequences_less_ud = 0
                    bad_sequence_id_list_lesser =[]
                    bad_sequences_greater_ud = 0
                    bad_sequence_id_list_greater = []
                    potentially_dangerous_seq = 0
                    potentially_dangerous_seq_list = []
                    for value in crRNA:
                        if Unique == value:
                            current_score = Score[count_sum]
                            if current_score <= user_defined_value:
                                bad_sequences_less_ud = bad_sequences_less_ud + 1
                                bad_sequence_id_list_lesser.append(SeqID[count_sum])
                            if current_score > user_defined_value: 
                                bad_sequences_greater_ud = bad_sequences_greater_ud + 1
                                bad_sequence_id_list_greater.append(SeqID[count_sum]) 
                            if current_score <= user_defined_value: 
                                potentially_dangerous_seq = potentially_dangerous_seq + 1
                                potentially_dangerous_seq_list.append(SeqID[count_sum]) 
                    count_sum = count_sum + 1
                    if potentially_dangerous_seq == 0:
                        Good_crRNA_List.append(Unique)
                        Uniquedetails.append(Unique)
                        string2 = str("Total bad seq " + str(bad_sequences_greater_ud) + "\n")        
                        Uniquedetails.append(string2)
                        Uniquedetails = list(dict.fromkeys(Uniquedetails))
                if len(Good_crRNA_List) >= 1:
                    check_for_unique = 1
                    primer_name = FPrimer + "_" + RPrimer + "_" + Amplicon
                    primer_list.append(primer_name) #creates list of primer pairs with crRNA targets on both Pos and Neg Strand
                    crRNA_value_string = ""
                    for entries in Good_crRNA_List:
                        crRNA_value_string = str(crRNA_value_string) +  str(entries) + ","
                    primer_string.append(crRNA_value_string)
                if output_check == 0:
                    with open("OffTargetAnalysis.csv", "w") as output_file:
                        output_file.write("Off Target crRNA analysis" + "," + "\n")
                        output_check = output_check + 1
                else:                    
                    with open("OffTargetAnalysis.csv", "a") as output_file:
                        output_file.write("Forward Primer" + ',' + file_details[0] + "," + "\n")
                        output_file.write("Reverse Primer" + ',' + file_details[1] + "," + "\n")
                        output_file.write("Amplicon Size"  + ',' + file_details[2] + "," + "\n" + ",")  #This should fix it     
                        for Unique_info_end in Uniquedetails:
                            output_file.write(str(Unique_info_end)  + ',')
                        output_file.write("\n")
    print("[Primed Sherlock] Removing temporary Background CSV Files")
    if path.endswith("Cas-OffInput") == True:
        print("Correct Working Directory")
    else:
        try:
            path = path + '\Cas-OffInput'
            os.chdir(path)
            files = os.listdir(path)
        except:
            print(path + "Fix working directory")
    
    ###############################################################################
    #This shortens the list down to useable crRNA pairs using the output data from 
    #performed analysis
    ###############################################################################
    import os
    try:
        os.mkdir("Combos")
        os.chdir("Combos")
    except:
        os.chdir("Combos")
        print("[Primed Sherlock] Unable to open Combos directory")
    counting_limit = 0
    for primers in primer_list:
        name = str(primers) + ".csv"
        if os.path.isfile(name) == False:    
            with open(name, "w") as output_file:
                output_file.write("Primer Pairs" + "," + "\n")
        with open(name, "a") as output_file:            
            stringfile = str(primer_string[counting_limit])
            #print("this is current " + stringfile + "\n")
            li = list(stringfile.split(","))
            lengthofli = len(li) - 1
            listcount = 0
            for values in li:
                if listcount < lengthofli:
                    listcount = listcount + 1
                    output_file.write(values + "\n")
        counting_limit = counting_limit + 1
            
    ###############################################################################
    #Makes on target Code
    ###############################################################################
    
    import os
    working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' #FIX
    CombosPath = working_dir + "/Cas-OffInput/Combos" #FIX
    os.chdir(working_dir)
    Cas_On = working_dir + "/Cas-On"
    try:
        os.mkdir(Cas_On)
        os.chdir(Cas_On)
        print("[Primed Sherlock] Making Cas On Target Directory")
    except:
        print("[Primed Sherlock] Switching to Cas On Target Directory")
    import os
    import pandas as pd
    import shutil as sh
    temp_on = working_dir + "/Cas-On-Temp"
    try:
        os.mkdir(temp_on)
    except:
        catch = 1
    os.chdir(CombosPath)
    path = os.getcwd()
    crRNAfiles = os.listdir(path)
    screening_list = []
    amount = 0
    #print(crRNAfiles)
    for value in crRNAfiles:
        #print(value)
        line_check = 0 #Checks amount of crRNA targets
        AAA_Check = 0
        TTT_Check = 0
        print(os.getcwd())
        with open(value, "r") as input_file:
            for line in input_file:
                line_check += 1
            if line_check >= 3:
                value_length = len(value)
                #print(value)
                amount += 1
                print(value)
                values_csv = pd.read_csv(value)
                #print(value_csv.head())
                values_csv.dropna()
                crRNA = values_csv[values_csv.columns[0]]
                crRNA = crRNA.tolist()
                crRNA = [x for x in crRNA if str(x) != 'nan']
                print(crRNA)
                os.chdir(temp_on)
                #crRNA.dropna()
                print(crRNA)
                for crRNA_ON in crRNA:
                    #print(crRNA_ON)
                    if crRNA_ON.startswith("TTT") == True:
                        if TTT_Check == 0:
                            with open("TTT_" + value + ".csv", "w") as output_file:
                                output_file.write("crRNA_Analysis" + "," + "\n")
                                TTT_Check += 1
                        with open("TTT_" + value + ".csv", "a") as output_file:
                            output_file.write(crRNA_ON + "," + "\n")
                    if crRNA_ON.endswith("AAA") == True:
                        if AAA_Check == 0:
                            with open("AAA_" + value + ".csv", "w") as output_file:
                                output_file.write("crRNA_Analysis" + "," + "\n")
                                AAA_Check += 1
                        with open("AAA_" + value + ".csv", "a") as output_file:
                            output_file.write(crRNA_ON + "," + "\n")
        os.chdir(path)
    ###############################################################################
    #Makes inputFiles
    ###############################################################################     
#%%
if __name__ == '__main__':        
        
    import os
    working_dir = 'C:/Users/James_Mann1/Desktop/Testbench'
    temp_on = working_dir + "/Cas-On-Temp"
    
    on_target_background_file_loc = working_dir #FIX
    Input_File = working_dir #FIX
    os.chdir(working_dir)
    os.chdir(temp_on)
    filepath = os.getcwd()
    crRNAfiles_on = os.listdir(filepath)
    for value in crRNAfiles_on:
        print(value)
        value_csv = pd.read_csv(value)
        value_csv.dropna()
        print(value_csv.head())
        if len(value_csv) >= 1:
            os.chdir(Input_File)
            try:
                os.mkdir("Cas-ONInput")
                os.chdir("Cas-ONInput")
            except:
                os.chdir("Cas-ONInput")
                catch = 1
            shortend_value = value.find(".csv")
            short = value[0:int(shortend_value)]
            input_file = open(str(short) + "input" + ".txt", "w")
            input_file.write(str(working_dir) + "/on" + "\n")
            rows = value_csv[value_csv.columns[0]]
            if value.startswith("TTT_"):
                input_file.write(str("TTTNNNNNNNNNNNNNNNNNNNNNNNNN") + "\n")
            if value.startswith("AAA_"):
                input_file.write(str("NNNNNNNNNNNNNNNNNNNNNNNNNAAA") + "\n")       
            currentRNA = []
            for crRNA in rows:
                currentRNA.append(crRNA)
            currentRNA = list(dict.fromkeys(currentRNA)) 
            for crRNA in currentRNA:
                input_file.write(str(crRNA) + " 10" + "\n")
        print(os.getcwd())
        os.chdir(temp_on)
        print(os.getcwd())
    ##############################################################################
    #Makes Cas-offinder for on target strains
    ###############################################################################
    import os
    import Bio
    from Bio import SeqIO
    #on_target_background_file_loc = config.get('main', 'Background_Blast_File') #FIX
    os.chdir(on_target_background_file_loc)
    items = os.listdir(".")
    newlist = []
    for names in items:
        try:
            if names == ("input.fna"):
                newlist.append(names)
                Background_File = newlist[0]
        except:
            print("[Primed Sherlock] No input.fna file found in " + str(on_target_background_file_loc))
        if len(newlist) >= 2:
            print("Ensure there is only one output file in the working directory, the first output file will be used otherwise")
    print(newlist)
    
    sequences={}
        # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(Background_File, "fasta"):
        # Take the current sequence
        sequence = str(seq_record.seq).upper()
        # Check if the current sequence is according to the user parameters
        if (len(sequence) >= 500 and
            (float(sequence.count("N"))/float(len(sequence)))*100 <= 10):
        # If the sequence passed in the test "is it clean?" and it isn't in the
        # hash table, the sequence and its id are going to be in the hash
            if sequence not in sequences:
                sequences[sequence] = seq_record.id
          # If it is already in the hash table, we're just gonna concatenate the ID
          # of the current sequence to another one that is already in the hash table
            #else:
                #sequences[sequence] += "_" + seq_record.id
    
    with open("clear_" + Background_File, "w+") as output_file:
            # Just read the hash table and write on the file as a fasta format
            for sequence in sequences:
                output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")
    
    
    Background_File = "clear_" + Background_File
    for record in SeqIO.parse(Background_File,"fasta"):
        file = str(record.id)
        try:
            os.mkdir("On")   
        except:
            catcha = []
        os.chdir("On") 
        with open(str(file), "w") as handle: 
            SeqIO.write(record, handle, "fasta")
        os.chdir(on_target_background_file_loc)
    #os.chdir(master_dir) #Fix
    
    ##############################################################################
    #Runs Cas off-finder for the on-target strains
    ##############################################################################
    # import subprocess
    # import shutil as sh
    # Input_File = config.get('main', 'Input.fna_Location')
    # try:
    #     os.chdir(Input_File)
    # except:
    #     catch = 1
    # try:
    #     sh.rmtree("Cas-On-Out")
    # except:
    #     catch = 1
    # try:
    #     os.mkdir("Cas-On-Out")
    # except:
    #     catch = 1
    # os.chdir("Cas-ONInput")
    # path = os.getcwd()
    # Input_Casoff = os.listdir(path)
    
    # good_check = []
    # for values in Input_Casoff:
    #     if values.endswith("input.txt"):
    #         good_check.append(values)
    # Input_Casoff = good_check
    # src = Input_File + "/cas-offinder.exe"
    # dst = Input_File + "/Cas-ONInput/cas-offinder.exe"
    # sh.copyfile(src, dst)
    # countsum = len(Input_Casoff)
    # for value in Input_Casoff:
    #    subprocess.call("cas-offinder " + value + " G0 " + " " + value + "out.txt")
    #    #os.system("cas-offinder " + value + " G0 " + " " + value + "out.txt")
    #    countsum = countsum - 1
    #    print("[Primed Sherlock] " + str(countsum) + " Remaining...")
       
    ##############################################################################
    #Runs Cas off-finder for the on-target strains
    ##############################################################################
    
# import subprocess
# import shutil as sh
# Input_File = config.get('main', 'Input.fna_Location')
# try:
#     os.chdir(Input_File)
# except:
#     catch = 1
# try:
#     sh.rmtree("Cas-On-Out")
# except:
#     catch = 1
# try:
#     os.mkdir("Cas-On-Out")
# except:
#     catch = 1
# os.chdir("Cas-ONInput")
# path = os.getcwd()
# Input_Casoff = os.listdir(path)
# good_check = []
# for values in Input_Casoff:
#     if values.endswith("input.txt"):
#         good_check.append(values)
# Input_Casoff = good_check
# src = Input_File + "/cas-offinder.exe"
# dst = Input_File + "/Cas-ONInput/cas-offinder.exe"
# sh.copyfile(src, dst)
# from multiprocessing import Pool
# def Cas_off(x):
#     subprocess.call("cas-offinder " + x + " G0 " + " " + x + "out.txt")
#     return x
#%%
import os
import shutil as sh
import subprocess
from multiprocessing import Pool

def cas_prep(Input_File, Output_Dir, Input_path):
    import subprocess
    import os
    import shutil as sh
    try:
        os.chdir(Input_File)
    except:
        catch = 1
    os.chdir(Output_Dir)
    Input_Casoff = os.listdir(Input_path)
    good_check = []
    for values in Input_Casoff:
        if values.endswith("input.txt"):
            good_check.append(values)
    Input_Casoff = good_check
    src = Input_File + "/cas-offinder.exe"
    dst = Output_Dir + "/" + "cas-offinder.exe"
    sh.copyfile(src, dst)
    return Input_Casoff

def Cas_off(x):
    subprocess.call("cas-offinder " + x + " G0 " + " " + x + "out.txt")
    #subprocess.call("cas-offinder " + x + " G0 " + " " + x + "out.txt")
    return x
#%%     
if __name__ == '__main__':
    Input_File ='c:/Users/James_Mann1/Desktop/TestBench'
    Output_Dir = Input_File + "/" + "Cas-ONInput"
    Input_path = Output_Dir
    Input_Casoff = cas_prep(Input_File, Output_Dir, Input_path, )
    with Pool(32) as p:
       print(p.map(Cas_off, Input_Casoff))    
    
    ##############################################################################
    #Goes through the sequences
    ##############################################################################
    
    import os
    working_dir = 'C:/Users/James_Mann1/Desktop/Testbench'
    os.chdir(working_dir)
    path = os.getcwd()
    print(path)
    if path.endswith("Cas-ONInput") == True:
        print("Correct Working Directory")
    else:
        try:
            path = path + '\Cas-ONInput'
            os.chdir(path)
            files = os.listdir(path)
        except:
            print(path + "Fix working directory")
    files = os.listdir(path)
    #Searches all the files
    print("[Primed Sherlock] Coverting Aquired Cas On-Target Data to CSV format")
    for file in files:
        if file.endswith('txtout.txt') == True:
            file_out = file + ".csv"
            with open(file) as infile, open(file_out, 'w') as outfile:
                outfile.write("Cas_off-Results" + "," + "SeqID" + "," + "Location" + "," + "Mismatch" + "," + "MismatchSeq" + "," + "Score" + "," + "\n")
                for line in infile:
                      outfile.write(" ".join(line.split()).replace(' ', ','))
                      outfile.write("," + "\n") # trailing comma shouldn't matter  
    
    print("DONE") 
#%%     
    ###############################################################################
    #Goes through ontarget sequences for data analysis
    ###############################################################################
    
    import os
    import csv
    import pandas as pd
    from itertools import zip_longest
    working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' #FIX
    os.chdir(working_dir)
    pathtest = os.getcwd()
    path = os.getcwd()
    if pathtest.endswith("Cas-OnInput") == True:
        print("Correct Working Directory")
    else:
        try:
            path = path + '\Cas-OnInput'
            os.chdir(path)
        except:
            print(pathtest + "Fix working directory")
    try:
        os.mkdir("Best_Pairs")
    except:
        print("Unable to create directory")
    current = os.getcwd()
    os.chdir("Best_Pairs")
    bestpairs = os.getcwd()
    os.chdir(current)
    try:
        os.mkdir("Temp")
    except:
        print("Temp Dir allready exists")
    os.chdir("Temp")
    temp_path = os.getcwd()
    os.chdir(path)
    files = os.listdir(path)
    primer_list = []
    primer_string = []
    primer = []
    #Sets up temp file
    current = os.getcwd()
    os.chdir(bestpairs)
    with open("MasterFile.csv", "w") as output_file:
        output_file.write("Primer" + "," + "crRNA" + "," + "Exact Match" + "," + "Matches within 3" + "," + "\n")
    os.chdir(current)
    for file in files:
        best_check = 0 
        Second_Best_check = 0
        best_check3 = 0
        Best_Unique = ""
        Second_Best_Unique = ""
        check_if_exists = 0
        if file.endswith('.csv'): 
           crRNA_List = []
           check_for_unique = 0
           Good_List = [] #List of good RPA Amplicons with crRNA and low off target
           if file.endswith('txt.csv'):
                file_data = pd.read_csv(file)
                #Grabs Primers + Amplicon Size
                Trim_loc = file.find("_") + 1
                length = len(file)
                Trim = file[Trim_loc:length]
                Trim_loc = Trim.find("input")
                Trimmed = Trim[0:Trim_loc]
                primer.append(Trimmed)
                file_name = Trimmed
                columns = list(file_data.columns.values.tolist()) 
                crRNA = file_data[file_data.columns[0]]
                SeqID = file_data[file_data.columns[1]]
                Loc = file_data[file_data.columns[2]]
                Score = file_data[file_data.columns[5]]
                #Grabs all Unique crRNA and puts into list
                for crRNA_Target in crRNA:
                    if crRNA_Target.startswith("TTTT") == False:
                        if crRNA_Target.endswith("AAAA") == False:
                            crRNA_List.append(crRNA_Target)
                Unique_crRNA_List = list(dict.fromkeys(crRNA_List))
                #print(Unique_crRNA_List)
                #NEED TO DELETE SEQS WHICH HAVE BAD INPUT
                Good_crRNA_List = []
                Uniquedetails = []
                usable_seq = []            
                
                #print(file_name)
                for Unique in Unique_crRNA_List:  
                    primer_string = file_name + "_" + Unique + ".csv"
                    count_sum = 0
                    seq_check = 0
                    usable_crRNA = 0
                    crRNA_matching_strains_0 = []
                    crRNA_matching_strains_3 = []
                    crRNA_matching_strains_greater = []
                    for value in crRNA:
                        if Unique == value:
                            current_score = Score[count_sum]
                            if current_score == 0:
                                crRNA_matching_strains_0.append(SeqID[count_sum])
                                #print(SeqID[count_sum])
                                usable_crRNA = 1
                            if current_score <=5 and current_score > 1:
                                usable_crRNA = 1
                                crRNA_matching_strains_3.append(SeqID[count_sum])
                                #print(SeqID[count_sum])
                            if current_score >= 5 and current_score <=7:
                                crRNA_matching_strains_greater.append(SeqID[count_sum])
                                #print(SeqID[count_sum])
                        count_sum = count_sum + 1
                    if usable_crRNA == 1:                    
                        os.chdir(temp_path)
                        with open(primer_string, "w") as output_file:
                            output_file.write("Matching Seq" + "," + "3 Mismatch" + "," + "More then 4 mismatch" + "," + "\n")
                        #Ranks Unique Seqs on matching etc
                        crRNA_matching_strains_0 = list(dict.fromkeys(crRNA_matching_strains_0))
                        crRNA_matching_strains_3 = list(dict.fromkeys(crRNA_matching_strains_3))
                        crRNA_matching_strains_greater = list(dict.fromkeys(crRNA_matching_strains_greater))
                        unique_rank_ontarget = len(crRNA_matching_strains_0)
                        unique_rank_mismatch = len(crRNA_matching_strains_3)
                        unique_rank_mismatch_bad = len(crRNA_matching_strains_greater)
                        Master_List = zip_longest(crRNA_matching_strains_0, crRNA_matching_strains_3, crRNA_matching_strains_greater)
                        with open(primer_string, "a", newline = "") as output_file:
                            for row in Master_List:
                                writer = csv.writer(output_file)
                                writer.writerow(row)
                                #print(row)
                        os.chdir(path)
                        #Grabs exact matches
    
    import os
    import pandas as pd
    temp_path = 'C:/Users/James_Mann1/Desktop/Testbench/Cas-ONInput/Temp' #FIX
    working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' #Fix
    os.chdir(temp_path)                   
    files = os.listdir(os.getcwd())
    RPA = ""
    RPA_STORE = []
    for file in files:
        if file.startswith("AAA_") == False: #Checks for erranious files
            if file.endswith(".csv") == True:
                    #Grabs Primers + Amplicon Size
                    #print(file)
                    Trim_loc = file.find("_") + 1
                    length = len(file)
                    FPrimer = file[0:Trim_loc - 1]
                    Trim = file[Trim_loc:length]
                    #print(FPrimer)
                    Trim_loc = Trim.find("_") + 1
                    length = len(Trim)
                    RPrimer = Trim[0:Trim_loc - 1]
                    Trim = Trim[Trim_loc:length]
                    #print(RPrimer)
                    Trim_loc = Trim.find("_") + 1
                    length = len(Trim)
                    AMP = Trim[0:Trim_loc - 1]
                    Trim = Trim[Trim_loc:length]
                    #print(AMP)                
                    Trim_loc = Trim.find(".csv")
                    crRNA = Trim[0:Trim_loc]
                    #print(crRNA)
                    RPASet = FPrimer + "_" + RPrimer
                    
                    RPA_STORE.append(RPASet)
    RPA_STORE = list(dict.fromkeys(RPA_STORE))
    current = os.getcwd()
    os.chdir(working_dir)
    with open("MasterAnalysis.csv", "w") as output_file:
        output_file.write("RPA PRIMER SET" + "," + "Amplicon Size" + "," + "crRNA Target" + "," + "On Target" + "," + "Mismatch 3" + "," + "Bad" + "," + "\n")
    os.chdir(current)
    for current_RPA in RPA_STORE:
        for file in files:
            if file.startswith(current_RPA):
                Trim_loc = file.find("_") + 1
                length = len(file)
                FPrimer = file[0:Trim_loc - 1]
                Trim = file[Trim_loc:length]
                #print(FPrimer)
                Trim_loc = Trim.find("_") + 1
                length = len(Trim)
                RPrimer = Trim[0:Trim_loc - 1]
                Trim = Trim[Trim_loc:length]
                #print(RPrimer)
                Trim_loc = Trim.find("_") + 1
                length = len(Trim)
                AMP = Trim[0:Trim_loc - 1]
                Trim = Trim[Trim_loc:length]
                #print(AMP)                
                Trim_loc = Trim.find(".csv")
                crRNA = Trim[0:Trim_loc]
                current_data = pd.read_csv(file)
                columns = list(current_data.columns.values.tolist()) 
                Match = current_data[current_data.columns[0]].count()
                Match3 = current_data[current_data.columns[1]].count()            
                MismatchBad = current_data[current_data.columns[2]].count()
                current_RPA = FPrimer + "_" + RPrimer 
                current_match = Match
                current_match3 = Match3
                current_MismatchBad = MismatchBad
                os.chdir(working_dir)
                with open("MasterAnalysis.csv", "a") as output_file:
                    output_file.write(current_RPA + "," + str(AMP) + "," + crRNA + "," + str(current_match) + "," +  str(current_match3) + "," + str(current_MismatchBad) + "," + "\n" )
                os.chdir(current)
    
    import os
    os.chdir(working_dir)
    file = "MasterAnalysis.csv"
    current_data = pd.read_csv(file)
    columns = list(current_data.columns.values.tolist()) 
    RPA_Primer = current_data[current_data.columns[0]]
    #print(RPA_Primer.head())
    Amplicon = current_data[current_data.columns[1]]
    crRNA = current_data[current_data.columns[2]]
    Match = current_data[current_data.columns[3]]
    Match3 = current_data[current_data.columns[4]]
    MismatchBad = current_data[current_data.columns[5]]
    Primer_List = []
    for primer in RPA_Primer:
        #print(primer)
        Primer_List.append(primer)
    Unique_Primer_List = list(dict.fromkeys(Primer_List))
    
    
    current = os.getcwd()
    os.chdir(working_dir)
    with open("MasterFile.csv", "w") as output_file:
        output_file.write("RPA PRIMER" + "," +  "AMPLICON" + "," + "Best crRNA" + "," + "Total Good" + "," + "MATCHING SEQ" + "," + "MISMATCH < 3" + "," + "MISMATCH > 7" + ","+ "\n")
    os.chdir(current)
    
    for Unique_RPA in Unique_Primer_List:
        catch = 0
        best_check = 0 
        Second_Best_check = 0
        best_check1 = 0
        Second_Best_check_1 = 0
        best_check3 = 0
        best_check4 = 0
        Second_Best_check_4 = 0
        Best_Unique = ""
        Second_Best_Unique = ""
        check_if_exists = 0
        Best_Amplicon = 0
        for primer in RPA_Primer:
            if Unique_RPA == primer:
                #Get
                current_RPA_Primer = RPA_Primer[catch]
                current_Amplicon = Amplicon[catch]
                current_crRNA = crRNA[catch]
                current_Match = Match[catch]
                current_Match3 = Match3[catch]
                current_MismatchBad = MismatchBad[catch]
                count1 = current_Match + current_Match3                 
                count2 = current_Match3
                count3 = current_Match
                count4 = current_MismatchBad
                if count1 > best_check:
                    if best_check > Second_Best_check:
                        Second_Best_check = best_check
                        Second_Best_AMP = Best_Amplicon
                        Second_Best_check_3 = best_check3
                        Second_Best_check_1 = best_check1
                        Second_Best_Unique = Best_Unique
                        Second_Best_check_4 = best_check4
                    best_check4 = count4    
                    best_check3 = count2
                    best_check1 = count3
                    best_check = count1
                    Best_Unique = current_crRNA
                    Best_Amplicon = current_Amplicon
                if count1 < best_check:
                    if count1 > Second_Best_check:
                        Second_Best_check_4 = count4
                        Second_Best_check = count1
                        Second_Best_check_3 = count2
                        Second_Best_check_1 = count3
                        Second_Best_Unique = current_crRNA
                        Second_Best_AMP = current_Amplicon
            catch = catch + 1 #if Unique == Primer:
        if best_check > 0:
            if Second_Best_check > 0:                
                #print(Unique)      
                #print(str(best_check) + " " + Best_Unique)
                #print(str(Second_Best_check) + " " + Second_Best_Unique)
                current = os.getcwd()
                os.chdir(working_dir)
                with open("MasterFile.csv", "a") as output_file:
                    output_file.write(str(Unique_RPA) + "," + str(Best_Amplicon) + "," + str(Best_Unique) + "," + str(best_check) + "," + str(best_check1) + "," + str(best_check3) + "," + str(best_check4) + "," + "\n")
                    output_file.write(str(Unique_RPA) + "," + str(Second_Best_AMP) + ","+  str(Second_Best_Unique) + "," + str(Second_Best_check) + "," + str(Second_Best_check_1) + "," + str(Second_Best_check_3) + "," + str(Second_Best_check_4) + "," + "\n")
                os.chdir(current)
   
#%%   
if __name__ == '__main__':    
    ##############################################################################
    #Sorts through output and finds best primer pairs based on three variables   #
    # 1. Amount of exact matches to input files                                  #
    # 2. Amount of mismatches less then 3                                        #
    # 3. Size of amplicon. Smaller amplicon == more efficent RPA                 #
    ##############################################################################
    
    working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' #Fix
    import os     
    import pandas as pd          
    print("[Primed Sherlock] Looking through output data to determine best RPA Primers + crRNA for maximized strain detection capability")
    os.chdir(working_dir)
    file = "MasterFile.csv"
    current_data = pd.read_csv(file)
    columns = list(current_data.columns.values.tolist()) 
    RPA_Primer = current_data[current_data.columns[0]]
    Amplicon = current_data[current_data.columns[1]]
    crRNA = current_data[current_data.columns[2]]
    Total_match = current_data[current_data.columns[3]]
    Match = current_data[current_data.columns[4]]
    Match3 = current_data[current_data.columns[5]]
    MismatchBad = current_data[current_data.columns[6]]
    Primer_List = []
    current = os.getcwd()
    temp_dir = working_dir + '/Cas-ONInput/Temp'
    for unique in RPA_Primer:
        Primer_List.append(unique)
    Unique_Primer_List = list(dict.fromkeys(Primer_List))
    #Gets info about matching pairs of crRNA
    
    best_seq_list = ""
    best_seq_amp = ""
    Best_crRNA_1 = ""
    Best_crRNA_2 = ""
    
    Unique_primer = []
    Unique_crRNA1 = []
    Unique_crRNA2 = []
    Unique_amp = []
    Unique_combined_match = []
    Unique_with_Mismatch = []
    Unique_combined_all = []
    Unique_combined_list = []
    
    for Unique in Unique_Primer_List:
        catch = 0
        catch2 = 0
        unique1 = 0
        matching_seq_id = []
        mismatch_seq_id = []
        combined_matching_seq_id = []
        combined_mismatching_seq_id = []
        for primer in RPA_Primer:
            if Unique == primer:
                file_name = str(RPA_Primer[catch]) + "_" + str(Amplicon[catch]) + "_" + str(crRNA[catch]) + ".csv"
                #print(file_name)
                os.chdir(temp_dir)
                current_file = pd.read_csv(file_name)
                columns2 = list(current_file.columns.values.tolist()) 
                #print(current_file.head())
                Exact = current_file[current_file.columns[0]]
                Mis3 = current_file[current_file.columns[1]]
                Bad = current_file[current_file.columns[2]]
                for value in Exact:
                    matching_seq_id.append(value)
                for value in Mis3:
                    mismatch_seq_id.append(value)  
                os.chdir(current)
                if catch2 == 0:
                    Unique_crRNA1.append(crRNA[catch])
                if catch2 == 1:
                    Unique_crRNA2.append(crRNA[catch])
                    Unique_primer.append(RPA_Primer[catch])
                    Unique_amp.append(Amplicon[catch])
                catch2 += 1
            catch += 1     
        combined_matching_seq_id = list(dict.fromkeys(matching_seq_id))
        combined_matching_seq_id = [x for x in combined_matching_seq_id if str(x) != 'nan']
        combined_mismatching_seq_id = list(dict.fromkeys(mismatch_seq_id))
        combined_mismatching_seq_id = [x for x in combined_mismatching_seq_id if str(x) != 'nan']
        combined_all_seq_id = combined_matching_seq_id + combined_mismatching_seq_id
        combined_all_seq_id = list(dict.fromkeys(combined_all_seq_id))
        Unique_combined_all.append(len(combined_all_seq_id))
        Unique_combined_match.append(combined_matching_seq_id)
        Unique_combined_list.append(list(combined_all_seq_id))
    best_all = 0 # value for best primer
    Current_Best_Amp = 0 # The current best amplicon size
    Second_Best_Amp = 0
    second_best = 0
    best_primer = ""
    second_best_primer = ""
    Current_Unique_combined_all = ""
    check = 0
    Best_list = []
    Second_best_list = []
    crRNA1 = ""
    crRNA2 = ""
    for values in Unique_primer:
        Current_Unique_crRNA1 = Unique_crRNA1[check]
        Current_Unique_crRNA2 = Unique_crRNA2[check]
        Current_Unique_crRNApair = Current_Unique_crRNA1 + "," +Current_Unique_crRNA2
        Current_Unique_amp = Unique_amp[check]
        Current_Unique_primer = values
        Current_Unique_combined_all = Unique_combined_all[check]
        Current_Combined_list = Unique_combined_list[check]
        Current_Unique_combined_match = str(Unique_combined_match[check])
    
        if Current_Unique_combined_all >= best_all:
            if Current_Unique_combined_all == best_all:
                if Current_Unique_amp < Current_Best_Amp:
                   current_best_primer = Current_Unique_primer
                   Current_Best_Amp = Current_Unique_amp
                   crRNA1 = Current_Unique_crRNApair
                   best_all = Current_Unique_combined_all
                   Best_List = Current_Combined_list
     
            if Current_Unique_combined_all > best_all:
                  current_best_primer = Current_Unique_primer
                  Current_Best_Amp = Current_Unique_amp 
                  crRNA1 = Current_Unique_crRNApair
                  best_all = Current_Unique_combined_all
                  Best_List = Current_Combined_list
            
        if Current_Unique_combined_all < best_all:
            if Current_Unique_combined_all >= second_best:
                if Current_Unique_amp < Second_Best_Amp:
                    second_best = Current_Unique_combined_all
                    Second_Best_Amp = Current_Best_Amp
                    crRNA2 = Current_Unique_crRNApair
                    second_best_primer = Current_Unique_primer
                    Second_best_list = Current_Combined_list
            
            if Current_Unique_combined_all > second_best:
                second_best = Current_Unique_combined_all
                Second_Best_Amp = Current_Best_Amp
                crRNA2 = Current_Unique_crRNApair
                second_best_primer = Current_Unique_primer
                Second_best_list = Current_Combined_list
        check = check + 1
    print(current_best_primer + " " + str(Current_Best_Amp) + " " + str(best_all))
    #print(Best_List)
    print(second_best_primer + " " + str(Second_Best_Amp) + " " + str(second_best))
    #print(Second_best_list)
    
    on_dir = working_dir + "/On"
    files = os.listdir(on_dir)
    Bad_strain = []
    strains = []
    for file in files:
        strains.append(file)
    for strain in strains:
        if strain not in Best_List:
            Bad_strain.append(strain)
    #print(len(Bad_strain))
    #print(Bad_strain)
    
    
    Bad_strain2 = []
    strains = []
    for file in files:
        strains.append(file)
    for strain in strains:
        if strain not in Second_best_list:
            Bad_strain2.append(strain)
    #print(len(Bad_strain2))
    #print(Bad_strain2)
    
    
    os.chdir(working_dir)
    with open("Final_Output.csv", "w") as output_file:
        output_file.write("PrimedSherlock Results" + "," + "\n")
        output_file.write("Best RPA Primer Set" + "," + current_best_primer + "," + "\n" + "Amplicon Size" + "," + str(Current_Best_Amp) + "," + "\n" + "Total Matches within 4 bases" + "," + str(best_all) + "," + "\n")
        output_file.write("Best crRNA" + "," + crRNA1 + "," + "\n")
        output_file.write("Input strains not covered" + "," + "\n")
        count = 0
        for strain in Bad_strain:
            if count <= 2:
                output_file.write(",") 
            if count >= 2:
                output_file.write(strain + ",")            
            if count == 19:
                output_file.write("\n")
                count = 0
            count = count + 1
        output_file.write("\n" + "\n" + "Second Best RPA Primer Set" + "," + second_best_primer + "," + "\n" + "Amplicon Size" + "," + str(Second_Best_Amp) + "," + "\n" + "Total Matches within 4 bases" + "," + str(second_best) + "," + "\n")
        output_file.write("Second Best crRNA" + "," + crRNA2 + "," + "\n")
        output_file.write("Input strains not covered" + "," + "\n")
        count = 0
        for strain in Bad_strain2:
    
            if count <= 2:
                output_file.write(",") 
            if count >= 2:
                output_file.write(strain + ",")            
            if count == 19:
                output_file.write("\n")
                count = 0
              
            count = count + 1
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
