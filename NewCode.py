# -*- coding: utf-8 -*-
#%%
if __name__ == '__main__':
    
    from timeit import default_timer as timer
    start = timer()
    import subprocess
    import sys
    import Bio
    import pandas as pd
    import os
    import shutil
    from configparser import ConfigParser
    #Starts by figuring out where the user is running the script from.
    #This is so we can setup the config file and move around the various script created folders
    import pathlib 
    ignition_dir = (pathlib.Path(__file__).parent.resolve()) 
    os.chdir(ignition_dir)
    #Now We Can get into creating the config and running the initial script.
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
    
    master_dir = os.getcwd() #Figures out the directory of the masterfile
    ##############################################################################
    #This section goes over the config file and sets up the file for user editing
    ##############################################################################
    print('[Primed Sherlock] Setting Working, Primed_RPA_Output and Input.fna directories' )
    print('[Primed Sherlock] Set as ' + '"' + str(master_dir) + '"')
    print('[Primed Sherlock] Do not edit the config to move from locations')
    
    #This looks for the config file. If it doesn't exist it promts the user to set it up.
    try:
        config = ConfigParser()
        config.read('config.ini')
        theta = config.get('main','Amplicon_Max_Size')
        print("[Primed Sherlock] Found valid Config File")
    except:
        print("[Primed Sherlock] First Time Execution, Unable to find config.ini, Fill out newly created config.ini. Instructions in Readme")
        config.read('config.ini')
        config.add_section('main')
        config.set('main', 'MismatchThreshold', '5')
        config.set('main', 'Use_Primed_RPA', 'Yes')
        config.set('main', 'Thread_Count', '23')
        config.set('main', 'Amplicon_Max_Size', '250')
        config.set('main', 'Amplicon_Min_Size', '150')
        config.set('main', 'Max_Background_Identitiy', '65')
        config.set('main', 'Off_Target_Mismatch_Score', '7')
        config.set('main', 'Contingency_crRNA_Utilization', 'Yes')
        config.set('main', 'Working_Directory', str(master_dir))
        config.set('main', 'Primed_RPA_Output', str(master_dir))
        config.set('main', 'Input.fna_Location',str(master_dir))
        config.set('main', 'Background_Blast_File', str(master_dir))
        config.add_section('Features')
        config.set('Features', 'Run_Cas-off', 'Yes')
        with open('config.ini', 'w') as f:
            config.write(f)
        sys.exit()
        
    #This pulls some of the config information for usage by the script.
    Multithread_Count = config.get('main', 'Thread_Count')
    working_dir = config.get('main', 'Working_Directory') 
    #You'll notice that I request config locations multiple times. 
    #This is because I (#%%) partition the code to run individual subsections in Spyder
    #Leaving it in doesn't do anything negative and allows for the end-user to troubleshoot
    
    
if __name__ == '__main__':
    #This program can take up 100GB of space depending on input and output files. 
    #We try to one by one delete old directories which are taking up harddrive space.
    #This isn't needed as it's done at the end, but is a redundant safeguard if the script crashes
    
    import shutil as sh #imported for rmtree function
    working_dir = config.get('main', 'Working_Directory')
    UpperLimitMismatch = config.get('main', 'MismatchThreshold')
    UpperLimitMismatch = int(UpperLimitMismatch)
    catch = 0
    #Lets setup the used directories and delete them if they still exist.
    dir_list = [ "ContingencycrRNA_Pairs", "PrimerVer", "PrimerCheck", "Cas-OffInput", "crRNA_Pairs", "PrimerComboFasta", "TempAlignDir", "TempDIR", 'On', 'Off', 'crRNA_Pairs', 'Cas-On-Temp', 'Cas-On-Out', 'Cas-ONInput', 'Cas-On']
    gen_file_list = ["MasterAnalysis.csv", "MasterFile.csv", "ContingencyMasterAnalysis.csv", "Final_Output.csv"]
    os.chdir(working_dir)
    directories = [f for f in os.listdir('.') if os.path.isdir(f)]
    directoryfiles = [f for f in os.listdir('.') if os.path.isfile(f)]
    for d in directories:
        if d in dir_list:
            try:
                sh.rmtree(str(d))
                print("[Primed Sherlock] Found existing " + d + " Directory, deleting directory")
            except:
                print("[Primed Sherlock] " + d + " directory could not be deleted please ensure you manually delete it")
    for g in directoryfiles:
        if g in gen_file_list:
            try:
                os.remove(str(g))
                print("[Primed Sherlock] Found existing " + g + " file, deleting file")
            except:
                print("[Primed Sherlock] " + g + " file could not be deleted please ensure you manually delete it")
                print("[Primed_Sherlock] If the file is open, save it in a different location")
    #Okay now we've established that there are no troublesome directories left or prompted the user to manually delete.
    catch = 0 #This value is reused elsewhere

if __name__ == '__main__':
    ###############################################################################
    #Extracting PrimedRPA File
    ###############################################################################    
    if str(config.get('Features', "Run_Cas-off")) == 'Yes':
        print("[Primed Sherlock] Prepping data for Cas-Offinder")
        PrimedOutputLoc = config.get('main','Primed_RPA_Output')
        os.chdir(PrimedOutputLoc)
        items = os.listdir(".")
        newlist = []
        
        
        
        Primed_RPA = config.get('main', 'Use_Primed_RPA') #Pulls user config to know if pulling an PrimedRPA Output Sets or Manual.
        for name in items:
        
       
            if(Primed_RPA != "Yes"):
                if name.endswith("Manual_Sets.csv"):
                    data = pd.read_csv("Manual_Sets.csv")
                    rpa_csv = open("Manual_Sets.csv")
            else:
               # for names in items:
                #    if names.endswith("Output_Sets.csv"):
                #        newlist.append(names)
              #  if len(newlist) >= 2:
               #     print("Ensure there is only one output file in the working directory, the first output file will be used otherwise")
                try:
                    data = pd.read_csv("Output_Sets.csv")
                    rpa_csv = open("Output_Sets.csv")
              #      data = pd.read_csv(newlist[0])
             #       rpa_csv = open(newlist[0], "r")
                except:
                    print("No Output file exists in directory: " + str(PrimedOutputLoc))
        #This section firstly pulls the user defined config parameters and implements them into the code such as amplicon size
        #and background identity, this is a variable pulled from PrimedRPA. For user defined primers, it is expected
        #that they understand if their primers are specific. I will implement a repeat value of 100 for user input. 
        #
        # Users will have to provide their own amplicon size as well, this may 
        
        theta = config.get('main','Amplicon_Max_Size')
        beta = config.get('main','Amplicon_Min_Size')
        maxback = config.get('main', 'Max_Background_Identitiy')
        #This next sectiond deals with how we put the Primer Data. Did we use PrimedRPA? 

        
        
        #If the user used their own primer sets. Here we assume they screened them.
        if(Primed_RPA != "Yes"):
            print("[Primed Sherlock] The script is configured to use user generated Primer Pairs. BEWARE! The program will not check for off-target amplification")
            #Pulls data from Output_Sets.csv
            columns = list(data.columns.values.tolist()) 
            Amp = data[data.columns[0]]
            FPrimer = data[data.columns[1]]
            RPrimer = data[data.columns[3]] 
            MaxBackground = 0
            FPrimerSite = data[data.columns[2]]
            RPrimerSite = data[data.columns[4]] 
            time_count = 0
            New_Amp = []
            New_FPrimer = []
            New_RPrimer = []
            New_MaxBackground = []
            NewFPrimerSite = []
            NewRPrimerSite = []
            for amplicons in Amp:
                if amplicons <= int(theta) and amplicons >= int(beta) :
                    if MaxBackground <=  int(maxback): #Again we assume they did the Background Check
                        New_Amp.append(Amp[time_count])
                        New_FPrimer.append(FPrimer[time_count])
                        New_RPrimer.append(RPrimer[time_count])
                        New_MaxBackground.append(MaxBackground)
                        NewFPrimerSite.append(FPrimerSite[time_count])
                        NewRPrimerSite.append(RPrimerSite[time_count])
                        #print(FPrimer[time_count]) #OCT9
                time_count = time_count + 1
                
        else:
            print("[Primed Sherlock] Importing Primed Sherlock Output.csv datasets")
            print("[Primed Sherlock] Using Output_Sets.csv as dataset source file")
            #Pulls data from Output_Sets.csv
            columns = list(data.columns.values.tolist()) 
            Amp = data[data.columns[0]]
            FPrimer = data[data.columns[3]]
            RPrimer = data[data.columns[12]] #10? #12
            MaxBackground = data[data.columns[6]] #6 This is supposed to be 6, I had it as 9???
            FPrimerSite = data[data.columns[1]]
            RPrimerSite = data[data.columns[10]] #10 #8?
            time_count = 0
            New_Amp = []
            New_FPrimer = []
            New_RPrimer = []
            New_MaxBackground = []
            NewFPrimerSite = []
            NewRPrimerSite = []
            BindingSiteNew = []
            BindingSiteLoc = []
            for amplicons in Amp:
                if amplicons <= int(theta) and amplicons >= int(beta) :
                   if MaxBackground[time_count] <=  int(maxback):
                    
                    BindingSiteNew.append(FPrimerSite[time_count])
                    BindingSiteLoc.append(FPrimer[time_count])
                    New_Amp.append(Amp[time_count])
                    New_FPrimer.append(FPrimer[time_count])
                    New_RPrimer.append(RPrimer[time_count])
                    New_MaxBackground.append(MaxBackground[time_count])
                    NewFPrimerSite.append(FPrimerSite[time_count])
                    NewRPrimerSite.append(RPrimerSite[time_count])
                    #print(NewRPrimerSite)
                time_count = time_count + 1   
    ##############################################################################
    #Creates Multifastas of each respective primersite
    # This section functions by getting the lenght of each primer, then essentially
    # cuts the amplicon from the region between the primers, excluding the primers
    # including them will lead to issues with crRNA design
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
        #print('Forward')
        #print(CCR)
        CutFix = CCR - CCF

        #print(CutFix)
        #print(CurrentReverse)
        #print(CurrentReverse + "this" + "RH")
        CutFix = int(CutFix) + len(CurrentReverse) #Cuts the primer to the length
        #print(CutFix)
        CurrentBackgroundScore = New_MaxBackground[n]
        FileName = str(CurrentForward + "_" + CurrentReverse + "_" + str(Amplicon_Size))
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
                #Fixes length
                # FPrimerLength = len(CurrentForward)
                # RPrimerLength = len(CurrentReverse)
                # totallength = len(finalsequence)
                # totallength = totallength - RPrimerLength
                # finalsequence = finalsequence[FPrimerLength:totallength]
                output_file.write(">" + seq_record.id + "\n" + finalsequence + "\n")
        #print(time_count , " primer pair amplicons extracted for input fasta." + "\n")
        seq_files.append(FileName)
        n = n + 1
    ##############################################################################
    #Makes Consensus of primer specific extracted input multifasta
    #The following getComplement definition was sourced from PrimedRPA
    #We aknowlege PrimedRPA and Mathew Higgins and utilize the code under agreements provided
    #in the open source license. 
    ##############################################################################

if __name__ == '__main__': 
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
    #Finds crRNA targets
    #TODO: Make it so 
    #############################################################################

if __name__ == '__main__':
    print("[Primed Sherlock] Finding best crRNA pairs to specifications in config")
    primer_delete_list = []
    delete_list = []
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
            #print(crRNA)
            crRNA2 = SeqUtils.nt_search(amplicon_sequence, pattern)     
            q = 1
            if q == 1:
            #try:
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
                        #UniquecrRNA.append(temp_crRNA) JAN5TH EDIT
                        if crRNASiteLocation[crRNA_count] >= 30 and crRNASiteLocation[crRNA_count] <= top_end:
                            UniquecrRNA.append(temp_crRNA)
                        merged_list = [(sequence_table[i], crRNASiteLocation[i]) for i in range (0, len(sequence_table))]
                        #UniquecrRNA = list(dict.fromkeys(merged_list))
                        #UniquecrRNA = list(dict.fromkeys(crRNAactual))
                    os.chdir(Input_File)
                    os.chdir("crRNA_Pairs")
                    UniquecrRNA = list(dict.fromkeys(UniquecrRNA))
                    for values in UniquecrRNA:
                        with open("TTT_" + primers + "_crRNA" + '.csv', "a") as output_file:
                            primer_delete_list.append("TTT_" + primers + "_crRNA" + '.csv')
                            #output_file = open("TTT_" + primers + "_crRNA" + '.csv', "a")
                            output_file.write(str(values + ","))
                            output_file.write(str("\n"))
                    #output_file.close()   
                    crRNA_count = crRNA_count + 1
                    os.getcwd()
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
                        #UniquecrRNA.append(temp_crRNA) EDIT JAN5
                        
                        if crRNASiteLocation[crRNA_count] >= 30 and crRNASiteLocation[crRNA_count] <= top_end:
                            UniquecrRNA.append(temp_crRNA)
                            #print(temp_crRNA)
                        merged_list = [(sequence_table[i], crRNASiteLocation[i]) for i in range (0, len(sequence_table))]
                        #UniquecrRNA = list(dict.fromkeys(merged_list))
                        #UniquecrRNA = list(dict.fromkeys(crRNAactual))
                    os.chdir(Input_File)
                    os.chdir("crRNA_Pairs")
                    UniquecrRNA = list(dict.fromkeys(UniquecrRNA))
                    for values in UniquecrRNA:
                        with open("AAA_" + primers + "_crRNA" + '.csv', "a") as output_file:
                            primer_delete_list.append("AAA_" + primers + "_crRNA" + '.csv')
                            #output_file = open("TTT_" + primers + "_crRNA" + '.csv', "a")
                            output_file.write(str(values + ","))
                            output_file.write(str("\n"))
                    #output_file.close()   
                    crRNA_count = crRNA_count + 1
                    os.getcwd()            
    primer_delete_list = list(dict.fromkeys(primer_delete_list))
    #print(str(primer_delete_list))
    dir_list = (os.listdir(os.getcwd()))
    
    for values in dir_list:
        if values not in primer_delete_list:
            delete_list.append(values)
    for values in delete_list:
        os.remove(values)
# =============================================================================
# 
#             #except: 
#                # print("\n")
#                # print('No crRNA site + ')  
#             try:
#                 crRNA2.pop(0)
#                 UniquecrRNA = []
#                 crRNA_count = 0
#                 for targetsite in crRNA2:
#                     crRNAlocation = targetsite
#                     crRNAlocation = crRNAlocation
#                     length = len(amplicon_sequence) #length of code
#                     Cut1Done = amplicon_sequence[(int(crRNAlocation)): (length)]
#                     #crRNAactual.append(Cut1Done[0:24])
#                     temp_crRNA = str(Cut1Done[0:24])
#                     
#                     if temp_crRNA not in sequence_table:
#                         crRNASiteLocation.append(targetsite)
#                         sequence_table.append(temp_crRNA)
#                         top_end = length
#                         if crRNASiteLocation[crRNA_count] >= 30 and crRNASiteLocation[crRNA_count] <= top_end:
#                             UniquecrRNA.append(temp_crRNA)
#                         merged_list = [(sequence_table[i], crRNASiteLocation[i]) for i in range (0, len(sequence_table))]
#                         #UniquecrRNA = list(dict.fromkeys(merged_list))
#                         #UniquecrRNA = list(dict.fromkeys(crRNAactual))
#                     os.chdir(Input_File)
#                     os.chdir("crRNA_Pairs")
#                     output_file = open("AAA_" + primers + "_crRNA" + '.csv', "a")
#                     UniquecrRNA = list(dict.fromkeys(UniquecrRNA))
#                     #print(UniquecrRNA)
#                     for values in UniquecrRNA:
#                         output_file.write(str(values + ","))
#                         output_file.write(str("\n"))
#                     crRNA_count = crRNA_count + 1
#             except: 
#                 print("\n")
#                 print('No crRNA site - ')  
# =============================================================================
    ###############################################################################
    #Grabs crRNA info
    #This following segment utilizes a third-party program CasOffinder
    #It first pulls all crRNAs discovered during the crRNA search process.
    #It then creates a seperate input file, containing an individaul crRNA
    #written in a specific format recognized by Cas-Offinder
    #This is then paired with the Cas-Offinder to search the genomic dataset for exact matches
    #
    #This is why multithreading is utilized, it creates an instance for each core which allows
    #multiple copies to be run at once.
    
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
        if len(value_csv) >= 1:
            os.chdir(Input_File)
            try:
                os.mkdir("Cas-OffInput")
                print("Made Cas-offInput in " + os.getcwd())
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
        #else:
            #print("[Primed Sherlock] No Values")     
        os.chdir(path)        
    ##############################################################################
    #Makes Cas-offinder off-target
    ##############################################################################
    
    
    
    
    ##############################################################################
    #Cleaning Background fasta files
    #If a fasta file contains to many missing nucleotides it is essentially useless
    #This finds genomic sequences with sufficient missing nucleotides and removes them
        
    ##############################################################################    
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
    #print(newlist) OCT10
    
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
    
    print("[Primed Sherlock] Running Cas-offinder using PrimedRPA results")
    print("[Primed Sherlock] This next section may take several hours or days to complete. See Multithreading in Config for settings to speed this up.")
    
#################################################################################
#Multiprocessing, Cas-Off
#Runs Cas-Off on background sequences
#################################################################################

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

# def Cas_off(x, crRNAdir, Original_files, cas_instance):
#     import os
#     import shutil as sh

#     #Checks if its the first instance of Cas "Off-Target" or second instance of Cas "On-Target".
#     if cas_instance == 1:    
#         number_files2 = Original_files
#         list3 = os.listdir(crRNAdir) # dir is your directory path
#         number_files = len(list3)
#         new_files = number_files - number_files2
#         left = number_files2 - new_files
#         print("[Primed Sherlock] Processed a total of " + str(new_files) + " CasOffinder Instances out of " + str(left))
#     if cas_instance == 2:
#         number_files2 = Original_files
#         list3 = os.listdir(crRNAdir) # dir is your directory path
#         number_files = len(list3)
#         new_files = number_files - number_files2
#         left = number_files2 - new_files
#         print("[Primed Sherlock] Processed a total of " + str(new_files) + " CasOffinder Instances out of " + str(left))
        
        
#     CREATE_NO_WINDOW = 0x08000000
#     subprocess.call("cas-offinder " + x + " G0 " + " " + x + "out.txt", creationflags=CREATE_NO_WINDOW )
#     #subprocess.call("cas-offinder " + x + " G0 " + " " + x + "out.txt")
#     return x

def Cas_off(x, y):
    import os
    import shutil as sh
    ################################################################################
    #Figures out total instances needed of Cas-Offinder
    ################################################################################
    path = os.getcwd()
    crRNAdir = path
    cas_instance = 1 #Setting the instance so we know we're handing the first instance see Cas_def for details. 
    #Checks if its the first instance of Cas "Off-Target" or second instance of Cas "On-Target".    

    number_files2 = y
    list3 = os.listdir(crRNAdir) # dir is your directory path
    number_files = len(list3)
    new_files = number_files - number_files2
    left = number_files2 - new_files
    print("[Primed Sherlock] Processed a total of " + str(new_files) + " CasOffinder Instances " + str(left) + " remaining")
        
        
    CREATE_NO_WINDOW = 0x08000000
    subprocess.call("cas-offinder " + x + " G0 " + " " + x + "out.txt", creationflags=CREATE_NO_WINDOW )
    #subprocess.call("cas-offinder " + x + " G0 " + " " + x + "out.txt"
    return x
###
#This next section sets the mutithreading for the Cas-Offinder
#Each thread is utilized to run a seperate copy of Cas-Offinder
#So for my 3900XT it utilizes 23 threads, 12 core/24threads 
#You want to ensure that this is less then your threadcount
###


if __name__ == '__main__':
    working_dir = config.get('main', 'Working_Directory')
    Input_File = working_dir
    #Input_File ='c:/Users/James_Mann1/Desktop/TestBench' #removed OCT8 DIRECTORYREPHASE
    Output_Dir = Input_File + "/" + "Cas-OffInput"
    
    list3 = os.listdir(Output_Dir) # dir is your directory path
    Original_files = len(list3)
    Input_path = Output_Dir
    Input_Casoff = cas_prep(Input_File, Output_Dir, Input_path)
    Iter_Value_List = []
    for values in Input_Casoff:
        Iter_Value_List.append(Original_files)
        zipped = zip(Input_Casoff, Iter_Value_List)
    #print(set(zipped))
    with Pool(int(Multithread_Count)) as p:
        #p.map(Cas_off, Input_Casoff) #OCT10
        #print(p.map(Cas_off, Input_Casoff))
        print(p.starmap(Cas_off, zipped))
        

        
    print("[Primed Sherlock] Done processing off-targets with Cas-Off Finder.")
##############################################################################
#After Primed RPA and Cas-OffFinder Is Ran, This checks for crRNA targets
#and finds those that have little to no matches to background sequences
##############################################################################

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
    #This code segment sets a threshold for off targets as well as discards crRNA that 
    #potentially could cause off-target activation of the enzyme, aka a false positive.
    #The threshold is a user defined value. It is provided this way, as the exact
    #mismatch threshold for causing a conformational change in Cas-12 / Cas-9 / Cas-13
    #enzymes is debatable. 
    #I find that 7 basepair mismatches hampers the ability to cause this conformational change
    #unpublished results, so i've subjectively set this as the default.
    #You can set this to any threshold.
    #As the code is currently written it does not take into account the distance of the mismatch to the seed region
    #the theoretical important region in which multiple mismatches signficantly hampers efficency. 
    #
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
        print("[Primed Sherlock] Found Correct Working Directory")
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
                #print(file_name)
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
    working_dir = config.get('main', 'Working_Directory')
    #working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' #FIX
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
        #print(os.getcwd())
        with open(value, "r") as input_file:
            for line in input_file:
                line_check += 1
            if line_check >= 3:
                value_length = len(value)
                #print(value)
                amount += 1
                #print(value)
                values_csv = pd.read_csv(value)
                #print(value_csv.head())
                values_csv.dropna()
                crRNA = values_csv[values_csv.columns[0]]
                crRNA = crRNA.tolist()
                crRNA = [x for x in crRNA if str(x) != 'nan']
                #print(crRNA)
                os.chdir(temp_on)
                #crRNA.dropna()
                #print(crRNA)
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
    
    import os
    working_dir = config.get('main', 'Working_Directory')
    #working_dir = 'C:/Users/James_Mann1/Desktop/Testbench'
    temp_on = working_dir + "/Cas-On-Temp"
    
    on_target_background_file_loc = working_dir #FIX
    Input_File = working_dir #FIX
    os.chdir(working_dir)
    os.chdir(temp_on)
    filepath = os.getcwd()
    crRNAfiles_on = os.listdir(filepath)
    for value in crRNAfiles_on:
        #print(value)
        value_csv = pd.read_csv(value)
        value_csv.dropna()
        #print(value_csv.head())
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
        #print(os.getcwd())
        os.chdir(temp_on)
        #print(os.getcwd())
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
    #print(newlist)
    
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
    ##############################################################################
    #Runs Cas off-finder for the on-target strains
    ##############################################################################
 
     #thi s 
if __name__ == '__main__':
    print("[PrimedSherlock] Running On-Target crRNA matching with Cas-Offinder")
    working_dir = config.get('main', 'Working_Directory')
    Input_File = working_dir
    #Input_File ='c:/Users/James_Mann1/Desktop/TestBench' OCT8
    Output_Dir = Input_File + str("/") + "Cas-ONInput"
    
    
    list3 = os.listdir(Output_Dir) # dir is your directory path
    Original_files = len(list3)
    Input_path = Output_Dir
    Input_Casoff = cas_prep(Input_File, Output_Dir, Input_path)
    Iter_Value_List = []
    for values in Input_Casoff:
        Iter_Value_List.append(Original_files)
        zipped = zip(Input_Casoff, Iter_Value_List)
    #print(set(zipped))
    with Pool(int(Multithread_Count)) as p:
        #p.map(Cas_off, Input_Casoff) #OCT10
        #print(p.map(Cas_off, Input_Casoff))
        print(p.starmap(Cas_off, zipped))
    ##############################################################################
    #Goes through the sequences
    ##############################################################################
    import os
    working_dir = config.get('main', 'Working_Directory')
    #working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' OCT8
    os.chdir(working_dir)
    path = os.getcwd()
    #print(path)
    if path.endswith("Cas-ONInput") == True:
        print("[Primed Sherlock] Correct Working Directory")
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
    print("[Primed Sherlock] Done Coverting Aquired Cas On-Target Data to CSV format")      
    ###############################################################################
    #Goes through ontarget sequences for data analysis
    ###############################################################################
if __name__ == '__main__': 
    print("[Primed Sherlock] Attempting to determine best pairs and write results to MasterAnalysis.csv")
    import os
    import csv
    import pandas as pd
    from itertools import zip_longest
    working_dir = config.get('main', 'Working_Directory')
    #working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' #FIX OCT8
    os.chdir(working_dir)
    pathtest = os.getcwd()
    path = os.getcwd()
    if pathtest.endswith("Cas-OnInput") == True:
        print("[Primed Sherlock] Code is halfway done with Contingency generation")
    else:
        try:
            path = path + '\Cas-OnInput'
            os.chdir(path)
        except:
            print(pathtest + "Fix working directory")
    try:
        os.mkdir("Best_Pairs")
    except:
        print("[Primed Sherlock] Temporary Directory Exists from Earlier Code, Changing Directory")
    current = os.getcwd()
    os.chdir("Best_Pairs")
    bestpairs = os.getcwd()
    os.chdir(current)
    try:
        os.mkdir("Temp")
    except:
        print("[Primed Sherlock] Temporary Directory Exists from Earlier Code, Changing Directory")
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
        
        #print(file_name)
        #This will quickly store mismatches 3 to 5 bp from crRNA for potential targeted amplicon approach. 
        crRNA_Mismatch_List =[] #Store crRNA
        crRNA_Mismatch_ListBoth = [] #Stores both     
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
                #NEED TO DELETE SEQS WHICH HAVE BAD INPUT
                Good_crRNA_List = []
                Uniquedetails = []
                usable_seq = []            
                for Unique in Unique_crRNA_List:  
                    primer_string = file_name + "_" + Unique + ".csv"
                    locfprimer = primer_string.split("_") #GRABS FPRIMER 
                    count_sum = 0
                    seq_check = 0
                    usable_crRNA = 0
                    crRNA_matching_strains_0 = []
                    crRNA_matching_strains_3 = []
                    crRNA_matching_strains_greater = []
                    crRNA_Mismatch3to5 = [] #Stores Strains
                    crRNA_Mismatch5to20 = []
                    locationcatch = 0
                    
                    currentbindingloc = [] #grabs locaiton info
                    currentbindingprimer = [] #grabs primer for location info
                    for FPRIMERLOC in BindingSiteLoc:
                        if FPRIMERLOC == locfprimer[0]:
                            currentbindingloc.append(locfprimer[0])
                            currentbindingprimer.append(BindingSiteNew[locationcatch])
                        locationcatch = locationcatch + 1
                    
                    currentbindinglocDict = currentbindingloc
                    BindingSiteNewDict = currentbindingprimer
                    for value in crRNA:                     
                        if Unique == value:
                            BindingMax = int(BindingSiteNewDict[0]) + 250
                            current_score = Score[count_sum]
                            if current_score == 0:
                                crRNA_matching_strains_0.append(SeqID[count_sum])
                                #print(os.getcwd())
                                #print(SeqID[count_sum] + " " + str(current_score))
                                usable_crRNA = 1
                            if current_score <= UpperLimitMismatch and current_score > 1: #SET USER VARIABLE #FIXME (was <=5 and > 1)
                                usable_crRNA = 1
                                crRNA_matching_strains_3.append(SeqID[count_sum])
                                #print(SeqID[count_sum])
                            if current_score >= 3 and current_score <=5: #Used for Contingency crRNA's
                                if Loc[count_sum] >= BindingSiteNewDict[0] and Loc[count_sum] <= BindingMax:
                                    #Gets nice list for 3rd crRNA
                                    crRNA_Mismatch3to5.append(SeqID[count_sum])                  
                            if current_score >= 5 and current_score <=7:
                                crRNA_matching_strains_greater.append(SeqID[count_sum])
                                #print(current_score)
                            if current_score >= 5: #Used for Contingency crRNA's
                                if Loc[count_sum] >= BindingSiteNewDict[0] and Loc[count_sum] <= BindingMax:
                                    crRNA_Mismatch5to20.append(SeqID[count_sum])
                                    #print(SeqID[count_sum] + " " + str(current_score))
                        count_sum = count_sum + 1
                    if usable_crRNA == 1:                    
                        os.chdir(temp_path)
                        ContingencyDir = os.getcwd()
                        with open(primer_string, "w") as output_file:
                            output_file.write("Matching Seq" + "," + "3 Mismatch" + "," + "More then 4 mismatch" + "," +  "Contingency crRNA" + "," +  "Contingency severe mismatch crRNA" + ",""\n")
                        #Ranks Unique Seqs on matching etc
                        crRNA_matching_strains_0 = list(dict.fromkeys(crRNA_matching_strains_0))
                        crRNA_matching_strains_3 = list(dict.fromkeys(crRNA_matching_strains_3))
                        crRNA_matching_strains_greater = list(dict.fromkeys(crRNA_matching_strains_greater))
                        crRNA_Mismatch3to5 = list(dict.fromkeys(crRNA_Mismatch3to5))
                        crRNA_Mismatch5to20 = list(dict.fromkeys(crRNA_Mismatch5to20))
                        unique_rank_ontarget = len(crRNA_matching_strains_0)
                        unique_rank_mismatch = len(crRNA_matching_strains_3)
                        unique_rank_mismatch_bad = len(crRNA_matching_strains_greater)
                        Master_List = zip_longest(crRNA_matching_strains_0, crRNA_matching_strains_3, crRNA_matching_strains_greater, crRNA_Mismatch3to5, crRNA_Mismatch5to20)
                        with open(primer_string, "a", newline = "") as output_file:
                            for row in Master_List:
                                writer = csv.writer(output_file)
                                writer.writerow(row)
                                #print(row)
                        os.chdir(path)
                        #Grabs exact matches  
    import os
    import pandas as pd
    

if __name__ == '__main__':
    working_dir = config.get('main', 'Working_Directory')
    temp_path = working_dir + "/Cas-ONInput/Temp"
    #temp_path = 'C:/Users/James_Mann1/Desktop/Testbench/Cas-ONInput/Temp' #FIX 
    #working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' #Fix OCT8  
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
    
    try:
        with open("MasterAnalysis.csv", "w") as output_file:
            output_file.write("RPA PRIMER SET" + "," + "Amplicon Size" + "," + "crRNA Target" + "," + "On Target" + "," + "Mismatch 3" + "," + "Bad" + "," + "\n")
        os.chdir(current)
        print("[Primed Sherlock] opened masteranalysis.csv")
    except:
        print("[Primed Sherlock] couldn't open masteranalysis.csv")
    try:
        with open("MasterAnalysis.csv", "w") as output_file:
            print("[Primed Sherlock] opened masteranalysis.csv")
            output_file.write("RPA PRIMER SET" + "," + "Amplicon Size" + "," + "crRNA Target" + "," + "On Target" + "," + "Mismatch 3" + "," + "Bad" + "," + "\n")
        os.chdir(current)
    except:
        print("[Primed Sherlock] couldn't open masteranalysis.csv")
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
    try:
        with open("MasterFile.csv", "w") as output_file:
            output_file.write("RPA PRIMER" + "," +  "AMPLICON" + "," + "Best crRNA" + "," + "Total Good" + "," + "MATCHING SEQ" + "," + "MISMATCH < 3" + "," + "MISMATCH > 7" + ","+ "\n")
    except:
        with open("MasterFile.csv", "w") as output_file:
            output_file.write("RPA PRIMER" + "," +  "AMPLICON" + "," + "Best crRNA" + "," + "Total Good" + "," + "MATCHING SEQ" + "," + "MISMATCH < 3" + "," + "MISMATCH > 7" + ","+ "\n")            
    os.chdir(current)
    print("[Primed Sherlock] Sometimes MasterFile.Csv doesnt want to open. It is working.")
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
                try:
                    with open("MasterFile.csv", "a") as output_file:
                        output_file.write(str(Unique_RPA) + "," + str(Best_Amplicon) + "," + str(Best_Unique) + "," + str(best_check) + "," + str(best_check1) + "," + str(best_check3) + "," + str(best_check4) + "," + "\n")
                        output_file.write(str(Unique_RPA) + "," + str(Second_Best_AMP) + ","+  str(Second_Best_Unique) + "," + str(Second_Best_check) + "," + str(Second_Best_check_1) + "," + str(Second_Best_check_3) + "," + str(Second_Best_check_4) + "," + "\n")
                except:
                    print("[Primed Sherlock] Unable to open Masterfile.csv on first try")
                    with open("MasterFile.csv", "a") as output_file:
                        output_file.write(str(Unique_RPA) + "," + str(Best_Amplicon) + "," + str(Best_Unique) + "," + str(best_check) + "," + str(best_check1) + "," + str(best_check3) + "," + str(best_check4) + "," + "\n")
                        output_file.write(str(Unique_RPA) + "," + str(Second_Best_AMP) + ","+  str(Second_Best_Unique) + "," + str(Second_Best_check) + "," + str(Second_Best_check_1) + "," + str(Second_Best_check_3) + "," + str(Second_Best_check_4) + "," + "\n")
                    
                
                os.chdir(current)
    

    ##############################################################################
    #Sorts through output and finds best primer pairs based on three variables   #
    # 1. Amount of exact matches to input files                                  #
    # 2. Amount of mismatches less then 3                                        #
    # 3. Size of amplicon. Smaller amplicon == more efficent RPA                 #
    ##############################################################################
if __name__ == '__main__':
    #working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' #Fix OCT8
    working_dir = config.get('main', 'Working_Directory')
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
    Best_List = ""
    current_best_primer = ""
    checksumfinalissue = 0
    for values in Unique_primer:
        #print(values)
        Current_Unique_crRNA1 = Unique_crRNA1[check]
        Current_Unique_crRNA2 = Unique_crRNA2[check]
        Current_Unique_crRNApair = Current_Unique_crRNA1 + "," + Current_Unique_crRNA2
        Current_Unique_amp = Unique_amp[check]
        Current_Unique_primer = values
        Current_Unique_combined_all = Unique_combined_all[check]
        Current_Combined_list = Unique_combined_list[check]
        Current_Unique_combined_match = str(Unique_combined_match[check])
        
        #print(Current_Unique_combined_all)
        #print(Current_Unique_amp)
        if Current_Unique_combined_all >= best_all:
            if Current_Unique_combined_all > best_all:
                  current_best_primer = Current_Unique_primer
                  Current_Best_Amp = Current_Unique_amp 
                  crRNA1 = Current_Unique_crRNApair
                  best_all = Current_Unique_combined_all
                  Best_List = Current_Combined_list
            if Current_Unique_combined_all == best_all:
                if Current_Unique_amp < Current_Best_Amp:
                    current_best_primer = Current_Unique_primer
                    Current_Best_Amp = Current_Unique_amp
                    crRNA1 = Current_Unique_crRNApair
                    best_all = Current_Unique_combined_all
                    Best_List = Current_Combined_list
            
            #To-do verify that this actually does what i intend it to do. This was made to fix,
            # an issue where, if the 1st value in values list was smaller, it wouldn't show up 
            # as the second best....
            if second_best == 0:
                    second_best = Current_Unique_combined_all
                    Second_Best_Amp = Current_Best_Amp
                    crRNA2 = Current_Unique_crRNApair
                    second_best_primer = Current_Unique_primer
                    Second_best_list = Current_Combined_list
        if Current_Unique_combined_all < best_all:
            if Current_Unique_combined_all == second_best:
                if Current_Unique_amp < Second_Best_Amp:
                    second_best = Current_Unique_combined_all
                    Second_Best_Amp = Current_Best_Amp
                    crRNA2 = Current_Unique_crRNApair
                    second_best_primer = Current_Unique_primer
                    Second_best_list = Current_Combined_list
            
            if Current_Unique_combined_all >= second_best:
                second_best = Current_Unique_combined_all
                Second_Best_Amp = Current_Unique_amp #Current_Best_Amp
                crRNA2 = Current_Unique_crRNApair
                second_best_primer = Current_Unique_primer
                Second_best_list = Current_Combined_list
        check = check + 1
    on_dir = working_dir + "/On"
    files = os.listdir(on_dir)
    Bad_strain = []
    strains = []
    for file in files:
        strains.append(file)
    for strain in strains:
        if strain not in Best_List:
            Bad_strain.append(strain)
    Bad_strain2 = []
    strains = []
    for file in files:
        strains.append(file)
    for strain in strains:
        if strain not in Second_best_list:
            Bad_strain2.append(strain)
    os.chdir(working_dir)
    try:
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
    except:
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
    ##############################################################################            
    # Activates crRNA generation for semi off-target strains                     #                                                                                                                                  #
    ##############################################################################
                    
    #ToDo Debug 1779 - Have to run segement twice to generate Cas-Offinder FILE
                    
                    
                    
#%%
if __name__ == '__main__':
    try:
        cycles = 0
        while cycles <=2:
            if cycles == 0:
                print("[Primed Sherlock] Generating output file for Cas-Offinder Input Run")
            if cycles == 2:
                print("[Primed Sherlock] Generating output file for Cas-OnFinder Input Run")
            config.read('config.ini') #Pulls config info
            Wintercontingency = config.get('main', 'contingency_crrna_utilization' ) #Do we generate a set of 2 crRNA for divergent strains?
            if Wintercontingency == "Yes":
                if cycles == 0:
                    print("[Primed Sherlock] Generating contingency crRNA targets.")
                Holding_Dir = os.getcwd() #Gets the current directory and holds it for a tiny bit.
                MasterMatch = []
                MasterContingency = []
                MasterContingencyOff = []
                ContingencyMatch = []
                ContingencycrRNA = []
                ContingencycrRNA2 = []
                UniqueContingencyMatch = []
                UniqueContingencycrRNA = []
                UniqueContingencycrRNA2 = []
                os.chdir(ContingencyDir) #Gets us to the directory where we can pull csv data for 5to7 and 5to20
                Contingency_files = os.listdir()
                for file in Contingency_files:
                    if file.startswith(str(current_best_primer)): 
                        crRNALIST = list(crRNA1.split(","))
                        if file.endswith(str(crRNALIST[0]) + ".csv"):
                            Contingency_Data = pd.read_csv(file)
                            columns = list(Contingency_Data.columns.values.tolist()) 
                            ContingencyMatch = Contingency_Data[Contingency_Data.columns[0]]
                            ContingencycrRNA = Contingency_Data[Contingency_Data.columns[3]]
                            ContingencycrRNA2 = Contingency_Data[Contingency_Data.columns[4]]             
                            for unique in ContingencyMatch:
                                UniqueContingencyMatch.append(unique)
                                UniqueContingencyMatch = list(dict.fromkeys(UniqueContingencyMatch))
                            for unique in ContingencycrRNA:
                                UniqueContingencycrRNA.append(unique)
                                UniqueContingencycrRNA = list(dict.fromkeys(UniqueContingencycrRNA))
                            for unique in ContingencycrRNA2:
                                UniqueContingencycrRNA2.append(unique)
                                UniqueContingencycrRNA2 = list(dict.fromkeys(UniqueContingencycrRNA2))
                            MasterMatch.append(UniqueContingencyMatch)
                            MasterContingency.append(UniqueContingencycrRNA)
                            MasterContingencyOff.append(UniqueContingencycrRNA2)
                    if file.startswith(str(current_best_primer)): 
                        crRNALIST = list(crRNA1.split(","))
                        if file.endswith(str(crRNALIST[1]) + ".csv"):
                            Contingency_Data = pd.read_csv(file)
                            Contingency_Data = Contingency_Data.dropna()
                            columns = list(Contingency_Data.columns.values.tolist()) 
                            ContingencyMatch = Contingency_Data[Contingency_Data.columns[0]]
                            ContingencycrRNA = Contingency_Data[Contingency_Data.columns[3]]
                            ContingencycrRNA2 = Contingency_Data[Contingency_Data.columns[4]]             
                            for unique in ContingencyMatch:
                                UniqueContingencyMatch.append(unique)
                                UniqueContingencyMatch = list(dict.fromkeys(UniqueContingencyMatch))
                                UniqueContingencyMatch = [x for x in UniqueContingencyMatch if pd.isnan(x) == False]
                            for unique in ContingencycrRNA:
                                UniqueContingencycrRNA.append(unique)
                                UniqueContingencycrRNA = list(dict.fromkeys(UniqueContingencycrRNA))
                                UniqueContingencycrRNA = [x for x in UniqueContingencycrRNA if pd.isnan(x) == False]
                            for unique in ContingencycrRNA2:
                                UniqueContingencycrRNA2.append(unique)
                                UniqueContingencycrRNA2 = list(dict.fromkeys(UniqueContingencycrRNA2))
                                UniqueContingencycrRNA2 = [x for x in UniqueContingencycrRNA2 if pd.isnan(x) == False]
                #Lets make a consensus amplicon from these.
                from Bio import SeqIO
                os.chdir(str(working_dir) + "/PrimerComboFasta" )
                fastafiles = os.listdir()
                contingency_fasta_seq = []
                output_file = "ContingencyCon.fna"
                primers = current_best_primer.split("_")
                FPRIMER_LEN = len(primers[0])
                RPRIMER_LEN = len(primers[1])
                for fasta in fastafiles:
                    if fasta.startswith(str(current_best_primer)):
                         for record in SeqIO.parse(fasta, 'fasta'):
                             for Contingency3crRNA in UniqueContingencycrRNA:
                                 if Contingency3crRNA == record.id:
                                    contingency_fasta_seq.append(record.seq)
                                    #Lets cut the amplicon not to include a full primer. Cas can detect from approximatly just 16 ssDNA basepairs...
                                    seqstring = str(record.seq)
                                    seqstringlength = len(seqstring)
                                    ReverseSite = int(seqstringlength) - int(RPRIMER_LEN) + 20
                                    final_seqstring = seqstring[FPRIMER_LEN:ReverseSite]
                                    with open(str(output_file), 'a') as inFile:
                                        inFile.write(">" + str(record.id) + "\n" + str(final_seqstring) + "\n")     
            #We gathered a list of amplicons that didn't share likely binding sites. So lets make a consensus of these
            seq_alignment = AlignIO.read(output_file, 'fasta')
            summary_align = AlignInfo.SummaryInfo(seq_alignment)
            summary_align.dumb_consensus(0.5)
            finalconsequence = summary_align.dumb_consensus(0.5)
            with open("ConsensusSeq.fasta", 'w') as inFile:
                inFile.write(">" + str(record.id) + "\n" + str(finalconsequence) + "\n")     
            crRNA = []
            crRNA2 = []
            leon = 0 
            files = os.listdir()
            for file in files:
                if file.startswith("ConsensusSeq"):
                    for seq_record in SeqIO.parse(file,'fasta'):        
                        output_file = open("TTT_" + current_best_primer + "_crRNA" + '.csv', "w")
                        output_file.write(str(current_best_primer + ","))
                        output_file.write("\n")     
                        output_file = open("AAA_" + current_best_primer + "_crRNA" + '.csv', "w")
                        output_file.write(str(current_best_primer + ","))
                        output_file.write("\n")
                        sequence_table = []
                        seq_id = seq_record.id
                        amplicon_sequence = str(seq_record.seq.ungap("-").upper())
                        amplicon_sequence_complement = getComplement(amplicon_sequence, False, 'N2N')
                        pattern = "TTTNNNNNNNNNNNNNNNNNNNNN"
                        crRNA = SeqUtils.nt_search(amplicon_sequence, pattern)
                        pattern = "NNNNNNNNNNNNNNNNNNNNNAAA"
                        #print(crRNA)
                        crRNA2 = SeqUtils.nt_search(amplicon_sequence, pattern)    
                        q = 1
                    if q == 1:
                    #try:
                        primers = current_best_primer
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
                                #UniquecrRNA.append(temp_crRNA) JAN5TH EDIT
                                if crRNASiteLocation[crRNA_count] >= 30 and crRNASiteLocation[crRNA_count] <= top_end:
                                    UniquecrRNA.append(temp_crRNA)
                                merged_list = [(sequence_table[i], crRNASiteLocation[i]) for i in range (0, len(sequence_table))]
                                #UniquecrRNA = list(dict.fromkeys(merged_list))
                                #UniquecrRNA = list(dict.fromkeys(crRNAactual))
                            os.chdir(Input_File)
                            try:
                                os.mkdir("ContingencycrRNA_Pairs")
                            except:
                                catch = 1
                            os.chdir("ContingencycrRNA_Pairs")
                            UniquecrRNA = list(dict.fromkeys(UniquecrRNA))
                            for values in UniquecrRNA:
                                with open("TTT_" + primers + "_" + str(seqstringlength) + '.csv', "a") as output_file:
                                    primer_delete_list.append("TTT_" + primers + "_crRNA" + '.csv')
                                    #output_file = open("TTT_" + primers + "_crRNA" + '.csv', "a")
                                    output_file.write(str(values + ","))
                                    output_file.write(str("\n"))
                            #output_file.close()   
                            crRNA_count = crRNA_count + 1
                            os.getcwd()
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
                                    #print(temp_crRNA)
                                merged_list = [(sequence_table[i], crRNASiteLocation[i]) for i in range (0, len(sequence_table))]
                                #UniquecrRNA = list(dict.fromkeys(merged_list))
                                #UniquecrRNA = list(dict.fromkeys(crRNAactual))
                            os.chdir(Input_File)
                            os.chdir("ContingencycrRNA_Pairs")
                            UniquecrRNA = list(dict.fromkeys(UniquecrRNA))
                            for values in UniquecrRNA:
                                with open("AAA_" + primers + "_" + str(seqstringlength) + '.csv', "a") as output_file:
                                    primer_delete_list.append("AAA_" + primers + "_crRNA" + '.csv')
                                    #output_file = open("TTT_" + primers + "_crRNA" + '.csv', "a")
                                    output_file.write(str(values + ","))
                                    output_file.write(str("\n"))
                            #output_file.close()   
                            crRNA_count = crRNA_count + 1
            path = os.getcwd()
            crRNAfiles = os.listdir(path)
            try:
                os.chdir(Input_File)
                sh.rmtree("Cas-OffInput")
                print("[Primed Sherlock] Removed Cas-Offinder Directory")
            except:
                catch = 1
            os.chdir(path)
            for value in crRNAfiles:
                value_csv = pd.read_csv(value)
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
            working_dir = config.get('main', 'Working_Directory')
            Input_File = working_dir
            #Input_File ='c:/Users/James_Mann1/Desktop/TestBench' #removed OCT8 DIRECTORYREPHASE
            Output_Dir = Input_File + "/" + "Cas-OffInput"
            os.chdir(Output_Dir)
            list3 = os.listdir(os.getcwd()) # dir is your directory path
            Original_files = len(list3)
            Input_path = Output_Dir
            Input_Casoff = cas_prep(Input_File, Output_Dir, Input_path)
            Iter_Value_List = []
            for values in Input_Casoff:
                Iter_Value_List.append(Original_files)
                zipped = zip(Input_Casoff, Iter_Value_List)
            #print(set(zipped))
            with Pool(int(Multithread_Count)) as p:
                #p.map(Cas_off, Input_Casoff) #OCT10
                #print(p.map(Cas_off, Input_Casoff))
                print(p.starmap(Cas_off, zipped))
                
            if cycles == 2:
                print("[Primed Sherlock] Done processing off-targets with Cas-Off Finder.")
                
            import pandas as pd
            #working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' (REMOVE)
            working_dir = config.get('main', 'Working_Directory')
            os.chdir(working_dir)
            path = os.getcwd() #Checks for sectional code piece
            if path.endswith("Cas-OffInput") == True:
                if cycles == 2:
                    print("[Primed Sherlock] Correct Working Directory Loaded")
            else:
                try:
                    path = path + '\Cas-OffInput'
                    os.chdir(path)
                    files = os.listdir(path)
                except:
                    catch = 1
            files = os.listdir(path) #Searches all the files
            if cycles == 2:
                print("[Primed Sherlock] Coverting Aquired Cas Off-Target Data to CSV format")
            for file in files:
                if file.endswith('txtout.txt') == True:
                    file_out = file + ".csv"
                    with open(file) as infile, open(file_out, 'w') as outfile:
                        outfile.write("Cas_off-Results" + "," + "SeqID" + "," + "Location" + "," + "Mismatch" + "," + "MismatchSeq" + "," + "Score" + "," + "\n")
                        for line in infile:
                              outfile.write(" ".join(line.split()).replace(' ', ','))
                              outfile.write("," + "\n") # trailing comma shouldn't matter 
            if cycles == 2:
                print("[Primed Sherlock] Program has finished converting Cas-Off results to CSV for analysis")
            
            catch = 0
            import os
            import pandas as pd
            user_defined_value = config.get('main', 'Off_Target_Mismatch_Score') #This is the value for Cas-Offinder's mismatch, goes after crRNA so e.g. crRNA 7
            user_defined_value = int(user_defined_value)
            #user_defined_value = 7 # defines value for specificity, less then or equal to. (7 normally, used for Cas-off output) (Remove if works)
            #exception_count = 0
            pathtest = os.getcwd()
            path = os.getcwd()
            if pathtest.endswith("Cas-OffInput") == True:
                if cycles == 2:
                    print("[Primed Sherlock] Found Correct Working Directory")
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
                        #print(file_name)
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
            if cycles == 2:
                print("[Primed Sherlock] Removing temporary Background CSV Files")
            # if path.endswith("Cas-OffInput") == True:
                
            #     print("Correct Working Directory")
            # else:
            #     try:
            #         path = path + '\Cas-OffInput'
            #         os.chdir(path)
            #         files = os.listdir(path)
            #     except:
            #         print(path + "Fix working directory")
            
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
            working_dir = config.get('main', 'Working_Directory')
            #working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' #FIX
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
                sh.rmtree(temp_on)
                os.mkdir(temp_on)
            except:
        
                catch = 1
            os.chdir(CombosPath)
            path = os.getcwd()
            crRNAfiles = os.listdir(path)
            screening_list = []
            amount = 0
            for value in crRNAfiles:
                #print(value)
                line_check = 0 #Checks amount of crRNA targets
                AAA_Check = 0
                TTT_Check = 0
                with open(value, "r") as input_file:
                    for line in input_file:
                        line_check += 1
                    if line_check >= 3:
                        value_length = len(value)
                        #print(value)
                        amount += 1
                        #print(value)
                        values_csv = pd.read_csv(value)
                        #print(value_csv.head())
                        values_csv.dropna()
                        crRNA = values_csv[values_csv.columns[0]]
                        crRNA = crRNA.tolist()
                        crRNA = [x for x in crRNA if str(x) != 'nan']
                        #print(crRNA)
                        os.chdir(temp_on)
                        #crRNA.dropna()
                        #print(crRNA)
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
            
            import os
            working_dir = config.get('main', 'Working_Directory')
            #working_dir = 'C:/Users/James_Mann1/Desktop/Testbench'
            temp_on = working_dir + "/Cas-On-Temp"
            on_target_background_file_loc = working_dir #FIX
            Input_File = working_dir #FIX
            os.chdir(working_dir)
            os.chdir(temp_on)
            filepath = os.getcwd()
            os.chdir(working_dir + "/" + "Cas-ONInput")
            crRNAfiles_on = os.listdir(filepath)
            ##Lets clean-up the previous Cas-ONInput
            originalcrRNAfiles = os.listdir()
            for leftover in originalcrRNAfiles:
                if leftover.endswith("Temp") != True:
                    if leftover.endswith("Best_Pairs") != True:
                        if leftover.endswith(".exe") != True:
                            os.remove(leftover)
            tempholder = os.getcwd()
            os.chdir("Temp")
            tempdirholder = os.listdir()
            for tempfile in tempdirholder:
                if tempfile.endswith("sis.csv") != True:
                    os.remove(tempfile)
            os.chdir(tempholder)
            os.chdir(filepath)
            for value in crRNAfiles_on:
                value_csv = pd.read_csv(value)
                value_csv.dropna()
                if len(value_csv) >= 1:
                    os.chdir(Input_File)
                    placeholderdir = os.getcwd()
                    try:
                        os.mkdir("Cas-ONInput")
                    except:
                        os.chdir("Cas-ONInput")
                        catch = 1
                    shortend_value = value.find(".csv")
                    short = value[0:int(shortend_value)]
                    filename_crRNA = str(short) + "input" + ".txt"
                    with open(filename_crRNA, 'w') as input_file:
                        input_file.write(str(working_dir) + "/on" + "\n")
                    rows = value_csv[value_csv.columns[0]]
                    if value.startswith("TTT_"):
                        with open(filename_crRNA, 'a') as input_file:
                            input_file = open(str(short) + "input" + ".txt", "a")
                            input_file.write(str("TTTNNNNNNNNNNNNNNNNNNNNNNNNN") + "\n")
                    if value.startswith("AAA_"):
                        with open(filename_crRNA, 'a') as input_file:
                            input_file = open(str(short) + "input" + ".txt", "a")              
                            input_file.write(str("NNNNNNNNNNNNNNNNNNNNNNNNNAAA") + "\n")
                    currentRNA = []
                    for crRNA in rows:
                        currentRNA.append(crRNA)
                    currentRNA = list(dict.fromkeys(currentRNA))
                    for crRNA in currentRNA:
                        with open(filename_crRNA, 'a') as input_file:
                            input_file.write(str(crRNA) + " 10" + "\n")
                os.chdir(temp_on)
            cycles +=1
            Contingency_Check = 0
    except:
        print("[Primed Sherlock] Was unable to generate contingency crRNA. This can happen readily if small amplicons don't contain binding site")
        Contingency_Check = 1 # Used to figure out if to run next segment or not
    ##############################################################################
    #Runs Cas off-finder for the on-target strains
    ##############################################################################
    #%%
if __name__ == '__main__':
    if Contingency_Check == 0:
        print("[Primed Sherlock] Was able to generate needed crRNA for on-target analysis")
        cycles = 0
        while cycles <= 2:
            if cycles == 0:
                print("[Primed Sherlock] This next code segement is designed to run in a 3x cycle. This due to sometimes a single run not generating an output file")
            if cycles == 2:
                print("[PrimedSherlock] Running On-Target crRNA matching with Cas-Offinder")
            working_dir = config.get('main', 'Working_Directory')
            Input_File = working_dir
            #Input_File ='c:/Users/James_Mann1/Desktop/TestBench' OCT8
            Output_Dir = Input_File + str("/") + "Cas-ONInput"
            Output_Dir = str(Output_Dir)
            txt_files =[]
            list3 = os.listdir(Output_Dir) # dir is your directory path
            for item3 in list3:
                if item3.endswith("txt") == True:
                    txt_files.append(item3)
            Original_files = len(txt_files)
            Input_path = Output_Dir
            Input_Casoff = cas_prep(Input_File, Output_Dir, Input_path)
            Iter_Value_List = []
            
            for values in Input_Casoff:
                Iter_Value_List.append(Original_files)
                zipped = zip(Input_Casoff, Iter_Value_List)
            with Pool(int(Multithread_Count)) as p:
                #p.map(Cas_off, Input_Casoff) #OCT10
                #print(p.map(Cas_off, Input_Casoff))
                print(p.starmap(Cas_off, zipped))
                
            ##############################################################################
            #Goes through the sequences
            ##############################################################################
            import os
            working_dir = config.get('main', 'Working_Directory')
            #working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' OCT8
            os.chdir(working_dir)
            path = os.getcwd()
            #print(path)
            if path.endswith("Cas-ONInput") == True:
                print("[Primed Sherlock] Correct Working Directory")
            else:
                try:
                    path = path + '\Cas-ONInput'
                    os.chdir(path)
                    files = os.listdir(path)
                except:
                    print(path + "Fix working directory")
            files = os.listdir(path)
            #Searches all the files
            if cycles == 2:
                print("[Primed Sherlock] Coverting Aquired Cas On-Target Data to CSV format")
            for file in files:
                if file.endswith('txtout.txt') == True:
                    file_out = file + ".csv"
                    with open(file) as infile, open(file_out, 'w') as outfile:
                        outfile.write("Cas_off-Results" + "," + "SeqID" + "," + "Location" + "," + "Mismatch" + "," + "MismatchSeq" + "," + "Score" + "," + "\n")
                        for line in infile:
                              outfile.write(" ".join(line.split()).replace(' ', ','))
                              outfile.write("," + "\n") # trailing comma shouldn't matter  
            
            if cycles == 2:
                    print("[Primed Sherlock] Done Coverting Aquired Cas On-Target Data to CSV format")      
            ###############################################################################
            #Goes through ontarget sequences for data analysis
            ###############################################################################
            if cycles == 2:
                print("[Primed Sherlock] Attempting to determine best pairs and write results to MasterAnalysis.csv")
            import os
            import csv
            import pandas as pd
            from itertools import zip_longest
            working_dir = config.get('main', 'Working_Directory')
            #working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' #FIX OCT8
            os.chdir(working_dir)
            pathtest = os.getcwd()
            path = os.getcwd()
            if pathtest.endswith("Cas-OnInput") == True:
                catch = 1
            else:
                try:
                    path = path + '\Cas-OnInput'
                    os.chdir(path)
                except:
                    catch = 1
            try:
                os.mkdir("Best_Pairs")
            except:
                catch = 1
            current = os.getcwd()
            os.chdir("Best_Pairs")
            bestpairs = os.getcwd()
            os.chdir(current)
            try:
                os.mkdir("Temp")
            except:
                catch = 1
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
            with open("ContingencyMasterFile.csv", "w") as output_file:
                output_file.write("Primer" + "," + "crRNA" + "," + "Exact Match" + "," + "Matches within 3" + "," + "\n")    
            os.chdir(current)
            for file in files:
                best_check = 0 
                Second_Best_check = 0
                best_check3 = 0
                Best_Unique = ""
                Second_Best_Unique = ""
                check_if_exists = 0
                #This will quickly store mismatches 3 to 5 bp from crRNA for potential targeted amplicon approach. 
                crRNA_Mismatch_List =[] #Store crRNA
                crRNA_Mismatch_ListBoth = [] #Stores both     
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
                        #NEED TO DELETE SEQS WHICH HAVE BAD INPUT
                        Good_crRNA_List = []
                        Uniquedetails = []
                        usable_seq = []            
                        for Unique in Unique_crRNA_List:  
                            primer_string = file_name + "_" + Unique + ".csv"
                            locfprimer = primer_string.split("_") #GRABS FPRIMER 
                            count_sum = 0
                            seq_check = 0
                            usable_crRNA = 0
                            crRNA_matching_strains_0 = []
                            crRNA_matching_strains_3 = []
                            crRNA_matching_strains_greater = []
                            crRNA_Mismatch3to5 = [] #Stores Strains
                            crRNA_Mismatch5to20 = []
                            locationcatch = 0
                            
                            currentbindingloc = [] #grabs locaiton info
                            currentbindingprimer = [] #grabs primer for location info
                            for FPRIMERLOC in BindingSiteLoc:
                                if FPRIMERLOC == locfprimer[0]:
                                    currentbindingloc.append(locfprimer[0])
                                    currentbindingprimer.append(BindingSiteNew[locationcatch])
                                locationcatch = locationcatch + 1
                            
                            currentbindinglocDict = currentbindingloc
                            BindingSiteNewDict = currentbindingprimer
                            for value in crRNA:                     
                                if Unique == value:
                                    BindingMax = int(BindingSiteNewDict[0]) + 250
                                    current_score = Score[count_sum]
                                    if current_score == 0:
                                        crRNA_matching_strains_0.append(SeqID[count_sum])
                                        #print(os.getcwd())
                                        #print(SeqID[count_sum] + " " + str(current_score))
                                        usable_crRNA = 1
                                    if current_score <= UpperLimitMismatch and current_score >= 1: #SET USER VARIABLE #FIXME (was <=5 and > 1)
                                        usable_crRNA = 1
                                        crRNA_matching_strains_3.append(SeqID[count_sum])
                                        #print(SeqID[count_sum])
                                    if current_score >= 3 and current_score <=5: #Used for Contingency crRNA's
                                        if Loc[count_sum] >= BindingSiteNewDict[0] and Loc[count_sum] <= BindingMax:
                                            #Gets nice list for 3rd crRNA
                                            crRNA_Mismatch3to5.append(SeqID[count_sum])                  
                                    if current_score >= 5 and current_score <=7:
                                        crRNA_matching_strains_greater.append(SeqID[count_sum])
                                        #print(current_score)
                                    if current_score >= 5: #Used for Contingency crRNA's
                                        if Loc[count_sum] >= BindingSiteNewDict[0] and Loc[count_sum] <= BindingMax:
                                            crRNA_Mismatch5to20.append(SeqID[count_sum])
                                            #print(SeqID[count_sum] + " " + str(current_score))
                                count_sum = count_sum + 1
                            if usable_crRNA == 1:                    
                                os.chdir(temp_path)
                                ContingencyDir = os.getcwd()
                                with open(primer_string, "a") as output_file:
                                    output_file.write("Matching Seq" + "," + "3 Mismatch" + "," + "More then 4 mismatch" + "," +  "Contingency crRNA" + "," +  "Contingency severe mismatch crRNA" + ",""\n")
                                #Ranks Unique Seqs on matching etc
                                crRNA_matching_strains_0 = list(dict.fromkeys(crRNA_matching_strains_0))
                                crRNA_matching_strains_3 = list(dict.fromkeys(crRNA_matching_strains_3))
                                crRNA_matching_strains_greater = list(dict.fromkeys(crRNA_matching_strains_greater))
                                crRNA_Mismatch3to5 = list(dict.fromkeys(crRNA_Mismatch3to5))
                                crRNA_Mismatch5to20 = list(dict.fromkeys(crRNA_Mismatch5to20))
                                unique_rank_ontarget = len(crRNA_matching_strains_0)
                                unique_rank_mismatch = len(crRNA_matching_strains_3)
                                unique_rank_mismatch_bad = len(crRNA_matching_strains_greater)
                                Master_List = zip_longest(crRNA_matching_strains_0, crRNA_matching_strains_3, crRNA_matching_strains_greater, crRNA_Mismatch3to5, crRNA_Mismatch5to20)
                                with open(primer_string, "a", newline = "") as output_file:
                                    for row in Master_List:
                                        writer = csv.writer(output_file)
                                        writer.writerow(row)
                                        #print(row)
                                os.chdir(path)
                                #Grabs exact matches     
            import os
            import pandas as pd
        
            working_dir = config.get('main', 'Working_Directory')
            temp_path = working_dir + "/Cas-ONInput/Temp"
            #temp_path = 'C:/Users/James_Mann1/Desktop/Testbench/Cas-ONInput/Temp' #FIX 
            #working_dir = 'C:/Users/James_Mann1/Desktop/Testbench' #Fix OCT8  
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
            
            try:
                with open("ContingencyMasterAnalysis.csv", "a") as output_file:
                    output_file.write("RPA PRIMER SET" + "," + "Amplicon Size" + "," + "crRNA Target" + "," + "On Target" + "," + "Mismatch 3" + "," + "Bad" + "," + "\n")
                os.chdir(current)
            except:
                catch = 1
            try:
                with open("ContingencyMasterAnalysis.csv", "a") as output_file:
                    output_file.write("RPA PRIMER SET" + "," + "Amplicon Size" + "," + "crRNA Target" + "," + "On Target" + "," + "Mismatch 3" + "," + "Bad" + "," + "\n")
                os.chdir(current)
            except:
                catch = 1
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
                        with open("ContingencyMasterAnalysis.csv", "a") as output_file:
                            output_file.write(current_RPA + "," + str(AMP) + "," + crRNA + "," + str(current_match) + "," +  str(current_match3) + "," + str(current_MismatchBad) + "," + "\n" )
                        os.chdir(current)
            cycles += 1
#%%                    
##############################################################################
#Verfies the primers
############################################################################
if __name__ == '__main__':
    import os
    #Testing
    #print(current_best_primer)
    #print(second_best_primer)
    #Delete before release
    Input_File = config.get('main', 'Input.fna_Location')
    working_dir = config.get('main', 'Working_Directory')
    os.chdir(working_dir)
    checksum = os.listdir(os.getcwd())
    #print(os.getcwd())
    checksum_value = 0 
    for value in checksum:
        if value.endswith("inal_Output.csv") == True:
            print("[Primed Sherlock] Found valid output csv file, checking primers.... ")
            print("[Primed Sherlock] This process requires that the input.fna was aligned, use Galaxy for quick and easy MAFFT")
            checksum_value = 1
    if checksum_value == 1:
        #Cuts first primer set for primed RPA analysis
        length = len(current_best_primer)
        gap = current_best_primer.find("_")
        BestF = current_best_primer[0:int(gap)]
        gap = gap + 1
        BestR = current_best_primer[int(gap):length]
        #Cuts second primer set for primed RPA analysis
        #Cuts first primer set for primed RPA analysis
        length = len(second_best_primer)
        gap = second_best_primer.find("_")
        TBestF = second_best_primer[0:int(gap)]
        gap = gap + 1
        TBestR = second_best_primer[int(gap):length]
        #Sets up dir, and dir switching for casoffinder
        loop_current_dir = os.getcwd()
        try:
            os.mkdir("PrimerVer")
        except:
            catch = 1
        os.chdir("PrimerVer")
        primer_dir = os.getcwd()
        primer_sets = []
        primer_sets.append(current_best_primer)
        primer_sets.append(second_best_primer)
        
        #Makes variable size casoffinder NNNNNNNNN inputs for different length primers....
        n_check_bf = len(BestF)
        n_check_br = len(BestR)
        n_check_sbf = len(TBestF)
        n_check_sbr = len(TBestR)
        n_best_f = ''
        n_best_r = ''
        n_sbest_f = ''
        n_sbest_r = ''     
        #Best Forward Loop
        loop_value = 0
        while loop_value < int(n_check_bf):
            loop_value += 1
            n_best_f = n_best_f + "N"
        #Best Reverse Loop
        loop_value = 0
        while loop_value < int(n_check_br):
            loop_value += 1
            n_best_r = n_best_r + "N"
        #Second Best Forward Loop
        loop_value = 0
        while loop_value < int(n_check_sbf):
            loop_value += 1
            n_sbest_f = n_sbest_f + "N"
        #Second Best Reverse Loop
        loop_value = 0
        while loop_value < int(n_check_sbr):
            loop_value += 1
            n_sbest_r = n_sbest_r + "N"
        for primers in primer_sets:
            if primers == current_best_primer:
                print("[Primed Sherlock] Processing Best Primer Set: F: " + BestF + " R: " + BestR)
                with open(str(BestF) + "_" + str(current_best_primer) + "input.txt", "w") as input_file:
                    input_file.write(str(working_dir) + "/on" + "\n")
                    input_file.write(str(n_best_f) + "\n")
                    input_file.write(str(BestF) + " 6" + "\n")
                with open(str(BestR) + "_" + str(current_best_primer) + "input.txt", "w") as input_file:
                    input_file.write(str(working_dir) + "/on" + "\n")
                    input_file.write(str(n_best_r) + "\n")
                    input_file.write(str(BestR) + " 6" + "\n")
            if primers == second_best_primer:
                print("[Primed Sherlock] Processing Second Best Primer Set: F: " + TBestF + " R: " + TBestR)
                with open(str(TBestF) + "_" + str(second_best_primer) + "input.txt", "w") as input_file:
                    input_file.write(str(working_dir) + "/on" + "\n")
                    input_file.write(str(n_best_f) + "\n")
                    input_file.write(str(TBestF) + " 6" + "\n")
                with open(str(TBestR) + "_" + str(second_best_primer) + "input.txt", "w") as input_file:
                    input_file.write(str(working_dir) + "/on" + "\n")
                    input_file.write(str(n_best_r) + "\n")
                    input_file.write(str(TBestR) + " 6" + "\n")
        print("[Primed Sherlock] Output files have been successfully created for Cas-Offinder.")
    #Runs CasOffinder on primer sets
    Input_File = config.get('main', 'Input.fna_Location')
    working_dir = config.get('main', 'Working_Directory')
    Output_Dir = Input_File + "/" + "PrimerVer"
    Input_path = Output_Dir
    # Input_Casoff = cas_prep(Input_File, Output_Dir, Input_path, )
    # with Pool(int(Multithread_Count)) as p:
    #     print(p.map(Cas_off, Input_Casoff))    #JAN 5


    Input_Casoff = cas_prep(Input_File, Output_Dir, Input_path)
    Iter_Value_List = []
    Original_files = len(Input_Casoff)
    for values in Input_Casoff:
        Iter_Value_List.append(Original_files)
        zipped = zip(Input_Casoff, Iter_Value_List)
    #print(set(zipped))
    with Pool(int(Multithread_Count)) as p:
        #p.map(Cas_off, Input_Casoff) #OCT10
        #print(p.map(Cas_off, Input_Casoff))
        print(p.starmap(Cas_off, zipped))


















###############################################################################
#Reads the CasOffinder Results for Primers
###############################################################################
#Grabs Primer Sites
    working_dir = config.get('main', 'Working_Directory')
    os.chdir(working_dir)
    checksum = os.listdir(os.getcwd())
    #print(os.getcwd())
    for value in checksum:
        if value.endswith("put_Sets.csv"):
            print("[Primed Sherlock] Found Output_Sets.csv file")
            try:
                data = pd.read_csv(value)
            except:
                print("[Primed Sherlock] No Valid Output_Sets.csv exists in directory: " + str(PrimedOutputLoc))
                system.os(exit)
            columns = list(data.columns.values.tolist()) 
            FPrimer = data[data.columns[3]]
            RPrimer = data[data.columns[12]]
            FPrimerSite = data[data.columns[1]]
            RPrimerSite = data[data.columns[10]]
            #Best Primers
            time_count = 0
            Primer_values= []
            Primer_Sites = []
            RPrimerlist = list(RPrimer)
            for primer in FPrimer:
                if primer == BestF:
                    if BestR == RPrimer[time_count]:
                        Primer_values.append(primer)
                        Primer_values.append(RPrimer[time_count])
                        Primer_Sites.append(FPrimerSite[time_count])
                        Primer_Sites.append(RPrimerSite[time_count])
                time_count +=1
            #Second Best Primers
            time_count = 0
            RPrimerlist = list(RPrimer)
            for primer in FPrimer:
                if primer == TBestF:
                    if TBestR == RPrimer[time_count]:
                        Primer_values.append(primer)
                        Primer_values.append(RPrimer[time_count])
                        Primer_Sites.append(FPrimerSite[time_count])
                        Primer_Sites.append(RPrimerSite[time_count])
                time_count +=1
            #print(Primer_values)
            #print(Primer_Sites)

    print("[Primed Sherlock] Successfully extracted Primer Site Data from Output_CSV")
    print("[Primed Sherlock] Comparing extracted Primer Site Data with Cas Results")
    ###########################################################################
    #Primer Fix
    ########################################################################### 

#Looks at the forward and reverse primers to determine best primer pairs. 
if __name__ == '__main__':
    import os
    import pandas as pd
    #Grabs Primer Sites
    working_dir = config.get('main', 'Working_Directory')
    os.chdir(working_dir)
    checksum = os.listdir(os.getcwd())
    #BestF
    #BestR
    #TBestF
    #TBestR
    PrimedOutputLoc = config.get('main','Primed_RPA_Output')
    os.chdir(PrimedOutputLoc)
    items = os.listdir(".")
    newlist = []
    for names in items:
        if names.endswith("Output_Sets.csv"):
            newlist.append(names)
    if len(newlist) >= 2:
        print("[Primed Sherlock] Ensure there is only one output file in the working directory, the first output file will be used otherwise")
    try:
        
        data = pd.read_csv(newlist[0])
        rpa_csv = open(newlist[0], "r")
    except:
        print("[Primed Sherlock] No Output file exists in directory: " + str(PrimedOutputLoc))
    
    theta = config.get('main','Amplicon_Max_Size')
    beta = config.get('main','Amplicon_Min_Size')
    maxback = config.get('main', 'Max_Background_Identitiy')

      
    columns = list(data.columns.values.tolist()) 
    Amp = data[data.columns[0]]
    F_Loc = data[data.columns[1]]
    FPrimer = data[data.columns[3]]
    R_Loc = data[data.columns[10]]
    RPrimer = data[data.columns[12]]
    MaxBackground = data[data.columns[6]]
    FPrimerSite = data[data.columns[1]]
    RPrimerSite = data[data.columns[10]]
    time_count = 0
    
    
    New_Amp = []
    New_FLoc = []
    New_RLoc = []
    New_FPrimer = []
    New_RPrimer = []
    New_MaxBackground = []
    NewFPrimerSite = []
    NewRPrimerSite = []

    for amplicons in Amp:
        New_Amp.append(Amp[time_count])
        New_FLoc.append(F_Loc[time_count])
        New_RLoc.append(R_Loc[time_count])
        New_FPrimer.append(FPrimer[time_count])
        New_RPrimer.append(RPrimer[time_count])
        New_MaxBackground.append(MaxBackground[time_count])
        NewFPrimerSite.append(FPrimerSite[time_count])
        NewRPrimerSite.append(RPrimerSite[time_count])
        time_count = time_count + 1
        #FPrimer
    FPrimer_Count = 0
    FPrimer_Loc = 0
    for FPrimers in New_FPrimer:
        if FPrimers == BestF:
            FPrimer_Loc = New_FLoc[FPrimer_Count]
           # print(str(FPrimer_Count) )
           # print(New_FPrimer[FPrimer_Count] + " aaa " + str(BestF))
        FPrimer_Count +=1
    FLoc_upper = FPrimer_Loc - 55
    
    FPrimer_Count = 0
    PossibleFPrimer = []
    for Current_FLoc in New_FLoc:
        if Current_FLoc >= FLoc_upper:
            if Current_FLoc < FPrimer_Loc:
                PossibleFPrimer.append(New_FPrimer[FPrimer_Count])
        FPrimer_Count +=1
    FPrimer_Candidates = list(dict.fromkeys(PossibleFPrimer))
    
    #Reverse Primer
    PossibleRPrimer = []
    FPrimer_Count = 0
    FPrimer_Loc = 0
    catch_low_sum = 0
    # for FPrimer in New_FPrimer:
    #     if FPrimer == BestF:
    #         if New_Amp[FPrimer_Count] <= 350:
    #             PossibleRPrimer.append(New_RPrimer[FPrimer_Count]) 
    #     FPrimer_Count +=1
    # print(PossibleRPrimer)
    for FPrimer in New_FPrimer:
        if FPrimer == BestF:
            if New_Amp[FPrimer_Count] <= 350:
                if New_RPrimer[FPrimer_Count] == BestR:
                    catch_low_sum = 1
                if catch_low_sum == 1:
                    PossibleRPrimer.append(New_RPrimer[FPrimer_Count]) 
        FPrimer_Count +=1
if __name__ == '__main__':
    import os
    import pandas as pd
    #Grabs Primer Sites
    working_dir = config.get('main', 'Working_Directory')
    os.chdir(working_dir)
    checksum = os.listdir(os.getcwd())
    
    try:
        shutil.rmtree("PrimerCheck")
    except:
        catch = 1
    try:
        os.mkdir("PrimerCheck")
    except:
        catch = 1
    current_dir = os.getcwd()
    os.chdir("PrimerCheck")
    for rprimer in PossibleRPrimer:
        rev_length = len(rprimer)
        ncheck = ""
        rev_length_loop = 0
        while rev_length_loop < rev_length:
            ncheck = ncheck + "n"
            rev_length_loop +=1
        with open("R_" + rprimer + "input.txt", "w") as input_file:
                    input_file.write(str(working_dir) + "/on" + "\n")
                    input_file.write(str(ncheck) + "\n")
                    input_file.write(str(rprimer) + " 6" + "\n")
    for fprimer in PossibleFPrimer:
        rev_length = len(fprimer)
        ncheck = ""
        rev_length_loop = 0
        while rev_length_loop < rev_length:
            ncheck = ncheck + "n"
            rev_length_loop +=1
        #print(len(ncheck))
        with open("F_" + fprimer + "input.txt", "w") as input_file:
                    input_file.write(str(working_dir) + "/on" + "\n")
                    input_file.write(str(ncheck) + "\n")
                    input_file.write(str(fprimer) + " 6" + "\n")
        
    rev_length = len(BestF)
    ncheck = ""
    rev_length_loop = 0
    while rev_length_loop < rev_length:
        ncheck = ncheck + "n"
        rev_length_loop +=1
    with open("BestF_" + BestF + "input.txt", "w") as input_file:
        input_file.write(str(working_dir) + "/on" + "\n")
        input_file.write(str(ncheck) + "\n")
        input_file.write(str(BestF) + " 6" + "\n")
        
    rev_length = len(BestR)
    ncheck = ""
    rev_length_loop = 0
    while rev_length_loop < rev_length:
        ncheck = ncheck + "n"
        rev_length_loop +=1
    with open("BestR_" + BestR + "input.txt", "w") as input_file:
        input_file.write(str(working_dir) + "/on" + "\n")
        input_file.write(str(ncheck) + "\n")
        input_file.write(str(BestR) + " 6" + "\n")
                
        
        
        
    os.chdir(current_dir)
if __name__ == '__main__':
    import os
    Input_File = config.get('main', 'Working_Directory')
    #Input_File ='c:/Users/James_Mann1/Desktop/TestBench' OCT8
    Output_Dir = Input_File + "/" + "PrimerCheck"
    Input_path = Output_Dir
    # Input_Casoff = cas_prep(Input_File, Output_Dir, Input_path, )
    # with Pool(int(Multithread_Count)) as p:
    #     print(p.map(Cas_off, Input_Casoff))    #JAN5
    
    Input_Casoff = cas_prep(Input_File, Output_Dir, Input_path)
    Iter_Value_List = []
    Original_files = len(Input_Casoff)
    for values in Input_Casoff:
        Iter_Value_List.append(Original_files)
        zipped = zip(Input_Casoff, Iter_Value_List)
    #print(set(zipped))
    print("[Primed Sherlock] Program is now running Cas-Offinder with the following files..." + "/n")
    with Pool(int(Multithread_Count)) as p:
        #p.map(Cas_off, Input_Casoff) #OCT10
        #print(p.map(Cas_off, Input_Casoff))
        print(p.starmap(Cas_off, zipped))

if __name__ == '__main__':
    import os
    import pandas as pd
    working_dir = config.get('main', 'Working_Directory')
    #print(os.getcwd())
    working_path = os.getcwd()
    files = os.listdir(working_path) #Searches all the files
    print("[Primed Sherlock] Coverting Cas Off-Target Output Data to CSV format")
    file_list = []
    for file in files:
        if file.endswith('txtout.txt') == True:
            file_out = file + ".csv"
            print(file_out)
            file_list.append(file_out)
            with open(file) as infile, open(file_out, 'w') as outfile:
                outfile.write("Cas_off-Results" + "," + "SeqID" + "," + "Location" + "," + "Mismatch" + "," + "MismatchSeq" + "," + "Score" + "," + "\n")
                for line in infile:
                      outfile.write(" ".join(line.split()).replace(' ', ','))
                      outfile.write("," + "\n") # trailing comma shouldn't matter  
    print("[Primed Sherlock] Finished converting, analysing results")
    print("  ")      
           
if __name__ == '__main__':
    files = os.listdir(os.getcwd())
    #print(os.getcwd())
    ranking = []
    score_total = []
    mismatch_seq = []
    good_seq = []
    Best = 0
    for file in files:
        filename = str(file)
        print(filename)
        if filename.endswith(".csv"):
            file_data = pd.read_csv(file)
            columns = list(file_data.columns.values.tolist()) 
            crRNA = file_data[file_data.columns[0]]
            SeqID = file_data[file_data.columns[1]]
            Loc = file_data[file_data.columns[2]]
            Mismatch = file_data[file_data.columns[3]]
            
            Score = file_data[file_data.columns[5]]
            count = 0
            score_count = 0
            mismatch_seq_id = []
            good_seq_id = []
            for scores in Score:
                if scores <= 2: #TODO MAKE USER EDITABLE
                    good_seq_id.append(SeqID[count])
                    score_count +=1
                if scores >= 3:
                    mismatch_seq_id.append(SeqID[count])
                count += 1
            ranking.append(file)
            score_total.append(score_count)     
            mismatch_seq.append(mismatch_seq_id)
            good_seq.append(good_seq_id)
            length = len(crRNA)
            length = str(length)
            length = int(length)
            if Best < length:
                Best = length
    count = 0 
    total = len(ranking)
    while count < total:
        if ranking[count].startswith("BestF") == True:
            bestFscore_new = mismatch_seq[count]
            bestFseq_1 = good_seq[count]
        count +=1 

if __name__ == '__main__':
    print("/n")
    print("/n")
    print("The current best F Primer " + BestF + " " + str(len(bestFseq_1)) + " has " + str(len(bestFscore_new)) + " strains with more then 3 mismatched primer pairs")
    print(bestFscore_new)
    
    scores = 0
    highscore = 0
    highest_seq = ""
    print("The following is a list of potential replacements followed by how many strains are within two basepairs")
    print(" ")
    
    count = 0
    for rank in ranking:
        if rank.startswith("F") == True:
            print(rank + "  " + str(score_total[scores]))
        scores +=1
    while count < total:
        if ranking[count].startswith("BestR") == True:
            BestRscore_new = mismatch_seq[count]
            BestRseq_1 = good_seq[count]
        count +=1 
    print("  ")
    print("\n")
    print("  ")
    print("The current best R Primer " + BestR + " " + str(len(BestRseq_1)) + " has " + str(len(BestRscore_new)) + " strains with more then 3 mismatched primer pairs")
    print(BestRscore_new)
    
    scores = 0
    highscore = 0
    highest_seq = ""
    print("The following is a list of potential replacements followed by how many strains are within two basepairs")    
    print(" ")
    for rank in ranking:
        if rank.startswith("R") == True:
            print(rank + "  " + str(score_total[scores]))
        scores += 1
    ##########################################################################
    #This code segment saves the original values, used at 2354 section
    ##########################################################################
    original_bestR_seq1 = BestRseq_1
    original_bestF_seq1 = bestFseq_1
    original_bestF_mismatch = bestFscore_new
    original_bestR_mismatch = BestRscore_new   
   

#Second best primer sets

if __name__ == '__main__':
    import os
    import pandas as pd
    #Grabs Primer Sites
    working_dir = config.get('main', 'Working_Directory')
    os.chdir(working_dir)
    checksum = os.listdir(os.getcwd())
    #BestF
    #BestR
    #TBestF
    #TBestR
    

    BestRM = BestR
    BestFM = BestF
    BestF = TBestF
    BestR = TBestR
    
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
    F_Loc = data[data.columns[1]]
    FPrimer = data[data.columns[3]]
    R_Loc = data[data.columns[10]]
    RPrimer = data[data.columns[12]]
    MaxBackground = data[data.columns[6]]
    FPrimerSite = data[data.columns[1]]
    RPrimerSite = data[data.columns[10]]
    time_count = 0
    
    
    New_Amp = []
    New_FLoc = []
    New_RLoc = []
    New_FPrimer = []
    New_RPrimer = []
    New_MaxBackground = []
    NewFPrimerSite = []
    NewRPrimerSite = []

    for amplicons in Amp:
        New_Amp.append(Amp[time_count])
        New_FLoc.append(F_Loc[time_count])
        New_RLoc.append(R_Loc[time_count])
        New_FPrimer.append(FPrimer[time_count])
        New_RPrimer.append(RPrimer[time_count])
        New_MaxBackground.append(MaxBackground[time_count])
        NewFPrimerSite.append(FPrimerSite[time_count])
        NewRPrimerSite.append(RPrimerSite[time_count])
        time_count = time_count + 1
        #FPrimer
    FPrimer_Count = 0
    FPrimer_Loc = 0
    for FPrimers in New_FPrimer:
        if FPrimers == BestF:
            FPrimer_Loc = New_FLoc[FPrimer_Count]
           # print(str(FPrimer_Count) )
           # print(New_FPrimer[FPrimer_Count] + " aaa " + str(BestF))
        FPrimer_Count +=1
    FLoc_upper = FPrimer_Loc - 55
    
    FPrimer_Count = 0
    PossibleFPrimer = []
    for Current_FLoc in New_FLoc:
        if Current_FLoc >= FLoc_upper:
            if Current_FLoc < FPrimer_Loc:
                PossibleFPrimer.append(New_FPrimer[FPrimer_Count])
        FPrimer_Count +=1
    FPrimer_Candidates = list(dict.fromkeys(PossibleFPrimer))
    
    #Reverse Primer
    PossibleRPrimer = []
    FPrimer_Count = 0
    FPrimer_Loc = 0
    catch_low_sum = 0
    for FPrimer in New_FPrimer:
        if FPrimer == BestF:
            if New_Amp[FPrimer_Count] <= 350:
                if New_RPrimer[FPrimer_Count] == BestR:
                    catch_low_sum = 1
                if catch_low_sum == 1:
                    PossibleRPrimer.append(New_RPrimer[FPrimer_Count]) 
        FPrimer_Count +=1
    
if __name__ == '__main__':
    import os
    import pandas as pd
    #Grabs Primer Sites
    working_dir = config.get('main', 'Working_Directory')
    os.chdir(working_dir)
    checksum = os.listdir(os.getcwd())
    
    try:
        shutil.rmtree("PrimerCheck")
        os.mkdir("PrimerCheck")
    except:
        catch = 1
    current_dir = os.getcwd()
    os.chdir("PrimerCheck")
    for rprimer in PossibleRPrimer:
        rev_length = len(rprimer)
        ncheck = ""
        rev_length_loop = 0
        while rev_length_loop < rev_length:
            ncheck = ncheck + "n"
            rev_length_loop +=1
        with open("R_" + rprimer + "input.txt", "w") as input_file:
                    input_file.write(str(working_dir) + "/on" + "\n")
                    input_file.write(str(ncheck) + "\n")
                    input_file.write(str(rprimer) + " 6" + "\n")
    for fprimer in PossibleFPrimer:
        rev_length = len(fprimer)
        ncheck = ""
        rev_length_loop = 0
        while rev_length_loop < rev_length:
            ncheck = ncheck + "n"
            rev_length_loop +=1
        with open("F_" + fprimer + "input.txt", "w") as input_file:
                    input_file.write(str(working_dir) + "/on" + "\n")
                    input_file.write(str(ncheck) + "\n")
                    input_file.write(str(fprimer) + " 6" + "\n")
        
    rev_length = len(BestF)
    ncheck = ""
    rev_length_loop = 0
    while rev_length_loop < rev_length:
        ncheck = ncheck + "n"
        rev_length_loop +=1
    with open("BestF_" + BestF + "input.txt", "w") as input_file:
        input_file.write(str(working_dir) + "/on" + "\n")
        input_file.write(str(ncheck) + "\n")
        input_file.write(str(BestF) + " 6" + "\n")
        
    rev_length = len(BestR)
    ncheck = ""
    rev_length_loop = 0
    while rev_length_loop < rev_length:
        ncheck = ncheck + "n"
        rev_length_loop +=1
    with open("BestR_" + BestR + "input.txt", "w") as input_file:
        input_file.write(str(working_dir) + "/on" + "\n")
        input_file.write(str(ncheck) + "\n")
        input_file.write(str(BestR) + " 6" + "\n")
                
        
        
        
    os.chdir(current_dir)
if __name__ == '__main__':
    import os
    Input_File = config.get('main', 'Working_Directory')
    #Input_File ='c:/Users/James_Mann1/Desktop/TestBench' OCT8
    Output_Dir = Input_File + "/" + "PrimerCheck"
    Input_path = Output_Dir
    # Input_Casoff = cas_prep(Input_File, Output_Dir, Input_path, )
    # with Pool(int(Multithread_Count)) as p:
    #     print(p.map(Cas_off, Input_Casoff))    #JAN5

    Input_Casoff = cas_prep(Input_File, Output_Dir, Input_path)
    Iter_Value_List = []
    Original_files = len(Input_Casoff)
    for values in Input_Casoff:
        Iter_Value_List.append(Original_files)
        zipped = zip(Input_Casoff, Iter_Value_List)
    #print(set(zipped))
    with Pool(int(Multithread_Count)) as p:
        #p.map(Cas_off, Input_Casoff) #OCT10
        #print(p.map(Cas_off, Input_Casoff))
        print(p.starmap(Cas_off, zipped))

if __name__ == '__main__':
    import os
    import pandas as pd
    working_dir = config.get('main', 'Working_Directory')
    working_path = os.getcwd()
    files = os.listdir(working_path) #Searches all the files
    print("[Primed Sherlock] Coverting Cas Off-Target Output Data to CSV format")
    file_list = []
    for file in files:
        if file.endswith('txtout.txt') == True:
            file_out = file + ".csv"
            file_list.append(file_out)
            with open(file) as infile, open(file_out, 'w') as outfile:
                outfile.write("Cas_off-Results" + "," + "SeqID" + "," + "Location" + "," + "Mismatch" + "," + "MismatchSeq" + "," + "Score" + "," + "\n")
                for line in infile:
                      outfile.write(" ".join(line.split()).replace(' ', ','))
                      outfile.write("," + "\n") # trailing comma shouldn't matter  
    print("[Primed Sherlock] Finished converting, analysing results")
    print("  ")      
if __name__ == '__main__':
    files = os.listdir(os.getcwd())
    ranking = []
    score_total = []
    mismatch_seq = []
    good_seq = []
    Best = 0
    for file in files:
        filename = str(file)
        if filename.endswith(".csv"):
            file_data = pd.read_csv(file)
            columns = list(file_data.columns.values.tolist()) 
            crRNA = file_data[file_data.columns[0]]
            SeqID = file_data[file_data.columns[1]]
            Loc = file_data[file_data.columns[2]]
            Mismatch = file_data[file_data.columns[3]]
            Score = file_data[file_data.columns[5]]
            count = 0
            score_count = 0
            mismatch_seq_id = []
            good_seq_id = []
            for scores in Score:
                if scores <= 2: #TODO MAKE USER EDITABLE
                    good_seq_id.append(SeqID[count])
                    score_count +=1
                if scores >= 3:
                    mismatch_seq_id.append(SeqID[count])
                count += 1
            ranking.append(file)
            score_total.append(score_count)     
            mismatch_seq.append(mismatch_seq_id)
            good_seq.append(good_seq_id)
            length = len(crRNA)
            length = str(length)
            length = int(length)
            if Best < length:
                Best = length
    count = 0 
    total = len(ranking)
    while count < total:
        if ranking[count].startswith("BestF") == True:
            bestFscore_new = mismatch_seq[count]
            bestFseq_1 = good_seq[count]
        count +=1 

    print("The current best F Primer " + BestF + " " + str(len(bestFseq_1)) + " has " + str(len(bestFscore_new)) + " strains with more then 3 mismatched primer pairs")
    print(bestFscore_new)
    
    scores = 0
    highscore = 0
    highest_seq = ""
    print("The following is a list of potential replacements followed by how many strains are within two basepairs")
    print(" ")
    count = 0
    for rank in ranking:
        if rank.startswith("F") == True:
            print(rank + "  " + str(score_total[scores]))
        scores +=1
    while count < total:
        if ranking[count].startswith("BestR") == True:
            BestRscore_new = mismatch_seq[count]
            BestRseq_1 = good_seq[count]
        count +=1 
    print("  ")
    print("\n")
    print("  ")
    print("The current best R Primer " + BestR + " " + str(len(BestRseq_1)) + " has " + str(len(BestRscore_new)) + " strains with more then 3 mismatched primer pairs")
    print(BestRscore_new)
    
    scores = 0
    highscore = 0
    highest_seq = ""
    print("The following is a list of potential replacements followed by how many strains are within two basepairs")    
    print(" ")
    for rank in ranking:
        if rank.startswith("R") == True:
            print(rank + "  " + str(score_total[scores]))
        scores += 1
        
        
 

    print(BestRscore_new)

    BestR = BestRM
    BestF = BestFM

    ###############################################################################
    #Finalizes Primer Findings in Output.csv file
    ############################################################################## 
if __name__ == '__main__':
    #It compares the files in the On target folder with the values to see if non compatiable strains exist".    
    os.chdir(working_dir)
    with open("Final_Output.csv", "a") as output_file:
        #PrimerSet1
        output_file.write("\n" + "," + "\n")
        output_file.write("Mismatch Analysis of Top Primer Candidates" + "," + "\n")
        output_file.write("Primer Set 1" + "," + BestF + " and " + BestR + "," + "\n")
        output_file.write("Forward Primer " + "," + BestF + "," + "\n")
        output_file.write("Total Seqs with 1 or less Mismatches" + "," + str(len(original_bestF_seq1)) + "," + "\n")
        seqs = ""
        count = 0
        for seq in original_bestF_seq1:
            if count == 0:
                output_file.write(",")
                output_file.write(",") 
            if count >= 0:
                output_file.write(seq + ",")            
            if count == 19:
                output_file.write("\n")
                count = -1
            count = count + 1 
        output_file.write("\n")
        output_file.write("Total Seqs which have more then 3 Mismatches" + "," + str(len(original_bestF_mismatch)) + "," + "\n")
        seqs = ""
        count = 0
        for seq in original_bestF_mismatch:
            if count == 0:
                output_file.write(",")
                output_file.write(",") 
            if count >= 0:
                output_file.write(seq + ",")            
            if count == 19:
                output_file.write("\n")
                count = -1
            count = count + 1  
        output_file.write("\n")
        output_file.write("\n")
        output_file.write("Reverse Primer " + "," + BestR + "," + "\n")
        output_file.write("Total Seqs with 1 or less  Mismatches" + "," + str(len(original_bestR_seq1)) + "," + "\n")
        seqs = ""
        count = 0
        for seq in original_bestR_seq1:
            if count == 0:
                output_file.write(",")
                output_file.write(",") 
            if count >= 0:
                output_file.write(seq + ",")            
            if count == 19:
                output_file.write("\n")
                count = -1
            count = count + 1 
        output_file.write("\n")
        output_file.write("Total Seqs which have more then 3 Mismatches" + "," + str(len(original_bestR_mismatch)) + "," + "\n")
        seqs = ""
        count = 0
        for seq in original_bestR_mismatch:
            if count == 0:
                output_file.write(",")
                output_file.write(",") 
            if count >= 0:
                output_file.write(seq + ",")            
            if count == 19:
                output_file.write("\n")
                count = -1
            count = count + 1 
       
        #PrimerSet2
        output_file.write("\n" + "," + "\n")
        output_file.write("\n" + "," + "\n")
        output_file.write("Mismatch Analysis of Top Primer Candidates" + "," + "\n")
        output_file.write("Primer Set 2" + "," + BestF + " and " + BestR + "," + "\n")
        output_file.write("Forward Primer " + "," + BestF + "," + "\n")
        output_file.write("Total Seqs with 1 or less Mismatches" + "," + str(len(bestFseq_1)) + "," + "\n")
        seqs = ""
        count = 0
        for seq in bestFseq_1:
            if count == 0:
                output_file.write(",")
                output_file.write(",") 
            if count >= 0:
                output_file.write(seq + ",")            
            if count == 19:
                output_file.write("\n")
                count = -1
            count = count + 1 
        output_file.write("\n")
        output_file.write("Total Seqs which have more then 3 Mismatches" + "," + str(len(bestFscore_new)) + "," + "\n")
        seqs = ""
        count = 0
        for seq in bestFscore_new:
            if count == 0:
                output_file.write(",")
                output_file.write(",") 
            if count >= 0:
                output_file.write(seq + ",")            
            if count == 19:
                output_file.write("\n")
                count = -1
            count = count + 1  
        output_file.write("\n")
        output_file.write("\n")
        output_file.write("Reverse Primer " + "," + BestR + "," + "\n")
        output_file.write("Total Seqs with 1 or less  Mismatches" + "," + str(len(BestRseq_1)) + "," + "\n")
        seqs = ""
        count = 0
        for seq in BestRseq_1:
            if count == 0:
                output_file.write(",")
                output_file.write(",") 
            if count >= 0:
                output_file.write(seq + ",")            
            if count == 19:
                output_file.write("\n")
                count = -1
            count = count + 1 
        output_file.write("\n")
        output_file.write("Total Seqs which have more then 3 Mismatches" + "," + str(len(BestRscore_new)) + "," + "\n")
        seqs = ""
        count = 0
        for seq in BestRscore_new:
            if count == 0:
                output_file.write(",")
                output_file.write(",") 
            if count >= 0:
                output_file.write(seq + ",")            
            if count == 19:
                output_file.write("\n")
                count = -1
            count = count + 1 
    
    ###############################################################################
    #Finalizes run by cleaning up the run files. Sometimes this is greater the 50gb
    #Later versions will lessen the amount of files created to lessen the write wear
    #on Solid State Drives.
    ############################################################################## 
   
    
    # working_dir = config.get('main', 'Working_Directory')
    # catch = 0
    # dir_list = [ "PrimerVer", "PrimerCheck", "Cas-OffInput", "crRNA_Pairs", "PrimerComboFasta", "TempAlignDir", "TempDIR", 'On', 'Off', 'crRNA_Pairs', 'Cas-On-Temp', 'Cas-On-Out', 'Cas-ONInput', 'Cas-On']
    # os.chdir(working_dir)
    # directories = [f for f in os.listdir('.') if os.path.isdir(f)]
    # for d in directories:
    #     if d in dir_list:
    #         try:
    #             sh.rmtree(str(d))
    #             print("[Primed Sherlock] Found existing " + d + " Directory, deleting directory")
    #         except:
    #             print("[Primed Sherlock] " + d + " directory could not be deleted these will be deleted upon program restart")
    # #Okay now we've established that there are no troublesome directories left or prompted the user to manually delete.
    # catch = 0 #This value is reused elsewhere

    end = timer()
    time = (end - start)
    time_minutes = int(time) / 60
    time_hours = int(time_minutes) / 60
    
    print("[Primed Sherlock] Thank you for using Primed Sherlock. Please cite this paper ")
    print("[Primed Sherlock] Completed generation of Primer & crRNA combos in " + str(time_minutes) + " minutes.")
    print("[Primed Sherlock] Completed generation of Primer & crRNA combos in " + str(time_hours) + " hours.")
    if time_hours >= 8:
        print("[Primed Sherlock] It is recommended to enable multithreading, or proceed with system upgrades to signficantly shorten time")
    
    
    
    
    
    
    
    
    
