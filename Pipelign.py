#!/usr/bin/env python3

'''
    Pipelign.py

    A python based program to align virus sequences.

    The program takes as input a single FASTA formatted file

    Returns a FASTA formatted alignment file
  
    Different flavours are present:
        - GNU Parallel for Linux and MacOS
        - Joblib for Linux, MacOS and Windows
'''

#*********************************************************************

import sys, os, shutil, subprocess, argparse, textwrap, stat
from Bio import SeqIO, Phylo
from Bio import AlignIO
from Bio.Align import AlignInfo
import tempfile
import time
import joblib
import math
from ete3 import Tree

#************************************
class Struct:
    '''
     - this class creates a structure of the parameter values needed for different functions.
     - parameters can be accessed as structName.itemName e.g. mArgs.gCode 
    '''
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self,k,v)
      
#*************************************      

class MyStruct(Struct):
    pass

#*****************************************

def lengthThreshold(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % (x,))
    return x

#***************************************
def simThreshold(x):
    x = float(x)
    if x < 0.4 or x > 1.0:
        raise argparse.ArgumentTypeError('%r not in range [0.4, 1.0]' % (x,))
    return x

#***************************************
def iterThreshold(x):
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError('%r must be a non-zero positive integer' % (x,))
    return x

#***************************************
def runType(x):
    x = x.upper()
    if x not in ['J','G']:
        raise argparse.ArgumentTypeError('%r must be either J or G' % (x,))
    return x

#***************************************
def mergeType(x):
    x = x.upper()
    if x not in ['F','P','C']:
        raise argparse.ArgumentTypeError('%r must be either F, P or C' % (x,))
    return x

#***************************************

def genCodeLimit(x):
    x = int(x)
    gCodes = list(1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25)
    if x not in gCodes:
        raise argparse.ArgumentTypeError('%r is not a valid genetic code' % (x,))
    else:
        return x

#****************************************

def cZip(cDir,tName,zName):
    '''
    creates a zip file of the temporary directory
    '''  
    os.chdir(cDir)
  
    #zName = 'pipelign.' + time.strftime('%Y-%m-%d-%H%M%S') 
    try:
        shutil.make_archive(zName,'zip',tName)
    except OSError as e: 
        sys.exit(e)
  
    print('\nArchive for all temporary files created in %s.zip\n' % zName)
    sys.exit()

#******************************************

def concatLogFiles(original,newLog):
    '''
    Concatenates two log files 
    '''
  
    lh = open(original,'a')
  
    with open(newLog,'rU') as f:
        for line in f:
            lh.write(line)
  
    lh.close()
  

#******************************************

def linsi(seqFile,alnFile,thread,mIterLong,log,cDir,tName,zName,p=None):
    '''
      - aligns sequences using MAFFT's L-INS-i method
    '''
  
    if p:
        lName = log + '.' + str(p)
    else:
        lName = log + '.temp'
  
    lh = open(lName,'a')
    ah = open(alnFile,'w')
  
    cl = ['mafft', '--localpair', '--thread', str(thread), '--maxiterate', str(mIterLong), '--preservecase', seqFile]

    try:
        subprocess.check_call(cl,stdout=ah,stderr=lh)
    except subprocess.CalledProcessError as e:
        lh.close()
        ah.close()
        print(e)
        cZip(cDir,tName,zName)
    
    lh.close()
    ah.close()
  
    if not p:
        concatLogFiles('pipelign.log',lName)  
        os.remove(lName)

#***********************************************************************

def fftnsi(seqFile,alnFile,thread,mIterLong,log,cDir,tName,zName,p=None):
    '''
      - aligns sequences using MAFFT's FFT-NS-i method
    '''

    if p:
        lName = log + '.' + str(p)
    else:
        lName = log + '.temp'
  
    lh = open(lName,'a')
    ah = open(alnFile,'w')
  
    cl = ['mafft', '--thread', str(thread), '--maxiterate', str(mIterLong), '--preservecase', seqFile]

    try:
        subprocess.check_call(cl, stdout=ah, stderr=lh)
    except subprocess.CalledProcessError as e:
        lh.close()
        ah.close()
        print(e)
        cZip(cDir,tName,zName)

    lh.close()
    ah.close()

    if not p:
        concatLogFiles('pipelign.log',lName)  
        os.remove(lName)

#***********************************************************************

def checkPresenceOfFile(fName):
    '''
      - Check whether file <fName> exists;
      - Terminate Pipelign run otherwise
    '''
   
    if not os.path.exists(fName) or os.stat(fName).st_size == 0:
        return False
  
    return True
  
#************************************************************************

def checkMPI():
    cmd = "mpirun"
    fh = open('mpi.log','a')
    mpi = True
    try:
        subprocess.call([cmd],stdout=fh,stderr=fh)
    except argparse.ArgumentTypeError as e:
        print('\nmpirun was not found. Running without MPI\n')
        mpi = False

    fh.close()
  
    return mpi

#************************************************************************  

def makeTempDir(tDir):
    '''
      - creates a temporary directory
      - creates inside user specified directory if path given
    '''
  
    if tDir is None: # no path provided for temp directory
        try:
            tempDir = tempfile.TemporaryDirectory() # create temporary directory to hold intermediary files
            tDirName = tempDir.name
        except OSError as e:
            sys.exit('\nError: system could not create temporary directory. Please try again')
  
    else: # user provided path using -d 
        if os.path.exists(tDir):
            tempDir = tDir + '/' + zName
            tDirName = tempDir
            try:
                os.mkdir(tempDir)
            except OSError as e:
                sys.exit('\nError: system could not create temporary directory. Please try again')
        else:
            sys.exit('\nError: Path for temporary directory does not exists. Please run again with correct path.')

    return tempDir, tDirName

#************************************************************************

def copyFile(sName,dName):
    '''
      - This function makes copy of a single file
      - file <sName> is copied to the destination <dName> 
    '''
  
    try:  
        shutil.copyfile(sName,dName)  
    except OSError as e:
        sys.exit(e)
    

#************************************************************************

def deAlign(iFile, dFile):
    '''
      - Removes gaps (if any) from the input sequence file
    ''' 
  
    #print("\nRemoving the gap characters '-'/'.' from the sequences")
  
    # Reading sequence file in fasta format
    seqs = list(SeqIO.parse(iFile,'fasta'))
  
    if len(seqs) > 1: # at least one sequence present | file is in FASTA format
  
        st = ''
  
        for seq in seqs:
          st += '>' + seq.id + '\n' + str(seq.seq).replace('-','').replace('.','') + '\n'
  
        fh = open(dFile,'w')
        fh.write(st)
        fh.close()
  
        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Gapless sequence file written in <%s>\n' % dFile
        print(msg)
  
    else: # no sequence present or wrong format
        msg = '\n\nError: Could not read the input sequence file.'
        msg += '\n       Make sure the file is in FASTA format'
        msg += '\n       and at least one sequnce present in file\n'
        sys.exit(msg) 
  
#*************************************************************************

def separateFullFragment(seqFile, lenThr, longName, fragName):
    '''
      Reads in the input sequence file.
          - finds length of the longest sequence
          - calculates minimum length required for full length sequences
          - writes full length sequences and fragments into two separate files 
    '''
  
    maxLen = 0
  
    # get maximum length of the input sequences
    handle = open(seqFile, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        if len(record.seq) > maxLen:
            maxLen = len(record.seq)
    handle.close()
  
    # calculate minimum length for long sequences
    minLengthLong = int(lenThr * maxLen)
  
    longS = []
    fragS = []
  
    # create separate lists for long sequences and fragments
    handle = open(seqFile, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        if len(record.seq) < minLengthLong:
            fragS.append(record)
        else:
            longS.append(record)
    handle.close()    
  
    # write long file
    SeqIO.write(longS,longName,'fasta') 
  
    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' <%s> ' % longName
  
    fEmpty = True
  
    if len(fragS) > 0:
        SeqIO.write(fragS, fragName,'fasta')
        msg += ' and <%s> ' % fragName
        fEmpty = False
    
    msg += 'created\n'
    print(msg)
  
    return fEmpty, len(fragS) 

#************************************************************************
def separateLongAmbigs(alphabet,ambigPer,cDir,tFileName,zName):
    '''
        Separates the long sequences that have ambiguous characters more than a limit
    '''
  
    # list the ambiguous characters
  
    if alphabet == 'aa':
        chars = ['X','!','*','?']
    else:
        chars = ['N','X','!','*','?']
  
    # count ambiguous characters for each sequences
    # if count is more than the threshold, separate them in another file
  
    gs = [] # contains good quality sequences
    bs = [] # contains bad quality sequences
  
    handle = open('long.fas','r')
  
    # counts number of ambigs for each sequences
    for seq in SeqIO.parse(handle,'fasta'):
        seqStr = str(seq.seq).upper()
        acount = 0
        for i in range(len(chars)):
            acount += seqStr.count(chars[i])
    
        propChars = int(len(seq.seq) * float(ambigPer))
    
        if acount > propChars:
            bs.append(seq)
        else:
            gs.append(seq)
    
    if len(gs) > 0:
        SeqIO.write(gs,'long.good.fas','fasta')
    else:
        print('\nAll long sequences have ambiguous characters more than allowed proportion %s' % ambigPer)
        print('\nPlease check whether correct alphabet (e.g. dna) was selected')
        print('\n\nPipelign is exiting\n')
        cZip(cDir,tFileName,zName)
    
    if len(bs) > 0:
        SeqIO.write(bs,'long.bad.fas','fasta')
    
    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Good quality long sequences written in <long.good.fas>\n' 
    print(msg)
  
    handle.close()

#************************************************************************

def runCDHIT(longName,alphabet,per,thread,cDir,tName,zName):
    '''
      CD-HIT is used to group similar sequences together in clusters for alignment
    '''
  
    # count number of sequences in long sequence file 
    seqCount = 0
    handle = open(longName,'rU')
    for record in SeqIO.parse(handle,'fasta'):  
        seqCount += 1
  
    '''
    # only one long sequence, one cluster
    if seqCount == 1:
        try:
            shutil.copy(longName,'grp')
        except OSError as e:
            print(e)
            cZip(cDir,tName,zName)
        print('Only one Long sequence present')
        return
    '''
    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
  
    # open log file
    lh = open('cdhit.log','a') 
    
    # choose word size for cd-hit
    if per >= 0.7:
        ws = 5
    elif per >= 0.6:
        ws = 4
    elif per >= 0.5:
        ws = 3
    else:
        ws = 2
        
    # argument string for CD-HIT
    if alphabet == 'dna' or alphabet == 'rna':
        cl = ['cd-hit-est','-c', str(per), '-n', str(ws), '-i', longName, '-o', 'grp', '-d', '0', '-T', str(thread)]
        msg += ' CD-HIT-EST started\n'

    elif alphabet == 'aa':
        cl = ['cd-hit','-c', str(per), '-n', str(ws), '-i', longName, '-o', 'grp', '-d', '0', '-T', str(thread)] 
        msg += ' CD-HIT started\n'

    print(msg)
    
    try:
        subprocess.check_call(cl, stdout=lh, stderr=lh)
    except subprocess.CalledProcessError as e:
        print(e)
        cZip(cDir,tName,zName)
  
    lh.close() # close log file   
  
    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' CD-HIT finished. Files created: <grp> and <grp.clstr>\n'
    print(msg)

    concatLogFiles('pipelign.log','cdhit.log')
    #os.remove('cdhit.log')
  
#*************************************************************************
def makeClusters(longName,cName):
    '''
      Separate files are created for each clusters
    '''

    if not os.path.exists(cName) or os.stat(cName).st_size == 0:
        msg = '\nError: the file <%s>  could not be found in the given path' % cName
        msg += '\n\nPipelign is exiting\n'
        sys.exit(msg)  
  
    # read in the long sequence file
    seqs = list(SeqIO.parse(longName,'fasta'))

    #clsSize = list() # will hold a list of cluster sizes
  
    # only one long sequence
    if len(seqs) == 1:
        copyFile(longName,'long.0.fas')
        #clsSize.append(1)
        fh = open('long.ClusterList.txt','w')
        fh.write('%s\t0' %seqs[0].id)
        fh.close()
        return 1 #, clsSize
  
    #cName = 'grp.clstr'
  
    lines = [line.strip('\n') for line in open(cName,'rU')] # read cluster file
  
    start = 0 # flag for the beginning of first cluster list
   
    cSeq = [] # hold sequences of a cluster
  
    cls = 0 # count clusters
  
    st = '' # clusterList string
  
    ids = [] # IDs for the full sequences
  
    for seq in seqs:
        ids.append(seq.id)
  
    # read cluster file and make separate cluster files 
    if 'Cluster' not in lines[0]:
        msg = '\n\nError: <%s> does not contain any cluster' % cName
        msg += '\nPlease try running the program again\n'
        sys.exit(msg)
  
    for i in range(1, len(lines)):
        if 'Cluster' in lines[i]: # start of a new cluster list
            gName = 'long.' + str(cls) + '.fas'
            cls = cls + 1
            if len(cSeq) > 0:
                SeqIO.write(cSeq,gName,'fasta')
                #clsSize.append(len(cSeq))
            cSeq = []
      
        else: # continue with the existing cluster
            seqID = lines[i].split()[2].replace('>','').replace('...','')
            st += seqID + '\t' + str(cls) + '\n' # updating clusterList file content
            sInd = ids.index(seqID)
            cSeq.append(seqs[sInd])
    
    gName = 'long.' + str(cls) + '.fas'  
    cls = cls + 1
    if len(cSeq) > 0:
        SeqIO.write(cSeq,gName,'fasta')
    
    fh = open('long.ClusterList.txt','w')
    fh.write(st)
    fh.close()
  
    if cls > 0:
        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' %d cluster file(s) created. File names: <long.x.fas>\n' % cls
        print(msg)
    else:
        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Could not create cluster files.'
        msg += '\nPipelign is exiting\n'
        sys.exit(msg)
    
    return cls #, clsSize   

#***********************************************************************
def addClusterNumberToReps(repName,lstFile,outFile):
    '''
        - Reads in the cluster representative FASTA file <grp> and the <long.ClusterList.txt> file
        - Adds cluster number and size to the sequence header e.g. >seq1_size_cluster   
        - Temporary file <clsReps.fas> is written
    '''
  
    # read in the <long.ClusterList.txt> file
    #cList = [line.strip() for line in open(lstFile,'r')]
  
    cID = [] # sequence ids 
    cNum = [] # cluster numbers
  
    with open(lstFile,'rU') as f:
        for line in f:
            words = line.split()
            cID.append(words[0])
            cNum.append(words[1])
    
    # read in the cluster reps file
    seqs = list(SeqIO.parse(repName,'fasta'))
  
    for seq in seqs:
        if seq.id in cID:
            ind = cID.index(seq.id)
      
            seq.id = seq.id + '_' + str(cNum.count(cNum[ind])) + '_' + cNum[ind] 
            seq.name = seq.id
            seq.description = seq.id
        else:
            msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
            msg += ' Error: %s was not found in <long.ClusterList.txt>.' % seq.id
            msg += '\nPipelign is exiting\n'
            sys.exit(msg)
  
    # write reps file
    SeqIO.write(seqs,outFile,'fasta')    

#***********************************************************************
def makeIQTree(alnFile,thread,cDir,tName,zName,alpha):
    '''
      - Constructs phylogenetic tree using IQ-TREE
    '''
  
    lh = open('tree.log','a')
  
    #print('\nCreating IQ-TREE from %s' % alnFile)
  
    if alpha == 'dna' or alpha == 'rna':
        cl = ['iqtree-omp', '-s', alnFile, '-m', 'GTR+R4', '-nt', str(thread)]
    elif alpha == 'aa':
        cl = ['iqtree-omp', '-s', alnFile, '-m', 'WAG', '-nt', str(thread)]
  
    try:
        subprocess.check_call(cl,stdout=lh,stderr=lh)
    except subprocess.CalledProcessError as e:
        print(e)
        lh.close()
        cZip(cDir,tName,zName)
  
    lh.close()
  
    concatLogFiles('pipelign.log','tree.log')
    os.remove('tree.log')
  
#***********************************************************************
def makeMidPointRootTree(treeFile):
    '''
     - Displays a dendogram of the tree generated from cluster representatives
    '''
  
    # Read the tree in newick format
    tree = Phylo.read(treeFile,'newick')
  
    # Root tree at midpoint
    tree.root_at_midpoint()
  
    # Write midpoint root tree
    Phylo.write(tree,'clsReps.aln.midpoint.treefile','newick')

#***********************************************************************

def drawAsciiTree(treeFile):
    '''
      - draws ASCII tree on the console using ete3
    '''
  
    # Read the tree in newick format
    #tree = Phylo.read(treeFile,'newick')

    #Phylo.draw_ascii(tree)
    tree = Tree(treeFile)
    #print(tree.get_ascii(show_internal=True))
    print(tree)
  
#***********************************************************************
def alnLongSequenceClustersGNUParallel(nClusters,thread,mIterL,cDir,tName,zName):
    '''
        aligns long sequence cluster files using GNU parallel
    '''

    log = 'align.log'
    fh = open(log,'a')
  
    numClusters = nClusters - 1 # for bash run
  
    #msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    #msg += ' Started aligning long sequence clusters\n'
    #print(msg)

    pStr = '#!/bin/bash\n\n'
    pStr += 'my_func(){\n\n'
    pStr += '\tinF="long.$1.fas"\n\n'
    pStr += '\tif [ -a $inF ]\n'
    pStr += '\tthen\n'
    pStr += '\t\tlen=$(grep -c ">" $inF)\n'
    pStr += '\tfi\n\n'
  
    pStr += '\tif [ "$len" -eq "1" ]\n'
    pStr += '\tthen\n'
    pStr += '\t\tcat $inF > long.$1.aln\n'
    pStr += '\t\techo "cluster $1 has only one sequence. <long.$1.aln> written"\n\n'
  
    pStr += '\telif [ "$len" -le "100" ]\n'
    pStr += '\tthen\n'
    pStr += '\t\tmafft --localpair --thread $2 --maxiterate $3 --preservecase $inF > long.$1.aln\n'
    pStr += '\t\techo "Alignment file <long.$1.aln> created with $len sequences"\n'    
    pStr += '\telse\n'
    pStr += '\t\tmafft --thread $2 --maxiterate $3 --preservecase $inF > long.$1.aln\n'
    pStr += '\t\techo "Alignment file <long.$1.aln> created with $len sequences"\n'    
    pStr += '\tfi\n\n'
    pStr += '}\n\n'
  
    pStr += 'export -f my_func\n\n'
    pStr += 'parallel -j %s my_func ::: $(seq 0 1 %d) ::: %s ::: %s' % (thread,numClusters,thread,mIterL)

    il = open('longAlign.sh','w')
    il.write(pStr)
    il.close()
  
    # assign executable permission to the file
    st = os.stat('longAlign.sh')
    os.chmod('longAlign.sh', st.st_mode | stat.S_IEXEC)
  
    # subprocess call to run mafft in parallel
    script = ['./longAlign.sh']
  
    try:
      subprocess.check_call(script,stderr=fh)
    except OSError as e:
      print(e)
      cZip(cDir,tName,zName)
  
#***********************************************************************
def alnLongSequenceClustersJoblibParallel(nClusters,thread,mIterL,cDir,tName,zName):
  '''
    - align long sequence cluster files in parallel using joblib 
  '''
  
  log = 'align.log'
  #lh = open('pipelign.log','a')
  
  msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
  msg += ' Started aligning long sequence clusters\n'
  print(msg)
  
  # Fork the worker processes to perform computation concurrently
  # Create parameter map
  aln = lambda i : ('long.' + str(i) + '.fas', 'long.' + str(i) + '.aln')

  to_run_tuples = list(map(aln, range(nClusters)))
  to_run_linsi = list(filter(lambda x : 1 < len(list(SeqIO.parse(str(x[0]),'fasta'))) < 101, to_run_tuples))
  to_run_fftnsi = list(filter(lambda x : len(list(SeqIO.parse(str(x[0]),'fasta'))) > 100, to_run_tuples))
  to_copy = list(filter(lambda x : len(list(SeqIO.parse(str(x[0]),'fasta'))) <= 1, to_run_tuples))

  if len(to_run_linsi):
    num_parallel_jobs = int(nClusters/thread) if nClusters > thread else thread
    num_threads_per_job = int(thread/nClusters) if thread > nClusters else 1
    joblib.Parallel(n_jobs=num_parallel_jobs)(joblib.delayed(linsi)(x[0],x[1],num_threads_per_job,mIterL,log,cDir,tName,zName,x[0].split('.')[1]) for x in to_run_linsi)

  if len(to_run_fftnsi):
    num_parallel_jobs = int(nClusters/thread) if nClusters > thread else thread
    num_threads_per_job = int(thread/nClusters) if thread > nClusters else 1
    joblib.Parallel(n_jobs=num_parallel_jobs)(joblib.delayed(fftnsi)(x[0],x[1],num_threads_per_job,mIterL,log,cDir,tName,zName,x[0].split('.')[1])  for x in to_run_fftnsi)

  if len(to_copy):
    num_parallel_jobs = math.ceil(len(to_copy)/thread) if nClusters < thread else thread
    joblib.Parallel(n_jobs=num_parallel_jobs)(joblib.delayed(shutil.copyfile)(x[0],x[1]) for x in to_copy)

  msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
  msg += ' Finished aligning long sequence clusters. Created files: <long.x.aln>\n'
  print(msg)
    
 
#***********************************************************************


def makeHMMdbParallel(nClusters,thread,alphabet,log,cDir,tName,zName):
    '''
        - builds HMM from long sequence clusters
        - builds the database from HMMs
    '''
  
    # Fork the worker processes to perform computation concurrently
    # Create parameter map
    hmm = lambda i : ('long.' + str(i) + '.aln', 'long.' + str(i) + '.hmm')
  
    to_run_tuples = list(map(hmm, range(nClusters)))  
  
    if len(to_run_tuples):
        msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Started building HMMs from long cluster alignments\n'
        print(msg)

        num_parallel_jobs = int(nClusters/thread) if nClusters > thread else thread
        num_threads_per_job = int(thread/nClusters) if thread > nClusters else 1
        joblib.Parallel(n_jobs=num_parallel_jobs)(joblib.delayed(buildHMM)(x[0],x[1],num_threads_per_job,alphabet,log,cDir,tName,zName,x[0].split('.')[1])  for x in to_run_tuples)
    
        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Building HMM database\n'
        print(msg)
    
        buildHMMdb(nClusters)

        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' HMM database created\n'
        print(msg)
    
#***********************************************************************
def buildHMM(aName,hName,thread,alphabet,log,cDir,tName,zName,p=None):
    '''
        - Builds HMM for an alignment
    '''
  
    if p:
        lName = log + '.' + str(p)
    else:
        lName = log + '.temp'
  
    lh = open(lName,'a')
  
    gName = aName + '.temp'
    # create a temporary file to hold the alignment to build HMM  
    try:
        shutil.copy(aName,gName)
    except OSError as e:
        print(e)
        cZip(cDir,tName,zName)

    #*****testing for failure of hmmbuild
    # hmmbuild fails when single sequence is present
    # also when alignment does not have any gaps
    # solution: make a copy of the single sequnce 
    # add a '-' at the end of each sequences
  
    # read in the alignment file
    aseqs = list(SeqIO.parse(gName,'fasta'))

    # for a single sequence cluster, add itself
    if len(aseqs) == 1: 
        aseqs.append(aseqs[0])
        aseqs[1].id = aseqs[0].id + '_1'
        aseqs[1].name = aseqs[1].id
        aseqs[1].description = aseqs[1].id

    # for cluster with multiple sequences
    if len(aseqs) > 1:
        sumGap = 0 # for counting gaps

        for seq in aseqs:
            sumGap = sumGap + seq.seq.count('-')
    
        # if the alignment has no gaps, add one
        if sumGap == 0:
            for seq in aseqs:
                seq.seq = seq.seq + '-'

            # Update 'temp.aln' file 
            try:
                os.remove(gName)
            except OSError as e:
                print(e)
                cZip(cDir,tName, zName)
        
            SeqIO.write(aseqs,gName,'fasta')  
      
  
    # create the command for hmmbuild
    if alphabet == 'dna':
        cl = ['hmmbuild', '--dna', '--cpu', str(thread), hName, gName]
  
    elif alphabet == 'rna':
        cl = ['hmmbuild', '--rna', '--cpu', str(thread), hName, gName]
  
    elif alphabet == 'aa':
        cl = ['hmmbuild', '--amino', '--cpu', str(thread), hName, gName]
    
    # run hmmbuild command
    try:
        subprocess.check_call(cl,stdout=lh,stderr=lh)
    except subprocess.CalledProcessError as e:
        #sys.exit(e)
        print(e)
        cZip(cDir,tName,zName)
    #print('\t<%s> created' % hName)
    
    # remove the temporary alignment file
    try:
        os.remove(gName)
    except OSError as e:
        print(e)
        cZip(cDir,tName,zName)
  
    lh.close()  
  
    if not p:
        concatLogFiles('pipelign.log',lName)
        os.remove(lName)
    
#***********************************************************************

def buildHMMdb(nClusters):
    '''
        - builds HMM database from cluster HMMs
    '''
  
    dbName = 'pipelign.hmm'
  
    # write all HMMs in one file
    fh = open(dbName,'w')
    for i in range(nClusters):
        hName = 'long.' + str(i) + '.hmm'
        shutil.copyfileobj(open(hName,'r'),fh)
    fh.close()

    # create HMM database
    cl = ['hmmpress', dbName]
  
    lh = open('buildhmm.log','a')
  
    # run command
    try:
        subprocess.check_call(cl, stdout=lh, stderr=lh)
    except subprocess.CalledProcessError as e:
        #sys.exit(e)
        print(e)
        cZip(cDir,tName)
  
    lh.close()
  
    concatLogFiles('pipelign.log','buildhmm.log')
    os.remove('buildhmm.log')
  
#***********************************************************************
def longSeqClusters(tempDirPath,cDir, zName, tDirName, tFileName):
    '''
    Copies the temporary working directory into current working directory
    The directory will contain long sequence cluster alignments, HMMs, log files
    '''

    msg = '\nPipelign has created long sequence cluster alignments and HMMs.\n'
    
    # save the pipelign run files into current directory
    
    # no path provided for temporary directory
    if tempDirPath is None:
        try:
            wName = cDir + '/' + zName
            shutil.copytree(tDirName,wName)
        except OSError as e:
            print(e)
            cZip(cDir,tFileName,zName)   
    
    msg += '\nAlignment files and HMMs can be found in <%s>\n' % zName
    msg += '\tLong sequence alignments have names <long.xx.aln>\n'
    msg += '\tHMM file written in <long.xx.hmm>\n'
    msg += '\tHMM database written in <pipelign.hmm>\n'
    print(msg)
    
    sys.exit('\nThank you for using Pipelign.\n')   

#***********************************************************************

def addNameToInternalNodes(tree):
    # annotate names to internal nodes
    edge = 0
    for node in tree.traverse():
        if not node.is_leaf():
            node.name = "NODE_%d" %edge
            edge += 1
    return tree

#************************************************************************

def changeLeafNamesToClusterNames(tree,tp):
    # annotate name of the leaf nodes to get cluster names
    for leaf in tree:
        lName = tp + '.' + leaf.name.split('_')[-1]
        leaf.name = lName
        #print(leaf.name)
    return tree

#************************************************************************

def mergeCherries_gnu_parallel(nCherries,cStr,iNodes,fPre,thread,mIterL,log,cDir,tName,zName):
    '''
     merge pair of alignments in parallel using GNU's parallel tool
    '''
    msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Merging pairs of clusters\n'
    print(msg)
    
    pStr = '#!/bin/bash\n\n'
    pStr += 'my_func(){\n\n'
    
    pStr += '\tind1=$(($1*2-1))\n'
    pStr += '\tind2=$(($1*2))\n\n'
    
    pStr += '\taln1=$(echo $4 | awk -v var="$ind1" \'{print $var}\')\n'
    pStr += '\taln2=$(echo $4 | awk -v var="$ind2" \'{print $var}\')\n'
    pStr += '\toName=$(echo $5 | awk -v var="$1" \'{print $var}\')\n\n'
    pStr += '\techo "Merging $aln1.aln and $aln2.aln into $oName.aln"\n\n'
    
    pStr += '\tl1=$(grep -c ">" $aln1.aln)\n'
    pStr += '\tl2=$(grep -c ">" $aln2.aln)\n\n'
    
    pStr += '\tif [ $l1 -eq "1" ] && [ $l2 -eq "1" ] \n'
    pStr += '\tthen\n'
    pStr += '\t\tcat $aln1.aln $aln2.aln > temp.$6.$1.fas\n'
    pStr += '\t\tmafft --localpair --preservecase --thread $2 --maxiterate $3 temp.$6.$1.fas > $oName.aln\n'
    pStr += '\telif [ $l1 -eq "1" ] \n'
    pStr += '\tthen\n'
    pStr += '\t\tmafft --preservecase --thread $2 --maxiterate $3 --addfragments $aln1.aln $aln2.aln > $oName.aln\n'
    pStr += '\telif [ $l2 -eq "1" ] \n'
    pStr += '\tthen\n'
    pStr += '\t\tmafft --preservecase --thread $2 --maxiterate $3 --addfragments $aln2.aln $aln1.aln > $oName.aln\n'
    pStr += '\telse\n'
    pStr += '\t\tcat $aln1.aln $aln2.aln > input.$6.$1.fas\n'
    pStr += '\t\tstr1=""\n'
    pStr += '\t\tcount=1\n'
    pStr += '\t\tgap=" "\n\n'
    
    pStr += '\t\tfor i in $(seq 1 1 $l1)\n'
    pStr += '\t\tdo\n'
    pStr += '\t\t\tstr1=$str1$count$gap\n'
    pStr += '\t\t\tcount=$(($count+1))\n'
    pStr += '\t\tdone\n\n'
    
    pStr += '\t\techo $str1 > msaTable.$6.$1\n'
    #pStr += '\t\techo "" >> msaTable.$6.$1\n\n'
    
    pStr += '\t\tstr2=""\n'
    pStr += '\t\tfor i in $(seq 1 1 $l2)\n'
    pStr += '\t\tdo\n'
    pStr += '\t\t\tstr2=$str2$count$gap\n'
    pStr += '\t\t\tcount=$(($count+1))\n'
    pStr += '\t\tdone\n\n'
    
    pStr += '\t\techo $str2 >> msaTable.$6.$1\n'
    
    pStr += '\t\tmafft --preservecase --thread $2 --maxiterate $3 --merge msaTable.$6.$1 input.$6.$1.fas > $oName.aln\n\n'
    
    pStr += '\tfi\n\n'
    
    pStr += '}\n\n'
    
    pStr += 'export -f my_func\n\n'
    
    pStr += 'parallel -j %s my_func ::: $(seq 1 1 %d)' % (thread,nCherries)
    pStr += ' ::: %s ::: %s ::: "%s"' % (thread,mIterL,cStr)
    pStr += ' ::: "%s" ::: "%s"\n' % (' '.join(iNodes), fPre)
    
    il = open('mergeParallel.sh','w')
    il.write(pStr)
    il.close()
    
    # assign executable permission to the bash script
    st = os.stat('mergeParallel.sh')
    os.chmod('mergeParallel.sh',st.st_mode | stat.S_IEXEC)
    
    # subprocess call to run merge in parallel
    script = ['./mergeParallel.sh']
    
    try:
        fh = open('merge.log','a')
        subprocess.check_call(script,stderr=fh)
    except OSError as e:
        print(e)
        fh.close()
        #cZip(cDir,tName,zName)
    
    fh.close()

#************************************************************************

def mergePair(aln1,aln2,alnFile,thread,mIterLong,log,cDir,tName,zName,p=None):
    '''
      - merge pair of alignments using MAFFT's --merge method
    '''
  
    if p:
        lName = log + '.' + str(p)
    else:
        lName = log + '.temp'
  
    lh = open(lName,'a')
    ah = open(alnFile,'w')
    
    seq1 = list(SeqIO.parse(aln1,'fasta'))
    seq2 = list(SeqIO.parse(aln2,'fasta'))
    
    lSeq1 = len(seq1)
    lSeq2 = len(seq2)
    
    nSeq = seq1 + seq2
    SeqIO.write(nSeq,'input','fasta')
    
    mStr = ''
    count = 1
    for i in range(lSeq1):
        mStr += str(count) + ' '
        count += 1
    mStr += '\n'
    for i in range(lSeq2):
        mStr += str(count) + ' '
        count += 1
        
    # write subMSAtable
    fh = open('subMSAtable','w')
    fh.write(mStr)
    fh.close()  
    
  
    cl = ['mafft', '--thread', str(thread), '--maxiterate', str(mIterLong), '--preservecase', '--merge', 'subMSAtable', 'input']

    try:
        subprocess.check_call(cl,stdout=ah,stderr=lh)
    except subprocess.CalledProcessError as e:
        lh.close()
        ah.close()
        print(e)
        cZip(cDir,tName,zName)
    
    lh.close()
    ah.close()
  
    if not p:
        concatLogFiles('pipelign.log',lName)  
        os.remove(lName)        
        
#******************************************

def merge_gnu_parallel(repTree,tp,thread,mIterL,log,cDir,tName,zName):
    # read the tree
    tree = Tree(repTree,format=1)

    # annotate internal node names
    tree = addNameToInternalNodes(tree)

    # display the tree         
    #print(tree.get_ascii(show_internal=True))

    # change leaf node names as cluster names
    if tp == 'long':
        tree = changeLeafNamesToClusterNames(tree,'long')
    elif tp == 'all':
        tree = changeLeafNamesToClusterNames(tree,'cls')
    print(tree.get_ascii(show_internal=True))
    
    # get cherries and update trees
    while(len(tree) > 1):
        cherries = [] # empty list for cherries
        iNodes = [] # empty list for internal nodes
        for leaf in tree:
            if len(leaf.up) == 2: # two leaf nodes
                #print(leaf.name)
                #print(leaf.up)
                cherries.append(leaf.name)
                if leaf.up.name not in iNodes:
                    iNodes.append(leaf.up.name)
    
        # call merge() 
        if len(cherries) % 2 != 0:
            sys.exit('incomplete cherries found: exiting')
        # even number of leafs in cherries list
        # get number of cherries
        nCherries = int(len(cherries) / 2)
        #print(nCherries)
        cStr = ' '.join(cherries)
        #print(cStr)
        mergeCherries_gnu_parallel(nCherries,cStr,iNodes,'merge',thread,mIterL,log,cDir,tName,zName)
        #break
    
        # get the pointers to the cherries leafs and remove them
        for cherry in cherries:
            l = tree.search_nodes(name=cherry)[0]
            #print(l)
            #print(tree.get_ascii(show_internal=True))
            l.detach()
            #print(tree.get_ascii(show_internal=True))
    
    fName = tree.name + '.aln'
    return fName
  

#************************************************************************
def addLongSeqs(sName,aName,oName,thread,mIterL,log,cDir,tName,zName):
    '''
        - add long sequences to the alignment using MAFFT's --addfull
    '''

    lName = log + '.temp'
  
    lh = open(lName,'a')
    ah = open(oName,'w')
  
    cl = ['mafft', '--preservecase', '--thread', str(thread), '--maxiterate', str(mIterL), '--addfull', sName, aName]
        
    try:
        subprocess.check_call(cl,stdout=ah,stderr=lh)
    except subprocess.CalledProcessError as e:
        lh.close()
        ah.close()
        print(e)
        cZip(cDir,tName,zName)
  
    lh.close()
    ah.close()

    concatLogFiles(log,lName)
    os.remove(lName)
  
#***********************************************************************

def mergeConsensus(numClusters,oName,tp,mArgs,log,cDir,tFileName,zName):
    '''
        This function merges clusters together by following:
            - first consensus sequence is created from each cluster alignments
            - these consensus sequences are aligned together using MAFFT
            - gaps are placed in the original alignments based on the gaps in the consensus alignment
        
        tp can be:
            - 'long' for long sequence alignment
            - 'all' for cluster alignments with fragments
        
        Final alignment is produced by sequentially adding each of the clusters in the output
    '''

    ## get the file names prefix
    if tp == 'all':
        tp = 'cls'
    
    ## create a string to hold the consensus sequences in a FASTA format
    st = ''
    
    ## create consensus from each of the cluster alignments and create a sequence file
    for i in range(numClusters):
        aName = tp + '.' + str(i) + '.aln'
        aln = AlignIO.read(open(aName),'fasta')
        con = AlignInfo.SummaryInfo(aln).dumb_consensus(threshold=0.5, ambiguous='N')
        
        st += '>' + aName + '\n'
        st += str(con) + '\n'
    
    ## write the consensus sequence file
    fh = open('consensus.fas','w')
    fh.write(st)
    fh.close()

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Aligning consensus sequences from each clusters\n'
    print(msg)
    
    ## align the consensus sequences using FFT-NS-i
    fftnsi('consensus.fas','consensus.aln',mArgs.thread,mArgs.mIterL,log,cDir,tFileName,zName)
    
    ## read in the consensus alignment
    ## get positions for gaps in each sequence
    ## add gaps to their respective positions in the cluster alignments
    fh = open(oName,'w') # create the output alignment file
    
    # read consensus alignment
    aln = AlignIO.read(open('consensus.aln'),'fasta')
        
    for i in range(numClusters):
        indices = [j for j, x in enumerate(aln[i].seq) if x == "-"]
        cName = aln[i].id
        
        # read in the cluster alignment file
        seqs = SeqIO.parse(open(cName),'fasta')
        
        
        for seq in seqs:
            # create a string object to write in the file
            st = '>' + seq.id + '\n'
            
            start = 0
            
            for j in range(len(indices)):
                end = indices[j] - j
                st += seq.seq[start:end] + '-'
                start = end
            
            st += seq.seq[start:] + '\n'
        
            fh.write(str(st))
    
    fh.close()
    
    
#***********************************************************************

def longSeqAlignmentConsensus(numClusters,mArgs,cDir,tFileName,zName):
    '''
    Creates alignment from only long sequence clusters
      - adds sequences from <long.bad.fas> if -b flag is used 
      
    '''
    
    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Merging long sequence clusters using consensus\n'
    print(msg)
    
    if numClusters == 1: # only one cluster present
        copyFile('long.0.aln','long.noBadSeq.aln')
    else:
        mergeConsensus(numClusters,'long.noBadSeq.aln','long',mArgs,'consensus.log',cDir,tFileName,zName)    

    # check if <long.bad.fas> needs to be added to the final alignment
    if mArgs.keepBadSeqs and checkPresenceOfFile('long.bad.fas'):
        addLongSeqs('long.bad.fas','long.noBadSeq.aln','final.all.aln',mArgs.thread,mArgs.mIterL,'addBadSeq.log',cDir,tFileName,zName) 
    else:
        copyFile('long.noBadSeq.aln','final.all.aln')

    '''
    elif numClusters == 2: # two clusters
        #mergePair('long.0.aln','long.1.aln','long.noBadSeq.aln',mArgs.thread,mArgs.mIterM,'merge.log',cDir,tFileName,zName)
        
    
    elif numClusters >= 3: # more than 2 clusters
        #resName = merge_gnu_parallel('clsReps.aln.midpoint.treefile','long',mArgs.thread,mArgs.mIterL,'merge.log',cDir,tFileName,zName)
        #copyFile(resName, 'long.noBadSeq.aln')
    '''
    '''
    # check if <long.bad.fas> needs to be added to the final alignment
    if mArgs.keepBadSeqs and checkPresenceOfFile('long.bad.fas'):
        addLongSeqs('long.bad.fas','long.noBadSeq.aln','final.aln',mArgs.thread,mArgs.mIterL,'addBadSeq.log',cDir,tFileName,zName) 
    else:
        copyFile('long.noBadSeq.aln','final.aln')
    
    if checkPresenceOfFile('final.aln'):
        copyFile('final.aln',mArgs.outFile)
    
    
    if checkPresenceOfFile(mArgs.outFile):
        msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Final alignment written in <%s>\n' % mArgs.outFile
        print(msg)
    else:
        msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Pipelign could not write the file <%s>\n' % mArgs.outFile
        msg += '\tSomething went wrong. Please check the zipped output directory.\n' 
        print(msg)
        
    '''
    #sys.exit('\nThank you for using Pipelign.\n')

#***********************************************************************
def longSeqAlignmentParallel(numClusters,mArgs,cDir,tFileName,zName):
    '''
    Creates alignment from only long sequence clusters
      - adds sequences from <long.bad.fas> if -b flag is used 
      
    '''
    
    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Merging long sequence clusters using parallel strategy\n'
    print(msg)
    
    if numClusters == 1: # only one cluster present
        copyFile('long.0.aln','long.noBadSeq.aln')
    
    elif numClusters == 2: # two clusters
        mergePair('long.0.aln','long.1.aln','long.noBadSeq.aln',mArgs.thread,mArgs.mIterM,'merge.log',cDir,tFileName,zName)
    
    elif numClusters >= 3: # more than 2 clusters
        resName = merge_gnu_parallel('clsReps.aln.midpoint.treefile','long',mArgs.thread,mArgs.mIterL,'merge.log',cDir,tFileName,zName)
        copyFile(resName, 'long.noBadSeq.aln')
    
    # check if <long.bad.fas> needs to be added to the final alignment
    if mArgs.keepBadSeqs and checkPresenceOfFile('long.bad.fas'):
        addLongSeqs('long.bad.fas','long.noBadSeq.aln','final.all.aln',mArgs.thread,mArgs.mIterL,'addBadSeq.log',cDir,tFileName,zName) 
    else:
        copyFile('long.noBadSeq.aln','final.all.aln')
    
                                          
#***********************************************************************

def runBlast(log):
    '''
        - creates BLAST database from cluster representatives
        - search fragments against BLAST database
        - assign clusters to fragments and writes into file <frags.ClusterList.txt> 
    '''
    # first create blast database of cluster reps
    bName = log + '.temp'
    bl = open(bName,'a')

    cStr = 'makeblastdb -in clsReps.fas -input_type fasta -title pipelign -dbtype '
    if mArgs.alphabet in ['dna','rna']:
        cStr += 'nucl '
    elif mArgs.alphabet == 'aa':
        cStr += 'prot '
      
    cStr += '-out pipelign.blastdb'   

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Creating BLAST database from Cluster representatives\n'
    print(msg)

    try:
        cl = cStr.split()
        subprocess.check_call(cl,stdout=bl,stderr=bl)
    except subprocess.CalledProcessError as e:
        print(e)
        cZip(cDir,tName,zName)
      
    # search fragments on the database to assign clusters
    cStr = ''
    if mArgs.alphabet in ['dna','rna']:
        cStr += 'blastn '
    elif mArgs.alphabet == 'aa':
        cStr += 'blastp '
      
    cStr += '-query frag.fas -db pipelign.blastdb -max_target_seqs 1 -outfmt 6 -evalue 0.05'
    #cStr += '-evalue 20 | sort -u -k1,2'
  
    bo = open('frag.blast.txt','w')    
  
    try:
        cl = cStr.split()
        subprocess.check_call(cl,stdout=bo,stderr=bl)
    except subprocess.CalledProcessError as e:
        print(e)
        cZip(cDir,tName,zName)
  
    # get uniq BLAST hits for each fragments
  
  
  
    # parse blast out file to create fragment clusterList file  
    lines = [line.strip() for line in open('frag.blast.txt','rU')]
  
    st = '' # holds string for frag.ClusterList.txt
  
    mids = [] # holds names of the matched fragments 
  
    for line in lines:
        words = line.split()
        if words[0] not in mids:
            st += words[0] + '\t' + words[1].split('_')[-1] + '\n'
            mids.append(words[0])
    
    fh = open('frag.ClusterList.txt','w')
    fh.write(st)
    fh.close()  
  
    bl.close()
  
    concatLogFiles('pipelign.log',bName)
    os.remove(bName)
  
    return mids  
  
#************************************************************************
def searchHMMdb(log,thread,alpha,res,cDir,tName,zName):
    '''
        HMM database is searched with the fragments to assign them a cliuster for alignment
    '''  
  
    lName = log + '.temp'
    lh = open(lName,'a')
  
    # generate the command for HMM search
    if alpha == 'dna' or alpha == 'rna':
        cStr = 'nhmmscan --cpu %d --tblout %s --noali  -E 0.05 pipelign.hmm frag.noBlast.fas' % (thread,res)
      
    elif alpha == 'aa':
        cStr = 'hmmscan --cpu %d --tblout %s --noali -E 0.05 pipelign.hmm frag.noBlast.fas' % (thread,res)
    
    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Searching HMM database to assign clusters to remaining fragments\n'
    print(msg)
  
    try:
        cl = cStr.split()
        subprocess.check_call(cl,stdout=lh,stderr=lh)
    except subprocess.CalledProcessError as e:
        print(e)
        cZip(cDir,tName,zName)
  
    lh.close()
  
    concatLogFiles('pipelign.log',lName)
    os.remove(lName)

#**********************************************************************
def createFragmentClusters(numClusters):
    '''
        Fragments are written into separate files based on their cluster assignments
            - Reads in <frag.fas>
            - Reads in <frag.ClusterList.txt>
            - Creates fragments clusters <frag.x.fas>
            - Creates <frag.noClusters.fas> if no match found
    '''

    # read in the fragment sequnece file
    fseqs = list(SeqIO.parse('frag.fas','fasta'))
  
    fids = [] # holds the ids of all fragments
  
    for f in fseqs:
        fids.append(f.id)
  
    # create an empty list to hold flags for fragments that are assigned clusters
    fragFlags = []
    for i in range(len(fids)):
        fragFlags.append(-1)  
  
  
    # read in the cluster assignment file
    clsList = [line.strip() for line in open('frag.ClusterList.txt','r')]
  
    ids = [] # holds IDs of the assigned fragments
    cnum= [] # holds the cluster numbers

    for line in clsList:
        words = line.split()
        ids.append(words[0])
        cnum.append(words[1])
    #print(cnum)
    for i in range(numClusters):
        tfrags = [] # temporary list of all the fragments of that cluster
    
        # get the indices with cluster 'i'
        indices = [j for j, x in enumerate(cnum) if int(x) == i]
        #print(indices)
    
        for ind in indices:
            fname = ids[ind] # get the sequence name/ID
            if fname in fids: # if id matches to the fragment list
                sind = fids.index(fname) # get the index in the fragment file
                tfrags.append(fseqs[sind]) # add fragment to the cluster list
                fragFlags[sind] = 1 # flag 
    
        if len(tfrags) > 0:
            fName = 'frag.' + str(i) + '.fas'
            SeqIO.write(tfrags,fName,'fasta')    # the the fragments to their cluster files
    
    # identify the orphan fragments and write them in a files
  
    orphans = []
  
    for i in range(len(fseqs)):
        if fragFlags[i] == -1:
            orphans.append(fseqs[i])
  
    if len(orphans) > 0:
        SeqIO.write(orphans,'frag.noClusters.fas','fasta')
        print("Unassigned fragments are written in frag.noClusters.fas\n")
        return 0
  
    else:
        return 1 # indicate that all fragments were assigned clusters  
  
#************************************************************************
def getHMMclusters():
    '''
        Parse results from <hmm.out> to assign clusters to fragments based on HMM search
            - Reads in <frag.hmm.fas> file to read fragments searched against HMM_DB
            - Reads in <hmm.out> to parse output
            - Adds assigned fragments into the file <frag.ClusterList.txt>
            - Unassigned fragments are written into <frag.noClusters.fas>
        
    '''
    # open <frag.ClusterList.txt> for appending new entries
    fh = open('frag.ClusterList.txt','a')
  
    # read in the hmm.out file  
    hmmOut = [line.strip() for line in open('hmm.out','rU')]
  
    qids = [] # list all matched query IDs from output file
    hids = [] # list the hmm that they matched with
  
    # make lists of queries and target HMMs
    for line in hmmOut:
        if not line.startswith('#'):
            words = line.split()
            query = words[2]
            target = words[0].split('.')[1]
            if query not in qids:
                qids.append(query)
                hids.append(target)
  
    # search each of the fragment ids in the qids list to find a match
    # get the first match as HMM search sorted the output based on best match
    handle = open('frag.noBlast.fas','rU')
    for seq in SeqIO.parse(handle,'fasta'):
        if seq.id in qids: # if the id matched
            ind = qids.index(seq.id)
            fh.write("%s\t%s\n" %(seq.id,hids[ind]))
  
    fh.close()   
    
    return len(qids) 
  
#**********************************************************************


def assignClusterstoFrags(numClusters,numFrag,mArgs,tDirName,cDir,tFileName,zName):
    '''
    Assigns clusters to fragments based on BLAST and HMM search
    Copies the content of the temporary working directory into current working directory
    
    '''
    
    # if no fragments present in the dataset
    if not checkPresenceOfFile('frag.fas'):
        print('\nFragment file <frag.fas> does not exist. No fragments found for assigning clusters.\n')

        msg = '\nPipelign has created only long sequence cluster alignments and HMMs.\n'
    
    else: # <frag.fas> file exists, fragments need to be assigned clusters
        # running BLAST on the fragments
        mIds = runBlast('pipelign.log') # running BLAST on the fragments      
        fragWithClusters = len(mIds) 
        #print(fragWithClusters)

        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' %d fragments were assigned clusters using BLAST search\n' % fragWithClusters
        print(msg)
  
        # not all the fragments were assigned clusters
        if fragWithClusters < int(numFrag):
            
            noMatchFrags = [] # hold fragments not matching in BLAST
      
            handle = open('frag.fas','rU')
            for seq in SeqIO.parse(handle,'fasta'):
                if seq.id not in mIds:
                    noMatchFrags.append(seq)
      
            if len(noMatchFrags) > 0:
                SeqIO.write(noMatchFrags,'frag.noBlast.fas','fasta')
          
            # search remaining fragments against HMM_DB
            searchHMMdb('pipelign.log', mArgs.thread, mArgs.alphabet, 'hmm.out', cDir,tFileName, zName)

            # get the clusters from HMM result
            numHMM = getHMMclusters()
            
            if numHMM > 0: # if fragments were assigned clusters based on HMM search
                msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
                msg += ' %d fragments were assigned clusters using HMM search\n' % numHMM
                print(msg)
                

        # create clusters of fragments, un-assigned ones are written in <frag.noCluster.fas>
        allAssigned = createFragmentClusters(numClusters) 
        
#***********************************************************************

def addFragmentToClustersGNUParallel(nClusters,log,mArgs,cDir,tName,zName):
    '''
        add fragments to clusters using GNU parallel
    '''
  
    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Adding fragments to long sequence clusters\n'
    print(msg)
  
    # get number of clusters for bash running
    numClusters = nClusters - 1
  
    # create bash script for running MAFFT --addfragment
    pStr = '#!/bin/bash\n\n'
    pStr += 'my_func(){\n\n'

    pStr += '\tinF="frag.$1.fas"\n\n'
    pStr += '\tif [ -a $inF ]\n' # file present
    pStr += '\tthen\n'
    pStr += '\t\tmafft --preservecase --thread %s --maxiterate %s --addfragments $inF long.$1.aln > cls.$1.aln\n' % (mArgs.thread,mArgs.mIterL)
    pStr += '\telse\n'
    pStr += '\t\tcat long.$1.aln > cls.$1.aln\n'
    pStr += '\tfi\n\n'
    pStr += '}\n\n'
  
    pStr += 'export -f my_func\n\n'
    pStr += 'parallel -j 2 my_func ::: $(seq 0 1 %d) ::: %s ::: %s' % (numClusters,mArgs.thread,mArgs.mIterL)

    il = open('fragAlign.sh','w')
    il.write(pStr)
    il.close()
  
    # assign executable permission to the file
    st = os.stat('fragAlign.sh')
    os.chmod('fragAlign.sh', st.st_mode | stat.S_IEXEC)
  
    # subprocess call to run mafft in parallel
    script = ['./fragAlign.sh']
  
    fh = open(log,'a')
  
    try:
        subprocess.check_call(script,stderr=fh,stdout=fh)
    except OSError as e:
        print(e)
        cZip(cDir,tName,zName)
   

#************************************************************************
def addFragments(fName,aName,oName,thread,mIterL,log,cDir,tName,zName,p=None):
    '''
        - add fragments to the alignment using MAFFT's --addfragments
    '''

    if p:
        lName = log + '.' + str(p)
    else:
        lName = log + '.temp'
  
    lh = open(lName,'a')
    ah = open(oName,'w')
  
    cl = ['mafft', '--preservecase', '--thread', str(thread), '--maxiterate', str(mIterL), '--addfragments', fName, aName]
        
    try:
        subprocess.check_call(cl,stdout=ah,stderr=lh)
    except subprocess.CalledProcessError as e:
        lh.close()
        ah.close()
        print(e)
        cZip(cDir,tName,zName)
  
    lh.close()
    ah.close()

    if not p:
        concatLogFiles('pipelign.log',lName)
        os.remove(lName)
    
#***********************************************************************

def addFragmentToClustersJoblibParallel(nClusters,log,thread,mIterL,cDir,tName,zName):
    '''
        This will add fragments to their cluster alignments using joblib
            - MAFFT's -addfragment will be used
            - will create cluster alignment files: <cls.x.aln>
    '''
  
    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Adding fragments to long sequence clusters\n'
    print(msg)

    # Fork the worker processes to perform computation concurrently
    # Create parameter map
    aln = lambda i : ('long.' + str(i) + '.aln', 'frag.' + str(i) + '.fas', 'cls.' + str(i) + '.aln')
  
    to_run_tuples = list(map(aln, range(nClusters)))
    to_run_addfragment = list(filter(lambda x : os.path.exists(x[1]) and os.stat(x[1]).st_size > 0, to_run_tuples)) # fragments present for that cluster
    to_copy = list(filter(lambda x : not(os.path.exists(x[1]) and os.stat(x[1]).st_size > 0), to_run_tuples)) # no fragment for that cluster
  
    # run in parallel
    if len(to_run_addfragment):
        num_parallel_jobs = int(nClusters/thread) if nClusters > thread else thread
        num_threads_per_job = int(thread/nClusters) if thread > nClusters else 1
        joblib.Parallel(n_jobs=num_parallel_jobs)(joblib.delayed(addFragments)(x[1],x[0],x[2],num_threads_per_job,mIterL,log,cDir,tName,zName,x[0].split('.')[1]) for x in to_run_addfragment)

    if len(to_copy):
        num_parallel_jobs = math.ceil(len(to_copy)/thread) if nClusters < thread else thread
        joblib.Parallel(n_jobs=num_parallel_jobs)(joblib.delayed(shutil.copyfile)(x[0],x[2]) for x in to_copy)

    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Finished adding fragments to long sequence clusters. Created files: <cls.x.aln>\n'
    print(msg)
      
#************************************************************************
def excludeClustersFromalignment(oAln,xAln):
    '''
    Excludes clusters from alignment 
    '''
    # draw midpoint root tree
    print('\nMidpoint rooted tree for cluster representatives is drawn below. Tips are labelled as <H_N_C>.')
    print('\tH: header of the cluster representative sequence')
    print('\tN: number of long sequences in that cluster')
    print('\tC: unique identifier of the cluster')
    drawAsciiTree('clsReps.aln.midpoint.treefile')
    print('\n')
       
    
    # get user input for clusters
    inStr = input('Enter clusters (space separated) to exclude from final alignment:')
    inStr = inStr.strip().split()
    #print(inStr)
    
    if len(inStr) < 1:
        print('No cluster selected for exclusion')
        copyFile(oAln,xAln)
        return
    
    exList = list() # contains list of sequences for exclusion
    
    # read in the long cluster list file to list sequences for exclusion
    longSeqs = [line.strip() for line in open('long.ClusterList.txt','r')]
    for lSeq in longSeqs:
        words = lSeq.split()
        if words[1] in inStr:
            exList.append(words[0])
    
    
    if checkPresenceOfFile('frag.ClusterList.txt'):
        # read in the fragment cluster list file to list fragments for exclusion
        fragSeqs = [line.strip() for line in open('frag.ClusterList.txt','r')]
        for fSeq in fragSeqs:
            words = fSeq.split()
            if words[1] in inStr:
                exList.append(words[0])
    
    msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Removing sequences from final alignment\n'
    print(msg)
    
    if len(exList) > 0:
        # check if original alignment file exists, read alignment
        if checkPresenceOfFile(oAln):
            aln = list(SeqIO.parse(oAln,'fasta'))
        else:
            print('\nOriginal alignment file not found. Pipelign is exiting') 
            cZip(cDir,tFileName,zName)
        
        # create a list of allowed sequences
        aSeq = list()
        for seq in aln:
            if seq.id not in exList:
                aSeq.append(seq)
                   
        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Writing new alignment file after removing sequences\n'
        print(msg)
        SeqIO.write(aSeq,xAln,'fasta')         

#************************************************************************

def getArguments():
    '''
        Parses all the command line arguments from the user
    '''
  
    parser = argparse.ArgumentParser(description="Pipelign: creates multiple sequence alignment from FASTA formatted sequence file", formatter_class=argparse.RawTextHelpFormatter)
    #formatter_class=argparse.RawDescriptionHelpFormatter)  
    parser.add_argument('-i', '--inFile', required=True, help="Input sequence file in FASTA format")
    parser.add_argument('-o', '--outFile', required=True, help="FASTA formatted output alignment file")
    parser.add_argument('-t', '--lenThr', type=lengthThreshold, help="Length threshold for full sequences (default: 0.9)", default=0.9)
    parser.add_argument('-c', '--code', type=int, help="Genetic code for translation (default: 1)",default=1, choices=[1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25])
    parser.add_argument('-a', '--alphabet', required=True, help='Input sequences can be dna/rna/aa (default: dna)', choices=['dna','aa','rna'], default='dna')
    parser.add_argument('-f', '--keepOrphans', help='Add fragments without clusters', action="store_true")
    parser.add_argument('-b', '--keepBadSeqs', help='Add long sequences with too many ambiguous residues', action="store_true")    
    parser.add_argument('-z', '--mZip', help='Create zipped temporary files', action="store_true")
    parser.add_argument('-p', '--simPer', type=simThreshold, help="Percent sequence similarity for clustering (default: 0.8)", default=0.8)
    parser.add_argument('-r', '--run', type=runType, help="Run either (J)oblib/(G)NU parallel version (default: G)", default='G')    
    parser.add_argument('-e', '--merge', type=mergeType, help="Merge using (P)arallel/(C)onsensus strategy  (default: P)", default='P')
    parser.add_argument('-q', '--thread', nargs='?', const=1, type=int, help="Number of CPU/threads to use (default: 1)", default=1)
    parser.add_argument('-s', '--mIterateLong', type=iterThreshold, help="Number of iterations to refine long alignments (default: 1)", default=1)
    parser.add_argument('-m', '--mIterateMerge', type=iterThreshold, help="Number of iterations to refine merged alignment (default: 1)", default=1)
    parser.add_argument('-d', '--tempDirPath', required=False, help="Path for temporary directory",default=None)
    #parser.add_argument('-l', '--longSeqsOnly', help='Only align long sequences', action="store_true")
    parser.add_argument('-w', '--ambigPer', type=lengthThreshold, help="Proportion of ambiguous characters allowed in the long sequences (default: 0.01)", default=0.01)  
    parser.add_argument('-n', '--stage', type=int,default=5, choices=[1,2,3,4,5,6],
        help=textwrap.dedent('''\
        1  Make cluster alignments and HMM of long sequences
        2  Align long sequences only
        3  Assign fragments to clusters
        4  Make cluster alignments with fragments
        5  Align all sequences
  	    '''))

    parser.add_argument('-x', '--excludeClusters', help='Exclude clusters from final alignment', action="store_true")  
  
    args = parser.parse_args()

    return args  
  
#****************************************************************************

if __name__=="__main__":
  
    # get all the command line arguments
    args = getArguments()
    
    mArgs = MyStruct(
        inFile = args.inFile,
        outFile = os.getcwd() + '/' + args.outFile,
        lenThr = args.lenThr,
        gCode = args.code,
        alphabet = args.alphabet,
        keepOrphans = args.keepOrphans,
        keepBadSeqs = args.keepBadSeqs,
        makeZip = args.mZip,
        simPer = args.simPer,
        run = args.run.upper(),
        thread = args.thread,
        mIterL = args.mIterateLong,
        mIterM = args.mIterateMerge,
        fragEmpty = 1,
        longName = 'long.fas',
        fragName = 'frag.fas',
        tempDirPath = args.tempDirPath,
        ambigPer = args.ambigPer,
        merge = args.merge.upper(),
        #longSeqsOnly = args.longSeqsOnly,
        stage = args.stage,
        excludeClusters = args.excludeClusters)

    # check whether input file exists
    if not checkPresenceOfFile(mArgs.inFile):
        msg = '\n==============================='
        msg += '\nInput file <%s> does not exist in the given path' % mArgs.inFile
        msg += '\nPlease try again with the correct file name and path'
        msg += '\n==============================='
        msg += '\n\nPipelign is exiting.\n'
        sys.exit(msg)

    # Pipeline process starts here
    cDir = os.getcwd() # save path for the current working directory
    tFileName = 'input.fas'  # input file in temporary directory
    deFileName = 'input.dealign.fas' # dealigned input file

  
    # get current time and generate name for the zipped temporary directory
    timeNow = time.strftime('%Y-%m-%d-%H%M%S')
    zName = 'pipelign.' + timeNow 

    # create temporary directory 
    tempDir, tDirName = makeTempDir(mArgs.tempDirPath)

    # copy input file inside the temporary directory
    dName = tDirName + '/' + tFileName
    copyFile(mArgs.inFile,dName)
  
    msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Pipelign started with sequence file <%s>\n' % mArgs.inFile 
    print(msg)
  
  
    # change current working directory to the temp
    os.chdir(tDirName)
  
    # removes any possible gaps from the sequence file
    deAlign(tFileName, deFileName) 

    # make separate sequence file for long sequences and fragments
    mArgs.fragEmpty, numFrag = separateFullFragment(deFileName, mArgs.lenThr, mArgs.longName, mArgs.fragName)

    # separate the long sequences containing ambiguous characters more than a given proportion
    separateLongAmbigs(mArgs.alphabet,mArgs.ambigPer,cDir,tFileName,zName)
    mArgs.longName = 'long.good.fas'

    # use CD-HIT to group sequences based on similarity
    runCDHIT(mArgs.longName, mArgs.alphabet, mArgs.simPer, mArgs.thread,cDir,tFileName,zName)  

    # create separate files for each cluster
    numClusters = makeClusters(mArgs.longName,'grp.clstr')

    # add cluster numbers and size to header
    addClusterNumberToReps('grp','long.ClusterList.txt','clsReps.fas')

    if numClusters < 3: # less than 3 clusters, tree will not be built
        msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Only %d cluster(s) created. Phylogenetic tree cannot be built\n' % numClusters
        print(msg)

    else: # at least 3 clusters
        # ********************************************************
        # align cluster reps; reconstruct phylogeny; draw midpoint rooted tree
        msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Aligning cluster representative sequences\n'
        print(msg)

        if numClusters <= 100:
            linsi('clsReps.fas','clsReps.aln',mArgs.thread,mArgs.mIterL,'pipelign.log',cDir,tFileName,zName)
        else:
            fftnsi('clsReps.fas','clsReps.aln',mArgs.thread,mArgs.mIterL,'pipelign.log',cDir,tFileName,zName)
    
        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Alignment of cluster representatives written in <clsReps.aln>\n' 
        print(msg)
    
        # generate phylogenetic tree using IQ-TREE  
        makeIQTree('clsReps.aln',mArgs.thread,cDir,tFileName,zName,mArgs.alphabet)

        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Phylogenetic tree of the cluster representatives written in <clsReps.aln.treefile>\n' 
        print(msg)
    
        # make midpoint root tree for cluster representatives
        makeMidPointRootTree('clsReps.aln.treefile')
    
        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Midpoint rooted tree written in <clsReps.aln.midpoint.treefile>\n' 
        print(msg)
    
        # draw midpoint root tree
        print('\nMidpoint rooted tree for cluster representatives is drawn below. Tips are labelled as <H_N_C>.')
        print('\tH: header of the cluster representative sequence')
        print('\tN: number of long sequences in that cluster')
        print('\tC: unique identifier of the cluster')
        drawAsciiTree('clsReps.aln.midpoint.treefile')
        print('\n')
        # ********************************************************
        
    # make individual long cluster alignments
    msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
    msg += ' Aligning long sequence clusters\n' 
    print(msg)
    
    # align long sequence clusters and create <long.xx.aln> files
    if mArgs.run == 'G': # run GNU parallel version
        alnLongSequenceClustersGNUParallel(numClusters,mArgs.thread,mArgs.mIterL,cDir,tFileName,zName) 
        
    elif mArgs.run == 'J':
        alnLongSequenceClustersJoblibParallel(numClusters,mArgs.thread,mArgs.mIterL,cDir,tFileName,zName)  
        
    
    # make HMMs and database from long sequence cluster alignments
    makeHMMdbParallel(numClusters,mArgs.thread,mArgs.alphabet,'hmmer.log',cDir,tFileName,zName)

    #***************
    # only create long sequence cluster alignments and HMMs and exit
    if mArgs.stage == 1:
        longSeqClusters(mArgs.tempDirPath,cDir,zName,tDirName,tFileName)

        # No path provided for temporary directory
        if mArgs.tempDirPath is None:
            try:
                wName = cDir + '/' + zName 
                shutil.copytree(tDirName,wName) 
                os.chmod(wName,0o777)
            except OSError as e:
                print(e)
                cZip(cDir,tFileName,zName)

        if mArgs.tempDirPath:
            os.chdir(cDir)
            if os.path.exists(zName):
                os.chmod(zName,0o777)    

        msg = '\nAlignment files and HMMs can be found in <%s>\n' % zName
        msg += '\tLong sequence alignments have names <long.xx.aln>\n'
        msg += '\tHMM file written in <long.xx.hmm>\n'
        msg += '\tHMM database written in <pipelign.hmm>\n'
        msg += '\nThank you for using Pipelign\n'
        sys.exit(msg)

    
    #***************
    # create alignment from only long sequences
    if mArgs.stage == 2:
        if mArgs.merge == 'P':
            longSeqAlignmentParallel(numClusters,mArgs,cDir,tFileName,zName)
        elif mArgs.merge == 'C':
            longSeqAlignmentConsensus(numClusters,mArgs,cDir,tFileName,zName)
        
        # remove clusters
        if mArgs.excludeClusters:
            excludeClustersFromalignment('final.all.aln','final.aln')
        else:
            copyFile('final.all.aln','final.aln')
            
        if checkPresenceOfFile('final.aln'):
            copyFile('final.aln',mArgs.outFile)
    
    
        if checkPresenceOfFile(mArgs.outFile):
            msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
            msg += ' Final alignment written in <%s>\n' % mArgs.outFile
            print(msg)
        else:
            msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
            msg += ' Pipelign could not write the file <%s>\n' % mArgs.outFile
            msg += '\tSomething went wrong. Please check the zipped output directory.\n' 
            print(msg)
        
        # give permission to output alignment to user; used mainly for Docker
        if checkPresenceOfFile(mArgs.outFile):
            os.chmod(mArgs.outFile, 0o777)
            
        if mArgs.tempDirPath:
            os.chdir(cDir)
            if os.path.exists(zName):
                os.chmod(zName,0o777)    
        
        sys.exit('\nThank you for using Pipelign.\n')

    #***************
    # Assigns clusters to fragments based on BLAST and HMM search
    if mArgs.stage == 3:
        assignClusterstoFrags(numClusters,numFrag,mArgs,tDirName,cDir,tFileName,zName)
        
        # No path provided for temporary directory
        if mArgs.tempDirPath is None:
            try:
                wName = cDir + '/' + zName 
                shutil.copytree(tDirName,wName) 
                os.chmod(wName,0o777)
            except OSError as e:
                print(e)
                cZip(cDir,tFileName,zName)
        
        if mArgs.tempDirPath:
            os.chdir(cDir)
            if os.path.exists(zName):
                os.chmod(zName,0o777)    

        msg = '\nAlignment files and HMMs can be found in <%s>\n' % zName
        msg += '\tLong sequence alignments have names <long.xx.aln>\n'
        msg += '\tHMM file written in <long.xx.hmm>\n'
        msg += '\tHMM database written in <pipelign.hmm>\n'
        if checkPresenceOfFile('frag.fas'):
            msg += '\tCluster files for fragments are written in <frag.xx.fas>\n'
        if checkPresenceOfFile('frag.noClusters.fas'):
            msg += '\tUnassigned fragments are written in <frag.noClusters.fas>\n'
        msg += '\nThank you for using Pipelign\n'
        sys.exit(msg)
        
          
    #***************
    # make cluster alignments with fragments
    if mArgs.stage == 4:
        assignClusterstoFrags(numClusters,numFrag,mArgs,tDirName,cDir,tFileName,zName)
        
        if mArgs.run == 'G':
            addFragmentToClustersGNUParallel(numClusters,'fragAlign.log',mArgs,cDir,tFileName,zName)
        elif mArgs.run == 'J':
            addFragmentToClustersJoblibParallel(numClusters,'fragAlign.log',mArgs.thread,mArgs.mIterL,cDir,tFileName,zName)

        # No path provided for temporary directory
        if mArgs.tempDirPath is None:
            try:
                wName = cDir + '/' + zName 
                shutil.copytree(tDirName,wName) 
            except OSError as e:
                print(e)
                cZip(cDir,tFileName,zName)

        msg = '\nAlignment files and HMMs can be found in <%s>\n' % zName
        msg += '\tLong sequence alignments have names <long.xx.aln>\n'
        msg += '\tHMM file written in <long.xx.hmm>\n'
        msg += '\tHMM database written in <pipelign.hmm>\n'
        if checkPresenceOfFile('frag.fas'):
            msg += '\tCluster files for fragments are written in <frag.xx.fas>\n'
        if checkPresenceOfFile('frag.noClusters.fas'):    
            msg += '\tUnassigned fragments are written in <frag.noClusters.fas>\n'
        msg += '\tCluster alignments are written in <cls.xx.aln>\n'        
        msg += '\nThank you for using Pipelign\n'
        sys.exit(msg)
        
        
    #************
    # complete alignment requested, -n 5
    if mArgs.stage == 5:
        assignClusterstoFrags(numClusters,numFrag,mArgs,tDirName,cDir,tFileName,zName)
        
        if mArgs.run == 'G':
            addFragmentToClustersGNUParallel(numClusters,'fragAlign.log',mArgs,cDir,tFileName,zName)
        elif mArgs.run == 'J':
            addFragmentToClustersJoblibParallel(numClusters,'fragAlign.log',mArgs.thread,mArgs.mIterL,cDir,tFileName,zName)

        if numClusters == 1:
            copyFile('cls.0.aln','final.noOrphans.aln')
        #elif numClusters > 1:
            #mergeConsensus(numClusters,'final.noOrphans.aln','all',mArgs,'consensus.log',cDir,tFileName,zName)
        
        # only two clusters, merge the pair
        elif numClusters == 2:
            if mArgs.merge == 'P':
                mergePair('cls.0.aln','cls.1.aln','final.noOrphans.aln',mArgs.thread,mArgs.mIterM,'merge.log',cDir,tFileName,zName)
            elif mArgs.merge == 'C':
                mergeConsensus(numClusters,'final.noOrphans.aln','all',mArgs,'consensus.log',cDir,tFileName,zName)
    
        # more than 2 clusters
        elif numClusters >= 3:
            if mArgs.merge == 'P':
                resName = merge_gnu_parallel('clsReps.aln.midpoint.treefile','all',mArgs.thread,mArgs.mIterL,'merge.log',cDir,tFileName,zName)
                copyFile(resName,'final.noOrphans.aln')
            elif mArgs.merge == 'C':
                mergeConsensus(numClusters,'final.noOrphans.aln','all',mArgs,'consensus.log',cDir,tFileName,zName)
        #'''
            
        # add unclustered fragments
        if mArgs.keepOrphans and checkPresenceOfFile('frag.noClusters.fas'):
                # create alignment by merging cluster alignments
                msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
                msg += ' Adding orphan fragments to final alignment\n' 
                print(msg)

                addFragments('frag.noClusters.fas','final.noOrphans.aln','final.noBadSeq.aln',mArgs.thread,mArgs.mIterL,'merge.log',cDir,tFileName,zName)
        else:
            copyFile('final.noOrphans.aln','final.noBadSeq.aln')
 
        msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' <final.noBadSeq.aln> file created\n' 
        print(msg)

        # check if <long.bad.fas> needs to be added to the final alignment
        if mArgs.keepBadSeqs and checkPresenceOfFile('long.bad.fas'):
            addLongSeqs('long.bad.fas','final.noBadSeq.aln','final.all.aln',mArgs.thread,mArgs.mIterL,'addBadSeq.log',cDir,tFileName,zName) 
        else:
            copyFile('final.noBadSeq.aln','final.all.aln')
        
        # if user wants to exclude clusters from the final alignment
        # first show the cluster representative tree
        # then asks for list of clusters separated by space
        # remove all the sequences from the <final.all.aln> that belong to the clusters
        if mArgs.excludeClusters:
            excludeClustersFromalignment('final.all.aln','final.aln')
        else:
            copyFile('final.all.aln','final.aln')    

        if checkPresenceOfFile('final.aln'):
            copyFile('final.aln',mArgs.outFile)
        
        if checkPresenceOfFile(mArgs.outFile):
            msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
            msg += ' Final alignment written in <%s>\n' % mArgs.outFile
            print(msg)
            
            # give permission to output file; for docker
            os.chmod(mArgs.outFile,0o777)
            
            # change to parent directory
            os.chdir(cDir)
            
            # give permission to temporary directory
            if mArgs.tempDirPath:
                if os.path.exists(zName):
                    os.chmod(zName,0o777)
                 
        else:
            msg = '\n[' + time.strftime('%d %b %H:%M:%S') + ']'
            msg += ' Pipelign could not write the file <%s>\n' % mArgs.outFile
            msg += '\tSomething went wrong. Please check the output directory.\n' 
            print(msg)
            cZip(cDir,tFileName,zName)
        
    
        sys.exit('\nThank you for using Pipelign.\n')

