#!/usr/bin/env python3

'''
  Pipelign.py

  A python based program to align virus sequences.

  The program takes as input a single FASTA formatted file

  Returns a FASTA formatted alignment file

'''
#*********************************************************************

import sys, os, shutil, subprocess, argparse, textwrap
from Bio import SeqIO, Phylo
import tempfile
import time
import joblib
import math

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

def iterThreshold(x):
  x = int(x)
  if x < 1:
    raise argparse.ArgumentTypeError('%r must be a non-zero positive integer' % (x,))
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

def linsi(seqFile,alnFile,thread,mIterLong,log,cDir,tName,zName,p=None):
  '''
    - aligns sequences using MAFFT's L-INS-i method
  '''
  
  if p:
    lName = log + '.' + str(p)
    lh = open(lName,'a')
  else:
    lh = open(log,'a')
  
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
  
#***********************************************************************

def fftnsi(seqFile,alnFile,thread,mIterLong,log,cDir,tName,zName,p=None):
  '''
    - aligns sequences using MAFFT's FFT-NS-i method
  '''
  if p:
    lName = log + '.' + str(p)
    lh = open(lName,'a')
  else:
    lh = open(log,'a')
  
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

#***********************************************************************

def checkPresenceOfFile(fName):
  '''
    - Check whether file <fName> exists;
    - Terminate Pipelign run otherwise
  '''
   
  if not os.path.exists(fName) or os.stat(fName).st_size == 0:
      msg = '\nError: the file <%s>  could not be found in the given path' % fName
      msg += '\n**Please check the name and the path of the file**'
      msg += '\n\nPipelign is exiting\n'
      sys.exit(msg)

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
  
    msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
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
  
  msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
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
  msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
  
  # open log file
  lh = open('pipelign.log','a') 
  
  # argument string for CD-HIT
  if alphabet == 'dna' or alphabet == 'rna':
    cl = ['cd-hit-est','-c', str(per), '-n', '5', '-i', longName, '-o', 'grp', '-d', '0', '-T', str(thread)]
    msg += ' CD-HIT-EST started\n'

  elif alphabet == 'aa':
    cl = ['cd-hit','-c', str(per), '-n', '5', '-i', longName, '-o', 'grp', '-d', '0', '-T', str(thread)] 
    msg += ' CD-HIT started\n'

  print(msg)
    
  try:
    subprocess.check_call(cl, stdout=lh, stderr=lh)
  except subprocess.CalledProcessError as e:
    print(e)
    cZip(cDir,tName,zName)
  
  lh.close() # close log file   
  
  msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
  msg += ' CD-HIT finished. Files created: <grp> and <grp.clstr>\n'
  print(msg)

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
  if len(seqs) < 2:
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
    msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
    msg += ' %d cluster file(s) created. File names <long.0.fas>, <long.1.fas>, ..., <long.xx.fas> \n' % cls
    print(msg)
  else:
    msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
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
      msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
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
  
  lh = open('pipelign.log','a')
  
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
    - draws ASCII tree on the console
  '''
  
  # Read the tree in newick format
  tree = Phylo.read(treeFile,'newick')

  Phylo.draw_ascii(tree)
  
#***********************************************************************

def alnLongSequenceClustersParallel(nClusters,thread,mIterL,cDir,tName,zName):
  '''
    - align long sequence cluster files in parallel using joblib 
  '''
  
  log = 'align.log'
  #lh = open('pipelign.log','a')
  
  msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
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

  msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
  msg += ' Finished aligning long sequence clusters. Created files: <long.x.aln>\n'
  print(msg)
    
#***********************************************************************

def buildHMM(aName,hName,thread,alphabet,log,cDir,tName,zName,p=None):
  '''
    - Builds HMM for an alignment
  '''
  
  if p:
    lName = log + '.' + str(p)
    lh = open(lName,'a')
  else:
    lh = open(log,'a')
  
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
  
  lh = open('pipelign.log','a')
  
  # run command
  try:
    subprocess.check_call(cl, stdout=lh, stderr=lh)
  except subprocess.CalledProcessError as e:
    #sys.exit(e)
    print(e)
    cZip(cDir,tName)

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
    msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
    msg += ' Started building HMMs from long cluster alignments\n'
    print(msg)

    num_parallel_jobs = int(nClusters/thread) if nClusters > thread else thread
    num_threads_per_job = int(thread/nClusters) if thread > nClusters else 1
    joblib.Parallel(n_jobs=num_parallel_jobs)(joblib.delayed(buildHMM)(x[0],x[1],num_threads_per_job,alphabet,log,cDir,tName,zName,x[0].split('.')[1])  for x in to_run_tuples)
    
    msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
    msg += ' Building HMM database\n'
    print(msg)
    
    buildHMMdb(nClusters)

    msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
    msg += ' HMM database created\n'
    print(msg)
    
#***********************************************************************

def addFragments(fName,aName,oName,thread,mIterL,log,cDir,tName,zName,p=None):
  '''
    - add fragments to the alignment using MAFFT's --addfragments
  '''

  if p:
    lName = log + '.' + str(p)
    lh = open(lName,'w')
  else:
    lh = open(log,'a')
  
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

#***********************************************************************

def addLongSeqs(sName,aName,oName,thread,mIterL,log,cDir,tName,zName):
  '''
    - add long sequences to the alignment using MAFFT's --addfull
  '''

  lh = open(log,'a')
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

  
#***********************************************************************

def mergeLongClusters(numClusters,outFile,thread,mIterM,cDir,tName,zName):
  '''
    - merges long sequence clusters <long.0.aln>, <long.1.aln>, ..., <long.xx.aln>
    - creates output alignment out.aln
  '''
  
  # only one cluster 
  if numClusters == 1:
    copyFile('long.0.aln',outFile)
    return
    
  elif numClusters > 1: # more than one cluster
    oLSeqs = [] # contains long sequences from clusters that have single members
    mFlag = 0 # is > 0 if at least one cluster contains multiple sequences
    #fFlag = False
    clsInd = None # holds the index of the alignment cluster 
    
    seqCount = 1
    catText = 'cat '
    mTab = ''
    
    for i in range(numClusters):
      cName = 'long.' + str(i) + '.aln' # name of cluster alignment
      seqs = list(SeqIO.parse(cName,'fasta'))
      if len(seqs) > 1:
        catText += cName + ' ' # add file name to cat string
        for k in range(len(seqs)):
          mTab += str(seqCount) + ' '
          seqCount += 1
        mTab += '\n'
        mFlag += 1
        clsInd = i
      
      else: # only one sequence in the cluster
        oLSeqs.append(seqs[0]) # add to singles sequence list
    
    # write single sequences into a file
    if len(oLSeqs) > 0: 
      SeqIO.write(oLSeqs,'longSingles.fas','fasta') 
      catText += 'longSingles.fas '     # add that file for merging
    
    
    #catText += '> merge'
    
    if mFlag == 0: # no cluster alignments, only single sequences in cluster
      if len(oLSeqs) > 1: # check list is not empty
        cl = ['mafft', '--preservecase', '--thread', str(thread), '--localpair', '--maxiterate', str(mIterM), 'longSingles.fas'] 
        
        lh = open('pipelign.log','a')
        ah = open('final.aln','w')
        
        try:
          subprocess.check_call(cl,stdout=ah,stderr=lh)
        except subprocess.CalledProcessError as e:
          print(e)
          cZip(cDir,tName,zName)
        
        lh.close()
        ah.close()    
               
      copyFile('final.aln',outFile)
      return  
    
    if mFlag == 1: # only one multiple sequences alignment cluster
      # don't merge, add unaligned sequences with --addFragments
      sName = 'long.' + str(clsInd) + '.aln'
      
      addLongSeqs('longSingles.fas',sName,'final.aln',thread,mIterM,'pipelign.log',cDir,tName,zName) 
      
      return
    
    if mFlag > 1: # at least two multiple sequence alignment 
      lh = open('pipelign.log','a')
      mh = open('merge','w')
      
      try: # create the input sequence file
        
        cl = catText.split() 
        subprocess.check_call(cl,stdout=mh,stderr=lh)
      except subprocess.CalledProcessError as e:
        print(e)
        cZip(cDir,tName,zName)
      
      lh.close()
      mh.close()
      
      # write subMSAtable
      fh = open('subMSAtable','w')
      fh.write(mTab)
      fh.close()  
      
      numSeqs = seqCount + len(oLSeqs)
      if numSeqs <= 100:
        cl = ['mafft', '--preservecase', '--thread', str(thread), '--localpair', '--maxiterate', str(mIterM), '--merge', 'subMSAtable', 'merge'] # L-INS-i
      
      else: # more than 100 sequences, use FFT-NS-i 
        cl = ['mafft', '--preservecase', '--thread', str(thread), '--maxiterate', str(mIterM), '--merge', 'subMSAtable', 'merge']
      
      lh = open('pipelign.log','a')
      ah = open('final.aln','w')
        
      try: # create the input sequence file
        subprocess.check_call(cl,stdout=ah,stderr=lh)
      except subprocess.CalledProcessError as e:
        print(e)
        cZip(cDir,tName,zName)
      
      lh.close()
      ah.close()
  
      copyFile('final.aln',outFile)
      return
   
#***********************************************************************

def runBlast(log):
  '''
    - creates BLAST database from cluster representatives
    - search fragments against BLAST database
    - assign clusters to fragments and writes into file <frags.ClusterList.txt> 
  '''
  # first create blast database of cluster reps
  bl = open(log,'a')

  cStr = 'makeblastdb -in clsReps.fas -input_type fasta -title pipelign -dbtype '
  if mArgs.alphabet in ['dna','rna']:
    cStr += 'nucl '
  elif mArgs.alphabet == 'aa':
    cStr += 'prot '
      
  cStr += '-out pipelign.blastdb'   

  msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
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
      
  cStr += '-query frag.fas -db pipelign.blastdb -max_target_seqs 1 -outfmt 6 -evalue 20'
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
  
  return mids  
  
#************************************************************************

def searchHMMdb(log,thread,alpha,res,cDir,tName,zName):
  '''
    HMM database is searched with the fragments to assign them a cliuster for alignment
  '''  
  
  lh = open(log,'a')
  
  # generate the command for HMM search
  if alpha == 'dna' or alpha == 'rna':
    cStr = 'nhmmscan --cpu %d --tblout %s --noali pipelign.hmm frag.hmm.fas' % (thread,res)
      
  elif alpha == 'aa':
    cStr = 'hmmscan --cpu %d --tblout %s --noali pipelign.hmm frag.hmm.fas' % (thread,res)
    
  msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
  msg += ' Searching HMM database to assign clusters to remaining fragments\n'
  print(msg)
  
  try:
    cl = cStr.split()
    subprocess.check_call(cl,stdout=lh,stderr=lh)
  except subprocess.CalledProcessError as e:
    print(e)
    cZip(cDir,tName,zName)
  
  lh.close()

#**********************************************************************

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
  handle = open('frag.hmm.fas','rU')
  for seq in SeqIO.parse(handle,'fasta'):
    if seq.id in qids: # if the id matched
      ind = qids.index(seq.id)
      fh.write("%s\t%s\n" %(seq.id,hids[ind]))
  
  fh.close()    
  
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

def addFragmentToClustersParallel(nClusters,log,thread,mIterL,cDir,tName,zName):
  '''
    This will add fragments to their cluster alignments using joblib
      - MAFFT's -addfragment will be used
      - will create cluster alignment files: <cls.x.aln>
  '''
  
  msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
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

  msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
  msg += ' Finished adding fragments to long sequence clusters. Created files: <cls.x.aln>\n'
  print(msg)
      

#************************************************************************

def mergeAllClusters(numClusters,thread,mIterM,cDir,tName,zName,allAssigned,outFile):
  '''
    This will merge all the cluster alignments with long sequecnes and fragments
  '''
  
  # a single cluster alignment
  if numClusters == 1:
    # all fragments are assigned clusters
    if allAssigned: 
      copyFile('cls.0.aln','final.aln')
    
    # not all fragments are assigned and some orphans need to be aligned
    else: 
      addFragments('frag.noClusters.fas','cls.0.aln','final.aln',thread,mIterM,'pipelign.log',cDir,tName, zName)
      
    copyFile('final.aln',outFile)
    return
  
  elif numClusters > 1: # more than one clusters: merge
    oSeqs = [] # contains single long sequences and orphan fragments
    mFlag = 0 # > 0 if at least one cluster contains multiple sequences
    clsInd = None # holds the last index of the cluster alignment
    
    seqCount = 1 # for subMSAtable content
    catText = 'cat '
    mTab = ''
    
    for i in range(numClusters):
      cName = 'cls.' + str(i) + '.aln' # name of cluster alignment
      seqs = list(SeqIO.parse(cName,'fasta'))
      if len(seqs) > 1:
        catText += cName + ' ' # add file name to cat string
        for k in range(len(seqs)):
          mTab += str(seqCount) + ' '
          seqCount += 1
        mTab += '\n'
        mFlag += 1
        clsInd = i
      
      else:
        oSeqs.append(seqs[0]) # add to long orphan file if only sequence in the cluster    
    
    if not allAssigned:
      # add orphan sequences to oSeqs list
      nSeqs = list(SeqIO.parse('frag.noClusters.fas','fasta'))
      for seq in nSeqs:
        oSeqs.append(seq)
    
    if len(oSeqs) > 0: # at least one single unaligned sequence present
      SeqIO.write(oSeqs,'singles.fas','fasta')
      catText += 'singles.fas '
    
    #catText += '> merge'
    
    if mFlag == 0: # no cluster alignments, only single sequences in clusters and orphan fragments
      if len(oSeqs) > 1: # check list is not empty
        cl = ['mafft', '--preservecase', '--thread', str(thread), '--localpair', '--maxiterate', str(mIterM), 'singles.fas'] 
        
        lh = open('pipelign.log','a')
        ah = open('final.aln','w')
        
        try:
          subprocess.check_call(cl,stdout=ah,stderr=lh)
        except subprocess.CalledProcessError as e:
          print(e)
          cZip(cDir,tName,zName)
        
        lh.close()
        ah.close()    
               
      return  


    if mFlag == 1: # only one multiple sequence alignment cluster
      # don't merge, add unaligned sequences with --addFragments
      sName = 'cls.' + str(clsInd) + '.aln'
      addFragments('singles.fas',sName,'final.aln',thread,mIterM,'pipelign.log',cDir,tName,zName)
      
      copyFile('final.aln',outFile)

      return

            
    elif mFlag > 1: # at least two multiple sequence alignment 
      lh = open('pipelign.log','a')
      mh = open('merge','w')
      
      try: # create the input text file
        cl = catText.split()
        subprocess.check_call(cl,stdout=mh,stderr=lh)
      except subprocess.CalledProcessError as e:
        print(e)
        cZip(cDir,tName,zName)
      
      lh.close()
      mh.close()
      
      # write subMSAtable
      fh = open('subMSAtable','w')
      fh.write(mTab)
      fh.close()  
      
      numSeqs = seqCount + len(oSeqs)
      if numSeqs <= 100:
        cStr = 'mafft --preservecase --thread %d --localpair --maxiterate %d --merge subMSAtable merge' % (thread,mIterM) # uses L-INS-i
      else: # too many sequences, choose fft-ns-i
        cStr = 'mafft --preservecase --thread %d --maxiterate %d --merge subMSAtable merge' % (thread,mIterM)
      
      lh = open('pipelign.log','a')
      ah = open('final.aln','w')
      
      try:
        cl = cStr.split()
        subprocess.check_call(cl,stdout=ah,stderr=lh)
      except subprocess.CalledProcessError as e:
        print(e)
        cZip(cDir,tName,zName)

      
      lh.close()
      ah.close()
      
      copyFile('final.aln',outFile)
      return
      
#************************************************************************
  
def getArguments():
  '''
    Parses all the command line arguments from the user
  
  '''
  parser = argparse.ArgumentParser(description="Pipelign: creates multiple sequence alignment from FASTA formatted sequence file", formatter_class=argparse.RawTextHelpFormatter) #formatter_class=argparse.RawDescriptionHelpFormatter)  
  parser.add_argument('-i', '--inFile', required=True, help="Input sequence file")
  parser.add_argument('-o', '--outFile', required=True, help="Output alignment file")
  parser.add_argument('-t', '--lenThr', type=lengthThreshold, help="Length threshold for full sequences", default=0.9)
  parser.add_argument('-c', '--code', type=int, help="Genetic code for translation",default=1, choices=[1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25])
  parser.add_argument('-a', '--alphabet', required=True, help='Input sequences can be DNA/Protein', choices=['dna','aa','rna'], default='dna')
  parser.add_argument('-f', '--keepOrphans', help='Add fragments without clusters', action="store_true")
  parser.add_argument('-z', '--mZip', help='Create zipped temporary files', action="store_true")
  parser.add_argument('-p', '--simPer', type=lengthThreshold, help="percent sequence similarity for clustering", default=0.8)
  parser.add_argument('-q', '--thread', nargs='?', const=1, type=int, help="Number of CPU to use for multithreads", default=1)
  parser.add_argument('-s', '--mIterateLong', type=iterThreshold, help="Number of iterations to refine long alignments", default=2)
  parser.add_argument('-m', '--mIterateMerge', type=iterThreshold, help="Number of iterations to refine merged alignment", default=2)
  parser.add_argument('-d', '--tempDirPath', required=False, help="Path for temporary directory",default=None)
  #parser.add_argument('-l', '--longSeqsOnly', help='Only align long sequences', action="store_true")
  parser.add_argument('-n', '--stage', type=int,default=4, choices=[1,2,3,4],
    help=textwrap.dedent('''\
    1  Make alignments and HMMs of long sequence clusters
    2  Align only long sequences
    3  Make cluster alignments with fragments
    4  Align all sequences
    5  Filter sequences from alignment
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
      makeZip = args.mZip,
      simPer = args.simPer,
      thread = args.thread,
      mIterL = args.mIterateLong,
      mIterM = args.mIterateMerge,
      fragEmpty = 1,
      longName = 'long.fas',
      fragName = 'frag.fas',
      tempDirPath = args.tempDirPath,
      #longSeqsOnly = args.longSeqsOnly,
      stage = args.stage,
      excludeClusters = args.excludeClusters)
  
  # Pipeline process starts here
  cDir = os.getcwd() # save current working directory
  tFileName = 'input.fas' # input file copied to temporary directory
  deFileName = 'input.dealign.fas' # dealigned input file

  
  # get current time and generate name for the zipped temporary directory
  timeNow = time.strftime('%Y-%m-%d-%H%M%S')
  zName = 'pipelign.' + timeNow 
  
  # check whether input file exists
  checkPresenceOfFile(mArgs.inFile)
  
  # check if mpi environment is present
  #mpi = checkMPI()
  
  # create temporary directory 
  tempDir, tDirName = makeTempDir(mArgs.tempDirPath)
  
  #print('%s %s' % (tempDir,tName))
  
  # copy input file inside the temporary directory
  dName = tDirName + '/' + tFileName
  copyFile(mArgs.inFile,dName)
  
  msg = '\n[' + time.strftime('%d %b %H:%m:%S') + ']'
  msg += ' Pipelign started with sequence file <%s>\n' % mArgs.inFile 
  print(msg)
  
  
  # change current working directory to the temp
  os.chdir(tDirName)
  
  # removes any possible gaps from the sequence file
  deAlign(tFileName, deFileName) 
  
  # make separate sequence file for long sequences and fragments
  mArgs.fragEmpty, numFrag = separateFullFragment(deFileName, mArgs.lenThr, mArgs.longName, mArgs.fragName)
  
  # use CD-HIT to group sequences based on similarity
  runCDHIT(mArgs.longName, mArgs.alphabet, mArgs.simPer, mArgs.thread,cDir,tFileName,zName)  
  
  # create separate files for each cluster
  numClusters = makeClusters(mArgs.longName,'grp.clstr')

  # add cluster numbers and size to header
  addClusterNumberToReps('grp','long.ClusterList.txt','clsReps.fas')
  
  if numClusters < 3: # less than 3 clusters, tree will not be built
    msg = '\n[' + time.strftime('%d %b %H:%m:%S') + ']'
    msg += ' Only %d cluster(s) created. Phylogenetic tree cannot be built\n' % numClusters
    print(msg)
  
  else: # at least 3 clusters
    # align cluster reps
    if numClusters <= 100:
      linsi('clsReps.fas','clsReps.aln',mArgs.thread,mArgs.mIterL,'pipelign.log',cDir,tFileName,zName)
    else:
      fftnsi('clsReps.fas','clsReps.aln',mArgs.thread,mArgs.mIterL,'pipelign.log',cDir,tFileName,zName)
    
    msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
    msg += ' Alignment of cluster representatives written in <clsReps.aln>\n' 
    print(msg)
    
    # generate phylogenetic tree using IQ-TREE  
    makeIQTree('clsReps.aln',mArgs.thread,cDir,tFileName,zName,mArgs.alphabet)

    msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
    msg += ' Phylogenetic tree of the cluster representatives written in <clsReps.aln.treefile>\n' 
    print(msg)
    
    # make midpoint root tree for cluster representatives
    makeMidPointRootTree('clsReps.aln.treefile')
    
    msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
    msg += ' Midpoint rooted tree written in <clsReps.aln.midpoint.treefile>\n' 
    print(msg)
    
    # draw midpoint root tree
    print('\nMidpoint rooted tree is drawn below:\n')
    drawAsciiTree('clsReps.aln.midpoint.treefile')
    print('\n')
  
  # make individual long cluster alignments
  alnLongSequenceClustersParallel(numClusters,mArgs.thread,mArgs.mIterL,cDir,tFileName,zName) 
  
  # make HMMs and database from cluster alignments
  makeHMMdbParallel(numClusters,mArgs.thread,mArgs.alphabet,'hmmer.log',cDir,tFileName,zName)
  
  # only create long sequence alignments and HMMs
  if mArgs.stage == 1:
    msg = '\nPipelign has created long sequence alignments and HMMs\n'
    
    # copy entire work directory inside current directory, name it <pipelign.run>
    
    # No path provided for temporary directory
    if mArgs.tempDirPath is None:
      try:
        wName = cDir + '/' + zName 
        shutil.copytree(tDirName,wName) 
      except OSError as e:
        print(e)
        cZip(cDir,tFileName,zName)
    
    msg += '\nAlignment files and HMMs can be found in <%s>\n' % zName
    print(msg)
    
    sys.exit('\nPipelign has finished working.\n')
    
  # create alignment for only long sequences
  if mArgs.stage == 2:
    msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
    msg += ' Merging long sequence clusters\n'
    print(msg)
 
    mergeLongClusters(numClusters,mArgs.outFile,mArgs.thread,mArgs.mIterM,cDir,tFileName,zName) 
    
    msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
    msg += ' Final alignment written in %s\n' % mArgs.outFile
    print(msg)
    

    sys.exit('\nPipelign has finished working.\n')
      
  # if no fragments present
  if mArgs.fragEmpty:
    print('\nNo fragments present in the input file')
      
    msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
    msg += ' Merging long sequence clusters\n'
    print(msg)

    mergeLongClusters(numClusters,mArgs.outFile,mArgs.thread,mArgs.mIterM,cDir,tFileName,zName) 
  
    msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
    msg += ' Final alignment written in %s\n' % mArgs.outFile
    print(msg)
    
  else: # at least one fragment present
    
    # fragments need to be assigned clusters
    # we create BLASTDB using cluster reps
    # Search fragments against BLASTDB
    # assign clusters based on blast search
      
    mIds = runBlast('pipelign.log') # running BLAST on the fragments      
    fragWithClusters = len(mIds) 
    #print(fragWithClusters)

    msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
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
        SeqIO.write(noMatchFrags,'frag.hmm.fas','fasta')
          
      # search remaining fragments against HMM_DB
      searchHMMdb('pipelign.log', mArgs.thread, mArgs.alphabet, 'hmm.out', cDir,tFileName, zName)
      
      # get the clusters from HMM result
      getHMMclusters()
        
    # create clusters of fragments, un-assigned ones are written in <frag.noCluster.fas>
    allAssigned = createFragmentClusters(numClusters)  
          
    addFragmentToClustersParallel(numClusters,'fragAlign.log',mArgs.thread,mArgs.mIterL,cDir,tFileName,zName)
    
    # only alignments with fragments were asked  
    if mArgs.stage == 3:
      # copy entire work directory inside current directory, name it <pipelign.run>
    
      # No path provided for temporary directory
      if mArgs.tempDirPath is None:
        try:
          wName = cDir + '/' + zName 
          shutil.copytree(tDirName,wName) 
        except OSError as e:
          print(e)
          cZip(cDir,tFileName,zName)
    
      msg += '\nCluster alignments with fragments can be found in <%s>\n' % zName
      print(msg)
      
      sys.exit('\nPipelign has finished working.\n')      
      
  # create final alignment by merging cluster alignments
  msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
  msg += ' Merging cluster alignments\n' 
  print(msg)

  mergeAllClusters(numClusters,mArgs.thread,mArgs.mIterM,cDir,tFileName,zName,allAssigned,mArgs.outFile)
  
  # create final alignment by merging cluster alignments
  msg = '[' + time.strftime('%d %b %H:%m:%S') + ']'
  msg += ' Pipelign finished merging cluster alignments\n' 
  print(msg)

  print('\n\nFinal alignment written in %s\n' % mArgs.outFile)
  sys.exit('\nPipelign has finished working.\n')
