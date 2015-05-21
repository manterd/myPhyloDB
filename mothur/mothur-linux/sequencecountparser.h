#ifndef Mothur_sequencecountparser_h
#define Mothur_sequencecountparser_h

//
//  sequencecountparser.h
//  Mothur
//
//  Created by Sarah Westcott on 8/7/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "mothur.h"
#include "mothurout.h"
#include "sequence.hpp"
#include "counttable.h"

/* This class reads a fasta and count file and parses the data by group. The countfile must contain group information.
 
 Note: The sum of all the groups unique sequences will be larger than the original number of unique sequences. 
 This is because when we parse the count file we make a unique for each group instead of 1 unique for all
 groups. 
 
 */

class SequenceCountParser {
	
public:
	
    SequenceCountParser(string, string);			//count, fasta - file mismatches will set m->control_pressed = true
    SequenceCountParser(string, CountTable&);		//fasta, counttable - file mismatches will set m->control_pressed = true
    ~SequenceCountParser();
    
    //general operations
    int getNumGroups();
    vector<string> getNamesOfGroups();	
    
    int getNumSeqs(string);		//returns the number of unique sequences in a specific group
    vector<Sequence> getSeqs(string); //returns unique sequences in a specific group
    map<string, int> getCountTable(string); //returns seqName -> numberOfRedundantSeqs for a specific group - the count file format, but each line is parsed by group.
    
    int getSeqs(string, string, bool); //prints unique sequences in a specific group to a file - group, filename, uchimeFormat=false
    int getCountTable(string, string); //print seqName -> numberRedundantSeqs for a specific group - group, filename
    
    map<string, string> getAllSeqsMap(){ return allSeqsMap; }  //returns map where the key=sequenceName and the value=representativeSequence - helps us remove duplicates after group by group processing
private:
	
    CountTable countTable;
    MothurOut* m;
	
    int numSeqs;
    map<string, string> allSeqsMap;
    map<string, vector<Sequence> > seqs; //a vector for each group
    map<string, map<string, int> > countTablePerGroup; //countTable for each group
    vector<string> namesOfGroups;
};



#endif
