/*
 *  groupmap.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 12/1/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "groupmap.h"

/************************************************************/

 GroupMap::GroupMap(string filename) {
	m = MothurOut::getInstance();
	groupFileName = filename;
	m->openInputFile(filename, fileHandle);
	index = 0;
}

/************************************************************/
 GroupMap::~GroupMap(){}
/************************************************************/
int GroupMap::readMap() {
    try {
        string seqName, seqGroup;
		int error = 0;
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
    
        while (!fileHandle.eof()) {
            if (m->control_pressed) { fileHandle.close();  return 1; }
        
            fileHandle.read(buffer, 4096);
            vector<string> pieces = m->splitWhiteSpace(rest, buffer, fileHandle.gcount());
        
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  seqName = pieces[i]; columnOne=false; }
                else  { seqGroup = pieces[i]; pairDone = true; columnOne=true; }
            
                if (pairDone) { 
                    setNamesOfGroups(seqGroup);
                    
                    if (m->debug) { m->mothurOut("[DEBUG]: name = '" + seqName + "', group = '" + seqGroup + "'\n"); }
                    m->checkName(seqName);
                    it = groupmap.find(seqName);
                    
                    if (it != groupmap.end()) { error = 1; m->mothurOut("Your groupfile contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine();  }
                    else {
                        groupmap[seqName] = seqGroup;	//store data in map
                        seqsPerGroup[seqGroup]++;  //increment number of seqs in that group
                    }
                    pairDone = false; 
                } 
            }
        }
		fileHandle.close();
        
        if (rest != "") {
            vector<string> pieces = m->splitWhiteSpace(rest);
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  seqName = pieces[i]; columnOne=false; }
                else  { seqGroup = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    setNamesOfGroups(seqGroup);
                    
                    if (m->debug) { m->mothurOut("[DEBUG]: name = '" + seqName + "', group = '" + seqGroup + "'\n"); }
                    m->checkName(seqName);
                    it = groupmap.find(seqName);
                    
                    if (it != groupmap.end()) { error = 1; m->mothurOut("Your groupfile contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine();  }
                    else {
                        groupmap[seqName] = seqGroup;	//store data in map
                        seqsPerGroup[seqGroup]++;  //increment number of seqs in that group
                    }
                    pairDone = false; 
                } 
            }
        }
        
		m->setAllGroups(namesOfGroups);
		return error;
    }
	catch(exception& e) {
		m->errorOut(e, "GroupMap", "readMap");
		exit(1);
	}
}
/************************************************************/
int GroupMap::readDesignMap() {
    try {
        string seqName, seqGroup;
		int error = 0;
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        
        while (!fileHandle.eof()) {
            if (m->control_pressed) { fileHandle.close();  return 1; }
            
            fileHandle.read(buffer, 4096);
            vector<string> pieces = m->splitWhiteSpace(rest, buffer, fileHandle.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  seqName = pieces[i]; columnOne=false; }
                else  { seqGroup = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    setNamesOfGroups(seqGroup);
                    
                    if (m->debug) { m->mothurOut("[DEBUG]: name = '" + seqName + "', group = '" + seqGroup + "'\n"); }
                    m->checkName(seqName);
                    it = groupmap.find(seqName);
                    
                    if (it != groupmap.end()) { error = 1; m->mothurOut("Your designfile contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine();  }
                    else {
                        groupmap[seqName] = seqGroup;	//store data in map
                        seqsPerGroup[seqGroup]++;  //increment number of seqs in that group
                    }
                    pairDone = false; 
                } 
            }
        }
		fileHandle.close();
        
        if (rest != "") {
            vector<string> pieces = m->splitWhiteSpace(rest);
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  seqName = pieces[i]; columnOne=false; }
                else  { seqGroup = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    setNamesOfGroups(seqGroup);
                    
                    if (m->debug) { m->mothurOut("[DEBUG]: name = '" + seqName + "', group = '" + seqGroup + "'\n"); }
                    m->checkName(seqName);
                    it = groupmap.find(seqName);
                    
                    if (it != groupmap.end()) { error = 1; m->mothurOut("Your designfile contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine();  }
                    else {
                        groupmap[seqName] = seqGroup;	//store data in map
                        seqsPerGroup[seqGroup]++;  //increment number of seqs in that group
                    }
                    pairDone = false; 
                } 
            }

        }
        
		m->setAllGroups(namesOfGroups);
		return error;
    }
	catch(exception& e) {
		m->errorOut(e, "GroupMap", "readDesignMap");
		exit(1);
	}
}
/************************************************************/
int GroupMap::readMap(string filename) {
    try {
        groupFileName = filename;
        m->openInputFile(filename, fileHandle);
        index = 0;
        string seqName, seqGroup;
		int error = 0;
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        
        while (!fileHandle.eof()) {
            if (m->control_pressed) { fileHandle.close();  return 1; }
            
            fileHandle.read(buffer, 4096);
            vector<string> pieces = m->splitWhiteSpace(rest, buffer, fileHandle.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  seqName = pieces[i]; columnOne=false; }
                else  { seqGroup = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    setNamesOfGroups(seqGroup);
                    
                    if (m->debug) { m->mothurOut("[DEBUG]: name = '" + seqName + "', group = '" + seqGroup + "'\n"); }
                    m->checkName(seqName);
                    it = groupmap.find(seqName);
                    
                    if (it != groupmap.end()) { error = 1; m->mothurOut("Your group file contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine();  }
                    else {
                        groupmap[seqName] = seqGroup;	//store data in map
                        seqsPerGroup[seqGroup]++;  //increment number of seqs in that group
                    }
                    pairDone = false; 
                } 
            }
        }
		fileHandle.close();
        
        if (rest != "") {
            vector<string> pieces = m->splitWhiteSpace(rest);
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  seqName = pieces[i]; columnOne=false; }
                else  { seqGroup = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    setNamesOfGroups(seqGroup);
                    
                    if (m->debug) { m->mothurOut("[DEBUG]: name = '" + seqName + "', group = '" + seqGroup + "'\n"); }
                    m->checkName(seqName);
                    it = groupmap.find(seqName);
                    
                    if (it != groupmap.end()) { error = 1; m->mothurOut("Your group file contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine();  }
                    else {
                        groupmap[seqName] = seqGroup;	//store data in map
                        seqsPerGroup[seqGroup]++;  //increment number of seqs in that group
                    }
                    pairDone = false; 
                } 
            }
        }
        
		m->setAllGroups(namesOfGroups);
		return error;
    }
	catch(exception& e) {
		m->errorOut(e, "GroupMap", "readMap");
		exit(1);
	}
}
/************************************************************/
int GroupMap::readDesignMap(string filename) {
    try {
        groupFileName = filename;
        m->openInputFile(filename, fileHandle);
        index = 0;
        string seqName, seqGroup;
		int error = 0;
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        
        while (!fileHandle.eof()) {
            if (m->control_pressed) { fileHandle.close();  return 1; }
            
            fileHandle.read(buffer, 4096);
            vector<string> pieces = m->splitWhiteSpace(rest, buffer, fileHandle.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  seqName = pieces[i]; columnOne=false; }
                else  { seqGroup = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    setNamesOfGroups(seqGroup);
                    
                    if (m->debug) { m->mothurOut("[DEBUG]: name = '" + seqName + "', group = '" + seqGroup + "'\n"); }
                    m->checkName(seqName);
                    it = groupmap.find(seqName);
                    
                    if (it != groupmap.end()) { error = 1; m->mothurOut("Your designfile contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine();  }
                    else {
                        groupmap[seqName] = seqGroup;	//store data in map
                        seqsPerGroup[seqGroup]++;  //increment number of seqs in that group
                    }
                    pairDone = false; 
                } 
            }
        }
		fileHandle.close();
        
        if (rest != "") {
            vector<string> pieces = m->splitWhiteSpace(rest);
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  seqName = pieces[i]; columnOne=false; }
                else  { seqGroup = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    setNamesOfGroups(seqGroup);
                    
                    if (m->debug) { m->mothurOut("[DEBUG]: name = '" + seqName + "', group = '" + seqGroup + "'\n"); }
                    m->checkName(seqName);
                    it = groupmap.find(seqName);
                    
                    if (it != groupmap.end()) { error = 1; m->mothurOut("Your designfile contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine();  }
                    else {
                        groupmap[seqName] = seqGroup;	//store data in map
                        seqsPerGroup[seqGroup]++;  //increment number of seqs in that group
                    }
                    pairDone = false; 
                } 
            }
        }
        
		m->setAllGroups(namesOfGroups);
		return error;
    }
	catch(exception& e) {
		m->errorOut(e, "GroupMap", "readDesignMap");
		exit(1);
	}
}
/************************************************************/
int GroupMap::getNumGroups() { return namesOfGroups.size();	}
/************************************************************/

string GroupMap::getGroup(string sequenceName) {
			
	it = groupmap.find(sequenceName);
	if (it != groupmap.end()) { //sequence name was in group file
		return it->second;	
	}else {
        //look for it in names of groups to see if the user accidently used the wrong file
        if (m->inUsersGroups(sequenceName, namesOfGroups)) {
            m->mothurOut("[WARNING]: Your group or design file contains a group named " + sequenceName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group."); m->mothurOutEndLine(); 
        }
		return "not found";
	}
}

/************************************************************/

void GroupMap::setGroup(string sequenceName, string groupN) {
	setNamesOfGroups(groupN);
	m->checkName(sequenceName);
	it = groupmap.find(sequenceName);
	
	if (it != groupmap.end()) {  m->mothurOut("Your groupfile contains more than 1 sequence named " + sequenceName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine();  }
	else {
		groupmap[sequenceName] = groupN;	//store data in map
		seqsPerGroup[groupN]++;  //increment number of seqs in that group
	}
}

/************************************************************/
void GroupMap::setNamesOfGroups(string seqGroup) {
	int i, count;
	count = 0;
	for (i=0; i<namesOfGroups.size(); i++) {
		if (namesOfGroups[i] != seqGroup) {
			count++; //you have not found this group
		}else {
			break; //you already have it
		}
	}
	if (count == namesOfGroups.size()) {
		namesOfGroups.push_back(seqGroup); //new group
		seqsPerGroup[seqGroup] = 0;
		groupIndex[seqGroup] = index;
		index++;
	}
}
/************************************************************/
bool GroupMap::isValidGroup(string groupname) {
	try {
		for (int i = 0; i < namesOfGroups.size(); i++) {
			if (groupname == namesOfGroups[i]) { return true; }
		}
		
		return false;
	}
	catch(exception& e) {
		m->errorOut(e, "GroupMap", "isValidGroup");
		exit(1);
	}
}
/************************************************************/
int GroupMap::getCopy(GroupMap* g) {
	try {
        vector<string> names = g->getNamesSeqs();
        for (int i = 0; i < names.size(); i++) {
            if (m->control_pressed) { break; }
            string group = g->getGroup(names[i]);
            setGroup(names[i], group);
        }
        return names.size();
	}
	catch(exception& e) {
		m->errorOut(e, "GroupMap", "getCopy");
		exit(1);
	}
}
/************************************************************/
int GroupMap::getNumSeqs(string group) {
	try {
		
		map<string, int>::iterator itNum;
		
		itNum = seqsPerGroup.find(group);
		
		if (itNum == seqsPerGroup.end()) { return 0; }
		
		return seqsPerGroup[group];
		
	}
	catch(exception& e) {
		m->errorOut(e, "GroupMap", "getNumSeqs");
		exit(1);
	}
}
/************************************************************/
int GroupMap::renameSeq(string oldName, string newName) {
	try {
		
		map<string, string>::iterator itName;
		
		itName = groupmap.find(oldName);
		
		if (itName == groupmap.end()) {
            m->mothurOut("[ERROR]: cannot find " + toString(oldName) + " in group file");
            m->control_pressed = true;
            return 0;
        }else {
            string group = itName->second;
            groupmap.erase(itName);
            groupmap[newName] = group;
        }
        
        return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "GroupMap", "renameSeq");
		exit(1);
	}
}
/************************************************************/
int GroupMap::print(ofstream& out) {
	try {
		
		for (map<string, string>::iterator itName = groupmap.begin(); itName != groupmap.end(); itName++) {
            out << itName->first << '\t' << itName->second << endl;
        }
             
        return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "GroupMap", "print");
		exit(1);
	}
}
/************************************************************/
int GroupMap::print(ofstream& out, vector<string> userGroups) {
	try {
		
		for (map<string, string>::iterator itName = groupmap.begin(); itName != groupmap.end(); itName++) {
            if (m->inUsersGroups(itName->second, userGroups)) {
                out << itName->first << '\t' << itName->second << endl;
            }
        }
        
        return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "GroupMap", "print");
		exit(1);
	}
}
/************************************************************/
vector<string> GroupMap::getNamesSeqs(){
	try {
	
		vector<string> names;
		
		for (it = groupmap.begin(); it != groupmap.end(); it++) {
			names.push_back(it->first);
		}
		
		return names;
	}
	catch(exception& e) {
		m->errorOut(e, "GroupMap", "getNamesSeqs");
		exit(1);
	}
}
/************************************************************/
vector<string> GroupMap::getNamesSeqs(vector<string> picked){
	try {
		
		vector<string> names;
		
		for (it = groupmap.begin(); it != groupmap.end(); it++) {
			//if you are belong to one the the groups in the picked vector add you
			if (m->inUsersGroups(it->second, picked)) {
				names.push_back(it->first);
			}
		}
		
		return names;
	}
	catch(exception& e) {
		m->errorOut(e, "GroupMap", "getNamesSeqs");
		exit(1);
	}
}

/************************************************************/

