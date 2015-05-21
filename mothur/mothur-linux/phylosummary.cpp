/*
 *  rawTrainingDataMaker.cpp
 *  Mothur
 *
 *  Created by westcott on 4/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "phylosummary.h"
#include "referencedb.h"
/**************************************************************************************************/

PhyloSummary::PhyloSummary(string refTfile, CountTable* c, bool r){
	try {
		m = MothurOut::getInstance();
		maxLevel = 0;
		ignore = false;
        numSeqs = 0;
        relabund = r;
		
		ct = c;
        groupmap = NULL;
        
		//check for necessary files
        if (refTfile == "saved") { ReferenceDB* rdb = ReferenceDB::getInstance(); refTfile = rdb->getSavedTaxonomy(); }
		string taxFileNameTest = m->getFullPathName((refTfile.substr(0,refTfile.find_last_of(".")+1) + "tree.sum"));
		ifstream FileTest(taxFileNameTest.c_str());
		
		if (!FileTest) { 
			m->mothurOut("Error: can't find " + taxFileNameTest + "."); m->mothurOutEndLine(); exit(1);
		}else{
			readTreeStruct(FileTest);
		}
		
		tree[0].rank = "0";
		assignRank(0);
        
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "PhyloSummary");
		exit(1);
	}
}

/**************************************************************************************************/

PhyloSummary::PhyloSummary(CountTable* c, bool r){
	try {
		m = MothurOut::getInstance();
		maxLevel = 0;
		ignore = true;
        numSeqs = 0;
        relabund = r;
		
		ct = c;
        groupmap = NULL;
		
		tree.push_back(rawTaxNode("Root"));
		tree[0].rank = "0";
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "PhyloSummary");
		exit(1);
	}
}
/**************************************************************************************************/
PhyloSummary::PhyloSummary(string refTfile, GroupMap* g, bool r){
	try {
		m = MothurOut::getInstance();
		maxLevel = 0;
		ignore = false;
        numSeqs = 0;
        relabund = r;
		
		groupmap = g;
        ct = NULL;
				
		//check for necessary files
        if (refTfile == "saved") { ReferenceDB* rdb = ReferenceDB::getInstance(); refTfile = rdb->getSavedTaxonomy(); }
		string taxFileNameTest = m->getFullPathName((refTfile.substr(0,refTfile.find_last_of(".")+1) + "tree.sum"));
		ifstream FileTest(taxFileNameTest.c_str());
		
		if (!FileTest) { 
			m->mothurOut("Error: can't find " + taxFileNameTest + "."); m->mothurOutEndLine(); exit(1);
		}else{
			readTreeStruct(FileTest);
		}
		
		tree[0].rank = "0";
		assignRank(0);

	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "PhyloSummary");
		exit(1);
	}
}

/**************************************************************************************************/

PhyloSummary::PhyloSummary(GroupMap* g, bool r){
	try {
		m = MothurOut::getInstance();
		maxLevel = 0;
		ignore = true;
        numSeqs = 0;
        relabund = r;
		
		groupmap = g;
        ct = NULL;
		
		tree.push_back(rawTaxNode("Root"));
		tree[0].rank = "0";
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "PhyloSummary");
		exit(1);
	}
}
/**************************************************************************************************/

int PhyloSummary::summarize(string userTfile){
	try {
		map<string, string> temp;
        m->readTax(userTfile, temp);
        
        for (map<string, string>::iterator itTemp = temp.begin(); itTemp != temp.end();) {
            addSeqToTree(itTemp->first, itTemp->second);
            temp.erase(itTemp++);
        }
        
        return numSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "summarize");
		exit(1);
	}
}

/**************************************************************************************************/

string PhyloSummary::getNextTaxon(string& heirarchy){
	try {
		string currentLevel = "";
		if(heirarchy != ""){
			int pos = heirarchy.find_first_of(';');
			currentLevel=heirarchy.substr(0,pos);
			if (pos != (heirarchy.length()-1)) {  heirarchy=heirarchy.substr(pos+1);  }
			else { heirarchy = ""; }
		}
		return currentLevel;
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "getNextTaxon");
		exit(1);
	}
}

/**************************************************************************************************/

int PhyloSummary::addSeqToTree(string seqName, string seqTaxonomy){
	try {
				
		numSeqs++;
		
		map<string, int>::iterator childPointer;
		
		int currentNode = 0;
		string taxon;
		
		int level = 0;
		
		//are there confidence scores, if so remove them
		if (seqTaxonomy.find_first_of('(') != -1) {  m->removeConfidences(seqTaxonomy);	}
		
		while (seqTaxonomy != "") {
			
			if (m->control_pressed) { return 0; }
			
			//somehow the parent is getting one too many accnos
			//use print to reassign the taxa id
			taxon = getNextTaxon(seqTaxonomy);
			
			childPointer = tree[currentNode].children.find(taxon);
			
			if(childPointer != tree[currentNode].children.end()){	//if the node already exists, update count and move on
				int thisCount = 1;
                
                if (groupmap != NULL) {
					//find out the sequences group
					string group = groupmap->getGroup(seqName);
					
					if (group == "not found") {  m->mothurOut("[WARNING]: " + seqName + " is not in your groupfile, and will be included in the overall total, but not any group total."); m->mothurOutEndLine();  }
					
					//do you have a count for this group?
					map<string, int>::iterator itGroup = tree[childPointer->second].groupCount.find(group);
					
					//if yes, increment it - there should not be a case where we can't find it since we load group in read
					if (itGroup != tree[childPointer->second].groupCount.end()) {
						tree[childPointer->second].groupCount[group]++;
					}
				}else if (ct != NULL) {
                    if (ct->hasGroupInfo()) {
                        vector<int> groupCounts = ct->getGroupCounts(seqName);
                        vector<string> groups = ct->getNamesOfGroups();
                        for (int i = 0; i < groups.size(); i++) {
                            
                            if (groupCounts[i] != 0) {
                                //do you have a count for this group?
                                map<string, int>::iterator itGroup = tree[childPointer->second].groupCount.find(groups[i]);
                                
                                //if yes, increment it - there should not be a case where we can't find it since we load group in read
                                if (itGroup != tree[childPointer->second].groupCount.end()) {
                                    tree[childPointer->second].groupCount[groups[i]] += groupCounts[i];
                                }
                            }
                        }
                    }
                    thisCount = ct->getNumSeqs(seqName);
                }
				
				tree[childPointer->second].total += thisCount;

				currentNode = childPointer->second;
			}else{	
				if (ignore) {
						
					tree.push_back(rawTaxNode(taxon));
					int index = tree.size() - 1;
				
					tree[index].parent = currentNode;
					tree[index].level = (level+1);
					tree[currentNode].children[taxon] = index;
                    int thisCount = 1;
					
					//initialize groupcounts
					if (groupmap != NULL) {
						vector<string> mGroups = groupmap->getNamesOfGroups();
						for (int j = 0; j < mGroups.size(); j++) {
							tree[index].groupCount[mGroups[j]] = 0;
						}
						
						//find out the sequences group
						string group = groupmap->getGroup(seqName);
						
						if (group == "not found") {  m->mothurOut("[WARNING]: " + seqName + " is not in your groupfile, and will be included in the overall total, but not any group total."); m->mothurOutEndLine();  }
						
						//do you have a count for this group?
						map<string, int>::iterator itGroup = tree[index].groupCount.find(group);
						
						//if yes, increment it - there should not be a case where we can't find it since we load group in read
						if (itGroup != tree[index].groupCount.end()) {
							tree[index].groupCount[group]++;
						}
					}else if (ct != NULL) {
                        if (ct->hasGroupInfo()) {
                            vector<string> mGroups = ct->getNamesOfGroups();
                            for (int j = 0; j < mGroups.size(); j++) {
                                tree[index].groupCount[mGroups[j]] = 0;
                            }
                            vector<int> groupCounts = ct->getGroupCounts(seqName);
                            vector<string> groups = ct->getNamesOfGroups();
                        
                            for (int i = 0; i < groups.size(); i++) {
                                if (groupCounts[i] != 0) {
                                   
                                    //do you have a count for this group?
                                    map<string, int>::iterator itGroup = tree[index].groupCount.find(groups[i]);
                                     
                                    //if yes, increment it - there should not be a case where we can't find it since we load group in read
                                    if (itGroup != tree[index].groupCount.end()) {
                                        tree[index].groupCount[groups[i]]+=groupCounts[i];
                                    }
                                }
                            }
                        }
                        thisCount = ct->getNumSeqs(seqName);
                    }
					
                    tree[index].total = thisCount;
					currentNode = index;
					
				}else{ //otherwise, error
					m->mothurOut("Warning: cannot find taxon " + taxon + " in reference taxonomy tree at level " + toString(tree[currentNode].level) + " for " + seqName + ". This may cause totals of daughter levels not to add up in summary file."); m->mothurOutEndLine();
					break;
				}
			}
			
			level++;
			
			if ((seqTaxonomy == "") && (level < maxLevel)) {  //if you think you are done and you are not.
				for (int k = level; k < maxLevel; k++) {  seqTaxonomy += "unclassified;";   }
			}
		}
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "addSeqToTree");
		exit(1);
	}
}
/**************************************************************************************************/

int PhyloSummary::addSeqToTree(string seqTaxonomy, map<string, bool> containsGroup){
	try {
		numSeqs++;
		
		map<string, int>::iterator childPointer;
		
		int currentNode = 0;
		string taxon;
		
		int level = 0;
		
		//are there confidence scores, if so remove them
		if (seqTaxonomy.find_first_of('(') != -1) {  m->removeConfidences(seqTaxonomy);	}
		
		while (seqTaxonomy != "") {
			
			if (m->control_pressed) { return 0; }
			
			//somehow the parent is getting one too many accnos
			//use print to reassign the taxa id
			taxon = getNextTaxon(seqTaxonomy);
			
			childPointer = tree[currentNode].children.find(taxon);
			
			if(childPointer != tree[currentNode].children.end()){	//if the node already exists, update count and move on
                for (map<string, bool>::iterator itGroup = containsGroup.begin(); itGroup != containsGroup.end(); itGroup++) {
                    if (itGroup->second == true) {
                        tree[childPointer->second].groupCount[itGroup->first]++;
                    }
                }
					
				tree[childPointer->second].total++;
				
				currentNode = childPointer->second;
			}else{	
				if (ignore) {
					
					tree.push_back(rawTaxNode(taxon));
					int index = tree.size() - 1;
					
					tree[index].parent = currentNode;
					tree[index].level = (level+1);
					tree[index].total = 1;
					tree[currentNode].children[taxon] = index;
						
                    for (map<string, bool>::iterator itGroup = containsGroup.begin(); itGroup != containsGroup.end(); itGroup++) {
                        if (itGroup->second == true) {
                            tree[index].groupCount[itGroup->first]++;
                        }
                    }
					
					currentNode = index;
					
				}else{ //otherwise, error
					m->mothurOut("Warning: cannot find taxon " + taxon + " in reference taxonomy tree at level " + toString(tree[currentNode].level) + ". This may cause totals of daughter levels not to add up in summary file."); m->mothurOutEndLine();
					break;
				}
			}
			
			level++;
			
			if ((seqTaxonomy == "") && (level < maxLevel)) {  //if you think you are done and you are not.
				for (int k = level; k < maxLevel; k++) {  seqTaxonomy += "unclassified;";   }
			}
		}
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "addSeqToTree");
		exit(1);
	}
}

/**************************************************************************************************/

void PhyloSummary::assignRank(int index){
	try {
		map<string,int>::iterator it;
		int counter = 1;
        
		for(it=tree[index].children.begin();it!=tree[index].children.end();it++){
			tree[it->second].rank = tree[index].rank + '.' + toString(counter);
			counter++;
			assignRank(it->second);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "assignRank");
		exit(1);
	}
}
/**************************************************************************************************/

void PhyloSummary::print(ofstream& out){
	try {
		
		if (ignore)     {  assignRank(0); }
        vector<string> mGroups;
		//print labels
		out << "taxlevel\t rankID\t taxon\t daughterlevels\t total\t";
		if (groupmap != NULL) {
			//so the labels match the counts below, since the map sorts them automatically...
			//sort(groupmap->namesOfGroups.begin(), groupmap->namesOfGroups.end());
            mGroups = groupmap->getNamesOfGroups();
			for (int i = 0; i < mGroups.size(); i++) {
				out << mGroups[i] << '\t';
			}
		}else if (ct != NULL) {
            if (ct->hasGroupInfo()) {
                mGroups = ct->getNamesOfGroups();
                for (int i = 0; i < mGroups.size(); i++) {
                    out << mGroups[i] << '\t';
                }
            }
        }
		out << endl;
		
		int totalChildrenInTree = 0;
		map<string, int>::iterator itGroup;
		map<string,int>::iterator it;
		for(it=tree[0].children.begin();it!=tree[0].children.end();it++){   
			if (tree[it->second].total != 0)  {   
				totalChildrenInTree++; 
				tree[0].total += tree[it->second].total;
				
				if (groupmap != NULL) {
					for (int i = 0; i < mGroups.size(); i++) { tree[0].groupCount[mGroups[i]] += tree[it->second].groupCount[mGroups[i]]; } 
				}else if ( ct != NULL) {
                    if (ct->hasGroupInfo()) { for (int i = 0; i < mGroups.size(); i++) { tree[0].groupCount[mGroups[i]] += tree[it->second].groupCount[mGroups[i]]; } }
                }
			}
		}
		
		//print root
        if (relabund) {
            out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
            out << tree[0].level << "\t" << tree[0].rank << "\t" << tree[0].name << "\t" << totalChildrenInTree << "\t" << (tree[0].total/(double) tree[0].total) << '\t';
            
            if (groupmap != NULL) {
                for (int i = 0; i < mGroups.size(); i++) {
                    double thisNum = tree[0].groupCount[mGroups[i]];
                    thisNum /= (double) groupmap->getNumSeqs(mGroups[i]);
                    out << thisNum << '\t';
                }
            }else if ( ct != NULL) {
                if (ct->hasGroupInfo()) {
                    for (int i = 0; i < mGroups.size(); i++) {
                        double thisNum = tree[0].groupCount[mGroups[i]];
                        thisNum /= (double) ct->getGroupCount(mGroups[i]);
                        out << thisNum << '\t';
                    }
                }
            }
            out << endl;
        }else {
            out << tree[0].level << "\t" << tree[0].rank << "\t" << tree[0].name << "\t" << totalChildrenInTree << "\t" << tree[0].total << "\t";
            if (groupmap != NULL) {
                for (int i = 0; i < mGroups.size(); i++) {  out << tree[0].groupCount[mGroups[i]] << '\t'; }
            }else if ( ct != NULL) {
                if (ct->hasGroupInfo()) { for (int i = 0; i < mGroups.size(); i++) {  out << tree[0].groupCount[mGroups[i]] << '\t'; } }
            }
            out << endl;
		}
		
		//print rest
		print(0, out);
		
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "print");
		exit(1);
	}
}
/**************************************************************************************************/

void PhyloSummary::print(ofstream& out, bool relabund){
	try {
		
		if (ignore) { assignRank(0); }
	
		int totalChildrenInTree = 0;
		map<string, int>::iterator itGroup;
		
		map<string,int>::iterator it;
		for(it=tree[0].children.begin();it!=tree[0].children.end();it++){
			if (tree[it->second].total != 0)  {
				totalChildrenInTree++;
				tree[0].total += tree[it->second].total;
				
				if (groupmap != NULL) {
                    vector<string> mGroups = groupmap->getNamesOfGroups();
					for (int i = 0; i < mGroups.size(); i++) { tree[0].groupCount[mGroups[i]] += tree[it->second].groupCount[mGroups[i]]; }
				}else if ( ct != NULL) {
                    vector<string> mGroups = ct->getNamesOfGroups();
                    if (ct->hasGroupInfo()) { for (int i = 0; i < mGroups.size(); i++) { tree[0].groupCount[mGroups[i]] += tree[it->second].groupCount[mGroups[i]]; } }
                }
			}
		}
		
		//print root
		out << tree[0].name << "\t" << "1.0000" << "\t"; //root relative abundance is 1, everyone classifies to root
		
		/*
		if (groupmap != NULL) {
			for (int i = 0; i < mGroups.size(); i++) {  out << tree[0].groupCount[mGroups[i]] << '\t'; }
        }else if ( ct != NULL) {
            if (ct->hasGroupInfo()) { for (int i = 0; i < mGroups.size(); i++) {  out << tree[0].groupCount[mGroups[i]] << '\t'; } }
        }*/
        
        if (groupmap != NULL) {
            vector<string> mGroups = groupmap->getNamesOfGroups();
			for (int i = 0; i < mGroups.size(); i++) {  out << "1.0000" << '\t'; }
        }else if ( ct != NULL) {
            vector<string> mGroups = ct->getNamesOfGroups();
            if (ct->hasGroupInfo()) { for (int i = 0; i < mGroups.size(); i++) {  out << "1.0000" << '\t'; } }
        }
        
		out << endl;
		
		//print rest
		print(0, out, relabund);
		
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "print");
		exit(1);
	}
}
/**************************************************************************************************/

void PhyloSummary::print(int i, ofstream& out){
	try {
		map<string,int>::iterator it;
		for(it=tree[i].children.begin();it!=tree[i].children.end();it++){
			
			if (tree[it->second].total != 0)  {
                
				int totalChildrenInTree = 0;
                
				map<string,int>::iterator it2;
				for(it2=tree[it->second].children.begin();it2!=tree[it->second].children.end();it2++){
					if (tree[it2->second].total != 0)  {   totalChildrenInTree++; }
				}
                
                if (relabund) {
                    out << tree[it->second].level << "\t" << tree[it->second].rank << "\t" << tree[it->second].name << "\t" << totalChildrenInTree << "\t" << (tree[it->second].total/(double) tree[0].total) << "\t";
                }else {
                    out << tree[it->second].level << "\t" << tree[it->second].rank << "\t" << tree[it->second].name << "\t" << totalChildrenInTree << "\t" << tree[it->second].total << "\t";
				}
                
                if (relabund) {
                    map<string, int>::iterator itGroup;
                    if (groupmap != NULL) {
                        vector<string> mGroups = groupmap->getNamesOfGroups();
                        for (int i = 0; i < mGroups.size(); i++) {  out << (tree[it->second].groupCount[mGroups[i]]/(double)groupmap->getNumSeqs(mGroups[i])) << '\t'; }
                    }else if (ct != NULL) {
                        if (ct->hasGroupInfo()) {
                            vector<string> mGroups = ct->getNamesOfGroups();
                            for (int i = 0; i < mGroups.size(); i++) {  out << (tree[it->second].groupCount[mGroups[i]]/(double)ct->getGroupCount(mGroups[i])) << '\t'; }
                        }
                    }
                }else {
                    map<string, int>::iterator itGroup;
                    if (groupmap != NULL) {
                        vector<string> mGroups = groupmap->getNamesOfGroups();
                        for (int i = 0; i < mGroups.size(); i++) {  out << tree[it->second].groupCount[mGroups[i]] << '\t'; }
                    }else if (ct != NULL) {
                        if (ct->hasGroupInfo()) {
                            vector<string> mGroups = ct->getNamesOfGroups();
                            for (int i = 0; i < mGroups.size(); i++) {  out << tree[it->second].groupCount[mGroups[i]] << '\t'; }
                        }
                    }
 
                }
                out << endl;
				
			}
			
			print(it->second, out);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "print");
		exit(1);
	}
}

/**************************************************************************************************/

void PhyloSummary::print(int i, ofstream& out, bool relabund){
	try {
		map<string,int>::iterator it;
		for(it=tree[i].children.begin();it!=tree[i].children.end();it++){
			
			if (tree[it->second].total != 0)  {
			
				int totalChildrenInTree = 0;
		
				map<string,int>::iterator it2;
				for(it2=tree[it->second].children.begin();it2!=tree[it->second].children.end();it2++){   
					if (tree[it2->second].total != 0)  {   totalChildrenInTree++; }
				}
                
                string nodeName = "";
                int thisNode = it->second;
                while (tree[thisNode].rank != "0") { //while you are not at top
                    if (m->control_pressed) { break; }
                    nodeName = tree[thisNode].name + "|" + nodeName;
                    thisNode = tree[thisNode].parent;
                }
                if (nodeName != "") { nodeName = nodeName.substr(0, nodeName.length()-1); }
                
				out << nodeName << "\t" << (tree[it->second].total / (float)tree[i].total) << "\t";
				
				map<string, int>::iterator itGroup;
				if (groupmap != NULL) {
					vector<string> mGroups = groupmap->getNamesOfGroups();
					for (int j = 0; j < mGroups.size(); j++) {
                        if (tree[i].groupCount[mGroups[j]] == 0) {
                            out << 0 << '\t';
                        }else { out << (tree[it->second].groupCount[mGroups[j]] / (float)tree[i].groupCount[mGroups[j]]) << '\t'; }
                    }
				}else if (ct != NULL) {
                    if (ct->hasGroupInfo()) {
                        vector<string> mGroups = ct->getNamesOfGroups();
                        for (int j = 0; j < mGroups.size(); j++) {
                            if (tree[i].groupCount[mGroups[j]] == 0) {
                                out << 0 << '\t';
                            }else { out << (tree[it->second].groupCount[mGroups[j]] / (float)tree[i].groupCount[mGroups[j]]) << '\t'; }
                        }
                    }
                }
				out << endl;
				
			}
			
			print(it->second, out, relabund);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "print");
		exit(1);
	}
}
/**************************************************************************************************/
void PhyloSummary::readTreeStruct(ifstream& in){
	try {
	
		//read version
		string line = m->getline(in); m->gobble(in);
		
		int num;
		
		in >> num; m->gobble(in);
		
		tree.resize(num);
		
		in >> maxLevel; m->gobble(in);
	
		//read the tree file
		for (int i = 0; i < tree.size(); i++) {
	
			in >> tree[i].level >> tree[i].name >> num; //num contains the number of children tree[i] has
			
			//set children
			string childName;
			int childIndex;
			for (int j = 0; j < num; j++) {
				in >> childName >> childIndex;
				tree[i].children[childName] = childIndex;
			}
			
			//initialize groupcounts
			if (groupmap != NULL) {
				for (int j = 0; j < (groupmap->getNamesOfGroups()).size(); j++) {
					tree[i].groupCount[(groupmap->getNamesOfGroups())[j]] = 0;
				}
			}else if (ct != NULL) {
                if (ct->hasGroupInfo()) {
                    for (int j = 0; j < (ct->getNamesOfGroups()).size(); j++) {
                        tree[i].groupCount[(ct->getNamesOfGroups())[j]] = 0;
                    }
                }
            }
			
			tree[i].total = 0;
			
			m->gobble(in);
			
			//if (tree[i].level > maxLevel) {  maxLevel = tree[i].level;  }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "readTreeStruct");
		exit(1);
	}
}
/**************************************************************************************************/


	
