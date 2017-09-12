import numpy as np
import pandas as pd
import re
import glob
import shutil
import multiprocessing
import subprocess
import time
import os

import Bio.PDB
import Bio.SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqUtils import seq3

# launch pymol

import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
from time import sleep
import pymol
from pymol import cmd

pymol.finish_launching()

nproc = 2

# paths

hla_fasta = '../data/mhc_alleles/hla/hla_prot.fasta'
corrected_pdb_path = '../sandbox/corrected_pdb'
output_dir = '../sandbox/mutated'

# tables

annotation = pd.read_csv('../sandbox/detailed_mhc', sep='\t')
alleles_vs_pdb = pd.read_csv('../sandbox/alleles_vs_pdb', sep='\t', index_col=0, low_memory=False)
peptides_vs_pdb = pd.read_csv('../sandbox/peptides_vs_pdb', sep='\t', index_col=0, low_memory=False)
ligand_assays = pd.read_csv('../data/epitopes/mhc_ligand_full.csv', skiprows=1, low_memory=False)

# table conversion

ligand_assays = ligand_assays[ligand_assays['MHC allele class'] == 'I']
affinity = ligand_assays[['Description', 'Allele Name', 'Qualitative Measure', 'Quantitative measurement']]
affinity['Length'] = [len(x) for x in affinity.Description]
affinity = affinity[(affinity.Length <= 13) & (affinity.Length >= 8)]
affinity = affinity[affinity['Allele Name'].str.contains('HLA-[ABC]\*\S+$', regex=True)]
affinity['Allele Name'] = [x[4:] for x in affinity['Allele Name']]
affinity.drop_duplicates(inplace = True)
affinity = affinity.reset_index(drop=True)

iclass = annotation[annotation.mhc_class == 'I']
iclass = iclass[iclass.a_chain_allele.str.contains('[ABC]\*', regex=True)]
iclass.a_chain_allele = [':'.join(x.split(':')[:2]) for x in iclass.a_chain_allele]
iclass = iclass.reset_index(drop=True)

alleledict = Bio.SeqIO.to_dict(Bio.SeqIO.parse(hla_fasta, 'fasta'))
for key in alleledict.keys():
    alleledict[key.split(';')[1]] = alleledict.pop(key)


def ungapped_align(a, b, matrix, peptide=False):
    score = -np.inf
    pos = []
    
    if not peptide:
        for i in range(len(a) + len(b)):
            a_beg, a_end = (max(0, i-len(b)), min(i, len(a)))
            b_beg, b_end = (max(0, len(b)-i), min(len(b), len(b)+len(a)-i))
            
            loc_score = sum([matrix[(k, t)] for k, t in zip(a[a_beg:a_end], b[b_beg:b_end])])
            if loc_score >= score:
                score = loc_score
                pos = [a_beg, b_beg]
            else:
                pos = [0, 0]
    
    a_aligned = a
    b_aligned = b
    
    if (pos[0] > pos[1]):
        b_aligned = '-'*(pos[0]-pos[1])+b_aligned
    else:
        a_aligned = '-'*(pos[1]-pos[0])+a_aligned
        
    if (len(a_aligned) > len(b_aligned)):
        b_aligned += '-'*(len(a_aligned) - len(b_aligned))
    else:
        a_aligned += '-'*(len(b_aligned) - len(a_aligned))
        
    return a_aligned, b_aligned

# s1 - mutated sequence
# s2 - reference sequence
def align(s1, s2, peptide=False):
    matrix = {}
    for key, val in matlist.blosum62.iteritems(): # do it outside the function
        matrix.update({(key[1], key[0]): val})
        matrix.update({(key[0], key[1]): val})
        
    top_aln = ungapped_align(s1, s2, matrix, peptide)
    
    print("\tOld: "+top_aln[0])
    print("\tNew: "+top_aln[1])
    
    # check that there is no gaps inside the alignment
    assert(re.search('[A-Z]-+[A-Z]', top_aln[0]) == None)
    assert(re.search('[A-Z]-+[A-Z]', top_aln[1]) == None)
    
    mutations = []
    gap = re.match('(-+)', top_aln[0])
    shift = 0
    if gap:
        shift = gap.end(1) - gap.start(1)
        
    for i, a1, a2 in zip(list(range(-shift, len(top_aln[0])-shift)), top_aln[0], top_aln[1]):
        if a1 != '-' and a2 != '-':
            if a1 != a2:
                mutations.append((i, a1, a2))
    return mutations

# mutate pdb structure in accordance to input sequences
# here we assume that all chains are continuous: CHECK THIS
def mutate(pdb, newseq, newpep, table):
    
    start_time = time.time()
    
    parser = Bio.PDB.PDBParser(QUIET=True)
    ppb=Bio.PDB.PPBuilder()
    
    name = pdb;
    amhc = table[table.pdb == pdb].iloc[0, 3]
    pept = table[table.pdb == pdb].iloc[0, 7]
    b2m  = table[table.pdb == pdb].iloc[0, 5]
    
    struct = parser.get_structure(id=name, file=corrected_pdb_path+'/'+name+'.pdb')[0]
    oldseq = ''.join([str(x.get_sequence()) for x in ppb.build_peptides(struct[amhc])])
    oldpep = table[table.pdb == pdb].iloc[0, 8]
    
    #print("\n0 time: " + str(time.time() - start_time)); start_time = time.time()
    
    assert(len(oldpep) == len(newpep))
    print("\n======= MHC business =======\n")
    mhc_mutations = align(oldseq, newseq)
    print("\nMHC mutations: %i\n" % len(mhc_mutations))
    
    #print("\n1 time: " + str(time.time() - start_time)); start_time = time.time()
    
    print("\n===== Peptide business =====\n")
    pep_mutations = align(oldpep, newpep, True)
    print("\nPeptide mutations: %i\n\n" % len(pep_mutations))
    
    #print("2 time: " + str(time.time() - start_time)); start_time = time.time()
    
    print("Mutagenesis: \n")
    # remove alternative atom locations
    cmd.remove("not alt ''+A")
    cmd.alter('all', "alt=''")
    
    # mutate
    cmd.wizard('mutagenesis')
    cmd.refresh_wizard()
    
    first_aa_id = ppb.build_peptides(struct[amhc])[0][0].get_id()[1]
    for i, old, new in mhc_mutations:
        print(i, old, new)
        cmd.get_wizard().do_select(amhc+'/'+str(i+first_aa_id)+'/')
        cmd.get_wizard().set_mode(seq3(new).upper())
        cmd.get_wizard().apply()
        
    first_aa_id = ppb.build_peptides(struct[pept])[0][0].get_id()[1]
    for i, old, new in pep_mutations:
        print(i, old, new)
        cmd.get_wizard().do_select(pept+'/'+str(i+first_aa_id)+'/')
        cmd.get_wizard().set_mode(seq3(new).upper())
        cmd.get_wizard().apply()
        
    cmd.set_wizard()
    
    return (len(mhc_mutations), len(pep_mutations))

def find_and_mutate(outpath, pepseq, allele, alleledict, pep2pdb, mhc2pdb, pdb_table, first_best=10):
    if pepseq not in pep2pdb.index:
        print(pepseq + " is not found in the table. Skipping.")
        return
        
    chunk = mhc2pdb.loc[mhc2pdb.index.str.startswith(allele), :]
    best_pdbs = (chunk.apply(min) + (pep2pdb.loc[pepseq, :] / 10.0)).sort_values().dropna()
    
    report = open(outpath + '/report', 'w')
    report.write('name\tpdb\ta_chain\told_allele\tnew_allele\tmhc_mutations\
    \tpeptide_chain\told_peptide\tnew_peptide\tpeptide_mutations\tcompatibility_score\n')
    
    cmd.reinitialize()
    
    counter = 0
    for i in xrange(len(best_pdbs)):
        try:
            print("\n=============================================")
            
            best_pdb   = best_pdbs.index[i]
            comp_score = best_pdbs[best_pdb]
            
            print("Operating on " + best_pdb + "\n")
            
            allele_list = chunk.loc[:, best_pdb]
            full_allele = allele_list[allele_list == min(allele_list)].index[0]
            print("New allele: "+full_allele)
            
            alleleseq = str(alleledict[full_allele].seq)
            print('Compability score: %f' % comp_score)
            
            cmd.load(corrected_pdb_path+'/'+best_pdb+'.pdb', best_pdb)
            
            mhc_mut_num, pep_mut_num = mutate(best_pdb, alleleseq, pepseq, pdb_table)
            line = pdb_table[pdb_table.pdb == best_pdb].iloc[0, :]
            
            cmd.set_name(best_pdb, "obj")
            
            mut_path = outpath + '/mutated%i.pdb' % counter
            hla_chain = line["a_chain_id"]
            b2m_chain = line["b_chain_id"]
            pep_chain = line["antigen_id"]
            
            # manipulations with segi are needed because some
            # chains are initially called A,B or C and it causes a bug
            cmd.alter('obj', 'segi="XZ"')
            cmd.alter('obj/XZ/%s//' % hla_chain, 'chain="A"; segi="XY"')
            cmd.alter('obj/XZ/%s//' % b2m_chain, 'chain="B"; segi="XY"')
            cmd.alter('obj/XZ/%s//' % pep_chain, 'chain="C"; segi="XY"')
            cmd.alter('obj', 'segi=""')
            cmd.save(mut_path, "obj")
            
            cmd.delete("all")
            
            subprocess.call(["pdbprep.pl", mut_path])
            subprocess.call(["pdbnmd.pl",  mut_path, "?"])
            mut_path_nmin = outpath + "/mutated%i_nmin.pdb" % counter
            
            cmd.set("pdb_retain_ids", 1)
            cmd.load(mut_path_nmin, "obj")
            
            cmd.select("min1", "obj//A+B//")
            cmd.select("chunk", "br. (obj//A+B// nto. 10 of obj//C//)")
            cmd.select("min2", "(obj//A+B// and (not chunk))")
            cmd.select("anchor_ca", "chunk and name CA")
            
            cmd.save("../sandbox/mutated/fixed%i_min1.pdb" % counter, "min1")
            cmd.save("../sandbox/mutated/fixed%i_min2.pdb" % counter, "min2")
            
            resi = []
            space = {'ar' : resi}
            cmd.iterate_state(0, "obj//C//CA", "ar.append(resi)", space=space)
            
            cmd.save("../sandbox/mutated/anchor%i_1.pdb" % counter, "obj//C/%s/CA" % resi[0])
            cmd.save("../sandbox/mutated/anchor%i_2.pdb" % counter, "obj//C/%s/CA" % resi[-1])
            
            ind = []
            space = {'ar' : ind}
            cmd.iterate_state(0, "anchor_ca", "ar.append(index)", space=space)
            for ix in ind:
                cmd.save("../sandbox/mutated/anchorCA%i_%s.pdb" % (counter, ix), "id %s" % ix)
            
            cmd.set("pdb_retain_ids", 0)
            cmd.delete("all")
            
            report.write('%s\t%s\t%s\t%s\t%s\t%i\t%s\t%s\t%s\t%i\t%.1f\n' %
                         ('mutated%i.pdb' % counter, best_pdb, line.a_chain_id, line.a_chain_allele,
                          allele, mhc_mut_num, line.antigen_id, line.antigen_seq, pepseq,
                          pep_mut_num, comp_score))
                          
        except KeyError:
            print('Key error: %s, %s, %s' % (pepseq, allele, best_pdb))
            continue
        
        except AssertionError:
            print('Assertion error: %s, %s, %s' % (pepseq, allele, best_pdb))
            continue
        
        else:
            if (counter + 1) >= first_best:
                break
            else:
                counter += 1
         
    report.close()
    return

def worker(pnum, pid):
    print("Worker %i launched\n" % pid)
    
    stdoutcopy = sys.stdout
    cmd.reinitialize()
    
    affinity_loc = affinity[affinity.index % pnum == pid]
    affinity_loc = affinity.iloc[:2, :]
    
    for ind, row in affinity_loc.iterrows():
        print('Worker %i working on %06i:' % (pid, ind))
        dirname = "%06i" % ind
        path = output_dir + "/" + dirname
        
        if os.path.exists(path): shutil.rmtree(path)
        os.mkdir(path)
        
        sys.stdout = open(path + "/log", "w")
        find_and_mutate(path, row["Description"], row["Allele Name"], alleledict, peptides_vs_pdb, alleles_vs_pdb, iclass, first_best=5)
        
        sys.stdout.close()
    
    sys.stdout = stdoutcopy
    
    print("Worker %i done!" % pid)
    return

if __name__ == '__main__':
    jobs = []
    for i in range(nproc):
        p = multiprocessing.Process(target=worker, args=(4,i,))
        jobs.append(p)
        p.start()
        p.join()
