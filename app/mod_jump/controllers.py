import pandas as pd
import numpy as np
import os
import shutil
from werkzeug.utils import secure_filename
from zipfile import ZipFile

from Bio.Restriction import BsaI, BsmBI
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt
from primers import primers
#import primer3 as pr

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(APP_ROOT, '../resources/')

### FILE-RELATED FUNCTIONS ###

def handle_file(input_file, format='tabular'):

	valid_tab = ['csv']
	valid_seq = ['fasta', 'zip']
	filename = input_file.filename
	ext = filename.rsplit('.', 1)[1].lower()

	if '.' in filename and ((format=='tabular' and ext in valid_tab) or (format=='sequence' and ext in valid_seq)):
		path_to_file = os.path.join(RESOURCES, secure_filename(filename))
		input_file.save(path_to_file)
		if ext == 'zip':
			with ZipFile(path_to_file, 'r') as zip_obj:
				plasmid_files = zip_obj.namelist()
				zip_obj.extractall(RESOURCES)
			return plasmid_files
		else:
			return [filename]
	return None

def write_fasta(assembly_result):

	shutil.rmtree(RESOURCES + 'fasta')
	os.mkdir(RESOURCES + 'fasta')
	for i, row in assembly_result.iterrows():
		record = SeqRecord(Seq(row['sequence']), id=row['name'], description='')
		with open(RESOURCES + 'fasta/{}.fasta'.format(row['name']), 'w') as handle:
			SeqIO.write(record, handle, 'fasta')


### DOMESTICATE PARTS ###

def find_templates(parts, plasmids):

	result = pd.DataFrame([(part[0], plasmid[0]) for _, part in parts.iterrows() \
                      for _, plasmid in plasmids.iterrows() if plasmid[3].find(part[3])!=-1], \
                      columns=['name', 'template_id']).groupby('name')['template_id'] \
                      .apply(lambda x: ', '.join(x)).reset_index()

	missing = list(set(parts['name'].tolist()).difference(set(result['name'].tolist())))

	return result, missing

def check_internal_sites(parts_ori, missing, forbidden_sites):

	parts = parts_ori.copy()
	#forbidden_sites = np.array([(enz, str(Seq(enz).reverse_complement())) \
    #                        for enz in [BsaI.site, BsmBI.site]]).ravel().tolist()

	#check any forbidden sites
	parts['num_forbidden_sites'] = parts['sequence'].apply(lambda x: np.sum([x.count(a) for a in forbidden_sites]))

	invalid_parts = parts[parts['num_forbidden_sites']>0]
	removed = invalid_parts['name'].tolist() + missing
	valid_parts = parts[~parts['name'].isin(removed)].reset_index(drop=True)

	return valid_parts, invalid_parts

def design_fragments(valid_parts, prefix, suffix):

	fragments = pd.merge(valid_parts, pd.read_csv(RESOURCES + '/overhang.csv'), on='sites', how='left') \
                [['name', 'sequence', 'left_site', 'right_site']]
	#forbidden_sites = np.array([(enz, str(Seq(enz).reverse_complement())) \
    #                        for enz in [BsaI.site, BsmBI.site]]).ravel().tolist()
	#prefix, suffix = forbidden_sites[2] + forbidden_sites[0] + 'A', 'T' + forbidden_sites[1] + forbidden_sites[3]
	fragments['left_overhang'] = prefix + fragments['left_site']
	fragments['right_overhang'] = fragments['right_site'] + suffix
	fragments['bases'] = fragments['left_overhang'] + fragments['sequence'] + fragments['right_overhang']
	#following benchling convention to start from index 1
	fragments['size'] = fragments['bases'].apply(lambda x: len(x)+1)

	return fragments

def design_primers(parts, prefix, suffix):
    
    primers_list = []
    #counter = 0
    for _, part in parts.iterrows():
        
        try:
            fp, rp = primers(part['sequence'], add_fwd=prefix + part['left_site'],
                             add_rev=str(Seq(suffix).reverse_complement()) + str(Seq(part['right_site']).reverse_complement()))
            
            #forward primers
            primers_list.append((#'(P{}J-RM)_{}.F'.format(str(counter+start).zfill(3), part['name']), \
                                 part['name'], fp.seq, 
                                 fp.tm, fp.tm_total, fp.gc, fp.dg, fp.fwd, fp.offtargets, fp.penalty))
            #counter += 1
            #reverse primers
            primers_list.append((#'(P{}J-RM)_{}.R'.format(str(counter+start).zfill(3), part['name']), \
                                 part['name'], rp.seq,
                                 rp.tm, rp.tm_total, rp.gc, rp.dg, rp.fwd, rp.offtargets, rp.penalty))
            #counter += 1
        
        except:
            print('Error on', part['name'], ': cannot find any feasible primer design, check your fragments!')
            continue
        
    raw_primers = pd.DataFrame(primers_list, columns=['part_name', 'sequence', 'tm', 'tm_total', 'gc', 'dg', 'fwd', 'offtargets', 'penalty'])
    
    start = 1
    primer_seq = raw_primers[['fwd', 'sequence']].drop_duplicates().reset_index(drop=True)
    primer_seq['primer_name'] = pd.Series(primer_seq.index).apply(lambda x: 'P' + str(x+start).zfill(3) + 'J')
    primer_seq.loc[primer_seq['fwd'], 'primer_name'] = primer_seq['primer_name'] + '.F'
    primer_seq.loc[~primer_seq['fwd'], 'primer_name'] = primer_seq['primer_name'] + '.R'
    final_primers = pd.merge(primer_seq[['primer_name', 'sequence']], raw_primers.drop('part_name', axis=1), \
							on='sequence', how='left').drop_duplicates().reset_index(drop=True)
    
    return final_primers, raw_primers, primer_seq

def get_annealed_part(a, b):

    return [a[i].replace(b[i], '') for i in range(len(a))]

def calculate_tm(seq):

	return np.median([mt.Tm_Wallace(seq), mt.Tm_GC(seq), mt.Tm_NN(seq)])

    #return pr.calcTm(seq, dna_conc=250, dntp_conc=10,
    #                  tm_method='santalucia', salt_corrections_method='santalucia')

def generate_pcr(fragments, raw_primers, primer_seq):

	pcr_primers = pd.merge(pd.merge(primer_seq, raw_primers[['part_name', 'sequence']], on='sequence', how='right'), \
                fragments[['name', 'bases', 'left_overhang', 'right_overhang', 'size']], \
                left_on='part_name', right_on='name', how='left')
	pcr_primers.loc[~pcr_primers['fwd'], 'sequence'] = pcr_primers['sequence']. \
				apply(lambda x: str(Seq(x).reverse_complement()))
	pcr_primers.loc[pcr_primers['fwd'], 'overhang_part'] = pcr_primers['left_overhang']
	pcr_primers.loc[~pcr_primers['fwd'], 'overhang_part'] = pcr_primers['right_overhang']
	pcr_primers['anneal_part'] = get_annealed_part(pcr_primers['sequence'], pcr_primers['overhang_part'])
	pcr_primers['tm'] = np.round(pcr_primers['anneal_part'].apply(calculate_tm), 1)
	pcr_primers['tm_product'] = np.round(pcr_primers['bases'].apply(calculate_tm), 1)
	pcr_primers['ext_time'] = pcr_primers['size'].apply(lambda x: max([60, x/1000 * 60]))
	fwd_primers = pcr_primers[pcr_primers['fwd']] \
				[['part_name', 'primer_name', 'tm', 'tm_product', 'size', 'ext_time']]
	rev_primers = pcr_primers[~pcr_primers['fwd']] \
				[['part_name', 'primer_name', 'tm', 'tm_product', 'size', 'ext_time']]

	pcr = pd.merge(fwd_primers, rev_primers, on=['part_name', 'tm_product', 'size', 'ext_time'], how='outer')
	pcr['ta'] = 0.3 * np.max([pcr['tm_x'], pcr['tm_y']], axis=0) + 0.7 * pcr['tm_product'] - 14.9
	pcr['primer_name_x'] = pcr['primer_name_x'].str.split('_', expand=True)[0] \
						.str.replace('(', '').str.replace(')', '')
	pcr['primer_name_y'] = pcr['primer_name_y'].str.split('_', expand=True)[0] \
						.str.replace('(', '').str.replace(')', '')
	pcr_rxn = pcr[['part_name', 'size', 'primer_name_x', 'primer_name_y', 'ta', 'ext_time']] \
						.sort_values('ta').reset_index(drop=True)
	return pcr_rxn

def domesticate_parts(parts_file, maps, raw_fastas):

	fastas = [f for f in raw_fastas if not f.startswith('__')]

	### plasmids ###
	plasmids = pd.DataFrame([(p.id, str(p.seq)) for fasta in fastas \
							for p in list(SeqIO.parse(RESOURCES + '{}'.format(fasta), 'fasta'))], \
							columns=['name', 'sequence'])
	plasmids['sequence'] = plasmids['sequence'].str.upper()

	if maps:
		plasmids = pd.merge(pd.read_csv(RESOURCES + '{}'.format(maps[0])), plasmids, \
                    left_on='full_name', right_on='name', how='right')[['id', 'short_name', 'name', 'sequence']]
	else:
		plasmids['id'] = ['p{}'.format(i) for i in np.arange(len(plasmids))]
		plasmids['short_name'] = plasmids['name']
		plasmids = plasmids[['id', 'short_name', 'name', 'sequence']]

	### parts ###
	parts = pd.read_csv(RESOURCES + '{}'.format(parts_file[0])) #always at index 0
	parts['sequence'] = parts['sequence'].str.upper()

	forbidden_sites = np.array([(enz, str(Seq(enz).reverse_complement())) \
                            for enz in [BsaI.site, BsmBI.site]]).ravel().tolist()
	prefix, suffix = forbidden_sites[2] + forbidden_sites[0] + 'A', 'T' + forbidden_sites[1] + forbidden_sites[3]

	result, missing = find_templates(parts, plasmids)
	valid_parts, invalid_parts = check_internal_sites(parts, missing, forbidden_sites)
	fragments = design_fragments(valid_parts, prefix, suffix)

	final_primers, raw_primers, primer_seq = design_primers(fragments, prefix, suffix)
	final_primers.to_csv(RESOURCES + 'primers.csv', index=False)

	pcr_rxn = generate_pcr(fragments, raw_primers, primer_seq)
	pcr_rxn.to_csv(RESOURCES + 'pcr_rxn.csv', index=False)

	return 'OK'

### SIMULATE ASSEMBLY ###

def reindex_ps1(plasmid):
    '''Reindex plasmids to start from the annealing region of PS1 so the fragment will be in the middle of the sequence.
    This will do nothing for a fragment from pcr products as they do not have PS1 region'''
    
    new_start = plasmid.find('AGGGCGGCGGATTTGTCC')
    return plasmid[new_start:] + plasmid[:new_start]

def len_amplicon(plasmid):
    '''Calculate the length of amplicon from PS1 to PS2, the plasmid always starts from PS1'''
    
    ps2 = 'GAACGCTCGGTTGCCGC' #reverse complement of PS2
    return len(plasmid[:plasmid.find(ps2)]) + len(ps2)

def get_sites(part, odd_level=True, vector=False):
    
    enz = BsaI if odd_level else BsmBI
    
    if vector:
        
        site = str(Seq(enz.site).reverse_complement())
        right_idx = part.find(site) - 5
        right_cut = part[right_idx: right_idx+4]

        site = enz.site
        left_idx = part.find(site) + len(site) + 1
        left_cut = part[left_idx: left_idx+4]
        
        fragment = part[left_idx:] + part[:right_idx+4]
        
    else:
        site = enz.site
        left_idx = part.find(site) + len(site) + 1
        left_cut = part[left_idx: left_idx+4]
        
        site = str(Seq(enz.site).reverse_complement())
        right_idx = part.find(site) - 5
        right_cut = part[right_idx: right_idx+4]
        
        fragment = part[left_idx: right_idx+4]
        
    return left_cut, right_cut, fragment

def assemble_lvl_0(fragments, uac):
    
    constructs = []
    for name in fragments:

        uac_sites = get_sites(uac, odd_level=False, vector=True)
        frag_sites = get_sites(fragments[name], odd_level=False, vector=False)
        
        if (uac_sites[1]==frag_sites[0] and uac_sites[0]==frag_sites[1]):
            assembly = reindex_ps1(frag_sites[2][:-4] + uac_sites[2][:-4])
            constructs.append((name, assembly))
            
    return pd.DataFrame(constructs, columns=['name', 'sequence'])

def assemble(assembly_plan, mapping, odd):

    constructs = []
    for i, entry in assembly_plan.iterrows():

        #print(entry)
        _id = entry[0]
        name = entry[1]
        vector = entry[-1]
        parts = entry[2:-1].tolist()

        fragments = []
        fragments.append(get_sites(mapping[vector], odd_level=odd, vector=True))
        for part in parts:
            fragments.append(get_sites(mapping[part], odd_level=odd, vector=False))
            
        sites = list(map(list, zip(*[fragment[:2] for fragment in fragments])))
        sites[0] = sites[0][1:] + [sites[0][0]]
        
        if (sites[0]==sites[1]):
            assembly = reindex_ps1(''.join([fragment[2][:-4] for fragment in fragments]))
            constructs.append((name, assembly))
        else:
            print('Error at', name)
        
    df_constructs = pd.DataFrame(constructs, columns=['name', 'sequence'])
    df_constructs['size'] = df_constructs['sequence'].apply(lambda x: len(x))
    df_constructs['amplicon'] = df_constructs['sequence'].apply(len_amplicon)
    return df_constructs

def simulate_assembly(assembly_plan_file, maps, raw_fastas, enzyme):

	#make sure it is correct files
	fastas = [f for f in raw_fastas if not f.startswith('__')]

	plasmids = pd.DataFrame([(p.id, str(p.seq), fasta) for fasta in fastas \
                        for p in list(SeqIO.parse(RESOURCES + '{}'.format(fasta), 'fasta'))], \
                        columns=['name', 'sequence', 'level'])
	plasmids['level'] = plasmids['level'].str.split('.', expand=True)[0]
	plasmids['sequence'] = plasmids['sequence'].str.upper()
	plasmids.loc[plasmids['level']!='fragments', 'sequence'] = plasmids['sequence'].apply(reindex_ps1)

	if maps:
		plasmids = pd.merge(pd.read_csv(RESOURCES + '{}'.format(maps[0])), plasmids, \
                    left_on='full_name', right_on='name', how='right')[['id', 'short_name', 'name', 'sequence', 'level']]
	else:
		plasmids['id'] = ['p{}'.format(i) for i in np.arange(len(plasmids))]
		plasmids['short_name'] = plasmids['name']
		plasmids = plasmids[['id', 'short_name', 'name', 'sequence', 'level']]

	mapping = {
		'UAC': 'pJUMP18-Uac',
		'1A': 'pJUMP29-1A(sfGFP)',
		'1B': 'pJUMP29-1B(sfGFP)',
		'1B*': 'pJUMP29-1B*(sfGFP)',
		'1C': 'pJUMP29-1C(sfGFP)',
		'1C*': 'pJUMP29-1C*(sfGFP)',
		'1D\'': 'pJUMP29-1D\'(sfGFP)',
		'1D': 'pJUMP29-1D(sfGFP)',
		'1Ep': 'pJUMP29-1E\'(sfGFP)',
		'1E': 'pJUMP29-1E(sfGFP)',
		'1F': 'pJUMP29-1F(sfGFP)',
		'2A': 'pJUMP49-2A(sfGFP)',
		'2B': 'pJUMP49-2B(sfGFP)',
		'2B*': 'pJUMP49-2B*(sfGFP)',
		'2C': 'pJUMP49-2C(sfGFP)',
		'2C*': 'pJUMP49-2C*(sfGFP)',
		'2Dp': 'pJUMP49-2D\'(sfGFP)',
		'2D': 'pJUMP49-2D(sfGFP)',
		'2E': 'pJUMP49-2E(sfGFP)',
		'B0033_RN': 'pJUMP18-B0033-MV_RN',
		'B0033_R': 'pJUMP18-B0033_R',
		'B0034_RN': 'pJUMP18-B0034-MV_RN',
		'B0034_R': 'pJUMP18-B0034_R',
		'mCherry_O': 'pJUMP18-mCherry_O',
		'sGFP_O': 'pJUMP18-sfGFP_O',
		'B0015_CT': 'pJUMP19-B0015_CT',
		'B0015_T': 'pJUMP19-B0015_T',
		'J100_P': 'pJUMP19-23100_P',
	}

	assembly_plan = pd.read_csv(RESOURCES + '{}'.format(assembly_plan_file[0]))

	if enzyme == 'Even / BsmbI' and assembly_plan.shape[1]==4:

		fragments = plasmids[plasmids['level']=='fragments']
		fragments_map = dict(zip(fragments['name'], fragments['sequence']))
		vectors = plasmids[plasmids['level']=='vectors']
		vectors_map = dict(zip(vectors['name'], vectors['sequence']))
		uac = vectors_map[mapping['UAC']]

		lvl0 = assemble_lvl_0(fragments_map, uac)
		new_promoters = lvl0[lvl0['name'].isin(['PBAD-RiboJ', 'PLuxB-RiboJ', 'PSalTTC-RiboJ', 'PBetI-RiboJ', 'PTac-RiboJ'])]
		new_promoters['name'] = 'pJ0-' + new_promoters['name'] + '_P'
		new_promoters['size'] = new_promoters['sequence'].apply(lambda x: len(x))
		#new_promoters.to_csv('datasets/jump/level-0-outputs.csv', index=False)

		print(new_promoters)

	else:

		#assembly_plan = pd.read_csv(RESOURCES + '{}'.format(assembly_plan_file[0]))
		assembly = assembly_plan.iloc[:, 2:].melt()
		assembly['name'] = assembly['value'].apply(lambda x: mapping[x] if x in mapping else '')
		assembly.loc[assembly['name']=='', 'name'] = 'pJ0-' + assembly['value'] if enzyme == 'Odd / BsaI' else 'pJ1-' + assembly['value']
		assembly = assembly.drop_duplicates().reset_index(drop=True)
		assembly = pd.merge(assembly, plasmids, on='name', how='left')
		assembly_plan_map = dict(zip(assembly['value'], assembly['sequence']))

		odd = True if enzyme == 'Odd / BsaI' else False
		assembly_result = assemble(assembly_plan, assembly_plan_map, odd)
		assembly_result = assembly_result[~assembly_result['name'].isin(plasmids['name'].tolist())]

		write_fasta(assembly_result)

		#else:

			#lvl2_assembly = assembly_plan.iloc[:, 2:].melt()
			#lvl2_assembly['name'] = lvl2_assembly['value'].apply(lambda x: mapping[x] if x in mapping else '')
			#lvl2_assembly.loc[lvl2_assembly['name']=='', 'name'] = lvl2_assembly['value']
			#lvl2_assembly = lvl2_assembly.drop_duplicates().reset_index(drop=True)
			#lvl2_assembly = pd.merge(lvl2_assembly, plasmids, on='name', how='left')
			#assembly_plan_map = dict(zip(lvl2_assembly['value'], lvl2_assembly['sequence']))

			#assemble_lvl_2_res = assemble(assembly_plan, assembly_plan_map, False)
			#assemble_lvl_2_res = assemble_lvl_2_res[~assemble_lvl_2_res['name'].isin(plasmids['name'].tolist())]



	


	return 'OK'
