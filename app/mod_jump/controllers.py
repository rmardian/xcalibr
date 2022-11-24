from http.client import FORBIDDEN
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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Wedge, FancyArrow, FancyBboxPatch, Rectangle, Circle
from matplotlib.lines import Line2D

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(APP_ROOT, '../resources/')
INPUT = os.path.join(RESOURCES, './input/')
OUTPUT = os.path.join(RESOURCES, './output/')

### FILE-RELATED FUNCTIONS ###

def handle_file(input_file, format='tabular'):

	valid_tab = ['csv']
	valid_seq = ['fa', 'fasta', 'zip']
	filename = input_file.filename
	ext = filename.rsplit('.', 1)[1].lower()

	if '.' in filename and ((format=='tabular' and ext in valid_tab) or (format=='sequence' and ext in valid_seq)):
		path_to_file = os.path.join(INPUT, secure_filename(filename))
		input_file.save(path_to_file)
		if ext == 'zip':
			with ZipFile(path_to_file, 'r') as zip_obj:
				plasmid_files = zip_obj.namelist()
				zip_obj.extractall(INPUT)
			return plasmid_files
		else:
			return [filename]
	return None

def zip_output(filenames, output_name, sequence=False, sbol=False):

    print(filenames)

    path_to_file = os.path.join(OUTPUT, secure_filename('{}.zip'.format(output_name)))
    with ZipFile(path_to_file, 'w') as zip_obj:
        for filename in filenames:
            zip_obj.write(OUTPUT + filename, filename)
        if sequence:
            for fasta in os.listdir(OUTPUT + '/fasta'):
                zip_obj.write(OUTPUT + '/fasta/' + fasta, '/fasta/' + fasta)
        if sbol:
            for img in os.listdir(OUTPUT + '/visualization'):
                zip_obj.write(OUTPUT + '/visualization/' + img, img)

def write_fasta(assembly_result):

    shutil.rmtree(OUTPUT + 'fasta')
    os.mkdir(OUTPUT + 'fasta')
    for i, row in assembly_result.iterrows():
        record = SeqRecord(Seq(row['sequence']), id=row['name'], description='')
        with open(OUTPUT + 'fasta/{}.fasta'.format(row['name']), 'w') as handle:
            SeqIO.write(record, handle, 'fasta')

### DOMESTICATE PARTS ###

def find_templates(parts, plasmids, mapping=None):

	templates = pd.DataFrame([(part['name'], plasmid['name']) for _, part in parts.iterrows() \
                      for _, plasmid in plasmids.iterrows() if plasmid['sequence'].find(part['sequence'])!=-1], \
                      columns=['name', 'template'])
	if mapping:
		templates['template'] = templates['template'].apply(lambda x: mapping[x] if x in mapping else x)
	templates = templates.groupby('name')['template'] \
                      .apply(lambda x: ', '.join(x)).reset_index()
	missing_parts = list(set(parts['name'].tolist()).difference(set(templates['name'].tolist())))

	return templates, missing_parts

def check_internal_sites(parts_ori, forbidden_sites):

	parts = parts_ori.copy()
	parts['num_forbidden_sites'] = parts['sequence'].apply(lambda x: np.sum([x.count(a) for a in forbidden_sites]))
	invalid_parts = parts[parts['num_forbidden_sites']>0].reset_index(drop=True)
	valid_parts = parts[parts['num_forbidden_sites']==0].reset_index(drop=True)

	return valid_parts, invalid_parts

def design_fragments(parts, prefix, suffix):

	fragments = pd.merge(parts, pd.read_csv(RESOURCES + '/overhang.csv'), on='overhang', how='left') \
                [['name', 'sequence', 'left_site', 'right_site']]
	fragments['left_overhang'] = prefix + fragments['left_site']
	fragments['right_overhang'] = fragments['right_site'] + suffix
	fragments['ext_sequence'] = fragments['left_overhang'] + fragments['sequence'] + fragments['right_overhang']
	fragments['size'] = fragments['ext_sequence'].apply(lambda x: len(x)+1) #following benchling convention to start from index 1

	return fragments

def design_primers(parts, prefix, suffix):

    print(parts)
    
    failed_parts = []
    primers_list = []
    for _, part in parts.iterrows():
        try:
            fp, rp = primers(part['sequence'], add_fwd=prefix + part['left_site'],
                             add_rev=str(Seq(suffix).reverse_complement()) + str(Seq(part['right_site']).reverse_complement()))
            #forward primers
            primers_list.append((part['name'], fp.seq, 
                                 fp.tm, fp.tm_total, fp.gc, fp.dg, fp.fwd, fp.offtargets, fp.penalty))
            #reverse primers
            primers_list.append((part['name'], rp.seq,
                                 rp.tm, rp.tm_total, rp.gc, rp.dg, rp.fwd, rp.offtargets, rp.penalty))
        except:
            print('Error on', part['name'], ': cannot find any feasible primer design, check your fragments!')
            failed_parts.append(part['name'])
            continue
        
    unnamed_primers = pd.DataFrame(primers_list, columns=['name', 'sequence', 'tm', 'tm_total', 'gc', 'dg', 'fwd', 'offtargets', 'penalty'])
	#assign names to primers, make sure to remove duplicated primers
    final_primers = unnamed_primers[['fwd', 'sequence']].drop_duplicates().reset_index(drop=True)
    final_primers['primer_name'] = pd.Series(final_primers.index).apply(lambda x: 'P' + str(x+1).zfill(len(str(final_primers.shape[0]))))
    final_primers.loc[final_primers['fwd'], 'primer_name'] = final_primers['primer_name'] + '.F'
    final_primers.loc[~final_primers['fwd'], 'primer_name'] = final_primers['primer_name'] + '.R'
    final_primers = pd.merge(final_primers[['primer_name', 'sequence']], unnamed_primers.drop('name', axis=1), \
							on='sequence', how='left').drop_duplicates().reset_index(drop=True)
	#assign primers to their respective parts
    parts_x_primers = pd.merge(unnamed_primers[['name', 'sequence', 'fwd']], final_primers[['primer_name', 'sequence']], on='sequence', how='left')
    
    return final_primers, parts_x_primers, failed_parts

def get_annealed_part(a, b):

    return [a[i].replace(b[i], '') for i in range(len(b))]

def calculate_tm(seq):

	return np.median([mt.Tm_Wallace(seq), mt.Tm_GC(seq), mt.Tm_NN(seq)])

def generate_pcr(fragments, parts_x_primers, templates):

	pcr_primers = pd.merge(parts_x_primers, \
                fragments[['name', 'ext_sequence', 'left_overhang', 'right_overhang', 'size']], \
                on='name', how='left')
	pcr_primers.loc[~pcr_primers['fwd'], 'sequence'] = pcr_primers['sequence']. \
				apply(lambda x: str(Seq(x).reverse_complement()))
	pcr_primers.loc[pcr_primers['fwd'], 'overhang_part'] = pcr_primers['left_overhang']
	pcr_primers.loc[~pcr_primers['fwd'], 'overhang_part'] = pcr_primers['right_overhang']
	pcr_primers['annealing_part'] = get_annealed_part(pcr_primers['sequence'], pcr_primers['overhang_part'])
	pcr_primers['tm'] = np.round(pcr_primers['annealing_part'].apply(calculate_tm), 1)
	pcr_primers['tm_product'] = np.round(pcr_primers['ext_sequence'].apply(calculate_tm), 1)
	pcr_primers['ext_time'] = pcr_primers['size'].apply(lambda x: max([60, x/1000 * 60]))
	fwd_primers = pcr_primers[pcr_primers['fwd']] \
				[['name', 'primer_name', 'tm', 'tm_product', 'size', 'ext_time']]
	rev_primers = pcr_primers[~pcr_primers['fwd']] \
				[['name', 'primer_name', 'tm', 'tm_product', 'size', 'ext_time']]
	pcr_rxn = pd.merge(fwd_primers, rev_primers, on=['name', 'tm_product', 'size', 'ext_time'], how='outer')
	pcr_rxn['ta'] = 0.3 * np.max([pcr_rxn['tm_x'], pcr_rxn['tm_y']], axis=0) + 0.7 * pcr_rxn['tm_product'] - 14.9
	pcr_rxn['primer_name_x'] = pcr_rxn['primer_name_x'].str.split('_', expand=True)[0] \
						.str.replace('(', '').str.replace(')', '')
	pcr_rxn['primer_name_y'] = pcr_rxn['primer_name_y'].str.split('_', expand=True)[0] \
						.str.replace('(', '').str.replace(')', '')
	pcr_rxn = pcr_rxn[['name', 'size', 'primer_name_x', 'primer_name_y', 'ta', 'ext_time']] \
	.sort_values('ta').reset_index(drop=True)
	pcr_rxn.columns = ['name', 'size', 'forward_primer', 'reversed_primer', 'annealing_temp', 'ext_time']
	pcr_rxn = pd.merge(pcr_rxn, templates, on='name', how='left')
	return pcr_rxn

def generate_lvl_0_map(fragments, mapping):

	uac = 'AGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGATGCCTTTAATTAAGGAGTTTTGCAGGTGCACCTGCTTTTCGCTGAATTCGCGGCCGCTTCTAGAGGGTCTGCGATGTTTGGTCTTGAGACGACTGTGACAAGGAGTTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGCGTACTTGTTTAACTTTAAGAAGGAGATATACAATGGTAGCCCGTAAAGGCGAAGAGCTGTTCACTGGTGTCGTCCCTATTCTGGTGGAACTGGATGGTGATGTCAACGGTCATAAGTTTTCCGTGCGTGGCGAGGGTGAAGGTGACGCAACTAATGGTAAACTGACGCTGAAGTTCATCTGTACTACTGGTAAACTGCCGGTACCTTGGCCGACTCTGGTAACGACGCTGACTTATGGTGTTCAGTGCTTTGCTCGTTATCCGGACCATATGAAGCAGCATGACTTCTTCAAGTCCGCCATGCCGGAAGGCTATGTGCAGGAACGCACGATTTCCTTTAAGGATGACGGCACGTACAAAACGCGTGCGGAAGTGAAATTTGAAGGCGATACCCTGGTAAACCGCATTGAGCTGAAAGGCATTGACTTTAAAGAAGATGGCAATATCCTGGGCCATAAGCTGGAATACAATTTTAACAGCCACAATGTTTACATCACCGCCGATAAACAAAAAAATGGCATTAAAGCGAATTTTAAAATTCGCCACAACGTGGAGGATGGCAGCGTGCAGCTGGCTGATCACTACCAGCAAAACACTCCAATCGGTGATGGTCCTGTTCTGCTGCCAGACAATCACTATCTGAGCACGCAAAGCGTTCTGTCTAAAGATCCGAACGAGAAACGCGATCATATGGTTCTGCTGGAGTTCGTAACCGCAGCGGGCATCACGCATGGTATGGATGAACTGTACAAAGGTTCGTAACTTCTGACTGAGTTGCACGCTCCTTGGTCAGCGTCTCAGACCTTTCATCGCGACCTACTAGTAGCGGCCGCTGCAGGGAGTTGTCTTCGAAGACTTCGCTCTAGTCTTGGACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATTTGTTCAGAACGCTCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAATCCAGGGGTCCCCAATAATTACGATTTAAATTAGTAGCCCGCCTAATGAGCGGGCTTTTTTTTAATTCCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGCATTCAGCATTTTCGTGTGGCGCTGATTCCGTTTTTTGCGGCGTTTTGCCTGCCGGTGTTTGCGCATCCGGAAACCCTGGTGAAAGTGAAAGATGCGGAAGATCAACTGGGTGCGCGCGTGGGCTATATTGAACTGGATCTGAACAGCGGCAAAATTCTGGAATCTTTTCGTCCGGAAGAACGTTTTCCGATGATGAGCACCTTTAAAGTGCTGCTGTGCGGTGCGGTTCTGAGCCGTGTGGATGCGGGCCAGGAACAACTGGGCCGTCGTATTCATTATAGCCAGAACGATCTGGTGGAATATAGCCCGGTGACCGAAAAACATCTGACCGATGGCATGACCGTGCGTGAACTGTGCAGCGCGGCGATTACCATGAGCGATAACACCGCGGCGAACCTGCTGCTGACGACCATTGGCGGTCCGAAAGAACTGACCGCGTTTCTGCATAACATGGGCGATCATGTGACCCGTCTGGATCGTTGGGAACCGGAACTGAACGAAGCGATTCCGAACGATGAACGTGATACCACCATGCCGGCAGCAATGGCGACCACCCTGCGTAAACTGCTGACGGGTGAGCTGCTGACCCTGGCAAGCCGCCAGCAACTGATTGATTGGATGGAAGCGGATAAAGTGGCGGGTCCGCTGCTGCGTAGCGCGCTGCCGGCTGGCTGGTTTATTGCGGATAAAAGCGGTGCGGGCGAACGTGGCAGCCGTGGCATTATTGCGGCGCTGGGCCCGGATGGTAAACCGAGCCGTATTGTGGTGATTTATACCACCGGCAGCCAGGCGACGATGGATGAACGTAACCGTCAGATTGCGGAAATTGGCGCGAGCCTGATTAAACATTGGTAAACCGATACAATTAAAGGCTCCTTTTGGAGCCTTTTTTTTTGGACGACCCTTGTCCTTTTCCGCTGCATAACCCTGCTTCGGGGTCATTATAGCGATTTTTTCGGTATATCCATCCTTTTTCGCACGATATACAGGATTTTGCCAAAGGGTTCGTGTAGACTTTCCTTGGTGTATCCAACGGCGTCAGCCGGGCAGGATAGGTGAAGTAGGCCCACCCGCGAGCGGGTGTTCCTTCTTCACTGTCCCTTATTCGCACCTGGCGGTGCTCAACGGGAATCCTGCTCTGCGAGGCTGGCCGTAGGCCGGCCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCAGGCGCGCCCAGCTGTCT'

	if mapping:
		fragments['name'] = fragments['name'].map(mapping)
	fragments_map = dict(zip(fragments['name'], fragments['ext_sequence']))
	return assemble_lvl_0(fragments_map, uac)

def domesticate_parts(parts, mapping, plasmids, output_name):

    mapping_full = dict(zip(mapping['full_name'], mapping['short_name'])) if mapping is not None else None
    mapping_short = dict(zip(mapping['short_name'], mapping['full_name'])) if mapping is not None else None

    recognition_sites = np.array([(enz, str(Seq(enz).reverse_complement())) for enz in [BsaI.site, BsmBI.site]]).ravel().tolist()
    prefix, suffix = recognition_sites[2] + recognition_sites[0] + 'A', 'T' + recognition_sites[1] + 'T' + recognition_sites[3]

    #generate templates and fragments, might need to print at some point?
    templates, missing_parts = find_templates(parts, plasmids, mapping_full)
    valid_parts, invalid_parts = check_internal_sites(parts[~parts['name'].isin(missing_parts)], recognition_sites)

    fragments = pd.DataFrame()
    failed_parts = pd.DataFrame()
    final_primers = pd.DataFrame()

    if len(valid_parts)>0:

        fragments = design_fragments(valid_parts, prefix, suffix)

        #generate primers
        final_primers, parts_x_primers, failed_parts = design_primers(fragments, prefix, suffix)
        final_primers.to_csv(OUTPUT + 'primers.csv', index=False)
        
        #generate pcr rxn
        pcr_rxn = generate_pcr(fragments, parts_x_primers, templates)
        pcr_rxn.to_csv(OUTPUT + 'pcr_rxn.csv', index=False)

        #generate plasmid maps
        level_0 = generate_lvl_0_map(fragments[['name', 'ext_sequence']], mapping_short)
        level_0.to_csv(OUTPUT + 'level-0-parts.csv', index=False)
        write_fasta(level_0)

        #generate visual SBOL
        # shutil.rmtree(RESOURCES + 'visualization')
        # os.mkdir(RESOURCES + 'visualization')
        # for _, sample in valid_parts[['name', 'overhang', 'color']].iterrows():
        #	generate_part_figure(sample, RESOURCES + 'visualization/{}.png'.format(sample[0]))

        #zip output
        zip_output(['primers.csv', 'pcr_rxn.csv', 'level-0-parts.csv'], output_name, sequence=True)

    else:
        zip_output([], output_name)

    #final output
    missing_msg = 'The following part(s) cannot be created due to missing template(s): <b>{}</b><br/><br/>'.format(missing_parts) if len(missing_parts)>0 else ''
    forbidden_msg = 'The following part(s) cannot be created due to internal restriction sites: <b>{}</b>. Future update will support further domestication to remove these sites.<br/><br/>'.format(invalid_parts['name'].tolist()) if invalid_parts.shape[0]>0 else ''
    failed_msg = 'The following part(s) cannot be created due to difficulty in finding the most optimum primers: <b>{}</b>.<br/><br/>'.format(failed_parts) if len(failed_parts) > 0 else ''
    success_msg = 'Generated <b>{}</b> primers for domesticating <b>{}</b> level-0 parts.<br/><br/>'.format(final_primers.shape[0], valid_parts.shape[0])

    return missing_msg + forbidden_msg + failed_msg + success_msg

### SIMULATE ASSEMBLY ###

def reindex_ps1(plasmid):
    '''Reindex plasmids to start from the annealing region of PS1 so the fragment will be in the middle of the sequence.
    This will do nothing for a fragment from pcr products as they do not have PS1 region'''
    
    new_start = plasmid.find('AGGGCGGCGGATTTGTCC')
    if new_start == -1:
        return plasmid
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
            
    df_constructs = pd.DataFrame(constructs, columns=['name', 'sequence'])
    df_constructs['size'] = df_constructs['sequence'].apply(lambda x: len(x))
    df_constructs['amplicon'] = df_constructs['sequence'].apply(len_amplicon)
    return df_constructs

def assemble(assembly_plan, mapping, odd):
	
    constructs = []
    incorrect_assembly = []
    for i, entry in assembly_plan.iterrows():

        _id = entry[0]
        name = entry[1]
        vector = entry[-1]
        parts = entry[2:-1].tolist()
        
        fragments = []
        if vector in mapping:
            fragments.append(get_sites(mapping[vector], odd_level=odd, vector=True))
        for part in parts:
            if part in mapping:
                fragments.append(get_sites(mapping[part], odd_level=odd, vector=False))
            
        sites = list(map(list, zip(*[fragment[:2] for fragment in fragments])))
        if len(sites) > 0:
            sites[0] = sites[0][1:] + [sites[0][0]]
            print(sites)
            if (sites[0]==sites[1]):
                assembly = reindex_ps1(''.join([fragment[2][:-4] for fragment in fragments]))
                constructs.append((name, assembly))
            else:
                incorrect_assembly.append(name)

    df_constructs = pd.DataFrame(constructs, columns=['name', 'sequence'])
    df_constructs['size'] = df_constructs['sequence'].apply(lambda x: len(x))
    df_constructs['amplicon'] = df_constructs['sequence'].apply(len_amplicon)

    return df_constructs, incorrect_assembly

def simulate_assembly(plan, mapping, plasmids, enzyme, output_name):

	mapping_full = dict(zip(mapping['full_name'], mapping['sample_id'])) if mapping is not None else None
	mapping_short = dict(zip(mapping['sample_id'], mapping['full_name'])) if mapping is not None else None

	#reindex plasmid to start from annealing region of PS1 primer
	plasmids['sequence'] = plasmids['sequence'].apply(reindex_ps1)

	assembly = plan.iloc[:, 2:].melt().rename(columns={'value': 'name'}).drop_duplicates().reset_index(drop=True)
	if mapping_full is not None:
		plasmids['name'] = plasmids['name'].apply(lambda x: mapping_full[x] if x in mapping_full else x)
	assembly = pd.merge(assembly, plasmids, on='name', how='left')

	print(plasmids, mapping_full, assembly)

	missing_plasmids = assembly[assembly['sequence'].isnull()]['name'].tolist()
	assembly = assembly.dropna().reset_index(drop=True)
	plan_mapping = dict(zip(assembly['name'], assembly['sequence']))

	odd = True if enzyme == 'Odd / BsaI' else False
	assembly_result, incorrect_assembly = assemble(plan, plan_mapping, odd)
	assembly_result = assembly_result[~assembly_result['name'].isin(plasmids['name'].tolist())]
	
	if mapping_short:
		assembly_result['name'] = assembly_result['name'].apply(lambda x: mapping_short[x] if x in mapping_short else x)
	assembly_result.to_csv(OUTPUT + 'assembled-constructs.csv', index=False)

	write_fasta(assembly_result)

	#zip output
	zip_output(['assembled-constructs.csv'], output_name, sequence=True)

	#final output
	success_msg = 'Assembled <b>{}</b> constructs out of <b>{}</b> assembly plans.<br/><br/>'.format(assembly_result.shape[0], plan.shape[0])
	missing_msg = 'Could not find the correct plasmid for the following part(s): <b>{}</b>. Please check the plasmids\' name.<br/><br/>'.format(missing_plasmids) if len(missing_plasmids) > 0 else ''
	failed_msg = 'Parts\' overhang did not match for the following plan(s): <b>{}</b>. Please make sure to use parts with correct overhangs.<br/><br/>'.format(incorrect_assembly) if len(incorrect_assembly) > 0 else ''

	return success_msg + missing_msg + failed_msg

### VISUALIZE CONSTRUCT ###

def create_promoter(ax, start, base, color='#000', text='promoter',
                    height=1, width=0.8, margin=0.25, lw=2, text_offset=1.2, z=-2):
    
    prom_dim = (height, width) #height, width
    arrow_dim = (0.5 * height, 0.5 * width)
    prom_len = (2 * margin) + (prom_dim[1] + arrow_dim[1]) # how many times longer than the width

    prom_start = start #plasmid_pos[0] + roundness
    prom_end = prom_start + prom_len
    prom_pos = prom_start + margin

    ax.add_line(Line2D([prom_pos, prom_pos], [base, base +prom_dim[0]], lw=lw, color=color, zorder=z))
    ax.add_line(Line2D([prom_pos, prom_pos+prom_dim[1]], [base+prom_dim[0], base+prom_dim[0]], lw=lw, color=color, zorder=z))
    ax.add_patch(Polygon([(prom_pos+prom_dim[1], base+prom_dim[0]-arrow_dim[0]/2),
                          (prom_pos+prom_dim[1]+arrow_dim[1], base+prom_dim[0]),
                          (prom_pos+prom_dim[1], base+prom_dim[0]+arrow_dim[0]/2)], color=color, zorder=z))
    ax.add_line(Line2D([prom_start, prom_end], [base, base], lw=lw, color='#000', zorder=-1))

    ax.text(prom_pos+((prom_dim[1]+arrow_dim[1])/2), base-text_offset, text, fontsize=12, horizontalalignment='center')
    
    return prom_end
    
def create_rbs(ax, start, base, color='#666', text='rbs',
               radius=0.5, margin=0.25, lw=2, text_offset=1):
    
    rbs_size = radius
    rbs_len = (2 * margin) + (2 * rbs_size)  # how many times longer than the radius

    rbs_start = start #prom_end
    rbs_end = rbs_start + rbs_len
    rbs_pos = rbs_start + margin + rbs_size

    ax.add_patch(Wedge((rbs_pos, base), rbs_size, 0, 180, lw=lw, color='gray', zorder=-2))
    ax.add_line(Line2D([rbs_start, rbs_end], [base, base], lw=lw, color='black', zorder=-1))
    
    ax.text(rbs_pos, base-text_offset, text, fontsize=12, horizontalalignment='center')
    
    return rbs_end

def create_cds(ax, start, base, color='#ff0000', text='cds',
               length=2, width=1, margin=0.25, lw=2, text_offset=1.5):
    
    cds_dim = (length, width) #length, width
    head_dim = (width, 1.6 * width) #length, width
    cds_len = (2 * margin) + (cds_dim[0] + head_dim[0]) #how many times longer than the width

    cds_start = start #rbs_end
    cds_end = cds_start + cds_len
    cds_pos = cds_start + margin
    ax.add_patch(FancyArrow(cds_pos, base, cds_dim[0], 0, width=cds_dim[1],
                            head_length=head_dim[0], head_width=head_dim[1], lw=lw,
                            edgecolor='#000', facecolor=color))
    ax.add_line(Line2D([cds_start, cds_end], [base, base], lw=lw, color='black', zorder=-1))
    
    ax.text(cds_pos+((cds_dim[1]+head_dim[1])/2), base-text_offset, text, fontsize=12, horizontalalignment='center')
    
    return cds_end

def create_terminator(ax, start, base, color='#000', text='terminator',
                      height=1, width=0.8, margin=0.25, lw=2, text_offset=1):
    
    term_dim = (height, width) #height, width
    term_len = (2*margin) + term_dim[1] #how many times longer than the width

    term_start = start #cds_end
    term_end = term_start + term_len
    term_pos = term_start + margin + (term_dim[1]/2)

    ax.add_line(Line2D([term_pos, term_pos], [base, base+term_dim[0]], lw=lw, color='black'))
    ax.add_line(Line2D([term_pos-term_dim[1]/2, term_pos+term_dim[1]/2],
                       [base+term_dim[0], base+term_dim[0]], lw=lw, color='black'))
    ax.add_line(Line2D([term_start, term_end], [base, base], lw=lw, color='black', zorder=-1))
    
    ax.text(term_pos, base-text_offset, text, fontsize=12, horizontalalignment='center')
    
    return term_end
    
def create_tu(ax, prom_start, base, cols=['#000', '#666', '#ff0000', '#000'], texts=['', '', '', ''], margin=0.25, lw=2):
    
    prom_end = create_promoter(ax, start=prom_start, base=base, color=cols[0], text=texts[0], margin=margin, lw=lw)
    rbs_end = create_rbs(ax, start=prom_end, base=base, color=cols[1], text=texts[1], margin=margin, lw=lw)
    cds_end = create_cds(ax, start=rbs_end, base=base, color=cols[2], text=texts[2], margin=margin, lw=lw)
    term_end = create_terminator(ax, start=cds_end, base=base, color=cols[3], text=texts[3], margin=margin, lw=lw)
    return term_end

def create_plasmid(ax, height, width, roundness=0.6, color='#fff', text='pBR322', lw=2, text_offset=1.5):
    
    plasmid_dim = (height, width)
    plasmid = FancyBboxPatch((roundness, roundness),
                           plasmid_dim[1]-2*roundness, plasmid_dim[0]-2*roundness,
                           boxstyle="round, pad={}".format(roundness), fill=None, lw=lw, zorder=-1)
    ax.add_patch(plasmid)
    ax.add_patch(Circle((plasmid_dim[1]/2, 0), 1.5, edgecolor='black', facecolor='white', lw=lw))
    ax.text(width/2, 0, text, fontsize=12, horizontalalignment='center')

def disassemble_parts(row, plan, plan_idx):
    
    output = []
    for i in row.iloc[2:-1].tolist():
        if i not in plan_idx:
            output.append(i)
        else:
            output.append(disassemble_parts(plan.iloc[plan_idx[i]], plan, plan_idx))
    return output

def generate_constructs(plan, parts):

	plan_idx = {a: b for a, b in zip(plan['name'], plan.index)}
	color_dict = {a: b for a, b in zip(parts['name'], parts['color'])}

	constructs = []
	for _, row in plan.iterrows():
		parts = np.array(disassemble_parts(row, plan, plan_idx)).flatten().tolist()
		colors = [color_dict[i] if i in color_dict else '#000000' for i in parts]
		constructs.append((row['name'], parts, colors))

	return constructs

def generate_part_figure(sample, filename, margin=0.4, w_offset=1, h_offset=4, lw=2):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # individual part
    if sample[1]=='P':
        point = create_promoter(ax, margin, 0, margin=margin, text=sample[0], color=sample[2])
    elif sample[1]=='R' or sample[1]=='RN':
        point = create_rbs(ax, margin, 0, margin=margin, text=sample[0], color=sample[2])
    if sample[1]=='O' or sample[1]=='NOC':
        point = create_cds(ax, margin, 0, margin=margin, text=sample[0], color=sample[2])
    elif sample[1]=='T' or sample[1]=='CT':
        point = create_terminator(ax, margin, 0, margin=margin, text=sample[0], color=sample[2])

    width = margin + point + w_offset
    height = h_offset
    
    fig.set_figwidth(width/5*4)
    fig.set_figheight(height/5*4)
    ax.set_xlim(-(w_offset/2), width-(w_offset/2))
    ax.set_ylim(-(h_offset/2), height-(h_offset/2))
    
    plt.axis('off')
    plt.savefig(filename, dpi=300)

def generate_tu_figure(sample, filename, show_plasmid=True, margin=0.18, plasmid_height=4, roundness=0.6, 
                   w_offset=1, h_offset=4, lw=2):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # transcriptional unit
    point = roundness
    num_part = 4
    for i in range(0, len(sample[1]), num_part):
        point = create_tu(ax, point, plasmid_height, sample[2][i:i+num_part], sample[1][i:i+num_part], margin=margin, lw=lw)
    
    # plasmid
    if show_plasmid:
        
        plasmid_width = roundness + point
        create_plasmid(ax, plasmid_height, plasmid_width, roundness, lw=lw)
        width = plasmid_width + w_offset
        height = plasmid_height + h_offset
    
    fig.set_figwidth(width/5*4)
    fig.set_figheight(height/5*4)
    ax.set_xlim(-(w_offset/2), width-(w_offset/2))
    ax.set_ylim(-(h_offset/2), height-(h_offset/2))
    
    plt.axis('off')
    plt.title(sample[0])
    plt.savefig(filename, dpi=300)

def visualize_construct(plan_file, parts_file, output_name):

	shutil.rmtree(OUTPUT + 'visualization')
	os.mkdir(OUTPUT + 'visualization')

	plan = pd.read_csv(INPUT + '{}'.format(plan_file[0]))
	parts = pd.read_csv(INPUT + '{}'.format(parts_file[0]))

	constructs = generate_constructs(plan, parts)

	#print(constructs)
	for sample in constructs:
		generate_tu_figure(sample, OUTPUT + 'visualization/{}.png'.format(sample[0]))

	#zip output
	zip_output([], output_name, sbol=True)