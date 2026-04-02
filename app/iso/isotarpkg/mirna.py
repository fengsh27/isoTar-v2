# -*- coding: utf-8 -*-

''' Libraries '''
import re
import os, errno
import os.path
import sys
import shutil
from operator import itemgetter
import csv
import json

class MiRNAsProcessing():

    def __init__(self, pre_mature_path, pre_mat_acc_path):
        ''' miRNA constructor '''
        self.pre_mature_path = pre_mature_path
	self.pre_mat_acc_path = pre_mat_acc_path

        if not os.path.isfile(self.pre_mature_path):
            self.ePrintClose("Error. No mature_pre_mirna.json file found.")

	if not os.path.isfile(self.pre_mat_acc_path):
	    self.ePrintClose("Error. No mat_pre_acc.json file found.")

        self.mirna_path = ''
        self.debug = False
        self.mirnas_list = []
        self.mirnas_split_list = []
        self.no_mirnas = 0
        self.longest_mirna_seq = 0
        self.no_cores = 0
        self.mirna_split_path_list = []
        self.mirna_split_details = []
        self.mirna_family_info = {}
        self.mirna_family_info_species_id = []
        self.mature_pre_mirnas_list = {}
	self.pre_mat_acc_list = {}

    def checkIntegerValue(self, val):
        try:
            int(val)
            return True
        except ValueError:
            return False

    def ePrint(self, message):
        print >> sys.stderr, message

    def ePrintClose(self, message):
        print >> sys.stderr, message
        sys.exit(1)

    def setDebugging(self, debug):
        self.debug = debug

    def createMaturePreList(self):
        if os.path.isfile(self.pre_mature_path):
            with open(self.pre_mature_path, 'r') as f:
                self.mature_pre_mirnas_list = json.load(f)

    def createMaturePreAccList(self):
        if os.path.isfile(self.pre_mat_acc_path):
            with open(self.pre_mat_acc_path, 'r') as f:
                self.pre_mat_acc_list = json.load(f)

    def removeDirectory(self, path):
        """ param <path> could either be relative or absolute. """
        if os.path.isfile(path):
            os.remove(path)  # remove the file
        elif os.path.isdir(path):
            try:
                shutil.rmtree(path) # remove dir and all contains
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
        else:
            raise ValueError("file {} is not a file or dir.".format(path))

    def createDirectory(self, path):
        # Create the directory
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    def getMirnaVariations(self, header):
        variations_check_list = []
        header_variations = {
           'substitution': [],
           'isomir': [],
           'error': []
        }

        match_obj = re.match(r'^>.*(hsa\-[^\s]+)\s*(hsa\-[^\s]+)\s*([^\s]*).*$',
                             header,
                             re.M|re.I)
	
        if match_obj:
            mirna_id = match_obj.group(1)
            # Get all specified variations removing duplicates
            variations = list(set(match_obj.group(3).split(';')))

            #>hsa-miR-XXX hsa-mir-XXX 0|+1;-1|+2;2:A|G;5:-|A;7:G|-
            for variation in variations:
                match_obj = re.match(r'^\s*([0-9]+)\s*:\s*([ACGTU\-]{1})\s*\|\s*([\-ACGTU]{1})\s*$',
                                     variation.upper(),
                                     re.M|re.I)
                # Nucleotide substitution
                if match_obj:
                    position = int(match_obj.group(1))
                    ref = match_obj.group(2)
                    mod = match_obj.group(3)
                    header_variations['substitution'].append({
                        'variation': variation,
                        'position': position,
                        'ref': ref.upper(),
                        'mod': mod.upper()
                    })
                else:
                    match_obj = re.match(r'^\s*([\+\-0-9]+)\s*\|\s*([\+\-0-9]+)\s*$',
                                         variation.upper(),
                                         re.M|re.I)
                    # isomiR
                    if match_obj:
                        if self.checkIntegerValue(match_obj.group(1)) and \
			    self.checkIntegerValue(match_obj.group(2)):
                            start = int(match_obj.group(1))
                            end = int(match_obj.group(2))
                            header_variations['isomir'].append({
                                'variation': variation,
                                'start': start,
                                'end': end,
                                '5p': match_obj.group(1),
                                '3p': match_obj.group(2)
                            })
                    else:
                        # Wrong variation
                        header_variations['error'].append({
                            'variation': variation,
                            'msg': 'Unrecognized variation'
                        })

        return header_variations

    def sequenceNucleotideSubstitution(self, sequence, position, ref, mod):
        mod_sequence = {
            'status': False,
            'status_msg': \
            'Unrecognized variation (Position: %d, Ref: %s, Mod %s)' \
            % (position, ref, mod),
            'sequence': ''
        }

        str_pos = position -1 # String position
        if str_pos >= 0 and str_pos < len(sequence):
            seq_ref = sequence[str_pos]

            # Nucleotide insertion
            if ref == '-' and mod in ['A', 'C', 'G', 'T', 'U']:
		s = sequence[0:str_pos] + mod + sequence[str_pos:]
                mod_sequence['status'] = True
                mod_sequence['status_msg'] = ''
                mod_sequence['sequence'] = s

            # Nucleotide deletion
            elif ref in ['A', 'C', 'G', 'T', 'U'] and mod == '-':
                s = sequence[0:str_pos] + sequence[str_pos+1:]
                mod_sequence['status'] = True
                mod_sequence['status_msg'] = ''
                mod_sequence['sequence'] = s
            else:
                # Nucleotide substitution
                if seq_ref != ref:
                    mod_sequence['status_msg'] = \
                        'Expected %s, specified %s (Position: %d, Ref: %s, Mod %s)' \
                        % (seq_ref, ref, position, ref, mod)
                elif seq_ref == ref and ref == mod:
                    mod_sequence['status_msg'] = \
                        'Specified the same nucleotide (Position: %d, Ref: %s, Mod %s)' \
                        % (position, ref, mod)
                elif seq_ref == ref and mod in ['A', 'C', 'G', 'T', 'U']:
                    s = list(sequence)
                    s[str_pos] = mod
                    mod_sequence['status'] = True
                    mod_sequence['status_msg'] = ''
                    mod_sequence['sequence'] = ''.join(s)

        return mod_sequence

    def getIsomirSequence(self, mirna_id, sequence, start, end, pre_mirna_id):
        mod_sequence = {
            'status': False,
            'status_msg': \
            'IsomiR outside the canonical pre-miRNA (miRNA ID: %s, Sequence: %s, Start: %d, End: %d)' \
            % (mirna_id, sequence, start, end),
            'sequence': '',
	    'pre_sequence': ''
        }

        if mirna_id in self.mature_pre_mirnas_list:
	    for pre in self.mature_pre_mirnas_list[mirna_id]:
		if pre['pre_id'] == pre_mirna_id:
		    pre_seq = pre['ext_pre_seq'] # pre_seq
		    mod_sequence['pre_sequence'] = pre['pre_seq']
		    # mature mirna local position
		    local_pos = pre_seq.find(sequence)
		    if local_pos != -1:
			isomir_start = local_pos - start
			isomir_end = (local_pos + len(sequence)) + end
			if isomir_start >= 0 and isomir_end <= len(pre_seq):
			    s = pre_seq[isomir_start:isomir_end]
			    mod_sequence['status'] = True
			    mod_sequence['status_msg'] = ''
			    mod_sequence['sequence'] = s
		    break

        return mod_sequence

    def getMirFamilyInfo(self, path):
        mir_map = {}
        if os.path.isfile(path):
            with open(path, 'r') as f:
                handler = csv.reader(f, delimiter='\t')
                next(handler, None)
                for line in handler:
                    if len(line) == 7:
                        seq = line[1].strip()
                        specie_id = line[2].strip()
                        if seq != '' and specie_id != '':

                            # Store all specie IDs
                            if specie_id not in self.mirna_family_info_species_id:
                                self.mirna_family_info_species_id.append(specie_id)

                            # Create the record in the dictionary
                            if seq not in mir_map:
                                mir_map[seq] = ['9606']

                            # Check if the specie ID exists within the seq record
                            if specie_id not in mir_map[seq]:
                                mir_map[seq].append(specie_id)

        self.mirna_family_info = mir_map

    def getMiRNAs(self, content):

        self.createMaturePreList()
	self.createMaturePreAccList()

        mirna_dir_indexes_list = {}
        tmp_mirnas_list = []

        mirna_list = []

        for i in xrange(0, len(content)):
            # mirna header
            m_header = content[i].replace('\n','')

            match_obj1 = re.match(r'^>.*(hsa\-[^\s]+)\s*(hsa\-[^\s]+)$',
                                 m_header,
                                 re.M|re.I)

            match_obj = re.match(r'^>.*(hsa\-[^\s]+)\s*(hsa\-[^\s]+)\s*(.+).*$',
                                 m_header,
                                 re.M|re.I)

            if match_obj1:
                # mirna ID
                m_id = match_obj1.group(1).strip()
		p_id = match_obj1.group(2).strip()
		m_dir = '%s__%s' % (m_id, p_id)

                # ==============================================================
                # Check for duplicate miRNAs
                # ==============================================================
                if m_dir in mirna_dir_indexes_list:
                    mirna_dir_indexes_list[m_dir]['idx'] += 1
                else:
                    mirna_dir_indexes_list[m_dir] = { 'idx': 1 }

		    mir_acc_info = self.pre_mat_acc_list[m_id] if m_id in self.pre_mat_acc_list else {}

		    if mir_acc_info and p_id in mir_acc_info['pre_mirnas']:
			mirna = {
			    'header': m_header,
			    'identifier': m_id,
			    'acc': mir_acc_info['mature_acc'] if 'mature_acc' in mir_acc_info else '',
			    'pre_mirna_id': p_id,
			    'pre_acc': mir_acc_info['pre_mirnas'][p_id]['pre_acc'] if p_id in mir_acc_info['pre_mirnas'] else '',
			    'pre_mirna_seq': '',
			    'pre_mirna_length': 0,
			    'duplicate': 0,
			    'seq': '',
			    'length': 0,
			    'output_path': '',
			    'variations': {
			        'substitution': [],
			        'isomir': [],
			        'error': []
			    }
			}

			# To get the mirna sequence
			i += 1
			tmp_seq = ''

			# Get the mirna sequence
			while i < len(content) and \
			      not re.match(r'^>.*$', content[i], re.M|re.I):
			    tmp_seq += content[i].replace('\n','')
			    i += 1

			# To get the next mirna header
			i -= 1

			# mirna sequence
			mirna['seq'] = tmp_seq.upper().strip()
			mirna['length'] = len(tmp_seq)

			# Get the longest mirna sequence
			if self.longest_mirna_seq < len(tmp_seq):
			    self.longest_mirna_seq = len(tmp_seq)

			# Store the mirna details
			if mirna['identifier'] != '':
			    # Update the number of duplicates
			    mirna['duplicate'] = mirna_dir_indexes_list[m_dir]['idx'] - 1
			    tmp_mirnas_list.append(mirna)

            elif match_obj:
                # mirna ID
                m_id = match_obj.group(1).strip()
		p_id = match_obj.group(2).strip()
		m_dir = '%s__%s' % (m_id, p_id)

                # ==============================================================
                # Check for duplicate miRNAs
                # ==============================================================
                if m_dir in mirna_dir_indexes_list:
                    mirna_dir_indexes_list[m_dir]['idx'] += 1
                else:
                    mirna_dir_indexes_list[m_dir] = { 'idx': 1 }

                    # Retrive all specified variations through the mirna header
                    variations = self.getMirnaVariations(m_header)

		    mir_acc_info = self.pre_mat_acc_list[m_id] if m_id in self.pre_mat_acc_list else {}

		    if mir_acc_info and p_id in mir_acc_info['pre_mirnas']:
			mirna = {
			    'header': m_header,
			    'identifier': m_id,
			    'acc': mir_acc_info['mature_acc'] if 'mature_acc' in mir_acc_info else '',
			    'pre_mirna_id': p_id,
			    'pre_acc': mir_acc_info['pre_mirnas'][p_id]['pre_acc'] if p_id in mir_acc_info['pre_mirnas'] else '',
			    'pre_mirna_seq': '',
			    'pre_mirna_length': 0,
			    'duplicate': 0,
			    'seq': '',
			    'length': 0,
			    'output_path': '',
			    'variations': {
			        'substitution': [],
			        'isomir': [],
			        'error': variations['error']
			    }
			}

			# To get the mirna sequence
			i += 1
			tmp_seq = ''

			# Get the mirna sequence
			while i < len(content) and \
			      not re.match(r'^>.*$', content[i], re.M|re.I):
			    tmp_seq += content[i].replace('\n','')
			    i += 1

			# To get the next mirna header
			i -= 1

			# mirna sequence
			mirna['seq'] = tmp_seq.upper().strip()
			mirna['length'] = len(tmp_seq)

			# ======================================================
			# Modify all sequences according to each specified
			# variation
			# ======================================================
			# Nt substitution
			for s in variations['substitution']:
			    s_seq = self.sequenceNucleotideSubstitution(tmp_seq.upper(),
				                                        s['position'],
				                                        s['ref'],
				                                        s['mod'])
			    if s_seq['status']:
				mirna['variations']['substitution'].append({
				    'position': s['position'],
				    'ref': s['ref'],
				    'mod': s['mod'],
				    'seq': s_seq['sequence'],
				    'output_path': '',
				    'identifier': '',
				    'length': len(s_seq['sequence'])
				})
			    else:
				mirna['variations']['error'].append({
				    'variation': s['variation'],
				    'msg': s_seq['status_msg']
				})

			# isomiR
			for i in variations['isomir']:
			    i_seq = self.getIsomirSequence(m_id,
				                           tmp_seq.upper(),
				                           i['start'],
				                           i['end'],
				                           p_id)
			    if i_seq['status']:
				mirna['variations']['isomir'].append({
				    'start': i['start'],
				    'end': i['end'],
				    '5p': i['5p'],
				    '3p': i['3p'],
				    'seq': i_seq['sequence'],
				    'output_path': '',
				    'identifier': '',
				    'length': len(i_seq['sequence'])
				})
				mirna['pre_mirna_seq'] = i_seq['pre_sequence']
				mirna['pre_mirna_length'] = len(i_seq['pre_sequence'])
			    else:
				mirna['variations']['error'].append({
				    'variation': i['variation'],
				    'msg': i_seq['status_msg']
				})

			# Get the longest mirna sequence
			if self.longest_mirna_seq < len(tmp_seq):
			    self.longest_mirna_seq = len(tmp_seq)

			# Store the mirna details
			if mirna['identifier'] != '':
			    # Update the number of duplicates
			    mirna['duplicate'] = mirna_dir_indexes_list[m_dir]['idx'] - 1
			    tmp_mirnas_list.append(mirna)

        # Descending sorting by sequence length
        self.mirnas_list = sorted(tmp_mirnas_list,
                                  key=itemgetter('length'),
                                  reverse=True)

        # No. of mirnas
        self.no_mirnas = len(self.mirnas_list)

    def loadMiRNAs(self, mirna_path):
        if os.path.isfile(mirna_path):
	    self.mirna_path = mirna_path
	    with open(self.mirna_path, 'r') as f:
		self.mirnas_list = json.load(f)
		self.no_mirnas = len(self.mirnas_list)

    def setMiRNAInputFiles(self, path, identifier, sequence):
        if os.path.isdir(path):
            # miRNA FASTA file
            mirna_fasta_path = '%s/mirna.fa' % path

            with open(mirna_fasta_path, 'w') as f:
                f.write('>%s\n' % identifier)
                f.write(sequence+'\n')

            # TargetScan miRNA file
            mirna_fasta_path = '%s/mirna_targetscan.txt' % path

            with open(mirna_fasta_path, 'w') as f:
                seed = sequence[1:8]

                # Check U-T nucleotides
                mirna_u = seed.find('U')

                # Replace miRNA's T with U
                if mirna_u == -1:
                    seed = seed.replace('T', 'U')
                seed = seed.upper()

                # Default
                species_id = '9606'

                # If the seed exists into miR Family Info
                if seed in self.mirna_family_info:
                    species_id = ';'.join(self.mirna_family_info[seed])

                line = '%s\t%s\t%s' % (identifier.replace('hsa-', ''),
                                       seed,
                                       species_id)
                f.write(line)

    def setMiRNAsOutputDirectories(self, path):
        if not os.path.isdir(path):
            self.ePrintClose('Error. No output directory found')

        for i in xrange(0, len(self.mirnas_list)):
            # miRNA output directory path
            m_mirna_dirpath = '%s/prediction/%s__%s' % (path,
	                                                self.mirnas_list[i]['identifier'],
	                                                self.mirnas_list[i]['pre_mirna_id'])
            # ==================================================================
            # ==================================================================
            # WildType
            # ==================================================================
            # ==================================================================
            if os.path.exists(m_mirna_dirpath):
                # Delete the miRNA existing output directory
                self.removeDirectory(m_mirna_dirpath)

            # Create the miRNA output directory
            self.createDirectory(m_mirna_dirpath)
            self.mirnas_list[i]['output_path'] = m_mirna_dirpath

            self.setMiRNAInputFiles(m_mirna_dirpath,
                                    self.mirnas_list[i]['identifier'],
                                    self.mirnas_list[i]['seq'])

            # ==================================================================
            # ==================================================================
            # Variations
            # ==================================================================
            # ==================================================================

            # Nucleotide substitution
            for j in xrange(0, len(self.mirnas_list[i]['variations']['substitution'])):
                mir = self.mirnas_list[i]['variations']['substitution'][j]
                variation_path = '%s__%s_%s_%s' % (self.mirnas_list[i]['identifier'],
                                                   str(mir['position']),
                                                   mir['ref'],
                                                   mir['mod'])

                mirna_dirpath = '%s/prediction/%s__%s/nt_sub/%s' % (path,
		                                                    self.mirnas_list[i]['identifier'],
		                                                    self.mirnas_list[i]['pre_mirna_id'],
		                                                    variation_path)

                # Create the miRNA variation output directory
                self.createDirectory(mirna_dirpath)

                self.setMiRNAInputFiles(mirna_dirpath,
                                        self.mirnas_list[i]['identifier'],
                                        mir['seq'])

                self.mirnas_list[i]['variations']['substitution'][j]['output_path'] = mirna_dirpath
                self.mirnas_list[i]['variations']['substitution'][j]['identifier'] = variation_path

            # isomiR
            for j in xrange(0, len(self.mirnas_list[i]['variations']['isomir'])):
                mir = self.mirnas_list[i]['variations']['isomir'][j]
                variation_path = '%s__%s_%s' % (self.mirnas_list[i]['identifier'],
                                                mir['5p'],
                                                mir['3p'])

                mirna_dirpath = '%s/prediction/%s__%s/isomir/%s' % (path,
		                                                    self.mirnas_list[i]['identifier'],
		                                                    self.mirnas_list[i]['pre_mirna_id'],
		                                                    variation_path)

                # Create the miRNA variation output directory
                self.createDirectory(mirna_dirpath)

                self.setMiRNAInputFiles(mirna_dirpath,
                                        self.mirnas_list[i]['identifier'],
                                        mir['seq'])

                self.mirnas_list[i]['variations']['isomir'][j]['output_path'] = mirna_dirpath
                self.mirnas_list[i]['variations']['isomir'][j]['identifier'] = variation_path

            # miRNA FASTA check
            mirna_fasta_parsing_path = '%s/mirna_check.json' % m_mirna_dirpath
            with open(mirna_fasta_parsing_path, 'w') as f:
                f.write(json.dumps(self.mirnas_list[i], indent=4))

    def printMiRNAsDetails(self):
        for t in self.mirnas_list:
            print t


