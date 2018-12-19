#!/usr/bin/env python3.6

import software.uniprotEntry as uE
import software.gene as Gene

__author__  = "Jens Loers"
__email__   = "jens.loers@uni-bielefeld.de"

class UnifiedFileHandler:
	"""
	Class to load a uniprot entry in GFF Format and store it in a dictionary
	of uniprot Entry objects
	"""
	
	def __init__(self):
		"""
		inititalize laoding of a GFF File. Returns an Dictonary of UniportEntry objcets
		"""
		
	def loadUniprotGFF(self, path_to_file):
		"""
		Function to load uniprot information in GFF Format
		"""
		self.allEntrys = {}
		with open(path_to_file, "r") as DataFile:
			entrys = DataFile.read().split('##sequence-region ')
			for entry in entrys[1:]:
				
				# evaluate first line of entry
				lines = entry.split('\n')
				uniID = lines[0].split(' ')[0]
				start = lines[0].split(' ')[1]
				end   = lines[0].split(' ')[2]

				## create a new uniprot object
				newEntry = uE.UniprotEntry(uniID, [int(start), int(end)])
				
				for line in lines[1:]:
					columns 	= line.split('\t')
					function	= []
					if len(columns) > 2:
						try:
								function.append(columns[2])

								function.append(int(columns[3]))
								function.append(int(columns[4]))
								function.append((columns[8]))
						except:
							#should add a logfile to store those conditions somewere (TODO)
							continue
					if len(function) != 0:						
						newEntry.functions.append(function)
				self.allEntrys[uniID] = newEntry
				
		return self.allEntrys
			
	def loadAnnotationGFF(self, path_to_file):
		"""
		Funtion to load Annotation in GFF Format. Might not work for 
		all different GFF files, because they are not that unified as 
		they should be
		"""
		self.allEntrys = {}
		with open(path_to_file, "r") as DataFile:
			entrys = DataFile.read().split('###')
			# go through GFF line by line
			current_entry 	= None
			current_variant	= None
			init			= True
			for blocks in entrys:
				block = blocks.split('\n')
				for entry in block:
					content = entry.split('\t')
					
					# create a new current_entry, if line contains 'gene'
					if len(content) > 1 and content[2] == 'gene':
						#print('1', content)
						# Attach current entry to allEntrys, if
						if current_entry != None:
							#print('3',current_entry.variants)
							self.allEntrys[current_entry.ID] = current_entry
							current_entry = None
							#exit()
						
						current_entry 				= Gene.Gene(content[8].split(';')[0].split('=')[1])
						#print('2',current_entry.ID)
						#exit()
						current_entry.chromosome 	= content[0] 
						current_entry.start 		= int(content[3])
						current_entry.end 			= int(content[4])
						current_entry.orientation 	= content[6]
						
					# create a new  current_variant, if line contains 'mRNA'
					if len(content) > 1 and content[2] == 'mRNA':
						#print(content)
						# Attach fromer current variant to the current entry
						#if current_variant != None:
						#	current_entry.variants[content[8].split(';')[0].split('=')[1]] = current_variant
						#	current_variant = None
							
						current_variant 			= Gene.Variant(content[8].split(';')[0].split('=')[1])
						current_variant.start		= int(content[3])
						current_variant.end			= int(content[4])
						current_entry.variants[content[8].split(';')[0].split('=')[1]] = current_variant
						
						
					# add  exons to current_variant
					if len(content) > 1 and content[2] == 'exon':
						current_variant.exons.append((int(content[3]), int(content[4])))
		
		return	self.allEntrys	
	
	def loadNavipOutputAlpha(self, path_to_file, genes, name_cultivar):
		"""
		Function to load output from NAVIP
		"""
		with open(path_to_file, "r") as DataFile:
			entrys = DataFile.read().split('\n')
			
			# go though all entrys lines
			current_gene = None
			current_gene_block = []
			
			for entry in entrys[6:]:
				content_all 	= entry.split('\t')
				
				# split up info part of VCF File, were NAVIP output is 
				# stored
				if len(content_all) > 6:
					content_NAVIP1	= content_all[7].split('|')
				
				# initialize current gene
				if current_gene == None:
					current_gene = content_NAVIP1[0].split('.')[0]
					
				# get all relevant entrys for the current gene:
				if current_gene == content_NAVIP1[0].split('.')[0]:
					#print(current_gene, 'IF', content_NAVIP1)
					# append and sample all VCF entrys for current gene
					if content_all not in current_gene_block and len(content_NAVIP1[2].split(',')) > 1:
						current_gene_block.append(content_all)
				else:
					#print(current_gene, 'Else', len(current_gene_block), content_NAVIP1)
					# exit for end of file
					if len(content_all) <= 1:
						continue
					# create new cultivar instance
					newCultivar = Gene.Cultivar(name_cultivar)
					# connect cultivar instance to corresponding gene instance
					genes[current_gene].cultivars.append(newCultivar)
					# identify redundant/different entrys (create ONE entry for all varaints that show same behavior
					nrE = getPositionEntries(current_gene_block)
					# create VariantPostionInstances and connect them to cultivar Instances
					for entry3 in nrE:
						positionEntries = []
						currentPosition = []
						a = 0
						for entry2 in entry3:
							content_NAVIP	= entry2[7].split('|')
							#print(content_NAVIP)
							# create nw variant Position Object
							newVariantPosition = Gene.VariationPosition(entry2[1])

							# write information to Variant Position Instance
							newVariantPosition.ref_amino_acid	= content_NAVIP[5].upper()
							newVariantPosition.alt_amino_acid	= content_NAVIP[7].upper()
							newVariantPosition.ref_base 		= entry2[3]
							newVariantPosition.alt_base			= entry2[4]
							newVariantPosition.pos_amino_acid 	= int(content_NAVIP[8])#int((int(content_NAVIP[8])/3.)+0.8)#int(((int(content_NAVIP[8]))/3.)+0.8)#int(round((int(content_NAVIP[8])+0.2)/3.,0))
							# append all types NAVIP put out (SUB, 
							# amino acid change, ...s
							[newVariantPosition.types.append(x) for x in content_NAVIP[2].split(',')]
							newVariantPosition.variant=content_NAVIP[0].split('.')[0]+'.'+content_NAVIP[0].split('.')[1].replace('A', '').replace('B', '').replace('C', '').replace('D', '').replace('E', '')
							newVariantPosition.variant_names.append(genes[content_NAVIP[0].split('.')[0]].variants[content_NAVIP[0].split('.')[0]+'.'+content_NAVIP[0].split('.')[1].replace('A', '').replace('B', '').replace('C', '').replace('D', '').replace('E', '')])
							positionEntries.append(newVariantPosition)
							currentPosition = entry2[1]
							del newVariantPosition

						newCultivar.variablePositions[int(currentPosition)] = positionEntries


					# reassign current gene and current gene Block
					current_gene 		= content_NAVIP1[0].split('.')[0]
					current_gene_block  =  []
					# assign to current gene block if criteria is 
					# fullfilled
					if len(content_NAVIP1[2].split(',')) > 1:
						current_gene_block.append(content_all)

	def loadNavipOutputBeta(self, path_to_file, genes, name_cultivar):
		"""
		Function to load output from NAVIP
		"""
		with open(path_to_file, "r") as DataFile:
			entrys = DataFile.read().split('\n')
			
			# go though all entrys lines
			current_gene = None
			current_gene_block = []
			
			for entry in entrys[6:]:
				content_all 	= entry.split('\t')
				#print(content_all)
				# split up info part of VCF File, were NAVIP output is 
				# stored
				if len(content_all) > 6:
					content_NAVIP1	= content_all[7].split('|')
				
				# initialize current gene
				if current_gene == None:
					current_gene = content_NAVIP1[0].split('.')[0]
					
				# get all relevant entrys for the current gene:
				if current_gene == content_NAVIP1[0].split('.')[0]:
					#print(current_gene, 'IF', content_NAVIP1)
					# append and sample all VCF entrys for current gene
					if content_all not in current_gene_block and len(content_NAVIP1[2].split(',')) > 1:
						current_gene_block.append(content_all)
				else:
					#print(current_gene, 'Else', len(current_gene_block), content_NAVIP1)
					# exit for end of file
					if len(content_all) <= 1:
						continue
					# create new cultivar instance
					newCultivar = Gene.Cultivar(name_cultivar)
					# connect cultivar instance to corresponding gene instance
					genes[current_gene].cultivars.append(newCultivar)
					# identify redundant/different entrys (create ONE entry for all varaints that show same behavior
					nrE = getPositionEntries(current_gene_block)
					# create VariantPostionInstances and connect them to cultivar Instances
					for entry3 in nrE:
						positionEntries = []
						currentPosition = []
						a = 0
						for entry2 in entry3:
							content_NAVIP	= entry2[7].split('|')
							#print(content_NAVIP)
							# create nw variant Position Object
							newVariantPosition = Gene.VariationPosition(entry2[1])
							
							# write information to Variant Position Instance
							newVariantPosition.ref_amino_acid	= content_NAVIP[5].upper()
							newVariantPosition.alt_amino_acid	= content_NAVIP[8].upper()
							newVariantPosition.ref_base 		= entry2[4]
							newVariantPosition.alt_base			= entry2[7]
							newVariantPosition.pos_amino_acid 	= int(content_NAVIP[6])
							
							# load genotype information if possible
							abreviation  = entry2[-2].split(':')
							tmpvalues	 = entry2[-1].split(':')
							#print(abreviation, entry2)
							
							genotype = {}
							for index in range(len(abreviation)-(len(abreviation)-len(tmpvalues))):
								try:
									temporal	= tmpvalues[index].split(',')
									temporal2	= [int(x) for x in temporal]
									genotype[abreviation[index]] = temporal2
								except:
									genotype[abreviation[index]] = tmpvalues[index].split(',')
							
							newVariantPosition.genotype = genotype
														
							# append all types NAVIP put out (SUB, amino acid change, ...s)
							[newVariantPosition.types.append(x) for x in content_NAVIP[2].split(',')]
							newVariantPosition.variant=content_NAVIP[0].split('.')[0]+'.'+content_NAVIP[0].split('.')[1].replace('A', '').replace('B', '').replace('C', '').replace('D', '').replace('E', '')
							newVariantPosition.variant_names.append(genes[content_NAVIP[0].split('.')[0]].variants[content_NAVIP[0].split('.')[0]+'.'+content_NAVIP[0].split('.')[1].replace('A', '').replace('B', '').replace('C', '').replace('D', '').replace('E', '')])
							positionEntries.append(newVariantPosition)
							currentPosition = entry2[1]
							del newVariantPosition

						newCultivar.variablePositions[int(currentPosition)] = positionEntries


					# reassign current gene and current gene Block
					current_gene 		= content_NAVIP1[0].split('.')[0]
					current_gene_block  =  []
					# assign to current gene block if criteria is 
					# fullfilled
					if len(content_NAVIP1[2].split(',')) > 1:
						current_gene_block.append(content_all)

	def loadSecondaryIDs(self, path_to_file, genes):
		"""
		Function to load secondary ids. This might be necessary to link
		the information from pdb database to the used gene idenfier
		"""
		with open(path_to_file, "r") as DataFile:
			entrys = DataFile.read().split('\n')
			for entry in entrys:
				tmp 	= entry.split('\t')
				tmp2	= tmp[0].split(',')
				for i in tmp2:
					try:
						genes[i].secondaryID.append(tmp[1])
					except:
						continue

	def loadAnnotationMapping(self, path_to_file, genes):
		"""
		Function to load an annotation mapping. This might be necessary
		to improve the quantity and quality of a given database.
		File to load should contain this format (tab sep.):
		
		ID				PosAnnot.		PosUniprot		Amount AnnPos mapped
		AT5G66150.1 	1,1047			1,1047			1.0
		AT5G66150.2 	1,474,475,1046	1,474,476,1047	1.0
		AT5G66850.1 	1,716			1,716			1.0
		AT5G66850.2 	1,658			1,658			0.982
		
		"""
		with open(path_to_file, "r") as DataFile:
			entrys = DataFile.read().split('\n')
			for entry in entrys:
				tmp 	= entry.split('\t')
				if len(tmp) < 2  or tmp[1] == '':
					continue
				geneID	= tmp[0].split('.')[0]
				Anot	= [int(x) for x in tmp[1].split(',')]
				Unip    = [int(x) for x in tmp[2].split(',')]
				
				tpR1, tpR2	= [], []
				
				for i in range(0,len(Anot), 2):
					tpR1.append([Anot[i], Anot[i+1]])
				for j in range(0,len(Unip), 2):
					tpR2.append([Unip[j], Unip[j+1]])	
				
				try:
					genes[geneID].variants[tmp[0]].AnnotationAaPos 		= tpR1
					genes[geneID].variants[tmp[0]].AnnotationDbPos 		= tpR2
					genes[geneID].variants[tmp[0]].AnnotationIdentiti	= round(float(tmp[3]),3)
				except:
					continue

	def loadSubsitutionMatrix(self, path_to_file):
		"""
		This module can load an amino acid substitution matrix of 
		following form:
			-	A	B	C	...
		-	0	0.1	0.1	0
		A	0	0.1	0.1	0
		B	0	0.1	0.1	0
		C	0	0.1	0.1	0
		...
		
		We store values in a dictionary with a string of length 2 as 
		key:
		{AB: 0.1, ...}
		AB indicates a change from A to B while 
		BA indicates a change from B to A.
		
		Return: Matrix as dictionary + name of Matrix (usefull to 
		connect it with the uniport entry)
		"""
		content = {}
		with open(path_to_file, "r") as DataFile:
			entries	= DataFile.read().split('\n')
			header	=  entries[0].split('\t')
			for entry in entries[1:]:
				tmp = entry.split('\t')				
				for i in range(1,len(tmp)):
					if tmp[0].replace(' ','') == header[i].replace(' ',''):
						continue
					content[(tmp[0]+header[i]).replace(' ','')] = float(tmp[i])
		return content, path_to_file.split('=')[1]	

	def loadProteinLength(self, path_to_file, genes):
		with open(path_to_file, "r") as DataFile:
			entrys = DataFile.read().split('\n')
			for entry in entrys:
				tmp 		= entry.split(' ')
				geneName	= tmp[0].split('.')[0]
				if len(tmp) > 1:
					try:
						genes[geneName].variants[tmp[0]].proteinLength = int(tmp[1])
					except:
						continue

def getPositionEntries(listOfEnties):
	"""
	Returns a list of non_redundant entrys and a list of all variants 
	which belong to the entry. This method was decoupled from 
	loadNavipOutput to make it more readable
	"""
	nrE			= []	# number of non redundant entries
	numbering	= []
	for i in listOfEnties:
			if i[1] not in numbering:
				nrE.append([i])
				numbering.append(i[1])
			else:
				nrE[numbering.index(i[1])].append(i)
	#print(nrE)
	return nrE
