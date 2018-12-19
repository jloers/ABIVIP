#!/usr/bin/env python3.6

__author__  = "Jens Loers"
__email__   = "jens.loers@uni-bielefeld.de"

class Analysis:
	"""
	Class used to store all analyis and evaluation routines to make the 
	code more readable 
	"""	
	def getAffectedFunctionalities(self, genes, database):
		"""
		Funktion to compare positions of interest in a gene  object 
		with information stored in an uniprot database object
		"""

		### analyse SNP genevise ###  
		for gene in genes:
			## preprocessing to get correct key for uniprot access ##
			
			# look up wether an db entry exists TODO: check also sec. IDs)
			functional_id = None
			for i in genes[gene].secondaryID:
				try:
					if  i in database:
						functional_id = i
						break
				except:
					continue
					
			# ignore entrys with insufficient mapping\no uniport information
			if functional_id == None:
				continue
				
			## start anlysis	##
			"""
			construct to interate through data strcuture
			"""
			
			#print(gene)
			# go through every cultivar
			for cultivar in genes[gene].cultivars:
				
				# access variantPosition objects in cultivar, store in VarPos
				for varPos in cultivar.variablePositions:
					# Go through all entrys for a postition (can be because of ploidity or variant)
					for i in cultivar.variablePositions[varPos]:
						
						## here we test whether the corresponding position of the cultivar has a valid genotype information.
						## also just ref coverages are skipped.
						if i.genotype['GT'] == ['./.'] or i.genotype['GT'] == ['0/0']:
							continue
						
						shift = 0
						# skip effectless positions
						if 'Hits before Start into CDS' in i.types:
							continue
						
						## create switch: Distinguish between 3 different Variant Types: INDEL, AASUBS, STOPGAIN/LOSS
						
						if 'Amino acid change' in i.types and ('DEL' in i.types or 'INS' in i.types):
							frameshift = [x for x in i.types if x.startswith('Frame')]
							if len(frameshift)> 0:
								shift += (int(frameshift[0].split('Frameshift')[1]))
								adaptedAaPosition = adaptAAPositionToMapping(i.pos_amino_acid+shift, i.variant_names[0].AnnotationAaPos, i.variant_names[0].AnnotationDbPos)
								if adaptedAaPosition == None:
									continue
								else:
									i.pos_amio_acid_adapted = adaptedAaPosition
								# go through DB functionalities	
								for function in database[functional_id].functions:
									if adaptedAaPosition <= function[1] and function[0] not in ['Non-terminal residue','Site','Chain','Repeat', 'Compositional bias', 'Topological domain', 'Initiator methionine', 'Sequence conflict', 'Natural variant', 'Alternative sequence', 'Peptide','Signal peptide','Transit peptide','Propeptide']:# and function[0] in ['Metal binding']:
										i.function.append(function)
									# covers the case for peptides which could have some restriction sites which go further that the end of the peptide site
									if adaptedAaPosition <= function[1]+3 and function[0] in ['Peptide','Signal peptide','Transit peptide','Propeptide']:
										i.function.append(function)
						
								continue
							else:
								# adapt affected amino acid position
								
								# Deal with insertions
								if 'INS' in i.types:
									adaptedAaPosition = adaptAAPositionToMapping(i.pos_amino_acid+shift, i.variant_names[0].AnnotationAaPos, i.variant_names[0].AnnotationDbPos)
									if adaptedAaPosition == None:
										continue
									else:
										i.pos_amio_acid_adapted = adaptedAaPosition
									for function in database[functional_id].functions:
										if adaptedAaPosition >= function[1] and adaptedAaPosition <= function[2] and function[0] not in ['Non-terminal residue','Site','Chain','Repeat', 'Compositional bias', 'Topological domain', 'Initiator methionine', 'Sequence conflict', 'Natural variant', 'Alternative sequence', 'Disulfide bond', 'Peptide','Signal peptide','Transit peptide','Propeptide']:
											i.function.append(function)
										if ((adaptedAaPosition >= function[1] and adaptedAaPosition <= function[1]) or (adaptedAaPosition >= function[2] and adaptedAaPosition <= function[2])) and function[0] == 'Disulfide bond':
											i.function.append(function)
										if adaptedAaPosition >= function[1] and adaptedAaPosition <= function[2]+3 and function[0] in ['Peptide','Signal peptide','Transit peptide','Propeptide']:
											i.function.append(function)
																						
								# Deal with deletions			
								if 'DEL' in i.types:
									aP1 = adaptAAPositionToMapping(i.pos_amino_acid+shift, i.variant_names[0].AnnotationAaPos, i.variant_names[0].AnnotationDbPos)
									if i.ref_amino_acid[0] == i.alt_amino_acid or i.ref_amino_acid[-1] == i.alt_amino_acid:
										aP2 = adaptAAPositionToMapping((i.pos_amino_acid+len(i.ref_base)-1+shift), i.variant_names[0].AnnotationAaPos, i.variant_names[0].AnnotationDbPos)
									else:
										aP2 = adaptAAPositionToMapping((i.pos_amino_acid+len(i.ref_base)+shift), i.variant_names[0].AnnotationAaPos, i.variant_names[0].AnnotationDbPos)
									#print(i.types, aP1, aP2, i.pos_amino_acid+shift, i.pos_amino_acid+len(i.ref_base)-1+shift, i.variant_names[0].AnnotationAaPos, i.variant_names[0].AnnotationDbPos, i.ref_base, i.alt_base)
									if aP1 == None or aP2 == None:
										continue
									else:
										i.delStart  = aP1
										i.delEnd	= aP2
										
									for function in database[functional_id].functions:
										if ((function[1] >= aP1 and function[1] <= aP2) or (function[2] >= aP1 and function[2] <= aP2) or (aP1 >= function[1] and aP2 <=function[1]) or (function[1] >= aP1 and function[2] <= aP2)) and function[0] not in ['Non-terminal residue','Site','Chain','Repeat', 'Compositional bias', 'Topological domain', 'Initiator methionine', 'Sequence conflict', 'Natural variant', 'Alternative sequence', 'Disulfide bond','Peptide','Signal peptide','Transit peptide','Propeptide']:
											i.function.append(function)
										if (((function[1] >= aP1 and function[1] <= aP2) or (function[1] >= aP1 and function[1] <= aP2) or (aP1 >= function[1] and aP2 <=function[1]) or (function[1] >= aP1 and function[1] <= aP2)) or ((function[2] >= aP1 and function[2] <= aP2) or (function[2] >= aP1 and function[2] <= aP2) or (aP1 >= function[2] and aP2 <=function[2]) or (function[2] >= aP1 and function[2] <= aP2))) and function[0] == 'Disulfide bond':
											i.function.append(function)
										if ((function[1] >= aP1 and function[1] <= aP2) or (function[2]+3 >= aP1 and function[2]+3 <= aP2) or (aP1 >= function[1] and aP2 <=function[1]) or (function[1] >= aP1 and function[2]+3 <= aP2)) and function[0] in ['Peptide','Signal peptide','Transit peptide','Propeptide']:
											i.function.append(function)
											
						if 'Amino acid change' in i.types and 'DEL' not in i.types and 'INS' not in i.types:

							# adapt affected amino acid position to numbering
							adaptedAaPosition = adaptAAPositionToMapping(i.pos_amino_acid+shift, i.variant_names[0].AnnotationAaPos, i.variant_names[0].AnnotationDbPos)
							if adaptedAaPosition == None:
								continue
							else:
								i.pos_amio_acid_adapted = adaptedAaPosition
							# go through DB functionalities	
							for function in database[functional_id].functions:
								if adaptedAaPosition >= function[1] and adaptedAaPosition <= function[2] and function[0] not in ['Non-terminal residue','Site','Chain','Repeat', 'Compositional bias', 'Topological domain', 'Initiator methionine', 'Sequence conflict', 'Natural variant', 'Alternative sequence', 'Disulfide bond','Peptide','Signal peptide','Transit peptide','Propeptide']:
									i.function.append(function)
								if ((adaptedAaPosition >= function[1] and adaptedAaPosition <= function[1]) or (adaptedAaPosition >= function[2] and adaptedAaPosition <= function[2])) and function[0] == 'Disulfide bond':
									i.function.append(function)
								if adaptedAaPosition >= function[1] and adaptedAaPosition <= function[2]+3 and function[0] in ['Peptide','Signal peptide','Transit peptide','Propeptide']:
									i.function.append(function)
						
						if 'Stop gained' in i.types:# or 'Stop lost' in i.types:
							adaptedAaPosition = adaptAAPositionToMapping(i.pos_amino_acid+shift, i.variant_names[0].AnnotationAaPos, i.variant_names[0].AnnotationDbPos)
							if adaptedAaPosition == None:
									continue
							else:
								i.pos_amio_acid_adapted = adaptedAaPosition
																
							# go through DB functionalities	
							for function in database[functional_id].functions:
								if adaptedAaPosition <= function[1] and function[0] not in ['Non-terminal residue','Site','Chain','Repeat', 'Compositional bias', 'Topological domain', 'Initiator methionine', 'Sequence conflict', 'Natural variant', 'Alternative sequence']:
									i.function.append(function)
			
	def compareCultivars(self, genes):
		"""
		This function compares the effects of SNPs/INDELS between all
		cultivars
		"""
		import software.stats as stats
		
		newStructure   		= {}
		genesDisrupted 		= {}
		genesFrameShifted	= {}
		struct = False
				
		# iterate through genes
		for gene in genes:
			statEntry = stats.GeneStat(gene, genes[gene])
			# interate through variants
			for variant in genes[gene].variants:
				#print('\tvariant:', variant)
			
			## get all different affected functions
				allAffectedFunctions = []
				stop_gained			 = False
			#iterate through cultivars
				for cultivar in genes[gene].cultivars:
				# iterate throug variant positions
					for VarPos in cultivar.variablePositions:
						for var in cultivar.variablePositions[VarPos]:
							if 'Stop gained' in var.types:
								
								stop_gained = True
								stop_pos	= var.ID
								genesDisrupted[gene] = genes[gene]
							if 'Frameshift+2' in var.types or 'Frameshift+1' in var.types or 'Frameshift-1' in var.types or 'Frameshift-2' in var.types:
								if 'Stop gained' not in var.types:
									#print(var.types)
									stop_gained = True
									stop_pos	= var.ID
									genesFrameShifted[gene] = genes[gene]	
							for funct in var.function:
								if funct not in allAffectedFunctions and len(var.function)>0:
									allAffectedFunctions.append(funct)

				# create GeneVar Object and connect it with VarStat Objects
				#if stop==True
				#print(allAffectedFunctions)
				if len(allAffectedFunctions):	
					statEntry.addVariantStat(variant, genes[gene].variants[variant])
					struct = True
				#if stop_gained == True:
				#	struct = True
					
				for affectedFunction in allAffectedFunctions:
					content  = []
					
					for cultivar in genes[gene].cultivars:
						# iterate throug variant positions
						for VarPos in cultivar.variablePositions:
							for var in cultivar.variablePositions[VarPos]:
								if  var.variant == variant and affectedFunction  in var.function:
									content.append([cultivar.ID, var])
					functionID = gene+'_'+str(affectedFunction[0])+'_'+str(affectedFunction[1])+'_'+str(affectedFunction[2])
					statEntry.variants[variant].addVariantStat(functionID, affectedFunction, content)
			
			if 	struct == True:
				newStructure[gene] =  statEntry
				struct =  False
		
		# cleanup
		for result in newStructure:
			for variant in newStructure[result].variants:
				variantO = (newStructure[result].variants[variant])
				newFunctions = []
				for function in variantO.functions:
					if len(function.varPos) != 0:
						newFunctions.append(function)
				variantO.functions = newFunctions
	
		return newStructure, genesDisrupted, genesFrameShifted
				
	def printComparence(self, genes):
		"""
		This function basically does the same as compareCultivars(),
		with the difference that it can print to console directely
		"""
		# cultivar list: name_cult, 
		
		# iterate through genes
		for gene in genes:
			print('gene:', gene)
			# interate through variants
			for variant in genes[gene].variants:
				print('\tvariant:', variant)
			
			## get all different affected functions
				allAffectedFunctions = []
			#iterate through cultivars
				for cultivar in genes[gene].cultivars:
				# iterate throug variant positions
					for VarPos in cultivar.variablePositions:
						for var in cultivar.variablePositions[VarPos]:
							for funct in var.function:
								if funct not in allAffectedFunctions and len(var.function)>0:
									allAffectedFunctions.append(funct)

					
				
				for affectedFunction in allAffectedFunctions:
					print('\t\taffected function: ',affectedFunction)
					for cultivar in genes[gene].cultivars:
						# iterate throug variant positions
						for VarPos in cultivar.variablePositions:
							for var in cultivar.variablePositions[VarPos]:
								#print(var.variant_names[0].ID, variant)
								#if var.variant_names[0].ID == variant and affectedFunction  in var.function:
								if  var.variant == variant and affectedFunction  in var.function:
									print('\t\t\t',cultivar.ID, VarPos,adaptAAPositionToMapping(var.pos_amino_acid, var.variant_names[0].AnnotationAaPos, var.variant_names[0].AnnotationDbPos), var.ref_amino_acid, var.alt_amino_acid, var.types, )
								if gene == 	'AT1G03300':
									print('NOQ:')
									print(cultivar.ID, VarPos, var.function, [x.ID for x in var.variant_names], var.variant)	

	def secondaryStructureAnalysis(self, comparenceResult, matrices, cultivars):
		"""
		Analyses the influence of SNP, Insertions and Deletions on the
		secondary structure elementes of a protein
		"""
		# import some mathematical opertations
		from math import log
		
		### region for all controllable variables
		# TODO
		stepLength				= 0.5
		stepStart				= 1		
		weightPrimaryTranscript	= 0.8
		
		###
		
		# get all structural element names out of matrices
		elements = [x for  x in matrices]
		
		# safe maxima of every matrix in dict
		maxima   = {}
		for element in elements:
			maxima[element] = getMaxValueOfMatrix(matrices[element])

		for result in comparenceResult:
			transcript_diff_scores	= {}
			transcript_max_scores	= {}
			
			for variant in comparenceResult[result].variants:
				variantO = (comparenceResult[result].variants[variant])
				# store function overlapping scorrings results per cultivar
				allFunctionsScores = {}
				for x in cultivars:
					allFunctionsScores[x] = []
				
				diffScore	=	0
				maxScore	=	0
				
				for function in variantO.functions:
					if function.properties[0] in elements:
						#print(function.properties[0])
						function_scores = {}		# samples all scores for variants in a function per cultivar
						# initialize dictionary to store culitvar information
						for x in cultivars:
							function_scores[x] = []
						# iterate throug variable positions
						for pos in function.varPos:
							
							# section for single amino acid exchanges. Different Approaches for different Variant Types are necessary
							if 'INS' not in pos[1].types and 'DEL' not in pos[1].types and 'Stop gained' not in pos[1].types:
								#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, round((1-(log(1+(matrices[function.properties[0]][pos[1].ref_amino_acid+pos[1].alt_amino_acid])))/log(maxima[function.properties[0]])),3))
								try:
									pos[1].secondaryScore = round((1-(log(1+(matrices[function.properties[0]][pos[1].ref_amino_acid+pos[1].alt_amino_acid])))/log(maxima[function.properties[0]])),3)
									function_scores[pos[0]].append(pos[1].secondaryScore)
								except:
									print('Error', pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties[0])
									continue
									
							if 'INS' in pos[1].types and 'Stop gained' not in pos[1].types and 'Frameshift+1' not in pos[1].types and 'Frameshift+2' not in pos[1].types and len(pos[1].alt_amino_acid)>1:
								insertion = pos[1].alt_amino_acid
								# corrrect insertion for 'replace site by itself'
								if insertion[0] == pos[1].ref_amino_acid:
									insertion = insertion[1:]
								elif insertion[-1] == pos[1].ref_amino_acid:
									insertion = insertion[:1]
								if len(insertion) <= 3:
									# calculate maximum influencing amino acid as score for the position	
									maxValue = 0	
									for a in insertion:
										tmpValue =  round((1-(log(1+(matrices[function.properties[0]]['-'+a])))/log(maxima[function.properties[0]])),3)
										if tmpValue > maxValue:
											maxValue = tmpValue
									pos[1].secondaryScore = maxValue
									function_scores[pos[0]].append(pos[1].secondaryScore)
									#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].secondaryScore )
								else:
									pos[1].secondaryScore = 1.
									function_scores[pos[0]].append(pos[1].secondaryScore)
									#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].secondaryScore )
							
							if 'DEL' in pos[1].types and 'Stop gained' not in pos[1].types and len(pos[1].ref_amino_acid)>1:
								if pos[1].delStart != None and pos[1].delEnd != None:
									deletion = pos[1].ref_amino_acid
									# corrrect insertion for 'replace site by itself'
									if deletion[0] == pos[1].alt_amino_acid:
										deletion = deletion[1:]
									elif deletion[-1] == pos[1].alt_amino_acid:
										deletion = deletion[:1]
										
									if len(deletion) <= 3:
										# calculate maximum influencing amino acid as score for the position	
										maxValue = 0	
										for a in deletion:
											tmpValue =  round((1-(log(1+(matrices[function.properties[0]][a+'-'])))/log(maxima[function.properties[0]])),3)
											if tmpValue > maxValue:
												maxValue = tmpValue
										pos[1].secondaryScore = maxValue
										function_scores[pos[0]].append(pos[1].secondaryScore)
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].secondaryScore )
									else:
										pos[1].secondaryScore = 1.
										function_scores[pos[0]].append(pos[1].secondaryScore)
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].secondaryScore )

									#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, pos[1].delStart, pos[1].delEnd, deletion )
								

						#print(sorted(function_scores))
						fs = {}
						for i in function_scores:
							fs[i] =  stepFunction(sorted(function_scores[i], reverse = True), stepStart, stepLength)
							allFunctionsScores[i].append(fs[i])
						function.secondaryScore = fs
				# evaluation part for the variants:
				#print(allFunctionsScores)
				max1 = []	
				for i in allFunctionsScores:
					max1.append(stepFunction(allFunctionsScores[i],stepStart,stepLength))				
				
				# calculate diffScore and maxScore
				diffScore = maxPairwiseDistance(max1)
				maxScore  = max(max1)	
				variantO.diffScore	=	diffScore
				transcript_diff_scores[variantO.name]	= diffScore
				variantO.maxScore	=	maxScore
				transcript_max_scores[variantO.name]	= maxScore
				
			# create Overall property score (weight representative Transcirpt with 0.8, rest distributed over 0.2)
			propertyScoreMax 	= 0
			propertyScoreDiff	= 0
			
			
			for i in transcript_diff_scores:
				if i.split('.')[1] == '1':
					if len(transcript_diff_scores) == 1:
						propertyScoreMax += 1*transcript_diff_scores[i]
					else:
						propertyScoreMax += 0.8*transcript_diff_scores[i]
				else:
					if len(transcript_diff_scores) > 1:
						propertyScoreMax += (0.2/(len(transcript_diff_scores)-1))*transcript_diff_scores[i]
			
			for i in transcript_max_scores:
				if i.split('.')[1] == '1':
					#print(i.split('.')[1])
					if len(transcript_max_scores) == 1:
						propertyScoreMax += 1*transcript_max_scores[i]
					else:
						propertyScoreMax += weightPrimaryTranscript*transcript_max_scores[i]
				else:
					if len(transcript_max_scores) > 1:
						propertyScoreMax += ((1-weightPrimaryTranscript)/(len(transcript_max_scores)-1))*transcript_max_scores[i]
			
			# store scores in comparenceResult
			comparenceResult[result].diffScoreSecondary	= round(propertyScoreDiff,3)
			comparenceResult[result].maxScoreSecondary	= round(propertyScoreMax,3)

	def membraneAnalysis(self, comparenceResult, matrices, cultivars):
		"""
		Analyses the influence of SNP, Insertions and Deletions on the
		membrane or membrane associated regions of a protein
		"""
		# import some mathematical opertations
		from math import log
		
		### region for all controllable variables
		# TODO
		stepLength				= 0.5
		stepStart				= 1		
		weightPrimaryTranscript	= 0.8
		
		###
		
		# get all structural element names out of matrices
		elements = [x for  x in matrices]

		# safe maxima of every matrix in dict
		maxima   = {}
		for element in elements:
			maxima[element] = getMaxValueOfMatrix(matrices[element])

		for result in comparenceResult:
			transcript_diff_scores	= {}
			transcript_max_scores	= {}
			
			for variant in comparenceResult[result].variants:
				variantO = (comparenceResult[result].variants[variant])
				# store function overlapping scorrings results per cultivar
				allFunctionsScores = {}
				for x in cultivars:
					allFunctionsScores[x] = []
				
				diffScore	=	0
				maxScore	=	0
				
				for function in variantO.functions:
					if function.properties[0] in ['Transmembrane','Intramembrane']:
						#print(function.properties[0])
						function_scores = {}		# samples all scores for variants in a function per cultivar
						# initialize dictionary to store culitvar information
						for x in cultivars:
							function_scores[x] = []
						# iterate throug variable positions
						for pos in function.varPos:
							
							# section for single amino acid exchanges. Different Approaches for different Variant Types are necessary
							if 'INS' not in pos[1].types and 'DEL' not in pos[1].types and 'Stop gained' not in pos[1].types:
								try:
									pos[1].membraneScore = round((1-(log(1+(matrices['membrane'][pos[1].ref_amino_acid+pos[1].alt_amino_acid])))/log(maxima['membrane'])),3)
									function_scores[pos[0]].append(pos[1].membraneScore)
									#print(pos[0], pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amino_acid, pos[1].pos_amio_acid_adapted,function.properties[1:3], pos[1].membraneScore)
								except:
									print('Error', pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties[0] )
									continue
							if 'INS' in pos[1].types and 'Stop gained' not in pos[1].types and len(pos[1].alt_amino_acid)>1:
								insertion = pos[1].alt_amino_acid
								# corrrect insertion for 'replace site by itself'
								if insertion[0] == pos[1].ref_amino_acid:
									insertion = insertion[1:]
								elif insertion[-1] == pos[1].ref_amino_acid:
									insertion = insertion[:1]
									
								if len(insertion) <= 3:
									# calculate maximum influencing amino acid as score for the position	
									maxValue = 0	
									for a in insertion:
										tmpValue =  round((1-(log(1+(matrices['membrane']['-'+a])))/log(maxima['membrane'])),3)
										if tmpValue > maxValue:
											maxValue = tmpValue
									pos[1].membraneScore = maxValue
									function_scores[pos[0]].append(pos[1].membraneScore)
									#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].membraneScore )
								else:
									pos[1].membraneScore = len(insertion)
									function_scores[pos[0]].append(pos[1].membraneScore)
										
						
							if 'DEL' in pos[1].types and 'Stop gained' not in pos[1].types and len(pos[1].ref_amino_acid)>1:
								if pos[1].delStart != None and pos[1].delEnd != None:
									#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, pos[1].delStart, pos[1].delEnd, deletion )
									deletion = pos[1].ref_amino_acid
									# corrrect insertion for 'replace site by itself'
									if deletion[0] == pos[1].alt_amino_acid:
										deletion = deletion[1:]
									elif deletion[-1] == pos[1].alt_amino_acid:
										deletion = deletion[:1]
										
									if len(deletion) <= 3:
										# calculate maximum influencing amino acid as score for the position	
										maxValue = 0	
										for a in deletion:
											tmpValue =  round((1-(log(1+(matrices['membrane'][a+'-'])))/log(maxima['membrane'])),3)
											if tmpValue > maxValue:
												maxValue = tmpValue
										pos[1].membraneScore = maxValue
										function_scores[pos[0]].append(pos[1].membraneScore)
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].peptideScore )
									else:
										pos[1].membraneScore = 1.
										function_scores[pos[0]].append(pos[1].membraneScore)
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].peptideScore )

									#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, pos[1].delStart, pos[1].delEnd, deletion )

						#print(sorted(function_scores))
						fs = {}
						for i in function_scores:
							fs[i] =  stepFunction(sorted(function_scores[i], reverse = True), stepStart, stepLength)
							allFunctionsScores[i].append(fs[i])
						function.membraneScore = fs
				# evaluation part for the variants:
				#print(allFunctionsScores)
				max1 = []	
				for i in allFunctionsScores:
					max1.append(stepFunction(allFunctionsScores[i],stepStart,stepLength))				
				
				# calculate diffScore and maxScore
				diffScore = maxPairwiseDistance(max1)
				maxScore  = max(max1)	
				variantO.diffScore	=	diffScore
				transcript_diff_scores[variantO.name]	= diffScore
				variantO.maxScore	=	maxScore
				transcript_max_scores[variantO.name]	= maxScore
				
			# create Overall property score (weight representative Transcirpt with 0.8, rest distributed over 0.2)
			propertyScoreMax 	= 0
			propertyScoreDiff	= 0
			
			
			for i in transcript_diff_scores:
				if i.split('.')[1] == '1':
					if len(transcript_diff_scores) == 1:
						propertyScoreMax += 1*transcript_diff_scores[i]
					else:
						propertyScoreMax += 0.8*transcript_diff_scores[i]
				else:
					if len(transcript_diff_scores) > 1:
						propertyScoreMax += (0.2/(len(transcript_diff_scores)-1))*transcript_diff_scores[i]
			
			for i in transcript_max_scores:
				if i.split('.')[1] == '1':
					#print(i.split('.')[1])
					if len(transcript_max_scores) == 1:
						propertyScoreMax += 1*transcript_max_scores[i]
					else:
						propertyScoreMax += weightPrimaryTranscript*transcript_max_scores[i]
				else:
					if len(transcript_max_scores) > 1:
						propertyScoreMax += ((1-weightPrimaryTranscript)/(len(transcript_max_scores)-1))*transcript_max_scores[i]
			
			# store scores in comparenceResult
			comparenceResult[result].diffScoreMembrane	= round(propertyScoreDiff,3)
			comparenceResult[result].maxScoreMembrane	= round(propertyScoreMax,3)
			
			#if round(propertyScoreMax,3) >0:
			#	print(round(propertyScoreMax,3))
			o=0
			
	def peptideAnalysis(self, comparenceResult, matrices, cultivars):
		"""
		Analyses the influence of SNP, Insertions and Deletions on the
		secondary structure elementes of a protein
		"""
		# import some mathematical opertations
		from math import log
			
		### region for all controllable variables
		# TODO
		stepLength				= 0.5
		stepStart				= 1		
		weightPrimaryTranscript	= 0.8
		
#		###
		
		# get all structural element names out of matrices
		elements = [x for  x in matrices]

		# safe maxima of every matrix in dict
		maxima   = {}
		for element in elements:
			maxima[element] = getMaxValueOfMatrix(matrices[element])
		
		# Tree Structure stores output information, which then can easily be printed.
		# Collector variables are indicated with the suffix Tree
		printTree = []
		
		for result in comparenceResult:
			
			transcript_diff_scores	= {}
			transcript_max_scores	= {}
			
			# gene, uniprotID, chromosome, number transkripts, transkriptTrees
			geneTree = [result, comparenceResult[result].gene.secondaryID[0], comparenceResult[result].gene.chromosome, len(comparenceResult[result].gene.variants), []]
			#print(geneTree)
			for variant in comparenceResult[result].variants:
				
				# Variant name, functionTree
				#transkriptTree = []
				
				variantO = (comparenceResult[result].variants[variant])
				# store function overlapping scorrings results per cultivar
				allFunctionsScores = {}
				for x in cultivars:
					allFunctionsScores[x] = []
			
				diffScore			=	0
				maxScore			=	0
				functionID 			= 	0
				
				# variant, function
				transkriptTrees	=	[variant, []]
				for function in variantO.functions:
					if function.properties[0] in ['Peptide','Signal peptide','Transit peptide','Propeptide']:
						
						# _Tree structure
						functionID += 1
						cultivarTree = {}							
						
						#print(function.properties)
						function_scores_restictionSides = {}		# samples all scores for variants in a function per cultivar
						function_scores_otherSides		= {}
						
						# initialize dictionary to store culitvar information
						for x in cultivars:
							function_scores_restictionSides[x]	= []
							function_scores_otherSides[x] 		= []
							cultivarTree[x] 					= []
							
						# iterate throug variable positions						
						for pos in function.varPos:
							#print(pos, pos[1].genotype)
							# position, types, ref_amino_acid, alt_amino_acid
							cultivarTree_tpm = [pos[1].ID, pos[1].types, pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].genotype]
							
							# section for single amino acid exchanges. Different Approaches for different Variant Types are necessary
							# append generalInfomration to cultivarTree
							
							if 'INS' not in pos[1].types and 'DEL' not in pos[1].types and 'Stop gained' not in pos[1].types:
								try:
									if pos[1].pos_amio_acid_adapted >= int(function.properties[2])-3 and pos[1].pos_amio_acid_adapted <= int(function.properties[2])+3:
										pos[1].peptideScore = round((1-(log(1+(matrices['all'][pos[1].ref_amino_acid+pos[1].alt_amino_acid])))/log(maxima['all'])),3)
										function_scores_restictionSides[pos[0]].append(pos[1].peptideScore)
										cultivarTree_tpm.extend([pos[1].ID, pos[1].ID, pos[1].pos_amio_acid_adapted, pos[1].pos_amio_acid_adapted, pos[1].peptideScore, 1])
										#print(pos[0], pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amino_acid, pos[1].pos_amio_acid_adapted)
									else:
										pos[1].peptideScore = round((1-(log(1+(matrices['all'][pos[1].ref_amino_acid+pos[1].alt_amino_acid])))/log(maxima['all'])),3)
										function_scores_otherSides[pos[0]].append(pos[1].peptideScore)
										cultivarTree_tpm.extend([pos[1].ID, pos[1].ID, pos[1].pos_amio_acid_adapted, pos[1].pos_amio_acid_adapted, pos[1].peptideScore, 0])
								except:
									continue
									print('problem with Key Error, '+pos[1].pos_amio_acid)

							if 'INS' in pos[1].types and 'Stop gained' not in pos[1].types and len(pos[1].alt_amino_acid)>1:
								insertion = pos[1].alt_amino_acid
								# corrrect insertion for 'replace site by itself'
								if insertion[0] == pos[1].ref_amino_acid:
									insertion = insertion[1:]
								elif insertion[-1] == pos[1].ref_amino_acid:
									insertion = insertion[:1]
								if pos[1].pos_amio_acid_adapted >= int(function.properties[2])-3 and pos[1].pos_amio_acid_adapted <= int(function.properties[2])+3:
									pos[1].peptideScore = 1
									function_scores_restictionSides[pos[0]].append(pos[1].peptideScore)
									cultivarTree_tpm.extend([pos[1].ID, pos[1].ID, pos[1].pos_amio_acid_adapted, pos[1].pos_amio_acid_adapted, pos[1].peptideScore, 0])
									#print(pos[0], pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amino_acid, pos[1].pos_amio_acid_adapted)
								else:
									if len(insertion) <= 3:
										# calculate maximum influencing amino acid as score for the position	
										maxValue = 0	
										for a in insertion:
											tmpValue =  round((1-(log(1+(matrices['all']['-'+a])))/log(maxima['all'])),3)
											if tmpValue > maxValue:
												maxValue = tmpValue
										pos[1].peptideScore = maxValue
										function_scores_otherSides[pos[0]].append(pos[1].peptideScore)
										cultivarTree_tpm.extend([pos[1].ID, pos[1].ID, pos[1].pos_amio_acid_adapted, pos[1].pos_amio_acid_adapted, pos[1].peptideScore, 0])
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].secondaryScore )
									else:
										pos[1].peptideScore = 1.
										function_scores_otherSides[pos[0]].append(pos[1].peptideScore)
										cultivarTree_tpm.extend([pos[1].ID, pos[1].ID, pos[1].pos_amio_acid_adapted, pos[1].pos_amio_acid_adapted, pos[1].peptideScore, 0])
										
						
							if 'DEL' in pos[1].types and 'Stop gained' not in pos[1].types and len(pos[1].ref_amino_acid)>1:
								if pos[1].delStart != None and pos[1].delEnd != None:
									#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, pos[1].delStart, pos[1].delEnd, deletion )
									deletion = pos[1].ref_amino_acid
									# corrrect insertion for 'replace site by itself'
									if deletion[0] == pos[1].alt_amino_acid:
										deletion = deletion[1:]
									elif deletion[-1] == pos[1].alt_amino_acid:
										deletion = deletion[:1]
									
									if (pos[1].delStart >= int(function.properties[2])-3 and pos[1].delStart <= int(function.properties[2])+3) or (pos[1].delEnd >= int(function.properties[2])-3 and pos[1].delEnd <= int(function.properties[2])+3) or (int(function.properties[2])-3 >= pos[1].delStart  and int(function.properties[2])+3 <= pos[1].delEnd):	
										if len(deletion) <= 3:
											# calculate maximum influencing amino acid as score for the position	
											maxValue = 0	
											for a in deletion:
												tmpValue =  round((1-(log(1+(matrices['all'][a+'-'])))/log(maxima['all'])),3)
												if tmpValue > maxValue:
													maxValue = tmpValue
											pos[1].peptideScore = maxValue
											function_scores_restictionSides[pos[0]].append(pos[1].peptideScore)
											cultivarTree_tpm.extend([pos[1].ID, pos[1].ID, pos[1].delStart, pos[1].delEnd, pos[1].peptideScore, 0])
											#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].peptideScore )
										else:
											pos[1].peptideScore = 1.
											function_scores_restictionSides[pos[0]].append(pos[1].peptideScore)
											cultivarTree_tpm.extend([pos[1].ID, pos[1].ID, pos[1].delStart, pos[1].delEnd, pos[1].peptideScore, 0])
											#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].peptideScore )
									else:
										if len(deletion) <= 3:
											# calculate maximum influencing amino acid as score for the position	
											maxValue = 0	
											for a in deletion:
												tmpValue =  round((1-(log(1+(matrices['all'][a+'-'])))/log(maxima['all'])),3)
												if tmpValue > maxValue:
													maxValue = tmpValue
											pos[1].peptideScore = maxValue
											function_scores_otherSides[pos[0]].append(pos[1].peptideScore)
											cultivarTree_tpm.extend([pos[1].ID, pos[1].ID, pos[1].delStart, pos[1].delEnd, pos[1].peptideScore, 0])
											#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].peptideScore )
										else:
											pos[1].peptideScore = 1.
											function_scores_otherSides[pos[0]].append(pos[1].peptideScore)
											cultivarTree_tpm.extend([pos[1].ID, pos[1].ID, pos[1].delStart, pos[1].delEnd, pos[1].peptideScore, 0])
											#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].peptideScore )
							
	
							cultivarTree[pos[0]].append(cultivarTree_tpm)
							#print(cultivarTree)
#						#print(sorted(function_scores))
						fs = {}
						#print(function_scores_otherSides, function_scores_restictionSides, variant)
						#print(function_scores_otherSides, function_scores_restictionSides)
				#		for i in function_scores_restictionSides:
				#			#print(i, function_scores_restictionSides[i], variant)
				#			fs[i] =  stepFunction(sorted(function_scores_restictionSides[i], reverse = True), stepStart, stepLength)
				#			try:
				#				fs[i] += sum(function_scores_otherSides[i])/float((int(function.properties[2])-int(function.properties[1])-2))
				#			except:
				#				fs[i] += 0
				#			#print(round(fs[i],3))
				#			allFunctionsScores[i].append(round(fs[i],3))

						for i in function_scores_restictionSides:
							fs[i] = 0
							#print(i, function_scores_restictionSides[i], variant)
							#print(float((int(function.properties[2])-int(function.properties[1]))), function.properties[2], function.properties[1])
							try:
								fs[i] += sum(function_scores_restictionSides[i])/float((int(function.properties[2])-int(function.properties[1])-1))
							except:
								fs[i] += 0
							try:
								fs[i] += sum(function_scores_otherSides[i])/float((int(function.properties[2])-int(function.properties[1])-1))
							except:
								fs[i] += 0
							fs[i] = round(fs[i],3)
							allFunctionsScores[i].append(round(fs[i],3))
						#print(fs)
						function.primaryScore = fs
						#print(function.primaryScore)
						# function_properties, positionInformation_protein, floatingIDForSort, CultivarInfomration
						#print(function, function.properties[0], function.functionID, functionID)
						transkriptTrees[1].append([function.properties[0], function.functionID, functionID, cultivarTree, function.primaryScore])
						
				# evaluation part for the variants:
				#print(allFunctionsScores)
				max1 = []	
				print(allFunctionsScores)
				for i in allFunctionsScores:
					max1.append(stepFunction(allFunctionsScores[i],stepStart,stepLength))				
				
				# calculate diffScore and maxScore
				diffScore = maxPairwiseDistance(max1)
				maxScore  = max(max1)	
				variantO.diffScore	=	diffScore
				print(variantO)
				transcript_diff_scores[variantO.name]	= diffScore
				variantO.maxScore	=	maxScore
				transcript_max_scores[variantO.name]	= maxScore
				
				geneTree[-1].append(transkriptTrees)

				
			# create Overall property score (weight representative Transcirpt with 0.8, rest distributed over 0.2)
			propertyScoreMax 	= 0
			propertyScoreDiff	= 0
			
			
			for i in transcript_diff_scores:
				if i.split('.')[1] == '1':
					if len(transcript_diff_scores) == 1:
						propertyScoreMax += 1*transcript_diff_scores[i]
					else:
						propertyScoreMax += 0.8*transcript_diff_scores[i]
				else:
					if len(transcript_diff_scores) > 1:
						propertyScoreMax += (0.2/(len(transcript_diff_scores)-1))*transcript_diff_scores[i]
			
			for i in transcript_max_scores:
				if i.split('.')[1] == '1':
					#print(i.split('.')[1])
					if len(transcript_max_scores) == 1:
						propertyScoreMax += 1*transcript_max_scores[i]
					else:
						propertyScoreMax += weightPrimaryTranscript*transcript_max_scores[i]
				else:
					if len(transcript_max_scores) > 1:
						propertyScoreMax += ((1-weightPrimaryTranscript)/(len(transcript_max_scores)-1))*transcript_max_scores[i]
			
			# store scores in comparenceResult
			comparenceResult[result].diffScorePeptide	= round(propertyScoreDiff,3)
			comparenceResult[result].maxScorePeptide	= round(propertyScoreMax,3)
			geneTree.append(round(propertyScoreMax,3))
			geneTree.append(round(propertyScoreDiff,3))
			#print(geneTree)
			printTree.append(geneTree)
		return printTree
		
	def externalModifications(self, comparenceResult, cultivars):
		"""
		External modifications differentiates three classes:
		Glycosylation:		Contains Glycolysations
		Lipidation:			Contains Lipidations
		Modified residue:	Contains all modifications like methylation, phosphorylation etc.
		Because it is difficult to know which external modifications will be chanced, we will 
		use an 'all or nothing' apporach, evaluating a change with 1 and no change with 0.
		"""
		
		### region for all controllable variables
		# TODO
		stepLength				= 0.5
		stepStart				= 1		
		weightPrimaryTranscript	= 0.8
		
		for result in comparenceResult:
			transcript_diff_scores	= {}
			transcript_max_scores	= {}
			
			for variant in comparenceResult[result].variants:
				variantO = (comparenceResult[result].variants[variant])
				# store function overlapping scorrings results per cultivar
				allFunctionsScores = {}
				for x in cultivars:
					allFunctionsScores[x] = []
			
				diffScore	=	0
				maxScore	=	0
				for function in variantO.functions:
					if function.properties[0] in ['Glycosylation','Lipidation','Modified residue']:
						#print(function.properties)
						function_scores = {}		# samples all scores for variants in a function per cultivar

						
						# initialize dictionary to store culitvar information
						for x in cultivars:
							function_scores[x]	= []
							
						# iterate throug variable positions						
						for pos in function.varPos:
							
							# section for single amino acid exchanges. Different Approaches for different Variant Types are necessary
							if 'INS' not in pos[1].types and 'DEL' not in pos[1].types and 'Stop gained' not in pos[1].types:
								if pos[1].pos_amio_acid_adapted >= int(function.properties[1]) and pos[1].pos_amio_acid_adapted <= int(function.properties[2]):
									pos[1].modificationScore = 1
									function_scores[pos[0]].append(pos[1].modificationScore)
									#print(pos[0], pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amino_acid, pos[1].pos_amio_acid_adapted)
													
							if 'INS' in pos[1].types and 'Stop gained' not in pos[1].types:
								if pos[1].pos_amio_acid_adapted >= int(function.properties[1]) and pos[1].pos_amio_acid_adapted <= int(function.properties[2]) and (pos[1].ref_amino_acid != pos[1].alt_amino_acid[0] and pos[1].ref_amino_acid != pos[1].alt_amino_acid[-1]):
									pos[1].modificationScore = 1
									function_scores[pos[0]].append(pos[1].modificationScore)
									#print(function.properties)
									#print(pos[0], pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amino_acid, pos[1].pos_amio_acid_adapted)

							if 'DEL' in pos[1].types and 'Stop gained' not in pos[1].types:
								if pos[1].delStart != None and pos[1].delEnd != None:
									if int(function.properties[1]) >= pos[1].delStart and int(function.properties[2] <=  pos[1].delEnd) and ((pos[1].ref_amino_acid != pos[1].alt_amino_acid[0] and pos[1].ref_amino_acid != pos[1].alt_amino_acid[-1]) or pos[1].delStart > function.properties[1])and len(pos[1].alt_amino_acid) > 1:
										pos[1].modificationScore = 1
										function_scores[pos[0]].append(pos[1].modificationScore)
										print(pos[0], pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amino_acid, pos[1].delStart, pos[1].delEnd,function.properties)
							
						#print(sorted(function_scores))
						fs = {}
						for i in function_scores:
							fs[i] =  stepFunction(sorted(function_scores[i], reverse = True), stepStart, stepLength)
							allFunctionsScores[i].append(round(fs[i],3))
	
											
						function.modificationScore = fs
				# evaluation part for the variants:
				#print(allFunctionsScores)
				max1 = []	
				for i in allFunctionsScores:
					max1.append(stepFunction(allFunctionsScores[i],stepStart,stepLength))				
				
				# calculate diffScore and maxScore
				diffScore = maxPairwiseDistance(max1)
				maxScore  = max(max1)	
				variantO.diffScore	=	diffScore
				transcript_diff_scores[variantO.name]	= diffScore
				variantO.maxScore	=	maxScore
				transcript_max_scores[variantO.name]	= maxScore
				
			# create Overall property score (weight representative Transcirpt with 0.8, rest distributed over 0.2)
			propertyScoreMax 	= 0
			propertyScoreDiff	= 0
			
			
			for i in transcript_diff_scores:
				if i.split('.')[1] == '1':
					if len(transcript_diff_scores) == 1:
						propertyScoreMax += 1*transcript_diff_scores[i]
					else:
						propertyScoreMax += 0.8*transcript_diff_scores[i]
				else:
					if len(transcript_diff_scores) > 1:
						propertyScoreMax += (0.2/(len(transcript_diff_scores)-1))*transcript_diff_scores[i]
			
			for i in transcript_max_scores:
				if i.split('.')[1] == '1':
					#print(i.split('.')[1])
					if len(transcript_max_scores) == 1:
						propertyScoreMax += 1*transcript_max_scores[i]
					else:
						propertyScoreMax += weightPrimaryTranscript*transcript_max_scores[i]
				else:
					if len(transcript_max_scores) > 1:
						propertyScoreMax += ((1-weightPrimaryTranscript)/(len(transcript_max_scores)-1))*transcript_max_scores[i]
			
			# store scores in comparenceResult
			comparenceResult[result].diffScoreModifications	= round(propertyScoreDiff,3)
			comparenceResult[result].maxScoreModifications	= round(propertyScoreMax,3)				

	def motifAnalysis(self, comparenceResult, cultivars):
		"""
		There are basically two classes of motifs in uniprot: The one we
		can extract a motif sequence from and the one we can not.
		We apply the following scoring sceme to account for that:
		If there is a known motif sequence; 	
			Full Effect if key side is exchanged, else 0.
		If there is not a known motif sequence;	
			amount of var / len(motif)**0.5, but with max 1
			-> this accounts for the fact that 5 out of 10 might be way more
			invasive than 0.5 would reflect (compared to 1 (instead of 1.58)
		Inserttions and Deletions are also counted as full effect
		(because motifs depend on their lenth)			
		"""
		
		### region for all controllable variables
		# TODO
		stepLength				= 0.5
		stepStart				= 1		
		weightPrimaryTranscript	= 0.8
		
		for result in comparenceResult:
			transcript_diff_scores	= {}
			transcript_max_scores	= {}
			
			for variant in comparenceResult[result].variants:
				variantO = (comparenceResult[result].variants[variant])
				# store function overlapping scorrings results per cultivar
				allFunctionsScores = {}
				for x in cultivars:
					allFunctionsScores[x] = []
			
				diffScore	=	0
				maxScore	=	0
				for function in variantO.functions:
					if function.properties[0] in ['Motif']:
						#print(function.properties)
						function_scores_motif = {}		# samples all scores for variants in a function per cultivar
						function_scores_notMotiv = {}
						
						# get sequence motif if possible
						motif		= 	None
						tmp			= 	(function.properties[3].split(';')[0].split('=')[1].split(' '))
						motifLength = 	function.properties[2]-function.properties[1]+1
						if 'motif' in tmp:
							if len(tmp[0]) == motifLength:
								motif = tmp[0]
						if 'region' in tmp:
							if len(tmp[0].replace('"', '')) == motifLength:
								motif = tmp[0].replace('"', '')
						if len(tmp) == 1 and  len(tmp[0].replace('"', '')) == motifLength:
							motif = tmp[0].replace('"', '')

						
						# initialize dictionary to store culitvar information
						for x in cultivars:
							function_scores_motif[x]	= []
							function_scores_notMotiv[x]	= []
							
						# iterate through variable positions						
						for pos in function.varPos:
							
							# section for single amino acid exchanges. Different Approaches for different Variant Types are necessary
							if 'INS' not in pos[1].types and 'DEL' not in pos[1].types and 'Stop gained' not in pos[1].types:
								if motif ==  None:
									if pos[1].pos_amio_acid_adapted >= int(function.properties[1]) and pos[1].pos_amio_acid_adapted <= int(function.properties[2]):
										pos[1].motifScore = round(1./(motifLength**0.5),3)
										function_scores_notMotiv[pos[0]].append(pos[1].motifScore)
										#if result =='AT5G45260':
										#	print(pos[0], pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amino_acid, pos[1].pos_amio_acid_adapted, motif, function.properties[1:3], pos[1].motifScore)
								else:
									if pos[1].pos_amio_acid_adapted >= int(function.properties[1]) and pos[1].pos_amio_acid_adapted <= int(function.properties[2]):
										posInMotiv = pos[1].pos_amio_acid_adapted - int(function.properties[1])
										if motif[posInMotiv] != 'X':
											pos[1].motifScore = 1
											function_scores_motif[pos[0]].append(pos[1].motifScore)

										#print(pos[0], pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amino_acid, pos[1].pos_amio_acid_adapted, motif, function.properties[1:3], pos[1].motifScore)
										

							if 'INS' in pos[1].types and 'Stop gained' not in pos[1].types:
								if motif ==  None:
									if pos[1].pos_amio_acid_adapted >= int(function.properties[1]) and pos[1].pos_amio_acid_adapted <= int(function.properties[2]):
										pos[1].motifScore = 1
										function_scores_notMotiv[pos[0]].append(pos[1].motifScore)
								else:
									if pos[1].pos_amio_acid_adapted >= int(function.properties[1]) and pos[1].pos_amio_acid_adapted <= int(function.properties[2]):
										pos[1].motifScore = 1
										function_scores_motif[pos[0]].append(pos[1].motifScore)
										
							if 'DEL' in pos[1].types and 'Stop gained' not in pos[1].types and len(pos[1].ref_amino_acid)>1:
								#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, pos[1].delStart, pos[1].delEnd)
								if pos[1].delStart != None and pos[1].delEnd != None:
									print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, pos[1].delStart, pos[1].delEnd, deletion )
									deletion = pos[1].ref_amino_acid
									# corrrect insertion for 'replace site by itself'
									if deletion[0] == pos[1].alt_amino_acid:
										deletion = deletion[1:]
									elif deletion[-1] == pos[1].alt_amino_acid:
										deletion = deletion[:1]
										
									if motif ==  None:							
										pos[1].motifScore = 1
										function_scores_notMotiv[pos[0]].append(pos[1].motifScore)
									else:
										pos[1].motifScore = 1
										function_scores_motif[pos[0]].append(pos[1].motifScore)

						fs = {}

						if motif  != None:
							for i in function_scores_motif:
								fs[i] =  stepFunction(sorted(function_scores_motif[i], reverse = True), stepStart, stepLength)
								allFunctionsScores[i].append(round(fs[i],3))
						else:
							for i in function_scores_notMotiv:
								fs[i] =  sum(function_scores_notMotiv[i])
								if fs[i] > 1:
									fs[i] = 1
								allFunctionsScores[i].append(round(fs[i],3))
						#print(allFunctionsScores)
				
						function.motifScore = fs
				# evaluation part for the variants:
				#print(allFunctionsScores)
				max1 = []	
				for i in allFunctionsScores:
					max1.append(stepFunction(allFunctionsScores[i],stepStart,stepLength))				
				
				# calculate diffScore and maxScore
				diffScore = maxPairwiseDistance(max1)
				maxScore  = max(max1)	
				variantO.diffScore	=	diffScore
				transcript_diff_scores[variantO.name]	= diffScore
				variantO.maxScore	=	maxScore
				transcript_max_scores[variantO.name]	= maxScore
				
			# create Overall property score (weight representative Transcirpt with 0.8, rest distributed over 0.2)
			propertyScoreMax 	= 0
			propertyScoreDiff	= 0			
			
			for i in transcript_diff_scores:
				if i.split('.')[1] == '1':
					if len(transcript_diff_scores) == 1:
						propertyScoreMax += 1*transcript_diff_scores[i]
					else:
						propertyScoreMax += 0.8*transcript_diff_scores[i]
				else:
					if len(transcript_diff_scores) > 1:
						propertyScoreMax += (0.2/(len(transcript_diff_scores)-1))*transcript_diff_scores[i]
			
			for i in transcript_max_scores:
				if i.split('.')[1] == '1':
					#print(i.split('.')[1])
					if len(transcript_max_scores) == 1:
						propertyScoreMax += 1*transcript_max_scores[i]
					else:
						propertyScoreMax += weightPrimaryTranscript*transcript_max_scores[i]
				else:
					if len(transcript_max_scores) > 1:
						propertyScoreMax += ((1-weightPrimaryTranscript)/(len(transcript_max_scores)-1))*transcript_max_scores[i]
			
			# store scores in comparenceResult
			comparenceResult[result].diffScoreMotif	= round(propertyScoreDiff,3)
			comparenceResult[result].maxScoreMotif	= round(propertyScoreMax,3)
			
	def bindingActiveSiteAnalysis(self, comparenceResult, matrices, cultivars):
		"""
		Analyses the influence of SNP, Insertions and Deletions on the
		binding sites and active sites of a protein
		"""
		# import some mathematical opertations
		from math import log
			
		### region for all controllable variables
		# TODO
		stepLength				= 0.5
		stepStart				= 1		
		weightPrimaryTranscript	= 0.8
		
#		###
		
		# get all structural element names out of matrices
		elements = [x for  x in matrices]

		# safe maxima of every matrix in dict
		maxima   = {}
		for element in elements:
			maxima[element] = getMaxValueOfMatrix(matrices[element])

		for result in comparenceResult:
			transcript_diff_scores	= {}
			transcript_max_scores	= {}
			
			for variant in comparenceResult[result].variants:
				variantO = (comparenceResult[result].variants[variant])
				# store function overlapping scorrings results per cultivar
				allFunctionsScores = {}
				for x in cultivars:
					allFunctionsScores[x] = []
			
				diffScore	=	0
				maxScore	=	0
				for function in variantO.functions:
					if function.properties[0] in ['Metal binding','Binding site','Active site']:
						#print(function.properties)
						function_scores = {}		# samples all scores for variants in a function per cultivar
						
						# initialize dictionary to store culitvar information
						for x in cultivars:
							function_scores[x]	= []

							
						# iterate throug variable positions						
						for pos in function.varPos:
							
							# section for single amino acid exchanges. Different Approaches for different Variant Types are necessary
							if 'INS' not in pos[1].types and 'DEL' not in pos[1].types and 'Stop gained' not in pos[1].types:
								if pos[1].pos_amio_acid_adapted >= int(function.properties[1])-3 and pos[1].pos_amio_acid_adapted <= int(function.properties[2])+3:
									pos[1].activeSiteScore = round((1-(log(1+(matrices['all'][pos[1].ref_amino_acid+pos[1].alt_amino_acid])))/log(maxima['all'])),3)
									function_scores[pos[0]].append(pos[1].activeSiteScore)
									#print(pos[0], pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amino_acid, pos[1].pos_amio_acid_adapted, function.properties, pos[1].activeSiteScore)
						

							if 'INS' in pos[1].types and 'Stop gained' not in pos[1].types:
								if pos[1].pos_amio_acid_adapted >= int(function.properties[1]) and pos[1].pos_amio_acid_adapted <= int(function.properties[2]) and (pos[1].ref_amino_acid != pos[1].alt_amino_acid[0] and pos[1].ref_amino_acid != pos[1].alt_amino_acid[-1]) and len(pos[1].alt_amino_acid) > 1:
									pos[1].activeSiteScore = 1
									function_scores[pos[0]].append(pos[1].activeSiteScore)
									#print(function.properties)
									#print(pos[0], pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amino_acid, pos[1].pos_amio_acid_adapted,function.properties, result)

							if 'DEL' in pos[1].types and 'Stop gained' not in pos[1].types:
								if pos[1].delStart != None and pos[1].delEnd != None:
									if int(function.properties[1]) >= pos[1].delStart and int(function.properties[2] <=  pos[1].delEnd) and ((pos[1].ref_amino_acid != pos[1].alt_amino_acid[0] and pos[1].ref_amino_acid != pos[1].alt_amino_acid[-1]) or pos[1].delStart > function.properties[1]):
										pos[1].activeSiteScore = 1
										function_scores[pos[0]].append(pos[1].activeSiteScore)
										#print(pos[0], pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amino_acid, pos[1].delStart, pos[1].delEnd,function.properties)

#						#print(sorted(function_scores))
						fs = {}
						for i in function_scores:
							fs[i] =  stepFunction(sorted(function_scores[i], reverse = True), stepStart, stepLength)
	
							#print(round(fs[i],3))
							allFunctionsScores[i].append(round(fs[i],3))
												
						function.primaryScore = fs
				# evaluation part for the variants:
				#print(allFunctionsScores)
				max1 = []	
				for i in allFunctionsScores:
					max1.append(stepFunction(allFunctionsScores[i],stepStart,stepLength))				
				
				# calculate diffScore and maxScore
				diffScore = maxPairwiseDistance(max1)
				maxScore  = max(max1)	
				variantO.diffScore	=	diffScore
				transcript_diff_scores[variantO.name]	= diffScore
				variantO.maxScore	=	maxScore
				transcript_max_scores[variantO.name]	= maxScore
				
			# create Overall property score (weight representative Transcirpt with 0.8, rest distributed over 0.2)
			propertyScoreMax 	= 0
			propertyScoreDiff	= 0
#			
#			
			for i in transcript_diff_scores:
				if i.split('.')[1] == '1':
					if len(transcript_diff_scores) == 1:
						propertyScoreMax += 1*transcript_diff_scores[i]
					else:
						propertyScoreMax += 0.8*transcript_diff_scores[i]
				else:
					if len(transcript_diff_scores) > 1:
						propertyScoreMax += (0.2/(len(transcript_diff_scores)-1))*transcript_diff_scores[i]
			
			for i in transcript_max_scores:
				if i.split('.')[1] == '1':
					#print(i.split('.')[1])
					if len(transcript_max_scores) == 1:
						propertyScoreMax += 1*transcript_max_scores[i]
					else:
						propertyScoreMax += weightPrimaryTranscript*transcript_max_scores[i]
				else:
					if len(transcript_max_scores) > 1:
						propertyScoreMax += ((1-weightPrimaryTranscript)/(len(transcript_max_scores)-1))*transcript_max_scores[i]
			
			# store scores in comparenceResult
			comparenceResult[result].diffScoreActiveSites	= round(propertyScoreDiff,3)
			comparenceResult[result].maxScoreActiveSites	= round(propertyScoreMax,3)	

	def mutagenesisAnalysis(self, comparenceResult, matrices, cultivars):
		"""
		Analyses the influence of SNP, Insertions and Deletions on the
		binding sites and active sites of a protein
		"""
		# import some mathematical opertations
		from math import log
			
		### region for all controllable variables
		# TODO
		stepLength				= 0.5
		stepStart				= 1		
		weightPrimaryTranscript	= 0.8
		
#		###
		
		# get all structural element names out of matrices
		elements = [x for  x in matrices]

		# safe maxima of every matrix in dict
		maxima   = {}
		for element in elements:
			maxima[element] = getMaxValueOfMatrix(matrices[element])

		for result in comparenceResult:
			transcript_diff_scores	= {}
			transcript_max_scores	= {}
			
			for variant in comparenceResult[result].variants:
				variantO = (comparenceResult[result].variants[variant])
				# store function overlapping scorrings results per cultivar
				allFunctionsScores = {}
				for x in cultivars:
					allFunctionsScores[x] = []
			
				diffScore	=	0
				maxScore	=	0
				for function in variantO.functions:
					if function.properties[0] in ['Mutagenesis']:
						#print(function.properties)
						function_scores = {}		# samples all scores for variants in a function per cultivar
						
						# initialize dictionary to store culitvar information
						for x in cultivars:
							function_scores[x]	= []

							
						# iterate throug variable positions						
						for pos in function.varPos:
							
							# section for single amino acid exchanges. Different Approaches for different Variant Types are necessary
							if 'INS' not in pos[1].types and 'DEL' not in pos[1].types and 'Stop gained' not in pos[1].types:
								if pos[1].pos_amio_acid_adapted >= int(function.properties[1]) and pos[1].pos_amio_acid_adapted <= int(function.properties[2]) and (int(function.properties[2])-int(function.properties[1])) <11:
									pos[1].mutagenesisScore = round((1-(log(1+(matrices['all'][pos[1].ref_amino_acid+pos[1].alt_amino_acid])))/log(maxima['all'])),3)
									function_scores[pos[0]].append(pos[1].mutagenesisScore)
									#print(pos[0], pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amino_acid, pos[1].pos_amio_acid_adapted, function.properties, pos[1].mutagenesisScore)
						

							if 'INS' in pos[1].types and 'Stop gained' not in pos[1].types and len(pos[1].alt_amino_acid)>1:
								insertion = pos[1].alt_amino_acid
								# corrrect insertion for 'replace site by itself'
								if insertion[0] == pos[1].ref_amino_acid:
									insertion = insertion[1:]
								elif insertion[-1] == pos[1].ref_amino_acid:
									insertion = insertion[:1]
									
									if len(insertion) <= 3:
										# calculate maximum influencing amino acid as score for the position	
										maxValue = 0	
										for a in insertion:
											tmpValue =  round((1-(log(1+(matrices['all']['-'+a])))/log(maxima['all'])),3)
											if tmpValue > maxValue:
												maxValue = tmpValue
										pos[1].mutagenesisScore = maxValue
										allFunctionsScores[pos[0]].append(pos[1].mutagenesisScore)
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].domainScore )
									else:
										pos[1].mutagenesisScore = len(insertion)
										allFunctionsScores[pos[0]].append(pos[1].mutagenesisScore)
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].domainScore )
						
							if 'DEL' in pos[1].types and 'Stop gained' not in pos[1].types and len(pos[1].ref_amino_acid)>1:
								if pos[1].delStart != None and pos[1].delEnd != None:
									deletion = pos[1].ref_amino_acid
									#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, pos[1].delStart, pos[1].delEnd, deletion )
									# corrrect insertion for 'replace site by itself'
									if deletion[0] == pos[1].alt_amino_acid:
										deletion = deletion[1:]
									elif deletion[-1] == pos[1].alt_amino_acid:
										deletion = deletion[:1]
										
									if len(deletion) <= 3:
										# calculate maximum influencing amino acid as score for the position	
										maxValue = 0	
										for a in deletion:
											tmpValue =  round((1-(log(1+(matrices['all'][a+'-'])))/log(maxima['all'])),3)
											if tmpValue > maxValue:
												maxValue = tmpValue
										pos[1].mutagenesisScore = maxValue
										function_scores[pos[0]].append(pos[1].mutagenesisScore)
										print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].mutagenesisScore )
									else:
										pos[1].mutagenesisScore = 1.
										function_scores[pos[0]].append(pos[1].mutagenesisScore)
										print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].mutagenesisScore )

									#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, pos[1].delStart, pos[1].delEnd, deletion )


						fs = {}
						for i in function_scores:
							fs[i] =  stepFunction(sorted(function_scores[i], reverse = True), stepStart, stepLength)
	
							#print(round(fs[i],3))
							allFunctionsScores[i].append(round(fs[i],3))
												
						function.primaryScore = fs
				# evaluation part for the variants:
				#print(allFunctionsScores)
				max1 = []	
				for i in allFunctionsScores:
					max1.append(stepFunction(allFunctionsScores[i],stepStart,stepLength))				
				
				# calculate diffScore and maxScore
				diffScore = maxPairwiseDistance(max1)
				maxScore  = max(max1)	
				variantO.diffScore	=	diffScore
				transcript_diff_scores[variantO.name]	= diffScore
				variantO.maxScore	=	maxScore
				transcript_max_scores[variantO.name]	= maxScore
				
			# create Overall property score (weight representative Transcirpt with 0.8, rest distributed over 0.2)
			propertyScoreMax 	= 0
			propertyScoreDiff	= 0
#			
#			
			for i in transcript_diff_scores:
				if i.split('.')[1] == '1':
					if len(transcript_diff_scores) == 1:
						propertyScoreMax += 1*transcript_diff_scores[i]
					else:
						propertyScoreMax += 0.8*transcript_diff_scores[i]
				else:
					if len(transcript_diff_scores) > 1:
						propertyScoreMax += (0.2/(len(transcript_diff_scores)-1))*transcript_diff_scores[i]
			
			for i in transcript_max_scores:
				if i.split('.')[1] == '1':
					#print(i.split('.')[1])
					if len(transcript_max_scores) == 1:
						propertyScoreMax += 1*transcript_max_scores[i]
					else:
						propertyScoreMax += weightPrimaryTranscript*transcript_max_scores[i]
				else:
					if len(transcript_max_scores) > 1:
						propertyScoreMax += ((1-weightPrimaryTranscript)/(len(transcript_max_scores)-1))*transcript_max_scores[i]
			
			# store scores in comparenceResult
			comparenceResult[result].diffScoreMutagenesis	= round(propertyScoreDiff,3)
			comparenceResult[result].maxScoreMutagenesis	= round(propertyScoreMax,3)
			
	def dnaNucleotideBindingAnalysis(self, comparenceResult, matrices, cultivars):
		"""
		Analyses the influence of SNP, Insertions and Deletions on the
		secondary structure elementes of a protein
		"""
		# import some mathematical opertations
		from math import log
			
		### region for all controllable variables
		# TODO
		stepLength				= 0.5
		stepStart				= 1		
		weightPrimaryTranscript	= 0.8
		
#		###
		
		# get all structural element names out of matrices
		elements = [x for  x in matrices]

		# safe maxima of every matrix in dict
		maxima   = {}
		for element in elements:
			maxima[element] = getMaxValueOfMatrix(matrices[element])

		for result in comparenceResult:
			transcript_diff_scores	= {}
			transcript_max_scores	= {}
			
			for variant in comparenceResult[result].variants:
				variantO = (comparenceResult[result].variants[variant])
				# store function overlapping scorrings results per cultivar
				allFunctionsScores = {}
				for x in cultivars:
					allFunctionsScores[x] = []
			
				diffScore	=	0
				maxScore	=	0
				for function in variantO.functions:
					if function.properties[0] in ['Calcium binding','Nucleotide binding','DNA binding','Zinc finger']:
					
						#print(function.properties)
						function_scores = {}		# samples all scores for variants in a function per cultivar
						
						# initialize dictionary to store culitvar information
						for x in cultivars:
							function_scores[x]	= []
							
						# iterate throug variable positions						
						for pos in function.varPos:
							
							# section for single amino acid exchanges. Different Approaches for different Variant Types are necessary
							if 'INS' not in pos[1].types and 'DEL' not in pos[1].types and 'Stop gained' not in pos[1].types:
								if pos[1].pos_amio_acid_adapted >= int(function.properties[1]) and pos[1].pos_amio_acid_adapted <= int(function.properties[2]):
									pos[1].dnaNucBindingScore = round((1-(log(1+(matrices['all'][pos[1].ref_amino_acid+pos[1].alt_amino_acid])))/log(maxima['all'])),3)
									function_scores[pos[0]].append(pos[1].dnaNucBindingScore)
									#print(pos[0], pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amino_acid, pos[1].pos_amio_acid_adapted, function.properties, pos[1].dnaNucBindingScore)
						

							if 'INS' in pos[1].types and 'Stop gained' not in pos[1].types and len(pos[1].alt_amino_acid)>1:
								insertion = pos[1].alt_amino_acid
								# corrrect insertion for 'replace site by itself'
								if insertion[0] == pos[1].ref_amino_acid:
									insertion = insertion[1:]
								elif insertion[-1] == pos[1].ref_amino_acid:
									insertion = insertion[:1]
									
									if len(insertion) <= 3:
										# calculate maximum influencing amino acid as score for the position	
										maxValue = 0	
										for a in insertion:
											tmpValue =  round((1-(log(1+(matrices['all']['-'+a])))/log(maxima['all'])),3)
											if tmpValue > maxValue:
												maxValue = tmpValue
										pos[1].dnaNucBindingScore = maxValue
										allFunctionsScores[pos[0]].append(pos[1].dnaNucBindingScore)
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].dnaNucBindingScore )
									else:
										pos[1].dnaNucBindingScore = len(insertion)
										allFunctionsScores[pos[0]].append(pos[1].dnaNucBindingScore)
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].dnaNucBindingScore )
						
							if 'DEL' in pos[1].types and 'Stop gained' not in pos[1].types and len(pos[1].ref_amino_acid)>1:
								if pos[1].delStart != None and pos[1].delEnd != None:
									deletion = pos[1].ref_amino_acid
									#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, pos[1].delStart, pos[1].delEnd, deletion )
									# corrrect insertion for 'replace site by itself'
									if deletion[0] == pos[1].alt_amino_acid:
										deletion = deletion[1:]
									elif deletion[-1] == pos[1].alt_amino_acid:
										deletion = deletion[:1]
										
									if len(deletion) <= 3:
										# calculate maximum influencing amino acid as score for the position	
										maxValue = 0	
										for a in deletion:
											tmpValue =  round((1-(log(1+(matrices['all'][a+'-'])))/log(maxima['all'])),3)
											if tmpValue > maxValue:
												maxValue = tmpValue
										pos[1].dnaNucBindingScore = maxValue
										function_scores[pos[0]].append(pos[1].dnaNucBindingScore)
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].dnaNucBindingScore )
									else:
										pos[1].domainScore = 1.
										function_scores[pos[0]].append(pos[1].dnaNucBindingScore)
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].dnaNucBindingScore )
								
#						#print(sorted(function_scores))
						fs = {}
						for i in function_scores:
							try:
								fs[i] = sum(function_scores[i])/float((int(function.properties[2])-int(function.properties[1]+1)))
							except:
								fs[i] = 0
							allFunctionsScores[i].append(round(fs[i],3))
												
						function.primaryScore = fs
				# evaluation part for the variants:
				#print(allFunctionsScores)
				max1 = []	
				for i in allFunctionsScores:
					max1.append(stepFunction(allFunctionsScores[i],stepStart,stepLength))				
				
				# calculate diffScore and maxScore
				diffScore = maxPairwiseDistance(max1)
				maxScore  = max(max1)	
				variantO.diffScore	=	diffScore
				transcript_diff_scores[variantO.name]	= diffScore
				variantO.maxScore	=	maxScore
				transcript_max_scores[variantO.name]	= maxScore
				
			# create Overall property score (weight representative Transcirpt with 0.8, rest distributed over 0.2)
			propertyScoreMax 	= 0
			propertyScoreDiff	= 0
#			
#			
			for i in transcript_diff_scores:
				if i.split('.')[1] == '1':
					if len(transcript_diff_scores) == 1:
						propertyScoreMax += 1*transcript_diff_scores[i]
					else:
						propertyScoreMax += 0.8*transcript_diff_scores[i]
				else:
					if len(transcript_diff_scores) > 1:
						propertyScoreMax += (0.2/(len(transcript_diff_scores)-1))*transcript_diff_scores[i]
			
			for i in transcript_max_scores:
				if i.split('.')[1] == '1':
					#print(i.split('.')[1])
					if len(transcript_max_scores) == 1:
						propertyScoreMax += 1*transcript_max_scores[i]
					else:
						propertyScoreMax += weightPrimaryTranscript*transcript_max_scores[i]
				else:
					if len(transcript_max_scores) > 1:
						propertyScoreMax += ((1-weightPrimaryTranscript)/(len(transcript_max_scores)-1))*transcript_max_scores[i]
			
			# store scores in comparenceResult
			comparenceResult[result].diffScorediffScoreDnaNucCal	= round(propertyScoreDiff,3)
			comparenceResult[result].maxScorediffScoreDnaNucCal		= round(propertyScoreMax,3)
			
	def domainAnalysis(self, comparenceResult, matrices, cultivars):
		"""
		Analyses the influence of SNP, Insertions and Deletions on the
		secondary structure elementes of a protein
		"""
		# import some mathematical opertations
		from math import log
			
		### region for all controllable variables
		# TODO
		stepLength				= 0.5
		stepStart				= 1		
		weightPrimaryTranscript	= 0.8
		
#		###
		
		# get all structural element names out of matrices
		elements = [x for  x in matrices]

		# safe maxima of every matrix in dict
		maxima   = {}
		for element in elements:
			maxima[element] = getMaxValueOfMatrix(matrices[element])

		for result in comparenceResult:
			transcript_diff_scores	= {}
			transcript_max_scores	= {}
			
			for variant in comparenceResult[result].variants:
				variantO = (comparenceResult[result].variants[variant])
				# store function overlapping scorrings results per cultivar
				allFunctionsScores = {}
				for x in cultivars:
					allFunctionsScores[x] = []
			
				diffScore	=	0
				maxScore	=	0
				for function in variantO.functions:
					if function.properties[0] in ['Domain','Region']:
					
						#print(function.properties)
						function_scores = {}		# samples all scores for variants in a function per cultivar
						
						# initialize dictionary to store culitvar information
						for x in cultivars:
							function_scores[x]	= []
							
						# iterate throug variable positions						
						for pos in function.varPos:
							
							# section for single amino acid exchanges. Different Approaches for different Variant Types are necessary
							if 'INS' not in pos[1].types and 'DEL' not in pos[1].types and 'Stop gained' not in pos[1].types:
								if pos[1].pos_amio_acid_adapted >= int(function.properties[1]) and pos[1].pos_amio_acid_adapted <= int(function.properties[2]):
									try:
										pos[1].domainScore = round((1-(log(1+(matrices['all'][pos[1].ref_amino_acid+pos[1].alt_amino_acid])))/log(maxima['all'])),3)
										function_scores[pos[0]].append(pos[1].domainScore)
									#print(pos[0], pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amino_acid, pos[1].pos_amio_acid_adapted, function.properties, pos[1].domainScore)
									except:
										print('Error', pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties[0] )
										continue
						

							if 'INS' in pos[1].types and 'Stop gained' not in pos[1].types and len(pos[1].alt_amino_acid)>1:
								insertion = pos[1].alt_amino_acid
								# corrrect insertion for 'replace site by itself'
								if insertion[0] == pos[1].ref_amino_acid:
									insertion = insertion[1:]
								elif insertion[-1] == pos[1].ref_amino_acid:
									insertion = insertion[:1]
									
									if len(insertion) <= 3:
										# calculate maximum influencing amino acid as score for the position	
										maxValue = 0	
										for a in insertion:
											tmpValue =  round((1-(log(1+(matrices['all']['-'+a])))/log(maxima['all'])),3)
											if tmpValue > maxValue:
												maxValue = tmpValue
										pos[1].domainScore = maxValue
										allFunctionsScores[pos[0]].append(pos[1].domainScore)
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].domainScore )
									else:
										pos[1].domainScore = len(insertion)
										allFunctionsScores[pos[0]].append(pos[1].domainScore)
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].domainScore )
						
							if 'DEL' in pos[1].types and 'Stop gained' not in pos[1].types and len(pos[1].ref_amino_acid)>1:
								if pos[1].delStart != None and pos[1].delEnd != None:
									deletion = pos[1].ref_amino_acid
									#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, pos[1].delStart, pos[1].delEnd, deletion )
									# corrrect insertion for 'replace site by itself'
									if deletion[0] == pos[1].alt_amino_acid:
										deletion = deletion[1:]
									elif deletion[-1] == pos[1].alt_amino_acid:
										deletion = deletion[:1]
										
									if len(deletion) <= 3:
										# calculate maximum influencing amino acid as score for the position	
										maxValue = 0	
										for a in deletion:
											tmpValue =  round((1-(log(1+(matrices['all'][a+'-'])))/log(maxima['all'])),3)
											if tmpValue > maxValue:
												maxValue = tmpValue
										pos[1].domainScore = maxValue
										function_scores[pos[0]].append(pos[1].domainScore)
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].domainScore )
									else:
										pos[1].domainScore = 1.
										function_scores[pos[0]].append(pos[1].domainScore)
										#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, function.properties, result, pos[1].domainScore )

									#print(pos[0], pos[1].ref_amino_acid,  pos[1].alt_amino_acid, pos[1].delStart, pos[1].delEnd, deletion )
								
#						#print(sorted(function_scores))
						fs = {}
						for i in function_scores:
							try:
								fs[i] = sum(function_scores[i])/float((int(function.properties[2])-int(function.properties[1]+1)))
							except:
								fs[i] = 0
							allFunctionsScores[i].append(round(fs[i],3))
												
						function.primaryScore = fs
				# evaluation part for the variants:
				#print(allFunctionsScores)
				max1 = []	
				for i in allFunctionsScores:
					max1.append(stepFunction(allFunctionsScores[i],stepStart,stepLength))				
				
				# calculate diffScore and maxScore
				diffScore = maxPairwiseDistance(max1)
				maxScore  = max(max1)	
				variantO.diffScore	=	diffScore
				transcript_diff_scores[variantO.name]	= diffScore
				variantO.maxScore	=	maxScore
				transcript_max_scores[variantO.name]	= maxScore
				
			# create Overall property score (weight representative Transcirpt with 0.8, rest distributed over 0.2)
			propertyScoreMax 	= 0
			propertyScoreDiff	= 0
			
			
			for i in transcript_diff_scores:
				if i.split('.')[1] == '1':
					if len(transcript_diff_scores) == 1:
						propertyScoreMax += 1*transcript_diff_scores[i]
					else:
						propertyScoreMax += 0.8*transcript_diff_scores[i]
				else:
					if len(transcript_diff_scores) > 1:
						propertyScoreMax += (0.2/(len(transcript_diff_scores)-1))*transcript_diff_scores[i]
			
			for i in transcript_max_scores:
				if i.split('.')[1] == '1':
					#print(i.split('.')[1])
					if len(transcript_max_scores) == 1:
						propertyScoreMax += 1*transcript_max_scores[i]
					else:
						propertyScoreMax += weightPrimaryTranscript*transcript_max_scores[i]
				else:
					if len(transcript_max_scores) > 1:
						propertyScoreMax += ((1-weightPrimaryTranscript)/(len(transcript_max_scores)-1))*transcript_max_scores[i]
			
			# store scores in comparenceResult
			comparenceResult[result].diffScorediffDomain	= round(propertyScoreDiff,3)
			comparenceResult[result].maxScorediffDomain		= round(propertyScoreMax,3)

#			if round(propertyScoreMax,3) > 0:
#				print(result, round(propertyScoreMax,3))
			n = 0

	def disruptiveStopEffects(self, comparenceResult, cultivars, genes):
		"""
		Function to reflect the disruptive effects of stop insertions and 
		frameshifts. For this, we account for loss of known domains as 
		well as for the amount of lost amio-acids. 
		"""
		# import some mathematical opertations
		from math import log, inf
		gene_scores = {}
		
		for gene in genes:
			currentGene = genes[gene]
			variant_scores 		   	= {}
			affectedVariantCounter 	= 0
			allVariants			 	= len(currentGene.variants)
			#print(currentGene.variants)
			positions				= {}
			averageVariantCoverage	= {}
			
			for variant in currentGene.variants:
				currentVariant = currentGene.variants[variant]
				#print(currentVariant.ID)
				cultivar_scores = {}
				for cultivar in currentGene.cultivars:
					# stores the scores for the cultivars
					# initalize start position for first stoips depending on orientation
					if currentGene.orientation == '+':
						earliestStop = [None, inf]
					if currentGene.orientation == '-':
						earliestStop = [None, 0]
					
					genomePosition = 0	
					for vP in cultivar.variablePositions:
						position = cultivar.variablePositions[vP]
						
						for pos in position:
							#print(pos.ID, pos.types)
							if 'Stop gained' not in pos.types:
								continue
							# get variant with the earliest Stop
							if currentGene.orientation == '+' and int(pos.ID) < earliestStop[1] and pos.variant == variant:
								genomePosition 					= int(pos.ID)
								positions[cultivar.ID] 			= pos
								earliestStop 					= [pos,int(pos.ID)]
								averageVariantCoverage[cultivar.ID]= pos.genotype['DP'][0]
							if currentGene.orientation == '-' and int(pos.ID) > earliestStop[1] and pos.variant == variant:
								genomePosition 					= int(pos.ID)
								earliestStop 					= [pos,int(pos.ID)]
								positions[cultivar.ID] 			= pos
								averageVariantCoverage[cultivar.ID]= pos.genotype['DP'][0]
							
					#allVariants 		   +=1
					if currentGene.orientation == '+' and earliestStop == [None, inf]:
						continue
					if currentGene.orientation == '-' and earliestStop == [None, 0]:
						continue

					# evaluate Stop
					ploididtyScore	= 0
					lostAAScore		= 0
					functionLoss	= 0
					lostFunctions	= []
					
					# get genotype influence
					if earliestStop[0].genotype['GT'] ==  ['1/1'] or earliestStop[0].genotype['GT'] ==  ['2/2']:
						# this case is for homozygot ploidity
						ploididtyScore = 1
					else:
						# this case is for heterozygot ploidity
						ploididtyScore = 0.5
						
					# getlostAA inlfuence
					adaptedAaPosition = adaptAAPositionToMapping(earliestStop[0].pos_amino_acid, earliestStop[0].variant_names[0].AnnotationAaPos, earliestStop[0].variant_names[0].AnnotationDbPos)
					try:
						lostAAScore = round(1-(float(adaptedAaPosition)/float(currentVariant.proteinLength)), 3)
					except:
						try:
							if earliestStop[0].pos_amino_acid/3.> currentVariant.proteinLength:
								lostAAScore = 0
							else:
								lostAAScore = round(1-(float(earliestStop[0].pos_amino_acid/3.)/float(currentVariant.proteinLength)), 3)
						except:
							lostAAScore = -1
					
					# check what functions are affected
					if gene in comparenceResult:
						if variant in comparenceResult[gene].variants:
							t = comparenceResult[gene].variants
							for functionality  in t[variant].functions:
								if [cultivar.ID, earliestStop[0]] in functionality.varPos:
									lostFunctions.append(functionality.functionID)
									functionLoss = 1
									cultivar_scores[cultivar] = [lostAAScore, functionLoss, ploididtyScore, round(((lostAAScore+functionLoss+ploididtyScore)/3.),3), 'STP' ,lostFunctions]
					else:
						functionLoss = lostAAScore
						cultivar_scores[cultivar] = [lostAAScore, functionLoss, ploididtyScore, round(((lostAAScore+functionLoss+ploididtyScore)/3.),3), 'STP' ,lostFunctions]
									
				#print('OUT', variant)
				if len(cultivar_scores) > 0:
					affectedVariantCounter += 1
					variant_scores[variant] = cultivar_scores
			if len(variant_scores) > 0:
				gene_scores[gene] = [variant_scores, affectedVariantCounter, allVariants, genomePosition, positions, averageVariantCoverage]
			
		return gene_scores
		
	def disruptiveFrameShiftEffects(self, comparenceResult, cultivars, genes):
		"""
		Function to reflect the disruptive effects of  
		frameshifts. For this, we account for loss of known domains as 
		well as for the amount of lost amio-acids. 
		"""
		# import some mathematical opertations
		from math import log, inf
		gene_scores = {}
		
		for gene in genes:
			currentGene 			= genes[gene]
			variant_scores 		   	= {}
			affectedVariantCounter 	= 0
			allVariants			 	= len(currentGene.variants)
			#print(currentGene.variants)
			positions				= {}
			averageVariantCoverage	= {}
			
			for variant in currentGene.variants:
				currentVariant = currentGene.variants[variant]
				#print(currentVariant.ID)
				cultivar_scores = {}
				for cultivar in currentGene.cultivars:
					# stores the scores for the cultivars
					# initalize start position for first stoips depending on orientation
					if currentGene.orientation == '+':
						earliestFS = [None, inf]
					if currentGene.orientation == '-':
						earliestFS = [None, 0]

					genomePosition = 0
					fsCount = 0
					for vP in cultivar.variablePositions:
						position = cultivar.variablePositions[vP]
						for pos in position:
							#print(pos.types)			
							if 'Frameshift+1' in pos.types or 'Frameshift+2' in pos.types or 'Frameshift-1' in pos.types or 'Frameshift-2' in pos.types:
								fsCount += 1
							else:
								continue
							#print(pos.ID, pos.types)
							# get variant with the earliest Stop
							if currentGene.orientation == '+' and int(pos.ID) < earliestFS[1] and pos.variant == variant:
								genomePosition 					= int(pos.ID)
								positions[cultivar.ID] 			= pos
								earliestFS 					= [pos,int(pos.ID)]
								averageVariantCoverage[cultivar.ID]= pos.genotype['DP'][0]
							if currentGene.orientation == '-' and int(pos.ID) > earliestFS[1] and pos.variant == variant:
								genomePosition 					= int(pos.ID)
								earliestFS 					= [pos,int(pos.ID)]
								positions[cultivar.ID] 			= pos
								averageVariantCoverage[cultivar.ID]= pos.genotype['DP'][0]

					#if fsCount >= 2:					
						#print(fsCount)
					#allVariants 		   +=1
					if currentGene.orientation == '+' and earliestFS == [None, inf]:
						continue
					if currentGene.orientation == '-' and earliestFS == [None, 0]:
						continue
						
					# get potential reversion of frameshift
					shift =  None
					for i in earliestFS[0].types:
						if i[0:10] == 'Frameshift':
							shift =  int(i[10:])
							#print(shift)
					
					reversedShift = [None, 0]
					for vP in cultivar.variablePositions:
						position = cultivar.variablePositions[vP]
						
						for pos in position:			
							if 'Frameshift+1' in pos.types or 'Frameshift+2' in pos.types or 'Frameshift-1' in pos.types or 'Frameshift-2' in pos.types:
								t=0
							else:
								continue
							shift2  = 0
							for i in pos.types:
								if i[0:10] == 'Frameshift':
									shift2 =  int(i[10:])
									#print(shift2)
							if pos != earliestFS[0]:	
								#print(shift, shift2)
								try:		
									if (shift2 + shift) == 0 and abs(int(pos.ID) - int(earliestFS[0].ID)) > int(reversedShift[1]):
										reversedShift = [pos, pos.ID]
								except:
									continue
											
					# evaluate Stop
					ploididtyScore	= 0
					lostAAScore		= 0
					functionLoss	= 0
					lostFunctions	= []
					
					# get genotype influence
					if earliestFS[0].genotype['GT'] ==  ['1/1'] or earliestFS[0].genotype['GT'] ==  ['2/2']:
						# this case is for homozygot ploidity
						ploididtyScore = 1
					else:
						# this case is for heterozygot ploidity
						ploididtyScore = 0.5
						
					# getlostAA inlfuence
					adaptedAaPosition1 = adaptAAPositionToMapping(earliestFS[0].pos_amino_acid, earliestFS[0].variant_names[0].AnnotationAaPos, earliestFS[0].variant_names[0].AnnotationDbPos)
					if reversedShift == [None, 0]: 
						try:
							lostAAScore = round(1-(float(adaptedAaPosition1)/float(currentVariant.proteinLength)), 3)
						except:
							try:
								if earliestFS[0].pos_amino_acid/3.> currentVariant.proteinLength:
									lostAAScore = 0
								else:
									lostAAScore = round(1-(float(earliestFS[0].pos_amino_acid/3.)/float(currentVariant.proteinLength)), 3)
							except:
								lostAAScore = -1
					else:
						adaptedAaPosition2 = adaptAAPositionToMapping(reversedShift[0].pos_amino_acid, reversedShift[0].variant_names[0].AnnotationAaPos, reversedShift[0].variant_names[0].AnnotationDbPos)
						try:
							lostAAScore = round((abs(float(adaptedAaPosition1)-float(adaptedAaPosition2))/float(currentVariant.proteinLength)), 3)
							#print(variant, adaptedAaPosition1, adaptedAaPosition2, currentVariant.proteinLength, lostAAScore)
						except:
							#print(variant, 'EXCEPT')
							try:
								if earliestFS[0].pos_amino_acid/3.> currentVariant.proteinLength or reversedShift[0].pos_amino_acid/3.> currentVariant.proteinLength:
									lostAAScore = 0
									#print('0')
								else:
									lostAAScore = round((abs(float(earliestFS[0].pos_amino_acid/3.)-(float(reversedShift[0].pos_amino_acid/3.)))/float(currentVariant.proteinLength)), 3)
									#print(variant, float(earliestFS[0].pos_amino_acid/3.), float(reversedShift[0].pos_amino_acid/3.), currentVariant.proteinLength, lostAAScore)
							except:
								lostAAScore = -1
								#print('-')
					


					# check what functions are affected
					if gene in comparenceResult:
						if variant in comparenceResult[gene].variants:
							t = comparenceResult[gene].variants
							for functionality  in t[variant].functions:
								#print([cultivar.ID, earliestFS[0]], functionality.varPos)
								if [cultivar.ID, earliestFS[0]] in functionality.varPos:
									if reversedShift == [None, 0]:
										lostFunctions.append(functionality.functionID)
										functionLoss = 1
										cultivar_scores[cultivar] = [lostAAScore, functionLoss, ploididtyScore, round(((lostAAScore+functionLoss+ploididtyScore)/3.),3), 'FS' ,lostFunctions, reversedShift]
									else:
										function_split = functionality.functionID.split('_')
										#print(function_split, reversedShift[0].pos_amino_acid/3., int(function_split[-2]))
										#print(int(reversedShift[0].pos_amino_acid/3.) >= int(function_split[-2]))
										if int(reversedShift[0].pos_amino_acid/3.) >= int(function_split[-2]):# or (int(reversedShift[0].pos_amino_acid/3.) >= int(function_split[-2]) and int(reversedShift[0].pos_amino_acid/3.) <= int(function_split[-1])):
											lostFunctions.append(functionality.functionID)
											functionLoss = 1
											cultivar_scores[cultivar] = [lostAAScore, functionLoss, ploididtyScore, round(((lostAAScore+functionLoss+ploididtyScore)/3.),3), 'PFS' ,lostFunctions,reversedShift]
										else:
											functionLoss = 0
											cultivar_scores[cultivar] = [lostAAScore, functionLoss, ploididtyScore, round(((lostAAScore+functionLoss+ploididtyScore)/3.),3), 'PFS' ,lostFunctions,reversedShift]
					else:
						functionLoss = 0
						if reversedShift == [None, 0]:
							cultivar_scores[cultivar] = [lostAAScore, functionLoss, ploididtyScore, round(((lostAAScore+functionLoss+ploididtyScore)/3.),3), 'FS' ,lostFunctions,reversedShift]
						else:
							cultivar_scores[cultivar] = [lostAAScore, functionLoss, ploididtyScore, round(((lostAAScore+functionLoss+ploididtyScore)/3.),3), 'PFS' ,lostFunctions,reversedShift]
									
				if len(cultivar_scores) > 0:
					affectedVariantCounter += 1
					variant_scores[variant] = cultivar_scores
			if len(variant_scores) > 0:
				gene_scores[gene] = [variant_scores, affectedVariantCounter, allVariants, genomePosition, positions, averageVariantCoverage]
			
		return gene_scores		
			
	def masterSumScore(self, compareResult, settings):
		"""
		This function generates  a simple master score by summing over 
		all individual scores for the properties weighted by their user 
		definded degree of importance.
		"""
		checksum = 0
		for result in  compareResult:
			for score in ['maxScoreMutagenesis','diffScoreMutagenesis','maxScoreDomain','diffScoreDomain','maxScoreDnaNucCal','diffScoreDnaNucCal','maxScoreActiveSites','diffScoreActiveSites','maxScoreMembrane','diffScoreMembrane','maxScoreMotifs','diffScoreMotifs','diffScoreSecondary', 'maxScoreSecondary', 'diffScorePeptide', 'maxScorePeptide', 'diffScoreModifications', 'maxScoreModifications']: 
				tmp = 'sw'+score[0].upper()+score[1:]
				if score[0] == 'd':
					#print(compareResult[result].masterScore#*settings[tmp]*settings['swMasterDiff'])
					compareResult[result].masterScore += getattr(compareResult[result], score,0)*settings[tmp]*settings['swMasterDiff']
				else:
					#print(compareResult[result].masterScore*settings[tmp]*settings['swMasterMax'])
					compareResult[result].masterScore += getattr(compareResult[result],score,0)*settings[tmp]*settings['swMasterMax']
			checksum += compareResult[result].masterScore
			#print(compareResult[result].masterScore)
		#print('checksum: '+str(checksum))
		i = 0
	
	def printScoreTable(self, compareResult):
		for i in ['ID','maxScoreMutagenesis','diffScoreMutagenesis','maxScoreDomain','diffScoreDomain','maxScoreDnaNucCal','diffScoreDnaNucCal','maxScoreActiveSites','diffScoreActiveSites','maxScoreMembrane','diffScoreMembrane','maxScoreMotifs','diffScoreMotifs','diffScoreSecondary', 'maxScoreSecondary', 'diffScorePeptide', 'maxScorePeptide', 'diffScoreModifications', 'maxScoreModifications', 'masterScore']:
			print(i, end=',')
		print()
		for result in compareResult:
			print(result, end=',')
			for j  in ['maxScoreMutagenesis','diffScoreMutagenesis','maxScoreDomain','diffScoreDomain','maxScoreDnaNucCal','diffScoreDnaNucCal','maxScoreActiveSites','diffScoreActiveSites','maxScoreMembrane','diffScoreMembrane','maxScoreMotifs','diffScoreMotifs','diffScoreSecondary', 'maxScoreSecondary', 'diffScorePeptide', 'maxScorePeptide', 'diffScoreModifications', 'maxScoreModifications', 'masterScore']:
				print(getattr(compareResult[result],j,-1), end=',')
			print()

def adaptAAPositionToMapping(position, mapping1, mapping2):
	"""
	Function to return an amino acid position adapted to a certain
	numbering scheme or mappig (e.g. Annotation -> uniport mapping)
	"""
	position2 =  int((int(position)/3.)+0.8)
	position  = position2
	if mapping1 != None:
		for i in range(len(mapping1)):
			if position >= mapping1[i][0] and position <= mapping1[i][1]:
				diff = mapping1[i][0]-mapping2[i][0]
				return abs(position+diff)
	else:
		return None

def getMaxValueOfMatrix(matrix):
	maximum  = 0
	for i in matrix:
		if matrix[i] > maximum:
			maximum = matrix[i]
	return maximum
	
def stepFunction(values, start, step):
	summ = 0
	if values == []:
		return 0
	for i in values:
		summ += i*start
		start = start/2.
	return round(summ,3)
	
def maxPairwiseDistance(value):
	tmp = sorted(value, reverse = True)
	return tmp[0]-tmp[-1]
