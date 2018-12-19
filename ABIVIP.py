#!/usr/bin/env python3.6

__author__  = "Jens Loers"
__email__   = "jens.loers@uni-bielefeld.de"

import software.unifiedFileHandler as uFH
import software.analysis as analyser



def test():
	# test1: Load the pdb database
	fH 		= uFH.UnifiedFileHandler()
	#values	= fH.loadUniprotGFF('/vol/cluster-data/jloers/BIFCoG/data/ath/Ath_gold_review.gff')
	genes   = fH.loadAnnotationGFF('/vol/cluster-data/jloers/BIFCoG/data/ath/Araport11_GFF3_genes_transposons.201606.gff')
	
	fH.loadNavipOutput('/vol/cluster-data/jloers/BIFCoG/data/ath/All_VCF.vcf', genes, 'Nd')
	
	for gene in genes:
		if len(genes[gene].cultivars) > 0 and gene == 'At1g71920':
			print(gene, len(genes[gene].cultivars[0].variablePositions) )
			for i in genes[gene].cultivars[0].variablePositions:
				for j in genes[gene].cultivars[0].variablePositions[i]:				
					print(i, j.ref_amino_acid, j.alt_amino_acid, j.ref_base, j.alt_base, j.pos_amino_acid, j.types, j.variant)
					
def simpleAnalysis(param):
	print('LoadData')
	#['uniprotDB','-u'], ['annotationGFF','-a'], ['secondaryIDs', '-s'], ['annotationMapping', '-m'], ['navipOutput', '-n'], ['substitutionMatrices', '-sub'], ['secondaryStructureMatrices', '-sMat'], ['generalizedSubstitutionMatrix', '-gMat'], ['membraneSubstiutionMatrix', '-mMat']
	fH 		= uFH.UnifiedFileHandler()
	
	genes   = fH.loadAnnotationGFF(param['annotationGFF'])
	uniprot	= fH.loadUniprotGFF(param['uniprotDB'])
	
	
	# load secondary IDs necessary for correct  link between genes and uniprot database
	fH.loadSecondaryIDs(param['secondaryIDs'], genes)
	
	# load annoation mapping for amnino acid position (annotation -> uniprot
	fH.loadAnnotationMapping(param['annotationMapping'], genes)
	
	# load NAvip Output
	fH.loadNavipOutput(param['navipOutput'], genes, 'Nd')
	fH.loadNavipOutput(param['navipOutput'], genes, 'TE')

	# create an analysis object
	analysis = analyser.Analysis()
	print('GetAffectedFunctions')
	# start a simple comparison
	analysis.getAffectedFunctionalities(genes, uniprot)
	#exit()
	print('CompareCultivars')
	# compare between the culativars
	comparenceResult = analysis.compareCultivars(genes)
	#displayComparanceResult(comparenceResult)
	# create secondary structure analysis
	matrices = {}
	for i in param['secondaryStructureMatrices']:
		m1, n1 = fH.loadSubsitutionMatrix(param['substitutionMatrices']+i)
		matrices[n1] =  m1
	print('Start secondaryAnalysis')
	analysis.secondaryStructureAnalysis(comparenceResult, matrices, ['Nd', 'TE'])	
	
	matrices = {}
	for i in [param['generalizedSubstitutionMatrix']]:
		m1, n1 = fH.loadSubsitutionMatrix(param['substitutionMatrices']+i)
		matrices[n1] =  m1
	print('Start peptideAnalysis')
	analysis.peptideAnalysis(comparenceResult, matrices, ['Nd', 'TE'])
	print('Start bindingActiveSiteAnalysis')
	analysis.bindingActiveSiteAnalysis(comparenceResult, matrices, ['Nd', 'TE'])
	print('Start dnaNucleotideBindingAnalysis')
	analysis.dnaNucleotideBindingAnalysis(comparenceResult, matrices, ['Nd', 'TE'])
	print('Start domainAnalysisAnalysis')
	analysis.domainAnalysis(comparenceResult, matrices, ['Nd', 'TE'])
	print('Start mutagenesisyAnalysis')
	analysis.mutagenesisAnalysis(comparenceResult, matrices, ['Nd', 'TE'])		
	print('Start externalModificationsAnalysis')
	analysis.externalModifications(comparenceResult,['Nd', 'TE'])
	
	analysis.motifAnalysis(comparenceResult, ['Nd', 'TE'])
	matrices = {}
	for i in [param['membraneSubstiutionMatrix']]:
		m1, n1 = fH.loadSubsitutionMatrix(param['substitutionMatrices']+i)
		matrices[n1] =  m1
	print('Start membraneAnalysis')
	analysis.membraneAnalysis(comparenceResult, matrices, ['Nd', 'TE'])
	print('Calculate masterscore')
	analysis.masterSumScore(comparenceResult, param)
	
	analysis.printScoreTable(comparenceResult)
	
def testDict():
	fH 				= uFH.UnifiedFileHandler()
	matrix, name 	= fH.loadSubsitutionMatrix('/vol/cluster-data/jloers/BIFCoG/data/substitutionMatrices/koshi_goldstein=Helix=.mat')
	matrix, name 	= fH.loadSubsitutionMatrix('/vol/cluster-data/jloers/BIFCoG/data/substitutionMatrices/koshi_goldstein=Beta strand=.mat')
	matrix 			= []

def displayComparanceResult(comparenceResult):

	for result in comparenceResult:
		if result in comparenceResult:
			for variant in comparenceResult[result].variants:
				variantO = (comparenceResult[result].variants[variant])
				for function in variantO.functions:
					for pos in function.varPos:
						#if 'Stop gained' in pos[1].types:
						print(result)
						print('\t', variant)
						print('\t\t', function.functionID, function.properties )
						print('\t\t\t', pos, pos[1].ref_amino_acid, pos[1].alt_amino_acid, pos[1].pos_amio_acid_adapted, pos[1].types)

def displayDisruptedGenes(disruptedGenes, analyzedCultivars, genes, outfolder, all_genes, mode):
	# create header
	if mode == 'STP':
		outfile = open(outfolder+'disrupted_genes.csv', 'w')
	else:
		outfile = open(outfolder+'frameshifted_genes.csv', 'w')
	#outfile.write('Transcript,Chr,GenomePos,AffectedTranscriptNr,AllTranscriptsNr,TranscriptRatio,')
	outfile.write('Transcript,Chr,AffectedTranscriptNr,AllTranscriptsNr,TranscriptRatio,')
	for name in analyzedCultivars:
		#outfile.write('LA_'+name+',LF_'+name+',PL_'+name+',AS_'+name+',TYPE_'+name+',')
		outfile.write('GenomePos_'+name+',AAStartPos_'+name+',AAEndPos_'+name+',Cov_'+name+',LA_'+name+',LF_'+name+',PL_'+name+',AS_'+name+',TYPE_'+name+',')
	outfile.write('maxScore, diffScore, PA_ratio, SecIDs, INFO\n')
	for gene in disruptedGenes:
		geneObject = disruptedGenes[gene]
		#print(geneObject)
		present		= 0
		notPresent	= 0
		for variant in geneObject[0]:
			#print(variant)
			maxScore 	= 0
			minScore 	= 10
			#print(geneObject[0][variant])
			outfile.write(variant+','+genes[gene].chromosome+','+str(geneObject[1])+','+str(geneObject[2])+','+str(round((geneObject[1])/float(geneObject[2]), 3)))
			variantObject = geneObject[0][variant]
			allFunctions = []
			#print(geneObject)
			for cultivar in analyzedCultivars:
				j = False
				for stats in variantObject:
					if stats.ID == cultivar:
						j = True
						present += 1
						if mode == 'STP':
							aaEndPos = all_genes[gene].variants[variant].proteinLength
						else:
							#print(variantObject[stats][6])
							if variantObject[stats][6] == [None, 0]:
								aaEndPos = all_genes[gene].variants[variant].proteinLength
							else:
								aaEndPos = int(variantObject[stats][6][0].pos_amino_acid/3.)
						if geneObject[4][cultivar].pos_amio_acid_adapted != None:
							aaPos	= geneObject[4][cultivar].pos_amio_acid_adapted
						else:
							aaPos	= int(geneObject[4][cultivar].pos_amino_acid/3.)
						outfile.write(','+geneObject[4][cultivar].ID+','+str(aaPos)+','+str(aaEndPos)+','+str(geneObject[4][cultivar].genotype['DP'][0])+','+str(variantObject[stats][0])+','+str(variantObject[stats][1])+','+str(variantObject[stats][2])+','+str(variantObject[stats][3])+','+str(variantObject[stats][4]))
						for i in variantObject[stats][5]:
							if i not in allFunctions:
								allFunctions.append(i)
					if variantObject[stats][3] > maxScore:
						maxScore = variantObject[stats][3]
					if variantObject[stats][3] < minScore:
						minScore = variantObject[stats][3]
				if j == False:
					coverage = None
					for i in all_genes[gene].cultivars:
						if i.ID == cultivar:
							coverage  = i.getAverageCoverage()
					if coverage == None:
						outfile.write(','+str('-')+','+str('-')+','+str('-')+','+str('-')+','+str(0)+','+str(0)+','+str(0)+','+str(0)+','+str('-'))
					else:
						outfile.write(','+str('-')+','+str('-')+','+str('-')+','+str(coverage)+','+str(0)+','+str(0)+','+str(0)+','+str(0)+','+str('-'))
					minScore = 0					
					notPresent += 1
			# Print DiffScore, MaxScore and P/A ratio
			outfile.write(','+str(maxScore)+','+str(maxScore-minScore)+','+str(round((present/(float(present)+float(notPresent))),2))+',')
			for i in genes[gene].secondaryID:		
				outfile.write('+'+i)
			outfile.write(',')	
			for i in allFunctions:		
				outfile.write('+'+i)
			outfile.write('\n')

def displayAnalysisresult(result, all_genes, cultivars, outfolder):
	'''
	Traverse the result tree and write the resut file for the corresponding analysis
	'''
	outfile = open(outfolder, 'w')
	
	outfile.write('GeneID\tUniprotID\tTranscript\tChrom\tAffectedTK\tAllTK\tTK_Ratio\tTag\tTagStart\tTagEnd\tTagID\tMax_Score\tDiff_Score\tPA_Ratio\t')
	for name in  cultivars:
		outfile.write(name+'_TYPE\t'+name+'_GenomePositionStart\t'+name+'_GenomePositionEnd\t'+name+'_AAPositionStart\t'+name+'_AAPositionEnd\t'+name+'_RefAA\t'+name+'_AltAA\t'+name+'_Cov\t'+name+'_Zyg\t'+name+'_Site_Score\t'+name+'_Function_Score\t'+name+'_Transkript_Score\t'+name+'_Cleav\t')
	outfile.write('\n')
	for gene in result:
		#print(gene)
		affectedTranskripts = 0
		for transcript in gene[4]:
			if transcript[1] != []:
				affectedTranskripts += 1
		for transcripts in gene[4]:
			if transcript[1] != []:
				for function in transcript[1]:
					allPositions = []
					for cultivar in function[3]:
						for i in function[3][cultivar]:
							#print(i)
							try:
								if i[0] not in allPositions:
									allPositions.append(i[0])
							except:
								continue
					#print(allPositions)
					#print(transcript[1][3])
					printString1 = ''
					
					
					for position in allPositions:
						cultivarCount = 0
						printString = ''
						#print(position)
						if gene[6] == 0:
							continue
						#outfile.write(gene[0]+'\t'+gene[1]+'\t'+transcripts[0]+'\t'+gene[2]+'\t'+str(affectedTranskripts)+'\t'+str(gene[3])+'\t'+str(round(affectedTranskripts/float(gene[3]),2))+'\t'+function[0]+'\t'+function[1].split('_')[3]+'\t'+function[1].split('_')[4]+'\t'+str(function[2])+'\t')
						printString1 = gene[0]+'\t'+gene[1]+'\t'+transcripts[0]+'\t'+gene[2]+'\t'+str(affectedTranskripts)+'\t'+str(gene[3])+'\t'+str(round(affectedTranskripts/float(gene[3]),2))+'\t'+function[0]+'\t'+function[1].split('_')[3]+'\t'+function[1].split('_')[4]+'\t'+str(function[2])+'\t'+str(gene[6])+'\t'+str(gene[7])+'\t'
						for cultivar in cultivars:#function[3]:
							#print(cultivar)
							if function[3][cultivar] == []:
								try:
									coverage = None
									for i in all_genes[gene[0]].cultivars:
										if i.ID == cultivar:
											coverage = i.getAverageCoverage()
									if coverage != None:
										printString += '-\t-\t-\t-\t-\t-\t-\t'+str(coverage)+'\t0\t-\t-\t-\t-\t'
									else:
										printString += '-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t'
								except:
										printString += '-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t'
								continue
							leap2 = 0	
							for i in function[3][cultivar]:
								#print('\t', i)
								if i[0] != position:
									leap2 = 0
									continue
								#print(i)
								genType = 0
								try: 
									#GT = function[3][cultivar][4]['GT']
									GT = i[4]['GT']
									if GT == ['1/1'] or GT == ['2/2']:
										genType = 1
									if GT == ['0/1'] or GT == ['0/2'] or GT == ['1/2']:
										genType = 0.5
								except:
									genType = '-'
								#print(genType)
								#print(function[3][cultivar])
								#print(transcript[1][0][3][cultivar])
								leap = 0
								try: 								
									if i[0] == position:
										functionIndex = 0
										#print(function[1])
										for j in range(len(transcript[1])):
											if transcript[1][j][1] == function[1]:
												functionIndex = j							
										#outfile.write(function[3][cultivar][1][0]+'\t'+function[3][cultivar][5]+'\t'+function[3][cultivar][6]+'\t'+str(function[3][cultivar][7])+'\t'+str(function[3][cultivar][8])+'\t'+str(function[3][cultivar][2])+'\t'+function[3][cultivar][3]+'\t'+str(function[3][cultivar][4]['DP'][0])+'\t'+str(genType)+'\t'+str(function[3][cultivar][9])+'\t'+str(round(transcript[1][0][4][cultivar],3))+'\t'+str(gene[5])+'\t'+str(gene[6])+'\t'+str(function[3][cultivar][10])+'\t')
										#printString += i[1][0]+'\t'+i[5]+'\t'+i[6]+'\t'+str(i[7])+'\t'+str(i[8])+'\t'+str(i[2])+'\t'+i[3]+'\t'+str(i[4]['DP'][0])+'\t'+str(genType)+'\t'+str(i[9])+'\t'+str(round(transcript[1][0][4][cultivar],3))+'\t'+str(gene[5][cultivar])+'\t'+str(i[10])+'\t'
										printString += i[1][0]+'\t'+i[5]+'\t'+i[6]+'\t'+str(i[7])+'\t'+str(i[8])+'\t'+str(i[2])+'\t'+i[3]+'\t'+str(i[4]['DP'][0])+'\t'+str(genType)+'\t'+str(i[9])+'\t'+str(round(transcript[1][functionIndex][4][cultivar],3))+'\t'+str(gene[5][cultivar])+'\t'+str(i[10])+'\t'										
										leap  = 1
										leap2 = 1
										cultivarCount += 1
										break

								except:
									#print(i[0], i)
									continue

							if leap2 == 0:
								printString += '-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t'
						if cultivarCount > 0:
							outfile.write(printString1+str(round(cultivarCount/len(cultivars),2))+'\t'+printString+'\n')
	
def traverseUniprotEntry():
	# test1: Load the pdb database
	fH 		= uFH.UnifiedFileHandler()
	uniprot	= fH.loadUniprotGFF('/vol/cluster-data/jloers/BIFCoG/data/ath/Ath_gold_review.gff')
	
	for entry in uniprot:
		for function in uniprot[entry].functions:
			if function[0] in ['Mutagenesis']:#['Glycosylation','Lipidation','Zinc finger','Modified residue']:##['Peptide','Signal peptide','Transit peptide','Propeptide']:
				print(function, entry)

def simpleAnalysisVitis(param):
	"""
	Function to analyse the several aspects of V. Vinifera
	"""
	
	# load Database and GFF file
	print('LoadData')
	fH 		= uFH.UnifiedFileHandler()
	genes   = fH.loadAnnotationGFF(param['annotationGFF'])
	uniprot	= fH.loadUniprotGFF(param['uniprotDB'])

	# load secondary IDs necessary for correct  link between genes and uniprot database
	fH.loadSecondaryIDs(param['secondaryIDs'], genes)
	
	# load annoation mapping for amnino acid position (annotation -> uniprot
	fH.loadAnnotationMapping(param['annotationMapping'], genes)
	
	# load ProteinLength	
	fH.loadProteinLength(param['proteinLengthFile'], genes)
	
	from os import listdir
	from os.path import isfile, join
	onlyfiles = [f.split('.')[0] for f in listdir(param['navipOutputFolder']) if isfile(join(param['navipOutputFolder'], f))]
	if '' in onlyfiles: onlyfiles.remove('')
	print(onlyfiles)
	for filename in  onlyfiles:
		fH.loadNavipOutputBeta(param['navipOutputFolder']+filename+'.vcf', genes, filename)
		print('Cultivar '+filename+' done!') 

	# create an analysis object
	analysis = analyser.Analysis()
	print('GetAffectedFunctions')
	# start a simple comparison
	analysis.getAffectedFunctionalities(genes, uniprot)

	# start a simple comparison
	analysis.getAffectedFunctionalities(genes, uniprot)

	print('CompareCultivars')
	# compare between the culativars
	comparenceResult, genesDisrupted, genesFrameShifted = analysis.compareCultivars(genes)

	matrices = {}
				

#	print('Start disruptive effect analysis')
	genesFrameShifted2 = {}
	for i in genesFrameShifted:
		if i not in genesDisrupted:
			genesFrameShifted2[i] = genesFrameShifted[i]
	#print(len(genesFrameShifted), len(genesFrameShifted2), len(genesDisrupted))

	print('Start secondaryAnalysis')	
#	disruptedStopGenes	= analysis.disruptiveStopEffects(comparenceResult, onlyfiles, genesDisrupted)
#	disruptedFSGenes	= analysis.disruptiveFrameShiftEffects(comparenceResult, onlyfiles, genesFrameShifted2)
#	displayDisruptedGenes(disruptedStopGenes, onlyfiles, genes, param['outfolder'], genes, 'STP')
#	displayDisruptedGenes(disruptedFSGenes, onlyfiles, genes, param['outfolder'], genes, 'FS')
	
	for i in param['secondaryStructureMatrices']:
		m1, n1 = fH.loadSubsitutionMatrix(param['substitutionMatrices']+i)
		matrices[n1] =  m1
#	print('Start secondaryAnalysis')
#	secondaryResult = analysis.secondaryStructureAnalysis(comparenceResult, matrices, onlyfiles)	
#	displayAnalysisresult(secondaryResult, genes, onlyfiles, param['outfolder']+'secondary_structure.csv')

	matrices = {}
	for i in [param['generalizedSubstitutionMatrix']]:
		m1, n1 = fH.loadSubsitutionMatrix(param['substitutionMatrices']+i)
		matrices[n1] =  m1
#	print('Start peptideAnalysis')
#	peptideAnalysisResult = analysis.peptideAnalysis(comparenceResult, matrices, onlyfiles)
#	displayAnalysisresult(peptideAnalysisResult, genes, onlyfiles, param['outfolder']+'peptide_effects.csv')
	
#	print('Start bindingActiveSiteAnalysis')
#	bindingSiteAnalysisResult= analysis.bindingActiveSiteAnalysis(comparenceResult, matrices, onlyfiles)
#	print('StartDisplay')
#	displayAnalysisresult(bindingSiteAnalysisResult, genes, onlyfiles, param['outfolder']+'binding_sites.csv')
	
	print('Start dnaNucleotideBindingAnalysis')
	dnaNucleotideBindingDomainsAnalysisResult = analysis.dnaNucleotideBindingAnalysis(comparenceResult, matrices, onlyfiles)
	displayAnalysisresult(dnaNucleotideBindingDomainsAnalysisResult, genes, onlyfiles, param['outfolder']+'dna_nucleotide_binding_domains.csv')

#	print('Start domainAnalysisAnalysis')
#	domainResult = analysis.domainAnalysis(comparenceResult, matrices, onlyfiles)
#	displayAnalysisresult(domainResult, genes, onlyfiles, param['outfolder']+'domains.csv')
		
#	print('Start mutagenesisyAnalysis')
#	mutagenesisResult = analysis.mutagenesisAnalysis(comparenceResult, matrices, onlyfiles)	
#	displayAnalysisresult(mutagenesisResult, genes, onlyfiles, param['outfolder']+'mutagenesis.csv')
	
#	print('Start externalModificationsAnalysis')
#	modificationResult = analysis.externalModifications(comparenceResult,onlyfiles)
#	displayAnalysisresult(modificationResult, genes, onlyfiles, param['outfolder']+'external_modifications.csv')

#	print('Start motifAnalysis')	
#	motifResult = analysis.motifAnalysis(comparenceResult, onlyfiles)
#	displayAnalysisresult(motifResult, genes, onlyfiles, param['outfolder']+'motifs.csv')
	
#	print('Start membraneAnalysis')	
	matrices = {}
	for i in [param['membraneSubstiutionMatrix']]:
		m1, n1 = fH.loadSubsitutionMatrix(param['substitutionMatrices']+i)
		matrices[n1] =  m1
#	membraneResult = analysis.membraneAnalysis(comparenceResult, matrices, onlyfiles)
#	displayAnalysisresult(membraneResult, genes, onlyfiles, param['outfolder']+'membranes.csv')
#	print('Calculate masterscore')
#	analysis.masterSumScore(comparenceResult, param)
	
	#analysis.printScoreTable(comparenceResult)

	#displayComparanceResult(comparenceResult)

def help():
	print('Flags:')
	print('--uniprotDB [path_to_file] or -u [path_to_file]')

def loadConfigFile(parameters):
	"""
	First function which should be called after starting the pogramm.
	Takes the input of an config file and translates them into variables
	Hirachie:
	Flag > Config File > Default Value > Warning and Exit 
	"""
	# get all flags and safe them in flags dict
	flags = {}

	if len(parameters)%2 != 0:
		print('Wrong console input: Call hast to be like: -flag [variable]], ...')
		exit()
	
	for i in range(0,len(parameters),2):
		flags[parameters[i]] = parameters[i+1]
	
	# get all entries in the configuŕation File
	config = {}
	if '-conf' not in flags:
		print('\n Could not load config File. Please provide a suited config file with -conf [config_file]\n')
		exit()
		
	with open(flags['-conf'], "r") as DataFile:
		entries	= DataFile.read().split('\n')
		for entry in entries:
			if len(entry) == 0 or entry[0] == '#':
				continue
			tmp = entry.split(':')
			if len(tmp) > 1:
				tmp[1] = tmp[1].strip(' ')
				tmp[1] = tmp[1].strip('\t')
				tmp[1] = tmp[1].strip(' ')
				tmp[1] = tmp[1].strip('\t')
				tmp[0] = tmp[0].strip(' ')
				tmp[0] = tmp[0].strip('\t')
				tmp[0] = tmp[0].strip(' ')
				tmp[0] = tmp[0].strip('\t')
				# test for array condition:
				tmp2 =  tmp[1].split(',')
				if len(tmp2) > 1:
					config[tmp[0]] =  tmp2
				else:
					config[tmp[0]] = tmp[1]
					
	## create dictionary with all variables
	settings = {}
	# set path variables 
	for key in [['uniprotDB','-u'], ['annotationGFF','-a'], ['secondaryIDs', '-s'], ['annotationMapping', '-m'], ['navipOutput', '-n'], ['navipOutputFolder', '-nf'], ['substitutionMatrices', '-sub'], ['secondaryStructureMatrices', '-sMat'], ['generalizedSubstitutionMatrix', '-gMat'], ['membraneSubstiutionMatrix', '-mMat'] , ['proteinLengthFile', 'pL'], ['outfolder', '-o']]:
		if '--'+key[0] in flags:
			settings[key[0]] = flags['--'+key[0]]
		elif key[1] in flags:
			settings[key[0]] = flags[key[1]]
		else:
			if key[0] in config:
				settings[key[0]] = config[key[0]]
			else:
				print('\nCould not find variable for '+key[0]+'. Provide a correct flag or config file or start with -h for further help.\n')
				exit()
				
	# set scorringWeights		
	for key in [['swDiffScoreSecondary', '-x1'], ['swMaxScoreSecondary', '-y1'],['swDiffScorePeptide', '-x2'], ['swMaxScorePeptide', '-y2'],['swDiffScoreModifications', '-x3'], ['swMaxScoreModifications', '-y3'],['swDiffScoreMotifs', '-x4'], ['swMaxScoreMotifs', '-y4'],['swDiffScoreMembrane', '-x5'], ['swMaxScoreMembrane', '-y5'],['swDiffScoreActiveSites', '-x6'], ['swMaxScoreActiveSites', '-y6']	,['swDiffScoreDnaNucCal', '-x7'], ['swMaxScoreDnaNucCal', '-y7']	,['swDiffScoreDomain', '-x8'], ['swMaxScoreDomain', '-y8']	,['swDiffScoreMutagenesis', '-x9'], ['swMaxScoreMutagenesis', '-y9'], ['swMasterDiff', '-ydiff'], ['swMasterMax', '-ymax'] ]:
		if '--'+key[0] in flags:
			settings[key[0]] = flags['--'+key[0]]
		elif key[1] in flags:
			settings[key[0]] = float(flags[key[1]])
		else:
			if key[0] in config:
				settings[key[0]] = float(config[key[0]])
			else:
				print('\nCould not find variable or flag for '+key[0]+'. Use default value of 1.\n')
				settings[key[0]] = float(1)
							
	return settings

if __name__ == '__main__':
	from sys import argv
	
	# load config file
	settings = loadConfigFile(argv[1:])
	#start vitis analysis
	simpleAnalysisVitis(settings)
	
