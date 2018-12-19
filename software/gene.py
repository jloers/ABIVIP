#!/usr/bin/env python3.6

__author__  = "Jens Loers"
__email__   = "jens.loers@uni-bielefeld.de"

class Gene:
	"""
	Datastructure to store all information and calculation results
	for a certain gene of interest. Also allows basic manipulation 
	of the data like asigning additional identifier etc.
	The datastructure has a special level (tree-) hirachie:
	1. Level:	Store general information about the gene itself
	2. Level:	Store information about the whole cultivar
	3. Level:	Store variant specific information (exhausting,
				better with representative transcript, the one which
				fits best to uniprot)
	4. Level:	Store variable position specific information
	
	TODO: It might be better for computational efficiency to just choose
	the representative Variant and analyse just this one
	"""
	
	def __init__(self, ID):
		"""
		inititalize and store all information of interest in a gene 
		entry. Might be extendend over developemental progress
		"""
		
		# information received from Annotatino file
		self.ID				=	ID		# stores the ID 
		self.chromosome		= 	0		# stores which chromosome
		self.start			= 	0		# stores start postion of gene
		self.end			= 	0		# stores end position of gene
		self.orientation	= 	''		# stores orientation of gene ( + or -) 
		
		# link to cultivar class
		self.cultivars		= 	[]		# List of cultivars (and below variable positions)
		self.variants		=	{}		# Dict of variants
		
		# additional gene parameters
		self.secondaryID	=	[]		# potentialMapping ID
		self.masterScore	=	0		# quantifies influences on biological effect
		
class Cultivar:
	"""
	inititalize and store all information of interest of a Cultivar
	"""
	def __init__(self, ID):
		self.ID					=	ID		# store the ID (Name of the Cultivar
		self.variablePositions	= 	{}		# List of variable postions
		self.coverage			= 	None
		
	def getAverageCoverage(self):
		coverage = 0
		counter  = 0
		for vP in self.variablePositions:
			for pos in self.variablePositions[vP]:
				try:
					coverage += pos.genotype['DP'][0]
					counter  += 1
				except:
					continue
		if counter > 0:
			self.coverage = round(coverage/float(counter))
			return self.coverage
		else:
			return None

class Variant: 		
	"""
	inititalize and store all information of interest of a Gene-Variant
	"""
	def __init__(self, ID):
		self.ID					=	ID		# store the ID
		self.canonical			= 	None	# True: cannonical, False: non-cannonical -> in this def. mapping of 100% with pdb entry
		self.AnnotationAaPos	=	None	# contains positions wich are covered by the sequences in the database (uniprot usually)
		self.AnnotationDbPos	=	None	# contains corresponding positions in the uniprot numbering
		self.AnnotationIdentity =	None	# refelcts how many amino acids mapped to the annotation
		self.start				=	0		# start of the mRNA Transcript for the variant (nucleotide position)
		self.end				= 	0		# start of the mRNA Transcript for the variant (nucleotide position)
		self.exons				= 	[]		# Exons of the variant
		self.proteinLength		= 	0		# Length of the protein transcript
		
class VariationPosition:
	"""
	Datastructure to store site specific information for every variant
	Position
	"""	
	def __init__(self, ID):
		self.ID						=	ID		# store the ID
		
		self.coverage				= 	0		# List of cultivars,
		
		self.genotype				= 	{}		# Dict, stores all information from the VCF 4 extra column: GT:AD:DP:GQ:PL	
		
		self.ref_amino_acid			=	None	# reference amino acid
		self.alt_amino_acid			=	None	# alternatice amio acid
		self.ref_base				= 	None	# reference base
		self.alt_base				=	None	# alternative base
		self.types					= 	[]		# kinds of effects
		self.pos_amino_acid			=	0		# position of aa with regard to effect predition tool
		self.pos_amio_acid_adapted 	= 	None
		self.delStart				=	None
		self.delEnd					=	None
		
		self.variant				= 	''		# name of variant
		self.variant_names			= 	[]		# pointer to variant object
		
		self.amio_acid_freq			=	[0.0, 0.0, 0.0, 0.0]

		self.function				=	[]		# contains all functions affected by this SNP or INDELS + description List 
		self.frameShiftVariants 	=	[]		# contains all alternative frame variants (e.g. if a vairant is spliced differently
		
		# scoring properties
		self.secondaryScore 		=	0
		self.motifScore				= 	0
		self.membraneScore			=	0


