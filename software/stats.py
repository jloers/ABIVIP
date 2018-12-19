#!/usr/bin/env python3.6

__author__  = "Jens Loers"
__email__   = "jens.loers@uni-bielefeld.de"

class GeneStat:
	"""
	This datastructure brings the data in a more suitable format after
	the comparency with a database and between cultivars took place. 
	The structure is tree like with followig hierarchy:
	
	GeneStat: (reference to gene)
		VariantStat: (reference to variant)
			FunctionStat: (reflects a function in the given gene)
				This Class contains a list of sites per cultivar, which
				makes the data 'with one look' comparable
				
	The datastructure is easy to traverse and brings all the information
	we need for scorring and comparence in one place  and right order 
	as well as in dependency of the affected functionality. Therefore
	it is convenient to indroduce another datastructure with features
	the already introduced gene class is not able to perform.
	"""
	
	def __init__(self, gene_ID, gene_pointer):
		"""
		Inititalizes the GeneStat Class and creates poiter to the
		corresponding gene
		"""
		self.ID			= gene_ID
		self.gene 		= gene_pointer
		self.variants	= {}
		
		# scores
		self.diffScoreSecondary		= 0		# Score for influence on secondary elements
		self.maxScoreSecondary		= 0		# Score for influence on secondary elements
		self.diffScorePeptide		= 0		# Score for influence on peptides
		self.maxScorePeptide		= 0		# Score for influence on peptides	
		self.diffScoreModifications	= 0		# Score for influence on modifications
		self.maxScoreModifications	= 0		# Score for influence on modifications	
		self.diffScoreMotifs		= 0		# Score for influence on motifs
		self.maxScoreMotifs			= 0		# Score for influence on motifs	
		self.diffScoreMembrane		= 0		# Score for influence on membranes
		self.maxScoreMembrane		= 0		# Score for influence on membranes	
		self.diffScoreActiveSites	= 0		# Score for influence on activeSites
		self.maxScoreActiveSites	= 0		# Score for influence on activeSites	
		self.diffScoreDnaNucCal		= 0		# Score for influence on dna, nucleotide and calcium binding domains
		self.maxScoreDnaNucCal		= 0		# Score for influence on dna, nucleotide and calcium binding domains	
		self.diffScoreDomain		= 0		# Score for influence on domains
		self.maxScoreDomain			= 0		# Score for influence on domains
		self.diffScoreMutagenesis	= 0		# Score for influence by mutagenesis experiments
		self.maxScoreMutagenesis	= 0		# Score for influence by mutagenesis experiments
		
		self.masterScore			= 0		# Score with reflects the 'interest' on an gene for an molecular biologist				
		
	def addVariantStat(self, variant_ID, variant_pointer):
		"""
		Create an variant instance of VariantStat and add it to self.variants
		"""
		varStat 			=  VariantStat(variant_ID)
		varStat.variant		=  variant_pointer
		self.variants[variant_ID] = varStat
	 
class VariantStat:
	"""
	Contains pointer to variant class as well as name of variants and 
	a list of FunctionStat objectes with are pointing towards the 
	certain instance
	"""
	
	def __init__(self, variant_ID):
		"""
		Inititalizes the VariantStat Class and creates pointer to the
		corresponding variant
		"""
		self.name 				= variant_ID		# Contains variant name
		self.variant			= None				# Shall contain pointer to variant
		self.functions			= []				# contains list of FunctionStat Objects
		
		# scores
		self.diffScore			= 0
		self.maxScore			= 0
		self.maxScorePeptide	= 0
		self.diffScorePeptide	= 0
		self.peptideScore		= {}
	
	def addVariantStatEmpty(self, function_ID, properties):
		"""
		Create a FunctionStat Instance and add it to self.functions,
		but without filling
		"""
		functionStat  			=  FunctionStat(function_ID)
		functionStat.properties =  properties 
		self.function.append(functionStat)
		
	def addVariantStat(self, function_ID, properties, content):
		"""
		Create a FunctionStat Instance and add it to self.functions
		"""
		functionStat  			= FunctionStat(function_ID)
		functionStat.properties = properties
		functionStat.varPos 	= content
		self.functions.append(functionStat)
	
class FunctionStat:
	"""
	Contains a list of cultivars together with  a pointer to the
	variant positions which affect this functionality 
	"""
	
	def __init__(self, function_ID):
		"""
		Inititalizes the FunctionStat  method, which stores cultivar and 
		VarPos pointer
		"""
		self.functionID = function_ID
		self.uniProtID	= None
		self.properties = None
		self.varPos 	= []				# [[cult1, VarPos1], [cult, VarPos2], ...]
		
		# scorring properties
		secondaryScore		= None				#
		peptideScore		= None				#	
		modificationScore	= None
		motifScore			= None
		membraneScore		= None
		activeSiteScore		= None
		dnaNucBindingScore	= None
		domainScore			= None
		mutagenesisScore	= None
