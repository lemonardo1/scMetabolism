"""
Gene Ontology (GO) analysis functionality for scMetabolism.
"""

import os
import requests
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Set, Union
import warnings
from concurrent.futures import ThreadPoolExecutor
import json
import gzip
from io import StringIO


class GOAnalysis:
    """
    Gene Ontology analysis for single-cell data.
    
    This class provides methods to download GO annotations, filter GO terms,
    and perform enrichment analysis using GO gene sets.
    """
    
    def __init__(self, organism: str = "human", cache_dir: str = None):
        """
        Initialize GO analysis.
        
        Parameters:
        -----------
        organism : str
            Organism name ("human", "mouse", "rat", etc.)
        cache_dir : str, optional
            Directory to cache downloaded GO files
        """
        self.organism = organism.lower()
        self.cache_dir = cache_dir or os.path.expanduser("~/.scmetabolism/go_cache")
        self.go_annotations = None
        self.go_terms = None
        
        # Create cache directory
        os.makedirs(self.cache_dir, exist_ok=True)
        
        # Organism-specific settings
        self.organism_settings = {
            "human": {
                "taxon_id": "9606",
                "annotation_url": "http://current.geneontology.org/annotations/goa_human.gaf.gz",
                "gene_info_url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
            },
            "mouse": {
                "taxon_id": "10090", 
                "annotation_url": "http://current.geneontology.org/annotations/mgi.gaf.gz",
                "gene_info_url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz"
            },
            "rat": {
                "taxon_id": "10116",
                "annotation_url": "http://current.geneontology.org/annotations/rgd.gaf.gz", 
                "gene_info_url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Rattus_norvegicus.gene_info.gz"
            }
        }
        
    def download_go_annotations(self, force_download: bool = False) -> pd.DataFrame:
        """
        Download GO annotations for the specified organism.
        
        Parameters:
        -----------
        force_download : bool
            Whether to force re-download even if cached file exists
            
        Returns:
        --------
        pd.DataFrame
            GO annotations dataframe
        """
        
        if self.organism not in self.organism_settings:
            raise ValueError(f"Organism '{self.organism}' not supported. "
                           f"Available: {list(self.organism_settings.keys())}")
        
        settings = self.organism_settings[self.organism]
        cache_file = os.path.join(self.cache_dir, f"{self.organism}_go_annotations.tsv")
        
        # Load from cache if exists and not forcing download
        if os.path.exists(cache_file) and not force_download:
            print(f"Loading cached GO annotations for {self.organism}...")
            self.go_annotations = pd.read_csv(cache_file, sep='\t')
            return self.go_annotations
        
        print(f"Downloading GO annotations for {self.organism}...")
        
        try:
            # Download GAF file
            response = requests.get(settings["annotation_url"], stream=True)
            response.raise_for_status()
            
            # Parse GAF file
            annotations = []
            
            with gzip.open(response.raw, 'rt') as f:
                for line in f:
                    if line.startswith('!'):  # Skip comments
                        continue
                        
                    fields = line.strip().split('\t')
                    if len(fields) < 15:
                        continue
                        
                    # Extract relevant fields
                    db_object_symbol = fields[2]  # Gene symbol
                    go_id = fields[4]  # GO ID
                    evidence_code = fields[6]  # Evidence code
                    aspect = fields[8]  # Aspect (P, F, C)
                    db_object_name = fields[9]  # Gene name
                    db_object_synonym = fields[10]  # Synonyms
                    
                    # Skip certain evidence codes if needed
                    if evidence_code in ['IEA']:  # Electronic annotations
                        continue
                        
                    annotations.append({
                        'gene_symbol': db_object_symbol,
                        'go_id': go_id,
                        'evidence_code': evidence_code,
                        'aspect': aspect,
                        'gene_name': db_object_name,
                        'synonyms': db_object_synonym
                    })
            
            # Create DataFrame
            self.go_annotations = pd.DataFrame(annotations)
            
            # Save to cache
            self.go_annotations.to_csv(cache_file, sep='\t', index=False)
            print(f"Downloaded {len(self.go_annotations)} GO annotations")
            
            return self.go_annotations
            
        except Exception as e:
            print(f"Error downloading GO annotations: {e}")
            raise
    
    def download_go_terms(self, force_download: bool = False) -> pd.DataFrame:
        """
        Download GO term definitions.
        
        Parameters:
        -----------
        force_download : bool
            Whether to force re-download
            
        Returns:
        --------
        pd.DataFrame
            GO terms dataframe
        """
        
        cache_file = os.path.join(self.cache_dir, "go_terms.tsv")
        
        # Load from cache if exists
        if os.path.exists(cache_file) and not force_download:
            print("Loading cached GO terms...")
            self.go_terms = pd.read_csv(cache_file, sep='\t')
            return self.go_terms
        
        print("Downloading GO terms...")
        
        try:
            # Download GO basic OBO file
            obo_url = "http://current.geneontology.org/ontology/go-basic.obo"
            response = requests.get(obo_url)
            response.raise_for_status()
            
            # Parse OBO file
            terms = []
            current_term = {}
            
            for line in response.text.split('\n'):
                line = line.strip()
                
                if line == '[Term]':
                    if current_term:
                        terms.append(current_term)
                    current_term = {}
                elif line.startswith('id: GO:'):
                    current_term['go_id'] = line.split('id: ')[1]
                elif line.startswith('name: '):
                    current_term['name'] = line.split('name: ')[1]
                elif line.startswith('namespace: '):
                    current_term['namespace'] = line.split('namespace: ')[1]
                elif line.startswith('def: '):
                    current_term['definition'] = line.split('def: ')[1].split('"')[1]
                elif line.startswith('is_obsolete: true'):
                    current_term['obsolete'] = True
            
            # Add last term
            if current_term:
                terms.append(current_term)
            
            # Create DataFrame
            self.go_terms = pd.DataFrame(terms)
            
            # Filter out obsolete terms
            if 'obsolete' in self.go_terms.columns:
                self.go_terms = self.go_terms[~self.go_terms['obsolete'].fillna(False)]
            
            # Save to cache
            self.go_terms.to_csv(cache_file, sep='\t', index=False)
            print(f"Downloaded {len(self.go_terms)} GO terms")
            
            return self.go_terms
            
        except Exception as e:
            print(f"Error downloading GO terms: {e}")
            raise
    
    def create_go_gene_sets(
        self,
        aspects: List[str] = ["biological_process", "molecular_function", "cellular_component"],
        min_genes: int = 5,
        max_genes: int = 500,
        evidence_codes: Optional[List[str]] = None
    ) -> Dict[str, List[str]]:
        """
        Create GO gene sets for pathway analysis.
        
        Parameters:
        -----------
        aspects : List[str]
            GO aspects to include ("biological_process", "molecular_function", "cellular_component")
        min_genes : int
            Minimum number of genes per GO term
        max_genes : int
            Maximum number of genes per GO term
        evidence_codes : List[str], optional
            Evidence codes to include (if None, excludes only IEA)
            
        Returns:
        --------
        Dict[str, List[str]]
            Dictionary mapping GO term names to gene lists
        """
        
        # Download data if not already done
        if self.go_annotations is None:
            self.download_go_annotations()
        if self.go_terms is None:
            self.download_go_terms()
        
        # Map aspect codes
        aspect_map = {
            "biological_process": "P",
            "molecular_function": "F", 
            "cellular_component": "C"
        }
        
        aspect_codes = [aspect_map[asp] for asp in aspects if asp in aspect_map]
        
        # Filter annotations
        filtered_annotations = self.go_annotations[
            self.go_annotations['aspect'].isin(aspect_codes)
        ].copy()
        
        # Filter by evidence codes if specified
        if evidence_codes is not None:
            filtered_annotations = filtered_annotations[
                filtered_annotations['evidence_code'].isin(evidence_codes)
            ]
        else:
            # Exclude electronic annotations by default
            filtered_annotations = filtered_annotations[
                filtered_annotations['evidence_code'] != 'IEA'
            ]
        
        # Merge with GO terms to get names
        merged = filtered_annotations.merge(
            self.go_terms[['go_id', 'name', 'namespace']], 
            on='go_id', 
            how='left'
        )
        
        # Group by GO term and collect genes
        go_gene_sets = {}
        
        for go_id, group in merged.groupby('go_id'):
            genes = group['gene_symbol'].unique().tolist()
            
            # Filter by gene count
            if len(genes) < min_genes or len(genes) > max_genes:
                continue
            
            # Get term name
            term_name = group['name'].iloc[0]
            namespace = group['namespace'].iloc[0]
            
            # Create descriptive name
            full_name = f"{term_name} ({namespace})"
            go_gene_sets[full_name] = genes
        
        print(f"Created {len(go_gene_sets)} GO gene sets")
        return go_gene_sets
    
    def get_metabolism_go_terms(self) -> Dict[str, List[str]]:
        """
        Get GO terms specifically related to metabolism.
        
        Returns:
        --------
        Dict[str, List[str]]
            Dictionary of metabolism-related GO gene sets
        """
        
        # Download data if needed
        if self.go_terms is None:
            self.download_go_terms()
        
        # Metabolism-related keywords
        metabolism_keywords = [
            'metabolic', 'metabolism', 'biosynthetic', 'biosynthesis',
            'catabolic', 'catabolism', 'glycolysis', 'gluconeogenesis',
            'oxidative phosphorylation', 'citric acid cycle', 'fatty acid',
            'amino acid', 'nucleotide', 'lipid', 'carbohydrate',
            'energy', 'ATP', 'NADH', 'respiratory chain'
        ]
        
        # Find metabolism-related GO terms
        metabolism_terms = []
        for _, term in self.go_terms.iterrows():
            term_name = term['name'].lower()
            if 'definition' in term and pd.notna(term['definition']):
                term_def = term['definition'].lower()
            else:
                term_def = ""
            
            # Check if any keyword is in name or definition
            if any(keyword in term_name or keyword in term_def 
                   for keyword in metabolism_keywords):
                metabolism_terms.append(term['go_id'])
        
        # Create gene sets for metabolism terms
        if self.go_annotations is None:
            self.download_go_annotations()
        
        metabolism_gene_sets = {}
        
        for go_id in metabolism_terms:
            # Get genes for this GO term
            term_annotations = self.go_annotations[
                (self.go_annotations['go_id'] == go_id) &
                (self.go_annotations['evidence_code'] != 'IEA')
            ]
            
            if len(term_annotations) == 0:
                continue
                
            genes = term_annotations['gene_symbol'].unique().tolist()
            
            # Filter by gene count (5-500 genes)
            if len(genes) < 5 or len(genes) > 500:
                continue
            
            # Get term name
            term_info = self.go_terms[self.go_terms['go_id'] == go_id]
            if len(term_info) == 0:
                continue
                
            term_name = term_info['name'].iloc[0]
            namespace = term_info['namespace'].iloc[0]
            
            full_name = f"{term_name} ({namespace})"
            metabolism_gene_sets[full_name] = genes
        
        print(f"Found {len(metabolism_gene_sets)} metabolism-related GO terms")
        return metabolism_gene_sets
    
    def save_gene_sets_as_gmt(
        self, 
        gene_sets: Dict[str, List[str]], 
        filename: str,
        description: str = "GO gene sets"
    ):
        """
        Save gene sets in GMT format.
        
        Parameters:
        -----------
        gene_sets : Dict[str, List[str]]
            Gene sets dictionary
        filename : str
            Output filename
        description : str
            Description for the gene sets
        """
        
        with open(filename, 'w') as f:
            for term_name, genes in gene_sets.items():
                # GMT format: name, description, genes...
                line = f"{term_name}\t{description}\t" + "\t".join(genes) + "\n"
                f.write(line)
        
        print(f"Saved {len(gene_sets)} gene sets to {filename}")


def download_go_gene_sets(
    organism: str = "human",
    gene_set_type: str = "metabolism",
    cache_dir: str = None
) -> Dict[str, List[str]]:
    """
    Convenience function to download GO gene sets.
    
    Parameters:
    -----------
    organism : str
        Organism name
    gene_set_type : str
        Type of gene sets ("metabolism", "all", "biological_process", etc.)
    cache_dir : str, optional
        Cache directory
        
    Returns:
    --------
    Dict[str, List[str]]
        GO gene sets
    """
    
    go_analysis = GOAnalysis(organism=organism, cache_dir=cache_dir)
    
    if gene_set_type == "metabolism":
        return go_analysis.get_metabolism_go_terms()
    elif gene_set_type == "all":
        return go_analysis.create_go_gene_sets()
    elif gene_set_type == "biological_process":
        return go_analysis.create_go_gene_sets(aspects=["biological_process"])
    elif gene_set_type == "molecular_function":
        return go_analysis.create_go_gene_sets(aspects=["molecular_function"])
    elif gene_set_type == "cellular_component":
        return go_analysis.create_go_gene_sets(aspects=["cellular_component"])
    else:
        raise ValueError(f"Unknown gene_set_type: {gene_set_type}")