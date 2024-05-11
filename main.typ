#import "@preview/unify:0.5.0": num
#import "@preview/big-todo:0.2.0": *

#set page(margin: 1.75in)
#set par(leading: 0.55em, justify: true)
#set text(font: "New Computer Modern")
#show raw: set text(font: "New Computer Modern Mono")
#show heading: set block(above: 1.4em, below: 1em, )
#set heading(numbering: "1.1")

// TODO: remove
#todo_outline
#pagebreak()

#align(center)[
  #text(size: 15pt)[
    Masaryk University
    #linebreak()
    Faculty of Informatics

    #linebreak()

    // color logo
    #image("img/fi-mu-logo-color.png", width: 40%)

    // BW logo
    // #image("img/FI-MU-logo-75-mm-1200-DPI-RGB.png", width: 40%)

    #linebreak()
    #linebreak()
    #linebreak()

    *Bachelor's Thesis*

    #linebreak()

    #text(size: 17pt, font: "Noto Sans", hyphenate: false)[
      #set par(justify: false)
      *Replication of PDB validation server functionality in the MetaCentrum environment*
    ]

    #linebreak()

    Martin Jediný

    Spring 2024
  ]
]

#pagebreak()

#heading(numbering: none, outlined: false)[
  Declaration 
]
Hereby I declare that this paper is my original authorial work, which I have
worked out on my own. All sources, references, and literature used or excerpted
during elaboration of this work are properly cited and listed in complete
reference to the due source.  


#align(right)[
  Martin Jediný 
]

#align(bottom)[
  *Advisor:* Mgr. Vladimír Horský, Ph.D.
]

#pagebreak()

#align(bottom)[
  #heading(numbering: none, outlined: false)[
    Acknowledgements 
  ]
  Computational resources were provided by the e-INFRA CZ project (ID:90254),
  supported by the Ministry of Education, Youth and Sports of the Czech Republic.
]

#pagebreak()

#outline(indent: auto)

#pagebreak()

#set page(numbering: "1")
#counter(page).update(1)

#heading(numbering: none)[
  Introduction
]

The quality of macromolecular structures has been the focus of biomolecular
chemists for as long as the methods determining the structure of these molecules
have been used.

The data from these experimental methods (X-ray crystallography, nuclear
magnetic resonance, and electron microscopy) are supplied with atomic
coordinates - electronic records containing the relative positions of atoms in
the structure. Because of the experimental nature of the methods, these
structural models can contain a wide range of errors. Therefore, determining the
quality of a structure is paramount in creating inferences from experiments. 

This reality is reflected by the fact that every structure deposited into the
_Protein Data Bank_ (the single global archive of three-dimensional structural
models @pdb; see more in @section-pdb) is validated using community-developed
tools @pdb-validation[p.  1917]. Based on the results of these tools, the
validation pipeline generates a report @pdb-validation that can be used for
further refining of the coordinate models.

However, the throughput of the validation pipeline provided by the _Protein Data
Bank_ is too low for use in some research projects (e.g., iterative validation
of continuously optimized structures or batch validation of up to hundreds of
millions of predicted simpler structures).

As a solution, in this thesis, I implement a scalable service which incorporates
the tools used by the Protein Data Bank validation pipeline and a Python library
for simple access. Thanks to the implemented queueing system, it is possible to
run batch validations of thousands of structures. The service is deployed via
_Ansible_ to the _Kubernetes_ cluster provided by the _MetaCentrum_ virtual
organization.

#todo[Explain sections of thesis.]

#pagebreak()

= Biomacromolecules
The IUPAC #footnote[International Union of Pure and Applied Chemistry] defines a
biopolymer or a biomacromolecule as a macromolecule #footnote[Large molecules 
consisting of many repeating subunits.] produced by living organisms
@iupac-glossary. These include proteins, nucleic acids and polysaccharrides
@iupac-glossary. In this chapter, we briefly introduce these three basic
biomacromolecules.

== Proteins
Proteins are polypeptides with a molecular mass of around 10,000 or more
@iupac-glossary-95[p. 1361]. They comprise one or more chains of $alpha$-amino
acids #footnote[There are over 500 different amino acids, but only 22 are
incorporated into proteins @amino-acids.] linked by peptide bonds #footnote[Covalent
bonds from the carbonyl carbon of one amino acid to the nitrogen atom of another
with loss of water.] @iupac-glossary-95[p. 1356].

#figure(
  image("img/AminoacidCondensation.svg"),
  caption: "The dehydration condensation of two amino acids to form a peptide bond (in red). Sourced from Wikimedia Commons."
) <peptide-bond>

Proteins perform a vast array of functions in organisms: catalyzing reactions,
providing structure to cells, replicating DNA, transporting molecules, and more
@mol-cell-bio[p. 59]. The sequence of amino acids determines the protein's
three-dimensional structure, which then dictates its function @mol-cell-bio[p. 60].

== Nucleic acids
Nucleic acids are polymers comprised of monomers (subunits) known as nucleotides
@mol-cell-bio[p. 40]. They are categorized into two classes: deoxyribonucleic
acid (DNA) and ribonucleic acid (RNA) @mol-cell-bio[p. 102].  All nucleotides
have a common structure: a phosphate group linked to a pentose #footnote[A
five-carbon sugar.] which in turn is linked to a _base_. The common bases include
_adenine_, _guanine_, and _cytosine_. _Thymine_ is exclusively found in DNA
while _uracil_ is exclusive to RNA molecules. The pentose is deoxyribose in DNA
and ribose in RNA.

#figure(
  image("img/DAMP_chemical_structure.svg", width: 80%), 
  caption: "Deoxyadenosine monophosphate, a nucleotide present in DNA. The phosphate group (blue) is linked to deoxyribose, which is in turn linked to an adenine base (red). Original image sourced from Wikimedia Commons and edited.", 
)

DNA molecules contain the information that dictates the sequences and, 
consequently, the structures of all proteins within a cell @mol-cell-bio[p.
101]. During protein synthesis, DNA is transcribed into ribonucleic acid (RNA).

#todo[Should I explain protein synthesis in short and types of RNA?]

== Polysacharrides

#todo[Fill in chapter. What is important?]

= Macromolecular Structural Data
Macromolecules can be represented in computers in one, two or three dimensions
@chemo-informatics.  In this chapter we first introduce these representations,
and then take a closer look at the three-dimensional representations, as they
are the most relevant to this thesis.

== One-dimensional structure
The simplest way of representing a molecule is by indicating the absolute counts
of atoms present (i.e. molecular formulae). To illustrate, the molecular formula
of _glucose_ is written as $C_6 H_12 O_6$, that is six carbon atoms, twelve
hydrogen atoms, and six oxygen atoms. 

However, this representation lacks information about the relative positions of
atoms and bonds, making it unable to differentiate molecules with identical
absolute atom counts.  For instance, molecules such as _fructose_ and
_galactose_ share the same formula as _glucose_, despite having different
structures.

Therefore, the one-dimensional representation has limited usage for polymers
containing thousands of atoms.

== Two-dimensional structure
A common way to represent structures in a two-dimensional manner is to use a
_molecular graph_ @chemo-informatics[p. 2]. A graph is a pair $G = (V, E)$,
where $V$ is a set of _vertices_ and $E$ is a set of pairs $E = {(v_1, v_2) |
v_1,v_2 in V}$, whose elements are called _edges_. Using a graph, it is possible
to capture the topology of a molecule. In a molecular graph, vertices represent
atoms and edges represent bonds between them @chemo-informatics[p.2]. Both
vertices and edges can hold additional information, such as bond orders for
edges or atomic numbers for vertices @chemo-informatics[p.2].

#todo[Add examples of molecular graphs?]

These molecular graphs can be encoded using various formats @chemical-formats,
with _line notation_ being one of the simpler methods. A line notation
represents a structure as a linear string of characters, making them relatively
simple to read and understand. Multiple standards are utilized, such as SLN
#footnote[Sybil Line Notation], WLN #footnote[Wiswesser Line Notation], and
ROSDAL #footnote[Representation of structure diagram arranged linearly]. 

One commonly used notation is SMILES #footnote[Simplified Molecular Input Line
Entry Specification] @smiles. Specifically, the OpenSMILES standard is widely adopted
and open-source @open-smiles.

SMILES enables the representation of molecular graphs using ASCII strings with
only a few rules. It works by doing a depth-first traversal of the graph and
printing appropriate symbols for vertices (atoms) and edges (bonds). Initially,
the graph is edited to remove hydrogen atoms and break cycles to create a
spanning tree (i.e. a subgraph that is a tree #footnote[A graph in which any two
vertices are connected by exactly one edge.] and contains all the original
vertices). @smiles-table illustrates a few examples of the SMILES format.

#figure(
  table(
    columns: (auto, auto, auto),
    align: center + horizon,
    table.header([*Molecule*], [*Chemical structure*], [*SMILES string*]),

    "Methane",
    $C H_4$,
    raw("C"),

    "Dinitrogen",
    $N≡N$,
    raw("N#N"),

    "Adenine",
    image("img/Adenine.svg", width: 60pt),
    raw("Nc1c2ncNc2ncn1"),

    "Glucose",
    image("img/Beta-D-Glucose.svg", width: 100pt),
    raw("OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)1"),

    "Nicotine",
    image("img/nicotine.svg", width: 60pt),
    raw("CN1CCC[C@H]1c2cccnc2"),
  ),
  caption: "Examples of SMILES representations of chemical structures."
) <smiles-table>

== Three-dimensional structure
Three-dimensional structures are captured in computer-readable form (and
human-readable to some extent) using a chemical file format. These exist in the
form of text files, which describe the locations of atoms in three-dimensional
space. Metadata about the represented structure may also be included. In this
chapter, the two formats relevant to this thesis are introduced.

=== PDB Format 
The Protein Data Bank format is the original format used by the Protein Data
Bank @section-pdb. It was originally developed in 1976 as a simple
human-readable format@pdb-history. 

Each line of the file contains a _record_ - information about some aspect of the
structure. The records can contain metadata (e.g. `AUTHOR`, `HEADER` or `TITLE`)
or data about the chemical strucutre of the molecule (e.g. `SEQRES` or `ATOM`).
Additionally, the wwPDB (worldwide Protein Data Bank) has used the `REMARK`
record type to extend the format to support new details about the experimental
methods used to obtain the macromolecular data @pdb-format-guides.

Unfortunately, the lines of the file are fixed-width as the format is based on
the original 80 column punched card @pdb-history. Because of this, limitations
exist on the stored structure:

- Maximum $num("100000")$ atoms in structure
- Maximum $num("10000")$ residues in one chain
- Maximum $num("62")$ chains in structure

One way to address these limitations is by dividing the structure into multiple
files. However, doing so may complicate tasks such as visualization or certain
validations of the structure. As a result, this makes the format less suitable
for handling very large structures.

Some attempts have been made to improve these limitations over the years (e.g.
the hybrid-36 counting system for atom and residue numbers @hybrid-36), but none
of them have been particularly prevalent, as it would be difficult to adapt
existing tools.
 
The PDB format has been deprecated by the Protein Data Bank in favor of the
PDBx/mmCIF format in 2014 @pdb-formats.

=== PDBx/mmCIF Format
The PDBx/mmCIF format was developed to address the limitations inherent in the
legacy PDB format @mmcif[p. 571]. With ever-increasing sizes of structures it
became clear that change was needed @mmcif-ecosystem[p. 3].

The format is based on the Crystallographic Information Framework (CIF), which
was adopted by the International Union of Crystallography (IUCr) in 1990. CIF
operates under the idea that all values within an ASCII text file are assigned
dictionary labels (keys). This concept was enabled by the use of a Dictionary
Definition Language (DDL), which is a versatile language allowing for the
description of dictionary data structures, controlled vocabularies, boundary  
conditions, and relationships between values.

Later, in 1997, the mmCIF dictionary was approved by the international Committee
for the Maintenance of the CIF Standard (COMCIFS) @mmcif-approval. It featured
expanded data types, which included support for protein and nucleic acid polymer
types, polymer chains, ligands, binding sites, macromolecular assemblies, amino
acid and nucleotide residues, atomic coordinates, and experimental data
@mmcif-ecosystem[p. 3].

While mmCIF already supported most of the structural and crystallographic
concepts present in the PDB format, additional key categories prefixed with
`pdbx_` were introduced to the dictionary and some existing categories have been
extended. This expansion aimed to guarantee complete compatibility and semantic
equivalence with the PDB format @crystallographic-data[p. 195].

In 2014, the PDBx/mmCIF format became the main format of the PDB.

= Tools and methods

== Protein Data Bank <section-pdb>
The Protein Data Bank (PDB) is the single global archive of three-dimensional
macromolecular structures. Established in 1971, its purpose is to serve as a
central repository for macromolecular data, ensuring their accessibility to all
@pdb-history.

#todo[add references]
Since 2003, it is managed by the wwPDB consortium @wwpdb, consisting of: 
- Research Collaboratory for Structural Bioinformatics (RCSB)
- Macromolecular Structure Database (MSD) at the European Bioinformatics Institute (EBI) 
- Protein Data Bank Japan (PDBj) at the Institute for Protein Research in Osaka University
- Electron Microscopy Data Bank (EMDB)

As of now (May 2024), it stores over two hundred and nineteen thousand
structures @pdb-entry-stats. Eighty-four percent of this data was obtained using
X-ray crystallography, nine percent using electron microscopy, and around six
percent by nuclear magnetic resonance @pdb-stats-summary.

As the number of large structures in the PDB continued to grow, the existing
infrastructure proved inadequate. That is why the wwPDB initiated the
development of OneDep, a novel system designed for deposition, biocuration, and
validation @onedep[p. 536]. 

During the deposition process, the deposited structure is validated using
community-made tools in OneDep's _validation pipeline_ @onedep[p. 539]. 
Additionally, a standalone server #footnote[https://validate.wwpdb.org] is
accessible to depositors, providing them with a platform for refining the
structure.

#todo[Mention why standalone server is not enough?]

=== OneDep validation pipeline
To implement the relevant validation tools, the wwPDB convened so-called
Validation Task Forces (VTFs) for each experimental method, which compiled
recommendations for the validation pipeline @pdb-validation[p. 1917].  Based on
these recommendations, the pipeline was integrated with the OneDep system
@pdb-validation[p. 1917].

// TODO: Think about explaining the categories a little bit.
The VTFs recommended that deposited structures undergo validation against
criteria falling into three categories. The first category encompasses criteria
for validating the resulting atomic model. The second category involves analysis
of the supplied experimental data. Lastly, the third category entails examining
the fit between the supplied experimental data and the atomic model
@pdb-validation[p. 1917].

The outcome of the validation pipeline is a validation report, which
incorporates the parsed outputs of the validation tools employed on the
structure. This report is presented in both a human-readable PDF format and a
machine-readable XML format, complemented by an XSD schema
#footnote[https://www.wwpdb.org/validation/schema/] @pdb-validation[p. 1917].

The pipeline is composed of component software tools that validate certain
aspects of the supplied data @pdb-validation[p. 1922]. Some of these tools are
publicly available, but access to some was only provided to the wwPDB. One of
the publicly available tools that was used in the implementation of this thesis
is _MolProbity_ @molprobity, which we will explore further in
@section-molprobity.

There are three ways to access the validation pipeline: an anonymous web user
interface, as part of the deposition and biocuration process, and as a web API
@pdb-validation[p. 1921]. For the purposes of validations of large numbers of
structures, the web API is most convenient, as it can be fully automated. It
offers both a CLI #footnote[Command Line Interface] application and a small
Python library.

== MolProbity <section-molprobity>

=== Validations
Describe the validations that molprobity validates.

== Ansible

== Kubernetes

= Design

== Requirements
Functional and non-functional requirements.

== Architecture
Describe how and why the architecture was chosen.

= Implementation

== Validation Service
Describe dockerification. MolProbity output parsing.

== Validation Client Library

= Deployment
Describe deployment architecture (in metacentrum k8s environment).

== Automation
Describe how to deploy SQC to a k8s namespace. 
Describe ansible setup (vaults, inventories).
Don't forget kubectl configuration.
Don't forget docker image build and push.

= Evaluation
Compare aspects of the service to the PDB standalone validation server.

== Throughput
Test throughput of variable-size structures against PDB with maximum k8s scaling.

== API 
Compare API access to PDB methods.

== Validation Results
Cry about different results and inaccessible reference data.

= Conclusion

== Future plans 
Describe stuff that's missing or can be improved.

#pagebreak()
  
#bibliography("bibliography.yaml", style: "ieee") 