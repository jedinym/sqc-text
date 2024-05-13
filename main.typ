#import "@preview/unify:0.5.0": num, unit
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

As a solution, in this thesis, I implement _SQC_ (Structure Quality Control), a
scalable service which incorporates the tools used by the Protein Data Bank
validation pipeline and a Python library for simple access. Thanks to the
implemented queueing system, it is possible to run batch validations of
thousands of structures. The service is deployed via _Ansible_ to the
_Kubernetes_ cluster provided by the _MetaCentrum_ virtual organization.

#todo[Explain sections of thesis.]

#pagebreak()

= Biomacromolecules
The IUPAC #footnote[International Union of Pure and Applied Chemistry] defines a
biopolymer or a biomacromolecule as a macromolecule #footnote[Large molecules 
consisting of many repeating subunits (monomers).] produced by living organisms
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

== Polysacharrides
Monosaccharides, or simple sugars, are the monomer units of polysaccharides
@iupac-glossary-95[p. 1360]. Polysaccharides are formed by linking
monosaccharides together through glycosidic linkages.

Some serve a storage function in a cell, preserving sugars for later use (e.g.
starch in plants, glycogen in animals). Others function as building material for
structures of cells @biology[p. 71].

#todo[Should I mention something more?]

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

#todo[Should I add examples of molecular graphs?]

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

=== OneDep validation pipeline
To implement the relevant validation tools, the wwPDB convened so-called
Validation Task Forces (VTFs) for each experimental method, which compiled
recommendations for the validation pipeline @pdb-validation[p. 1917].  Based on
these recommendations, the pipeline was integrated with the OneDep system
@pdb-validation[p. 1917].

The VTFs recommended that deposited structures undergo validation against
criteria falling into three categories. The first category encompasses criteria
for validating the resulting atomic model; that is only the three-dimensional
computer representation. The second category involves analysis of the supplied
experimental data. Lastly, the third category entails examining the fit between
the supplied experimental data and the atomic model @pdb-validation[p. 1917].

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
// TODO: explain why it was chosen
_MolProbity_ is the single validation tool used in this thesis. It provides
evaluation of atomic model quality for biomacromolecules @molprobity[p. 12]. The
current maintainers are Richardson Laboratory of Duke University. Its source
code can be found in their GitHub repository
#footnote[https://github.com/rlabduke/MolProbity].

The software package contains a web interface for simple use but also a command
line interface for bulk validations.

=== Validations
In this section, we briefly introduce the various validations that _MolProbity_
performs on atomic models.

==== All-atom contact analysis <section-clashes>
This validation option checks for overlaps of van der Waals surfaces of
nonbonded atoms. Overlaps over $0.4 angstrom$ #footnote[$1 angstrom = 0.1
unit("nano meter")$] are reported as _clashes_ @molprobity[p. 14].

==== Covalent-geometry analyses
_MolProbity_ additionally assesses outliers in bond lengths #footnote[The
average distance between nuclei of two bonded atoms.] and bond angles
#footnote[The geometric angle between two adjacent bonds.] of backbones
#footnote[The main chain of the polymer.] @molprobity[p. 15]. 

Based on derived parameters @molprobity[p. 15], the instances where the bond
lengths or bond angles deviate from the ideal value by at least four standard
deviations are identified and listed as bond length and bond angle outliers,
respectively @molprobity[p.  15].

==== Ramachandran and rotamer analyses
#todo[Not really understanding these right now.]

=== Interface
Multiple command-line programs can be used to run _MolProbity's_ validation
suite on structures. They can be found in the `/cmdline` directory in the source
tree. The two programs used in the implementation of this thesis are
_residue-analysis_ and _clashscore_. 

Firstly, _residue-analysis_ runs all available validations and outputs outliers
in each residue in the CSV #footnote[Comma Separated Values] format.

Secondly, the _clashscore_ program checks for clashes (@section-clashes) and
outputs detected clashes in a one clash per line format.

#todo[Work on following subchapters when implementation chapters are done to see what needs to be mentioned.]

== Kubernetes
// Kubernetes, originally designed by Google, is an open-source system for
// deployment automation. It enables simple execution of containerized workloads in
// any compatible cloud environment.

// What is a Cluster
// How does the k8s API work (objects)
// 

== Ansible
// TODO: Ansible is an open-source automation tool primarily maintained by Red Hat.

== MetaCentrum

== RabbitMQ <section-rabbitmq>

== Minio <section-minio>

= Requirements

#todo[Think about these some more.]

== Functional requirements
+ The service must support running validations on atomic models.
+ The service must support both the legacy PDB and PDBx/mmCIF formats as inputs.
+ The service must provide a machine-readable output containing results of validations.
+ The service must provide a schema for the output.
+ The service must support many concurrent validations using a queueing system.
+ The service must be easily deployable to the MetaCentrum Kubernetes cluster.
+ The service must be easily scalable to accommodate higher loads.
+ The service must provide a web-accessible API.
+ The service must provide a Python library for accessing the API.

== Non-functional requirements
+ The service must be implemented using the Python programming language.
+ The service must be containerized using Docker.
+ The service must be deployed to the MetaCentrum Kubernetes cluster.
+ The service must use Ansible to automate deployment to the MetaCentrum cluster.

= Implementation
The result of this thesis is a service named SQC (Structure Quality Control)
alongside a corresponding Python library, SQClib. This chapter begins with an
introduction to SQC's high-level architecture and design decisions, followed by
a closer examination of the implementation.

== Architecture
Since supporting many concurrent validations is crucial, we needed to devise a
solution that can handle them effectively. To guarantee that even during high
load, the system can accept new validation requests, we use a Producer Consumer
system design pattern.

=== Producer Consumer pattern
#todo[Does this chapter belong to tools and methods?]

In the Producer Consumer system design pattern, producers generate data or
events and put them in a shared queue. Consumers then retrieve the data or
events from the queue and process them (@figure-producer-consumer). The main
benefit of using the pattern is the decoupling of producers and consumers.

While a consumer is occupied processing prior data or events, producers can
continue adding more data to the queue. This additional data awaits consumption
and processing once at least one consumer completes all its tasks. Adding data
to the queue is fast compared to processing, ensuring producers encounter no
delays due to consumers.

During periods of increased load, when the queue accumulates more data and
consumers are slow in processing it, new consumers accessing the shared queue
can be created.

#figure(
  caption: "Diagram of the producer consumer pattern. Producers produce new data in parallel and consumers consume and process data in parallel.",
  image("img/producer_consumer.svg"),
) <figure-producer-consumer>

=== Producer Consumer pattern in SQC
To implement the Producer Consumer pattern in SQC, we utilize RabbitMQ
(@section-rabbitmq) and MinIO (@section-minio).

When a user submits a validation request through SQClib (refer to
@section-sqclib), the library code creates a new object in a MinIO bucket. This
object contains the atomic model along with associated metadata. The library
subsequently returns the pending request's ID to the user. The user uses this ID
to invoke another SQClib function, which awaits the completion of the request.

Upon the creation of the object, Minio publishes a bucket notification event to
the `requests` exchange in RabbitMQ. This notification contains the identifier
of the new request.

The SQC server listens to the `requests` exchange. Upon receiving a notification
from MinIO, it downloads the associated atomic model of the request to temporary
storage and runs validations on it.

Once the atomic model is validated, SQC uploads the results to another MinIO
bucket. When the result is submitted, the library returns the results of the
validation to the user.

In our scenario, users submitting validation requests act as producers,
RabbitMQ's exchange acts as a shared queue, and the SQC server acts as a
consumer. During periods of high load, additional instances of the SQC server
can be created, accessing the same shared queue. This approach boosts the
system's throughput.

#figure(
  caption: "UML Sequence Diagram of an SQC validation request.", 
  image("img/sequence.png"),
)

== Validation Service
Describe dockerification. MolProbity output parsing.

== Validation Client Library <section-sqclib>

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