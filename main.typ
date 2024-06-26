#import "@preview/unify:0.5.0": num, unit
#import "@preview/big-todo:0.2.0": *
#import "@preview/codly:0.2.1": *

#set page(margin: 1.75in)
#set par(leading: 0.55em, justify: true)
#set text(font: "New Computer Modern")
#show raw: set text(font: "Noto Sans Mono")
#show heading: set block(above: 1.4em, below: 1em, )
#set heading(numbering: "1.1")
#show figure: set block(breakable: true)

#show: codly-init.with()

#let icon(codepoint) = {
  box(
    height: 0.8em,
    baseline: 0.05em,
    image(codepoint)
  )
  h(0.1em)
}

#codly(
  languages: (
    python: (name: "Python", icon: icon("img/brand-python.svg"), color: rgb("#3380ff"))
  ),
  enable-numbers: true,
)

#show "SQC": emph 
#show "SQCLib": emph 
#show "Kubernetes": emph 
#show "RabbitMQ": emph
#show "MinIO": emph
#show "MetaCentrum": emph

#align(center)[
  #text(size: 15pt)[
    Masaryk University
    #linebreak()
    Faculty of Informatics

    #linebreak()

    // color logo
    #image("img/fi-mu-logo-color.png", width: 40%)

    #linebreak()
    #linebreak()
    #linebreak()

    *Bachelor's Thesis*

    #linebreak()

    #text(size: 17pt, font: "Noto Sans", hyphenate: false)[
      #set par(justify: false)
      *Replication of PDB validation server functionality in the _MetaCentrum_ environment*
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
Hereby, I declare that this thesis is my original authorial work, which I have
worked out on my own. All sources, references, and literature used or excerpted
during the elaboration of this work are properly cited and listed in complete
reference to the due source. 

During the preparation of this thesis, I used the following AI tools:
- Grammarly for checking grammar,
- ChatGPT to improve my writing style.
I declare that I used these tools in accordance with the principles of academic
integrity. I checked the content and take full responsibility for it.

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
  
  I would like to thank my advisor, Vladimír for his huge support and time
  investment. I also thank my dearest friends for showing me what's possible.
  Lastly, I must thank my family for making sure I don't starve.

  Computational resources were provided by the e-INFRA CZ project (ID:90254),
  supported by the Ministry of Education, Youth and Sports of the Czech Republic.
]

#pagebreak()

#heading(numbering: none, outlined: false)[
  Abstract
]

Assessing the quality of biomacromolecular structural models (models of proteins, nucleic acids, and polysaccharides) acquired using experimental methods is a vital step in making inferences in structural biology.  That is also why the validation of structures is an obligatory step in the Protein Data Bank's (PDB) deposition process. To allow structural biologists to refine their structures, the PDB provides a standalone validation service.  However, the throughput of this service is too low for research projects that need to validate millions of structures. In this thesis, we implement Structure Quality Control (SQC), a scalable, containerized, and easily deployable service, that incorporates community-made validation tools. The implemented tool allows for easy-to-use structure validations via both a Python and a shell API. Deployed in the MetaCentrum virtual organization, the service is now ready for use in research.

// Kontrola kvality biomakromolekulárnych štruktúr (modely bielkovín, nukleových kyselín a polysacharidov) získaných použitím experimentálnych metód je dôležitým krokom pri vyvodzovaní záverov v štrukturálnej biológii. Je to tiež jedným z dôvodov prečo validácia štruktúr je nutnou súčasťou depozičného procesu _Protein Data Bank_ (PDB). PDB poskytuje aj validačnú službu, ktorá je využívaná na iteratívnu kontrolu štruktúr. Priepustnosť tejto služby je však príliš nízka na využitie v niektorých projektoch. V tejto práci implementujeme _Structure Quality Control_ (SQC), škáľovateľnú, kontajnerizovanú a jednoducho nasaditeľnú službu, ktorá používa validačné nástroje vyvinuté komunitou štrukturálnych biológov. Implementovaný nástroj sprístupňuje validácie štruktúr pomocou Python rozhrania alebo rozhrania systémovej príkazovej riadky. Finálne riešenie je nasadené pomocou virtuálnej organizácie MetaCentrum a je pripravené na použitie vo výskume.

#align(bottom)[
  #heading(numbering: none, outlined: false)[
    Keywords
  ]

  Protein Data Bank, PDB, biomacromolecules, MolProbity, proteins, nucleic
  acids, _Kubernetes_, bioinformatics, validation, _MetaCentrum_
]


#pagebreak()

#outline(indent: auto)

#pagebreak()

#set page(numbering: "1")
#counter(page).update(1)

#heading(numbering: none)[
  Introduction
]

Understanding the structure of biomacromolecules is fundamental to structural
biology research. The data acquired using experimental methods are used to
assemble structures, i.e., computer models of large biomolecules. These models
are crucial in understanding the molecular basis of biological processes,
leading to advances in drug discovery or disease treatment.

Because of the experimental nature of the methods used, these structural models
can contain a wide range of errors.  Therefore, validating biomacromolecular
structural data is a critical step in ensuring their accuracy and reliability in
biological research. 

This reality is reflected by the fact that every structure deposited into the
_Protein Data Bank_ (the single global archive of three-dimensional structural
models @pdb) is validated using community-developed tools @pdb-validation[p.
1917]. Based on the results of these tools, the validation pipeline generates a
report @pdb-validation that can be used by structure authors for further
refining of the coordinate models, as well as by users when searching for the
best structure possible for their research.

However, the throughput of the validation pipeline provided by the _Protein Data
Bank_ is too low for use in some research projects (e.g., iterative validation
of continuously optimized structures or batch validation of up to hundreds of
millions of simpler structures predicted using neural networks @alphafoldDB).

As a solution, a new service will be implemented that incorporates the tools
used by the Protein Data Bank validation pipeline. Additionally, the service
will include a queueing system, enabling batch validations of large numbers of
structures. The service will utilize the computing resources of the
_MetaCentrum_ virtual organization.

#pagebreak()

= Biomacromolecules
The IUPAC #footnote[International Union of Pure and Applied Chemistry] defines a
biopolymer or a biomacromolecule as a macromolecule #footnote[Large molecules 
consisting of many repeating subunits (monomers).] produced by living organisms
@iupac-glossary. These include proteins, nucleic acids, and polysaccharides
@iupac-glossary. In this section, we briefly introduce these three basic
biomacromolecules.

== Proteins
Proteins are chains of amino acids with a molecular mass of around 10,000 or more
@iupac-glossary-95[p. 1361]. They comprise one or more chains of amino acids
#footnote[There are over 500 different amino acids, but only 22 are incorporated
into proteins @amino-acids.] linked by peptide bonds $dash.en$ covalent bonds
from the carbonyl carbon of one amino acid to the nitrogen atom of another with
loss of water @iupac-glossary-95[p. 1356]. An example of such a reaction can be
seen in @figure-peptide-bond.

#figure(
  image("img/AminoacidCondensation.svg"),
  caption: "The dehydration condensation of two amino acids to form a peptide bond (in red). Sourced from Wikimedia Commons."
) <figure-peptide-bond>

Proteins perform a vast array of functions in organisms: catalyzing reactions,
providing structure to cells, replicating DNA, transporting molecules, and more
@mol-cell-bio[p. 59]. The sequence of amino acids determines the protein's
three-dimensional structure, which then dictates its function @mol-cell-bio[p. 60].

== Nucleic acids
Nucleic acids are polymers comprised of monomers (subunits) known as nucleotides
@mol-cell-bio[p. 40]. They are categorized into two classes: deoxyribonucleic
acid (DNA) and ribonucleic acid (RNA) @mol-cell-bio[p. 102].  All nucleotides
have a common structure: a phosphate group linked to a pentose, #footnote[A
five-carbon sugar.] which in turn is linked to a _base_. The common bases
include _adenine_, _guanine_, and _cytosine_. _Thymine_ is exclusively found in
DNA, while _uracil_ is exclusive to RNA molecules. The pentose is deoxyribose in
DNA and ribose in RNA.

#figure(
  image("img/DAMP_chemical_structure.svg", width: 80%), 
  caption: "Deoxyadenosine monophosphate, a nucleotide present in DNA. The phosphate group (blue) is linked to deoxyribose (black), which is in turn linked to an adenine base (red). Original image sourced from Wikimedia Commons and edited.", 
)

DNA molecules contain the information that dictates the sequences and, 
consequently, the structures of all proteins within a cell @mol-cell-bio[p.
101]. During protein synthesis, DNA is transcribed into ribonucleic acid (RNA),
which is then translated into a protein @mol-cell-bio[p. 101].

== Polysacharides
Monosaccharides, or simple sugars, are the monomer units of polysaccharides
@iupac-glossary-95[p. 1360]. Polysaccharides are formed by linking
monosaccharides together through glycosidic linkages.

Some serve a storage function in a cell, preserving sugars for later use (e.g.,
starch in plants, glycogen in animals). Others function as building material for
structures of cells @biology[p. 71].

#pagebreak()
= Macromolecular structural data
Macromolecules can be represented in computers in one, two or three dimensions
@chemo-informatics.  In this section, we first introduce these representations,
and then take a closer look at the three-dimensional representations as they
are the most relevant to this thesis.

== One-dimensional structure
The simplest way of representing a molecule is by indicating the total number 
of atoms present (using a molecular formula). To illustrate, the molecular formula
of _glucose_ is written as $C_6 H_12 O_6$, that is, six carbon atoms, twelve
hydrogen atoms, and six oxygen atoms. 

However, this representation lacks information about the relative positions of
atoms and bonds, making it unable to differentiate molecules with identical
total atom counts.  For instance, molecules such as _fructose_ and
_galactose_ share the same formula as _glucose_ despite having different
structures.

Therefore, the one-dimensional representation has limited usage for polymers
containing thousands of atoms.

== Two-dimensional structure
A common way to represent structures in a two-dimensional manner is to use a
_molecular graph_ @chemo-informatics[p. 2]. A graph is a pair $G = (V, E)$,
where $V$ is a set of _vertices_, and $E$ is a set of pairs $E = {(v_1, v_2) |
v_1,v_2 in V}$, elements called _edges_. Using a graph, it is
possible to capture the topology of a molecule. In a molecular graph, vertices
represent atoms, and edges represent bonds between them @chemo-informatics[p.2].
Both vertices and edges can hold additional information, such as bond orders for
edges or atomic numbers for vertices @chemo-informatics[p.2].

These molecular graphs can be encoded using various formats @chemical-formats,
with _line notation_ one of the simpler methods. A line notation represents a
structure as a linear string of characters, making them relatively simple to
read and understand.

One commonly used notation is SMILES #footnote[Simplified Molecular Input Line
Entry Specification] @smiles. Specifically, the OpenSMILES standard is widely
adopted and open-source @open-smiles.

SMILES enables the representation of molecular graphs using ASCII strings with
only a few rules. It works by doing a depth-first traversal of the graph and
printing appropriate symbols for vertices (atoms) and edges (bonds). Initially,
the graph is edited to remove hydrogen atoms and break cycles to create a
spanning tree (i.e. a subgraph that is a tree #footnote[A graph in which any two
vertices are connected by exactly one edge.] and contains all the original
vertices) @smiles-algorithm. @smiles-table illustrates a few examples of the
SMILES format.

#figure(
  table(
    columns: (auto, auto, auto),
    align: center + horizon,
    table.header([*Molecule*], [*Chemical formula*], [*SMILES string*]),

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
  caption: "Examples of SMILES representations of chemical molecules."
) <smiles-table>

== Three-dimensional structure
Three-dimensional structures are represented in computer-readable form (and
human-readable to some extent) using a chemical file format. These exist in the
form of text files, which describe the locations of atoms in three-dimensional
space. Metadata about the represented structure may also be included. In this
section, the two formats relevant to this thesis are introduced.

=== PDB format
The Protein Data Bank format is the original format used by the Protein Data
Bank. It was originally developed in 1976 as a simple
human-readable format @pdb-history. 

Each line of the file contains a _record_ - information about some aspect of the
structure. The records can contain metadata (e.g., `AUTHOR`, `HEADER` or `TITLE`)
or data about the molecule's chemical structure (e.g. `SEQRES` or `ATOM`).
Additionally, the `REMARK` record type has been used to extend the format to
support new details about the experimental methods used to obtain the
macromolecular data @pdb-format-guides.

Unfortunately, the lines of the file are fixed-width as the format is based on
the original 80-column punched card @pdb-history. Because of this, limitations
exist on the stored structure:

- Maximum $num("100000")$ atoms in structure
- Maximum $num("10000")$ residues in one chain
- Maximum $num("62")$ chains in structure

One way to address these limitations is by dividing the structure into multiple
files. However, doing so may complicate tasks such as visualization or certain
structure validations. As a result, this makes the format less suitable
for handling very large structures.

Some attempts have been made to improve these limitations over the years (e.g.,
the hybrid-36 counting system for atom and residue numbers @hybrid-36), but none
of them have been particularly prevalent, as it would be difficult to adapt
existing tools.
 
The PDB format has been deprecated by the _Protein Data Bank_ in favor of the
PDBx/mmCIF format in 2014 @pdb-formats. However, the PDB still offers structures
in the PDB format; these are only best-effort conversions of the PDBx/mmCIF
format.

=== PDBx/mmCIF Format
The PDBx/mmCIF format was developed to address the limitations inherent in the
legacy PDB format @mmcif[p. 571]. With ever-increasing sizes of structures, it
became clear that change was needed @mmcif-ecosystem[p. 3].

The format is based on the Crystallographic Information Framework (CIF), which
was adopted by the International Union of Crystallography (IUCr) in 1990. CIF
operates under the idea that all values within an ASCII text file are assigned
dictionary labels (keys). This concept was enabled using a Dictionary
Definition Language (DDL), which is a versatile language allowing for the
description of dictionary data structures, controlled vocabularies, boundary  
conditions, and relationships between values @mmcif-ecosystem[p. 2].

Later, in 1997, the mmCIF dictionary was approved by the international Committee
for the Maintenance of the CIF Standard (COMCIFS) @mmcif-approval. It featured
expanded data types, which included support for protein and nucleic acid polymer
types, polymer chains, ligands, binding sites, macromolecular assemblies, amino
acid and nucleotide residues, atomic coordinates, and experimental data
@mmcif-ecosystem[p. 3].

While mmCIF already supported most of the structural and crystallographic
concepts present in the PDB format, additional key categories prefixed with
`pdbx_` were introduced to the dictionary, and some existing categories have been
extended. This expansion aimed to guarantee complete compatibility and semantic
equivalence with the PDB format @crystallographic-data[p. 195].

In 2014, the PDBx/mmCIF format became the primary format of the PDB
@mmcif-ecosystem[p. 3].

#pagebreak()
= Tools and methods
In this section, we introduce the tools employed in the implementation. First,
we discuss the _Protein Data Bank_ and its validation pipeline. Next, we examine
MolProbity, a validation tool used in the implementation. Finally, we explore
the components of the new system.

== Protein Data Bank <section-pdb>
The Protein Data Bank (PDB) is the single global archive of three-dimensional
macromolecular structures. Established in 1971, its purpose is to serve as a
central repository for macromolecular data, ensuring their accessibility to all
@pdb-history.

Since 2003, it is managed by the wwPDB consortium @wwpdb, consisting of: 
- Research Collaboratory for Structural Bioinformatics (RCSB) @pdb[p. D520]
- Macromolecular Structure Database (MSD) at the European Bioinformatics
  Institute (EBI) @pdb[p. D520]
- Protein Data Bank Japan (PDBj) at the Institute for Protein Research in Osaka
  University @pdb[p. D520]
- Electron Microscopy Data Bank (EMDB) @emdb

Currently (May 2024), it stores over two hundred and nineteen thousand
structures @pdb-entry-stats. Eighty-four percent of this data was obtained using
X-ray crystallography, nine percent using electron microscopy, and around six
percent by nuclear magnetic resonance @pdb-stats-summary.

As the number of large and complex structures in the PDB continued to grow, the
existing deposition and validation infrastructure proved inadequate. That is why
the wwPDB initiated the development of OneDep, a novel system designed for
deposition, biocuration, and validation @onedep[p. 536]. 

During deposition, the deposited structure is validated using
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
for validating the resulting atomic model, i.e., is only the three-dimensional
computer representation. The second category involves analysis of the supplied
experimental data. Lastly, the third category examines the fit between
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
is MolProbity @molprobity, which we will explore further in
@section-molprobity.

There are three ways to access the validation pipeline: via an anonymous web
user interface, as part of the deposition and biocuration process, and through a
web API @pdb-validation[p. 1921]. The web API is most convenient for validating
large numbers of structures, as it can be fully automated. It offers both a CLI
#footnote[Command Line Interface] application and a small Python library.

== MolProbity <section-molprobity>
MolProbity is the single validation tool used in this thesis. It provides
evaluation of atomic model quality for biomacromolecules @molprobity[p. 12]. The
current maintainers are the Richardson Laboratory of Duke University. Its source
code can be found in their GitHub repository
#footnote[https://github.com/rlabduke/MolProbity]. The repository also contains
a short guide on how to install MolProbity, however, the guide and
installation scripts are rather outdated and require some refactoring for use on
a modern system.

MolProbity was chosen because, unlike other tools in the PDB validation
pipeline, it is freely available and provides multiple validations in a single
package. The main focus of the validations is the geometry of the structure.
Validations such as too-close contacts, bond lenghts and angles and torsion
angles.

The software package contains a web interface for simple use but also a command
line interface for bulk validations, which is more useful for automated use.

// === Validations
// In this section, we briefly introduce the various validations that MolProbity
// performs on atomic models.
// 
// ==== All-atom contact analysis <section-clashes>
// This validation option checks for overlaps of van der Waals surfaces of
// nonbonded atoms. Overlaps over $0.4 angstrom$ #footnote[$1 angstrom = 0.1
// unit("nano meter")$] are reported as _clashes_ @molprobity[p. 14].
// 
// ==== Covalent-geometry analyses
// MolProbity additionally assesses outliers in bond lengths #footnote[The
// average distance between nuclei of two bonded atoms.] and bond angles
// #footnote[The geometric angle between two adjacent bonds.] of backbones
// #footnote[The main chain of the polymer.] @molprobity[p. 15]. 
// 
// Based on ideal parameters derived ahead of time @molprobity[p. 15]
// @struct-quality-data @backbone-conformation, the instances where the bond
// lengths or bond angles deviate from the ideal value by at least four standard
// deviations are identified and listed as bond length and angle outliers,
// respectively @molprobity[p. 15].
// 
// ==== Ramachandran and rotamer analyses
// #todo[Not really understanding these right now.]

=== Interface
Multiple command-line programs can be used to run _MolProbity's_ validation
suite on structures. They can be found in the `/cmdline` directory in the source
tree. The two programs used in the implementation of this thesis are
_residue-analysis_ and _clashscore_. 

Firstly, _residue-analysis_ runs all available validations and outputs outliers
in each residue in the CSV #footnote[Comma Separated Values] format.

Secondly, the _clashscore_ program checks for too-close contacts and outputs
detected clashes in a one-clash-per-line format.

== _Kubernetes_
Kubernetes is an open-source platform designed to automate the deployment,
scaling, and operation of containerized applications. It organizes these
applications into clusters, which are sets of nodes (machines) that run
containerized applications. A Kubernetes cluster typically includes a control
plane, which manages the overall state and lifecycle of the applications, and
worker nodes, which run the containers.

Central to Kubernetes is its API, which serves as the primary interface for
interaction with the cluster. Users and automated systems communicate with the
API to deploy, manage, and monitor applications. The API accepts requests in a
declarative manner, meaning users specify the desired state of the system, and
Kubernetes takes the necessary actions to achieve and maintain that state.

Kubernetes objects are persistent entities within the system that represent
the state and configuration of various cluster components. Key objects for this
thesis include:

+ *Pods*: The smallest deployable unit that can be created, managed, and
  run. It encapsulates one or more containers that share the same network
  namespace and storage. Pods are designed to host a single instance of
  a running process in a cluster, making them the basic building blocks of a
  Kubernetes application.

+ *Replicas*: Refer to multiple instances of a Pod running concurrently to
  ensure high availability and reliability. Kubernetes uses ReplicaSet objects,
  often defined within Deployments, to manage the desired number of replicas. This
  ensures that the specified number of pod replicas are always running, allowing
  for automatic scaling, self-healing, and load balancing across the cluster.

+ *Deployments*: Manage the deployment and scaling of application replicas. They
  ensure that the desired number of pod replicas are running and can update pods
  in a controlled manner.

+ *Persistent Volume Claims (PVCs)*: Abstract the storage resources that a pod can
  use. A PVC requests storage resources from the cluster, which are then
  provisioned by Persistent Volumes (PVs).

+ *Ingresses*: Manage external access to the services in a cluster, typically HTTP.
  An Ingress controller watches the Ingress resources and manages the routing of
  external traffic to the appropriate services.

+ *Services*: Define a logical set of pods and a policy by which to access them.
  Services enable communication within the cluster and provide stable IP addresses
  and DNS names for pods, meaning pods can communicate easily.

+ *Secrets*: Store sensitive information such as passwords, tokens, and keys.
  Secrets allow sensitive data to be used in pods without exposing it directly in
  pod specifications.

These objects are defined in YAML #footnote[YAML Ain't Markup Language] or JSON
#footnote[JavaScript Object Notation] format and submitted to the Kubernetes
API, which manages their lifecycle.

== Ansible
_Ansible_ is an open-source automation tool widely used for IT orchestration,
configuration management, and application deployment @ansible. It operates
through playbooks, which are scripts written in the YAML format that define a
series of tasks to be executed on managed nodes. Each task specifies an action,
such as installing software or configuring network settings, and is applied to
hosts listed in an inventory. Inventories categorize hosts into groups, allowing
for targeted management and execution of tasks across different subsets of
infrastructure.

A fundamental aspect of _Ansible_ is its focus on idempotency. This ensures that
tasks, when executed, bring the system to a desired state without causing
additional changes if the system is already in that state. This property is
crucial for maintaining consistency and preventing unintended side effects
during repeated executions.

_Ansible_ also supports the integration of custom modules to extend its
capabilities. These custom modules allow for the automation of specialized
tasks. In this thesis, we utilize the `kubernetes.core.k8s` module
#footnote[https://docs.ansible.com/ansible/latest/collections/kubernetes/core/k8s_module.html],
for deployment automation.

== _MetaCentrum_
The MetaCentrum virtual organization provides computing resources to all
students and employees of academic institutions in the Czech Republic @metacentrum.
Membership is free, with the requirement that members acknowledge MetaCentrum
in their publications.

It offers many different platforms for computation across the Czech Republic,
but crucially for this thesis, it also offers a Kubernetes cluster @metacentrum-k8s
via a _Rancher_ #footnote[https://www.rancher.com/] instance.

== _RabbitMQ_ <section-rabbitmq>
RabbitMQ is a messaging and streaming broker that supports several standard
messaging protocols @rabbitmq. It is used as a mediator between producers and
consumers of messages.

Publishers (producers) publish a message to an exchange. The message is then
routed to a queue. If the queue has any active consumers, the message is
delivered to them. If no consumers are active, the message is cached on disk and
delivered at the next opportunity.

== _MinIO_ Object Store <section-minio>
// TODO: Why MinIO and not a traditional database?
// To store atomic structures and validation reports, a storage solution is required. 

MinIO is an object store inspired by and compatible with Amazon's S3 service
#footnote[https://aws.amazon.com/s3/] @minio. MinIO is simple to deploy, as
it's intended for use in the cloud, including using Kubernetes.

MinIO offers high-performance storage of data by storing them as objects in
_buckets_. A bucket is simply a container for objects, where each object has a
key that uniquely identifies it in a _bucket_. Each object can also contain
additional text metadata.

Access to MinIO is mediated via an HTTP API designed initially for Amazon's
S3 service. However, MinIO also offers software development kits (SDKs) for
multiple programming languages, including Python. Another part of the MinIO
deployment is the managment server $dash.en$ a web application used for
configuration and management.

Crucially, for this thesis, MinIO can monitor events associated with a
_bucket_ and publish notifications via multiple protocols
#footnote[https://min.io/docs/minio/kubernetes/upstream/administration/monitoring.html].
One of the supported protocols is _AMQP 0-9-1_, which is the protocol used by
RabbitMQ. A few examples of such events are: 
- `s3:ObjectCreated:Put` which occurs when a new object is added to a bucket.  
- `s3:ObjectRemoved:Delete` which occurs when an object is deleted from a
  bucket.  
- `s3:ObjectCreated:CompleteMultipartUpload` which occurs when an object larger
  than 5MiB is added to a bucket. 

// == Miscellaneous
// In this section, we introduce the tools and libraries that are not the core of
// the solution.
// 

#pagebreak()
= Implementation
The result of this thesis is a validation service called _Structure Quality
Control_ (SQC) alongside a corresponding Python library, SQCLib. The source code
is available in the _sb-ncbr_ (Structural bioinformatics research group at
National Centre for Biomolecular Research) organization on GitHub for both SQC
#footnote[https://github.com/sb-ncbr/sqc] and SQCLib
#footnote[https://github.com/sb-ncbr/sqclib]. 

This section begins with an introduction to _SQC's_ high-level architecture and
design decisions, followed by a closer examination of the implementation.

== Requirements
In this section, we outline the functional and non-functional requirements of
the SQC system.

=== Functional requirements
+ The service must support running validations on atomic models.
+ The service must support both the legacy PDB and PDBx/mmCIF formats as inputs.
+ The service must provide a machine-readable output containing results of validations.
+ The service must provide a schema for the output.
+ The service must support many concurrent validations using a queueing system.
+ The service must be easily deployable to the MetaCentrum Kubernetes cluster.
+ The service must be easily scalable to accommodate higher loads.
+ The service must provide a web-accessible API.
+ The service must provide a Python library for accessing the API.

=== Non-functional requirements
+ The service must be implemented using the Python programming language.
+ The service must be containerized using Docker.
+ The service must be deployed to the MetaCentrum Kubernetes cluster.
+ The service must use Ansible to automate deployment to the MetaCentrum cluster.

== Architecture
Since supporting many concurrent validations is crucial, we needed to devise a
solution that could handle them effectively. To guarantee that even during high
load, the system can accept new validation requests, we use a Producer Consumer
system design pattern.

=== Producer Consumer pattern
In the Producer Consumer system design pattern, producers generate data or
events and put them in a shared queue. Consumers then retrieve the data or
events from the queue and process them (@figure-producer-consumer). The main
benefit of using the pattern is the decoupling of producers and consumers.

While a consumer is occupied processing prior data or events, producers can
continue adding more data to the queue. This additional data awaits consumption
and processing once at least one consumer completes all its tasks. Adding data
to the queue is faster than processing, ensuring producers encounter no
delays due to consumers.

During periods of increased load, when the queue accumulates more data and
consumers are slow in processing it; new consumers accessing the shared queue
can be created.

#figure(
  caption: "Diagram of the Producer Consumer pattern. Producers produce new data in parallel and consumers consume and process data in parallel.",
  image("img/producer_consumer.svg"),
) <figure-producer-consumer>

=== Producer Consumer pattern in SQC
To implement the Producer Consumer pattern in SQC, we utilize RabbitMQ
(@section-rabbitmq) and MinIO (@section-minio).

When a user submits a validation request through SQCLib (refer to
@section-sqclib), the library code creates a new object in a MinIO bucket. This
object contains the atomic model along with associated metadata. The library
subsequently returns the ID of the pending request to the user. The user uses
this ID to invoke another SQCLib function, which awaits the completion of the
request.

Upon the creation of the object, MinIO publishes a bucket notification event to
the `requests` exchange in RabbitMQ. This notification contains the identifier
of the new request.

The SQC server listens to the `requests` exchange. Upon receiving a notification
from MinIO, it downloads the associated atomic model of the request to temporary
storage and runs validations on it.

Once the atomic model is validated, SQC uploads the results to another MinIO
bucket. When the result is submitted, the library returns the results of the
validation to the user.

In our scenario, users submitting validation requests act as producers,
_RabbitMQ's_ exchange acts as a shared queue, and the SQC server acts as a
consumer. During periods of high load, additional instances of the SQC server
can be created, accessing the same shared queue. This approach boosts the
system's throughput.

A sequence diagram showing a user's interaction with the system when validating
a single structure can be seen in @figure-sequence.

#figure(
  caption: "UML Sequence Diagram of an SQC validation request.", 
  image("img/sequence.png"),
) <figure-sequence>

== Validation service
In this section, we introduce the validation service and its key components.
First, we discuss the service's source code. Next, we examine containerization.
Then, we explore how the MolProbity suite is utilized. Finally, we inspect the
service's output format.

=== Source code overview
In this section, we explore _SQC's_ source code repository and provide an overview
of its elements.

The `ansible/` directory contains _Ansible_ playbooks and inventories, used to
deploy SQC to MetaCentrum. @section-automation goes further into the
deployment setup.

The `compose.yaml` file contains the configuration file for _Docker Compose_
#footnote[https://docs.docker.com/compose/]. Using it, developers can start a
local instance of SQC for simpler development and testing. _Docker Compose_
starts the RabbitMQ, MinIO, and SQC containers, and exposes the default MinIO
port `9000`. The local instance can then be accessed for testing using the
SQCLib library.

The `Dockerfile` contains instructions for _Docker_ to build _SQC's_ image.
See more in @section-containerization.

The `pyproject.toml` contains definitions pertaining to the SQC Python project.
It defines the required Python version to run the project and the project's
production and development dependencies. To manage dependencies, we utilize
_PDM_ #footnote[https://pdm-project.org/en/latest/], a modern python dependency
manager. _PDM_ simplifies various development tasks, including dependency
resolution and updates, virtual environment support, and management of user
scripts.

The `pdm.lock` file is created by _PDM_ and contains exact versions of packages
required to run SQC. The advantage of using a lockfile is that during the Docker
image-building process, dependencies are not resolved again, but extracted from
the lockfile, which guarantees the reproducibility of all builds.

The `sqc/` directory contains the source code of the SQC validation service.

=== Containerization <section-containerization>
With each validation tool potentially having its own dependencies, a key
challenge in implementing the validation service is ensuring its simple setup.

We chose to simplify the setup by employing containerization, specifically
Docker. Containerizing the validation service offers additional advantages,
including straightforward deployment, isolation from the host system, and
consistent results, as all the dependencies of the tools (both software
dependencies and reference data) can be controlled.

The _Dockerfile_ utilizes a multi-stage build. The first stage installs the
necessary packages, fetches _MolProbity's_ reference data and source code from
GitHub and builds _MolProbity's_ binaries. In the second stage, SQC's
dependencies are installed. This method was selected to reduce the Docker
image's build time.  Building MolProbity can be time-consuming, taking up to
an hour on a recent laptop. Docker can cache unchanged parts of the image
between builds, improving build times.

Since _MolProbity's_ results are dependent on reference data, the repositories
containing it were forked into the _sb-ncbr_ organization on GitHub. These are
then cloned during the image build.

== Authentication
In order to restrict access to the service, we utilize MinIO's built-in
authentication mechanism. To access resources in a MinIO bucket, the user must
provide an access key and a secret key. These two values are generated by the
administrator using the MinIO management service and then provided to a new
user.

=== Running validations
The core of the validation service is utilizing the MolProbity suite to
validate atomic models and parsing its outputs into a validation report.

The SQC server listens for messages in the `requests` exchange in RabbitMQ. When
a request arrives, the respective atomic model is fetched from MinIO and
stored in the `/tmp` directory.

Since MolProbity exclusively supports validating models stored in the legacy
PDB format, it's necessary to initially convert the PDBx/mmCIF atomic models to
this legacy format. This is done using the _Gemmi_
#footnote[https://gemmi.readthedocs.io/en/latest/] library for structural
biology. _Gemmi_ has multiple features useful in macromolecular biology but
supports format conversions using the `gemmi convert` subprogram. After the
structure is fetched from MinIO, _Gemmi_ converts it to the legacy format.

Another constraint of MolProbity is that it can only operate on PDB files
containing only a single model. It is possible for PDB files to contain an
ensemble of models, something that is typical for structures determined by
nuclear magnetic resonance. We utilize the _BioPython_
#footnote[https://biopython.org/] library to split the original PDB file into
multiple single-model files. Once the files are prepared, MolProbity's programs
are executed on each one. The results of the validations are then returned on a
per-model basis.

After the input files are preprocessed, _MolProbity's_ _residue-analysis_
program is executed with the path to the structure as the first argument. The
process is started and managed using the built-in _subprocess_ library. The
output of the process is captured and then parsed using the _csv_ built-in
library into the format of the validation report.

The second program to be executed is _clashscore_. The output of _clashscore_ is
not in the CSV format and therefore requires quite intricate parsing. The
output is parsed into the validation report as well.

Before converting the validation report to JSON and sending it to MinIO, we
add the versions of the reference data used for validation of the model. These
versions are represented as URLs of the git repositories and commit hashes of
the latest commits.

=== Validation reports
We decided to use the JSON format for the validation reports, as it's a simple
human-readable format that is often used in the industry. To ensure that the
reports are convenient to use, we utilize _JSON Schema_
#footnote[https://json-schema.org/], a vocabulary used to annotate JSON
documents. By providing a _JSON Schema_, we allow users to easily parse the
validation reports.

To facilitate the generation of _JSON Schema_, we use the _Pydantic_
#footnote[https://docs.pydantic.dev/latest/] Python library. _Pydantic_ is
primarily used for validation of external data, but can also be
used to emit _JSON Schema_.

With _Pydantic_ a _model_ of the output is created as a Python object. From this
model, it is possible to generate a _JSON Schema_. Once populated with results
from MolProbity, the model is converted into JSON and returned to the user in
that format.

This approach offers the advantage that as the model expands in the future, the
schema is automatically adjusted accordingly. @figure-pydantic shows an example
Python object with its associated _JSON Schema_.

#show figure: set block(breakable: true)
#figure(
  grid(
    rows: 2,
    gutter: 2em,
    ```python
    import pydantic

    class SubObject(pydantic.BaseModel):
        integer: int
        opt_string: str | None

    class ExampleObject(pydantic.BaseModel):
        fields: list[str]
        sub: SubObject
    ```,
    ```JSON
    {
      "$defs": {
        "SubObject": {
          "properties": {
            "integer": {
              "title": "Integer",
              "type": "integer"
            },
            "opt_string": {
              "anyOf": [
                {
                  "type": "string"
                },
                {
                  "type": "null"
                }
              ],
              "title": "Opt String"
            }
          },
          "required": [
            "integer",
            "opt_string"
          ],
          "title": "SubObject",
          "type": "object"
        }
      },
      "properties": {
        "fields": {
          "items": {
            "type": "string"
          },
          "title": "Fields",
          "type": "array"
        },
        "sub": {
          "$ref": "#/$defs/SubObject"
        }
      },
      "required": [
        "fields",
        "sub"
      ],
      "title": "ExampleObject",
      "type": "object"
    }
    ```,
  ),
  caption: figure.caption(position: bottom)[An example of a small Pydantic model and its generated JSON Schema. The schema describes the types of every field of the model.] 
) <figure-pydantic>
#show figure: set block(breakable: false)

== Validation client library <section-sqclib>
In order to send a validation request to the SQC server, it is necessary to use
the SQCLib Python package. This package includes a library API for use in Python
code and a small command-line application for uploading structures from the
system shell. The library is easily installable and distributed as a Python
package.

The implementation wraps the Python API provided by MinIO
#footnote[https://min.io/docs/minio/linux/developers/python/API.html], while
providing a few additional validations.

=== Python API
The main part of the library consists of a simple Python API. Users utilize the
`SQCClient` object to interact with SQC instances. To create the object, users
must pass their credentials to the client, which will then use them to
communicate with SQC. Optionally, the users can specify a different URL of the
SQC server, which can be useful when running SQC locally or having a separately
deployed instance. The client offers the following methods:

+ The `submit(path)` method submits a structure to an SQC instance and returns
  an ID of the pending request. The structure is specified via a filesystem
  path. This method can be useful when validating a lot of structures, as they
  can be all uploaded first and all awaited later.

+ The `get_result(request_id, timeout_s)` method is used in combination with the
  `submit(path)` method. It takes an ID of a pending request as an argument and
  blocks until the SQC server returns a validation report.

+ The `validate(path, timeout_s)` is a combination of the `submit` and
  `get_result` methods. It submits a request to the SQC server and blocks until a
  report is returned. This method offers the simplest API and is useful when
  validating a single structure.

An example script utilizing the Python API can be seen in @figure-sqclib-python.

#figure(
  ```python
  from os import environ
  from sqclib import SQCClient

  client = SQCClient(
    access_key=environ.get("SQC_ACCESS_KEY"),
    secret_key=environ.get("SQC_SECRET_KEY"),
  )

  requests = []
  for path in ["183d.pdb", "2pde.mmcif"]:
    request_id = client.submit(path)
    requests.append((path, request_id))

  reports = []
  for path, request_id in requests:
    report = client.get_result(request_id)
    reports.append((path, report))
  ```,
  caption: "An example Python script showing SQCLib usage. Once the client is initialized, two structures are submitted to the SQC instance and their reports awaited."
) <figure-sqclib-python>

=== Shell API
The library also contains the program `sqc` that can be used to submit
structures for validation directly from the system shell. To use SQC, the
`SQC_ACCESS_KEY` and `SQC_SECRET_KEY` environment variables must be set with
before-provided credentials. @figure-sqclib-shell shows how to upload a
structure to the main instance of SQC in the _bash_ shell.

#figure(
  ```sh
  export SQC_ACCESS_KEY=access
  export SQC_SECRET_KEY=secret

  sqc structure.pdb
  sqc structure.mmcif
  ```,
  caption: "An example of submitting a structure for validation using the sqc program. Credentials must be provided by the SQC administrators."
) <figure-sqclib-shell>

By default, the `sqc` program submits validations to the public SQC instance. To
use a locally running instance of SQC, it is possible to use the `--url` and
`--insecure` parameters. @figure-sqclib-local shows how to utilize these
parameters.

#figure(
  ```sh
  # "minioadmin" is the default key of the local instance
  export SQC_ACCESS_KEY=minioadmin
  export SQC_SECRET_KEY=minioadmin

  # following two calls are equivalent
  sqc --insecure --url localhost:9000 structure.mmcif
  sqc -k -u localhost:9000 structure.mmcif
  ```,
  caption: "An example of submitting a structure for validation to a local instance of SQC. The local instance uses default \"minioadmin\" credentials."
) <figure-sqclib-local>

=== Documentation
Since the SQCLib library is the sole entry point to the SQC system, thorough
documentation is essential. To generate documentation from the code, we use
_Sphinx_, a documentation generator designed for Python that supports various
output formats #footnote[https://www.sphinx-doc.org/en/master/].

To automatically generate documentation from Python docstrings, we use the
_autodoc_ extension to _Sphinx_. _Autodoc_ parses Google-style docstrings
#footnote[https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings]
and generates a reference for a Python module.

Furthermore, it is important that the documentation is easily accessible. For
this end, we use _GitHub Pages_ #footnote[https://pages.github.com/], a system
that allows free hosting of websites directly from a GitHub repository.

To automate the process, we utilize _GitHub Actions_
#footnote[https://github.com/features/actions], a CI/CD solution from GitHub,
that allows automation to run on GitHub repository events. When a new commit
is added to the `main` branch in the SQCLib GitHub repository, the documentation
is built and uploaded to _GitHub Pages_. This process guarantees that the
accessible documentation is always up to date. The documentation is then
available in the `sb-ncbr` subdomain of the github.io domain
#footnote[https://sb-ncbr.github.io/sqclib/].

#pagebreak()
= Deployment
To deploy the main instance of SQC, we chose the Rancher Kubernetes
distribution provided by MetaCentrum
#footnote[https://rancher.cloud.e-infra.cz]. Since we had already decided to
utilize containerization and needed simple scaling, Kubernetes was the obvious
choice. 

In this section we first explore the array of objects used to deploy SQC
to Kubernetes, and then we describe how the process of deployment is automated
using _Ansible_.

== _Kubernetes_ project setup
To deploy SQC to a Kubernetes project, we take advantage of various
Kubernetes objects.

The credentials for the MinIO administrator account and RabbitMQ user are
stored using a Secret. These secrets are then copied to an environment variable
when creating respective Deployments.

The MinIO object store requires storage for objects, so we create a Persistent
Volume Claim (PVC) that is later mounted into a MinIO Deployment. This
guarantees that the data stored by MinIO are persistent even across Pod
restarts. Similarly, RabbitMQ requires some storage for buffering of messages
in queues. Again, we create a PVC and later mount it into the RabbitMQ
Deployment.

To deploy MinIO, RabbitMQ, and SQC, we create three Deployments. The SQC
Deployment is unique in the fact that the number replicas in the ReplicaSet is
configurable. By specifying a higher number of replicas, the throughput of the
service will be increased.

Exposing MinIO, RabbitMQ, and SQC is achieved by utilizing a Kubernetes
Service. Each of the three Deployments is exposed by a Service, allowing Pods to
communicate over the network.

Lastly, we create two Ingresses, exposing the main MinIO port and the MinIO
management port. We utilize special annotations provided by MetaCentrum
#footnote[https://docs.cerit.io/docs/kubectl-expose.html], which create a
subdomain under the `dyn.cloud.e-infra.cz` domain and automatically enable TLS
to use HTTPS.

== Automation <section-automation>
To automate the deployment of SQC, we use two _Ansible_ playbooks:
`deploy-sqc.yaml`, which creates all Kubernetes objects required for a minimal
setup of SQC and `teardown-sqc.yaml`, which deletes all objects pertaining to
SQC.

Both playbooks utilize an _Ansible_ vault, which contains secrets needed to
deploy SQC, like the token to _MetaCentrum's_ Kubernetes API. When running the
playbooks, a vault password is entered, which then decrypts the vault and makes
the secrets available to _Ansible_ as inventory variables.

The inventory only contains one host: `metacentrum`, because it is the only
deployment target at the moment. However, new deployment targets can be easily added
in the future. Another Kubernetes cluster can be added with minimal effort.

The inventory for the host `metacentrum`, contains some important values for the
deployment:
+ `k8s_sqc_namespace` specifies the Kubernetes project that SQC will be
  deployed to.
+ `k8s_sqc_replicas` specifies the base number of replicas of SQC pods.
+ `k8s_sqc_image` specifies the repository, name and tag of the SQC Docker image.
+ `k8s_sqc_request_cpu` and `k8s_sqc_request_mem` specify the requested memory
  for the SQC container.

When deploying a new version of SQC, it is necessary to rebuild the image and
push it to a docker repository, so Kubernetes can then pull it during the
deployment.

One important aspect of the deployment is fitting SQC's replicas, CPU and memory
requests to the Kubernetes project quota. The default project quota in the
MetaCentrum cluster is 20 CPUs request, 32 CPUs limit, 40GB memory requests,
and 64GB memory limits. Adjusting the requested CPUs, memory, and the number of
replicas is important to get the highest performance. If validating many small
structures, it is faster to have more replicas with fewer requested CPUs. When
validating larger structures, it might be more beneficial to have fewer replicas
with more requested CPU and memory. @figure-resources shows the default
configuration of resources.

#[
  #set par(justify: false)

  #figure(
    table(
      columns: (auto, auto, auto, auto, auto),
      align: center + horizon,
      table.header([*Deployment*], [*Replicas*], [*CPU Requests*], [*Memory Requests*], [*Memory Limits*]),

      "SQC",
      "20",
      "0.9",
      [$1.5$ GiB], 
      [$3$ GiB],

      "RabbitMQ",
      "1",
      "1",
      [$512$ MiB],
      [$512$ MiB],

      "MinIO",
      "1",
      "1",
      [$512$ MiB],
      [$512$ MiB],
     ),
    caption: "Default resource requests and limits of the SQC Deployments."
  ) <figure-resources>
]

#pagebreak()
= Evaluation
In this section, we address two aspects of the newly implemented solution.
First, we examine SQC's scaling performance, and second, we review the
validation results.

== Scaling
To test SQC's scaling capabilities, we validated thirty structures with an
increasing number of replicas. We chose thirty structures because it is the
highest number of replicas that can be provisioned in the Kubernetes project 
while fitting into quota memory limits. @figure-scaling shows that with the
increasing number of replicas, the time spent on multiple structure validations
drops rapidly. Scaling replicas further is possible; the limiting factor is the
resource quota in the Kubernetes cluster.

#figure(
  image("img/replica-perf.png"),
  caption: "Time spent on validating 30 183d structures (PDB id 183d) plotted against the number of SQC validation service replicas."
) <figure-scaling>

== Validation results
Certain validations performed by MolProbity (namely bond lengths, bond angles,
and torsion angles) are heavily influenced by reference data. Unfortunately, the
datasets used by PDB are not publicly available. After an email inquiry, PDB
referred us to the MolProbity paper from 2017 @molprobity-more-data. According
to this paper, MolProbity in the PDB pipeline should be using the data provided
by the Richardson Laboratory at Duke University. However, even after using their
datasets from GitHub, some validations produced slightly different results.
Without access to the specific files used by PDB, it is difficult to determine
the cause of these discrepancies.

#pagebreak()
= Conclusion
In this thesis, we aimed to reimplement parts of the PDB validation pipeline to
improve the throughput of structure validation and allow local deployment. We
developed a solution that can be simply accessed using a Python library or a
commandline program. According to performed tests, the implementation scales
well with increasing computational resources. Additionally, the service defines
an output schema, facilitating its use. 

The implementation leverages Kubernetes for deployment and scaling, and Ansible
for deployment automation. Other cloud-native technologies were also used in the
implementation, such as RabbitMQ and the MinIO object store.

The final solution is ready for upcoming projects in the Structural
bioinformatics research group at the National Centre for Biomolecular Research
such as validations of hundreds of millions of predicted structures from
_AlphaFoldDB_ @alphafoldDB or fast iterative validations during the process of structure
optimization by means of computational chemistry.

== Future plans
Even though the implemented solution offers practical validation of structures,
there are several ares that could be improved in the future. Notably, the
service only utilizes one validation tool. Adding more validation tools has the
potential to make the tool the one-stop shop for high throughput validations,
even though they would have to be implemented from scratch, as they are not
available.  Additionally, allowing users to specify versions of reference data
to use would significantly enhance the replicability of validations. 

#pagebreak()
  
#bibliography("bibliography.yaml", style: "ieee") 

#pagebreak()

= Appendix
The attached zip file contains the source code of _SQC_ and _SQCLib_ in two
directories.  Please refer to the `README.md` files in both directories for
instructions on running and using the SQC validation service.