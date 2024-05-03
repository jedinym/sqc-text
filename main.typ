#set page(margin: 1.75in)
#set par(leading: 0.55em, justify: true)
#set text(font: "New Computer Modern")
#show raw: set text(font: "New Computer Modern Mono")
#show heading: set block(above: 1.4em, below: 1em, )
#set heading(numbering: "1.1")



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
models @pdb; see more in @pdb-section) is validated using community-developed
tools @pdb-validation[p.  1917]. Based on the results of these tools, the
validation pipeline generates a report @pdb-validation that can be used for
further refining the coordinate models.

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

#pagebreak()

= Biomacromolecules

// alt title: Computer representation of structural data
= Macromolecular Structural Data
The composition of macromolecules is captured in computer-readable form (and
human-readable to some extent) using a chemical file format. These exist in the
form of text files, which describe the locations of atoms in three-dimensional
space. Metadata about the represented structure may also be included. In this
chapter, the two formats relevant to this thesis are introduced.

== PDB Format
The Protein Data Bank format is the original format used by the PDB. It was
originally developed in 1976 as a simple human-readable format@pdb-history. 

Each line of the file contains a _record_ - information about some aspect of the
structure. The records can contain metadata (e.g. `AUTHOR`, `HEADER` or `TITLE`)
or data about the chemical strucutre of the molecule (e.g. `SEQRES` or `ATOM`).
Additionally, the wwPDB (worldwide Protein Data Bank) has used the `REMARK`
record type to extend the format to support new details about the experimental
methods used to obtain the macromolecular data @pdb-format-guides.

// TODO: explain what is a residue in PDB
// TODO: maybe mention that large structures can be split into several files
Unfortunately, the lines of the file are fixed-width as the format is based on
the original 80 column punched card @pdb-history. Because of this, limitations
exist on the stored structure:

- Maximum 100000 atoms in structure
- Maximum 10000 residues in one chain
- Maximum 62 chains in structure

This renders the format less suitable for handling very large structures.

Some attempts have been made to improve these limitations over the years (e.g.
the hybrid-36 counting system for atom and residue numbers @hybrid-36), but none
of them have been particularly prevalent, as it would be difficult to adapt
existing tools.
 
The PDB format has been deprecated by the Protein Data Bank in favor of the
PDBx/mmCIF format in 2014 @pdb-formats.

== PDBx/mmCIF Format
The PDBx/mmCIF format was developed to address the limitations inherent in
the legacy PDB format @mmcif[p. 571].

// TODO: what do here?

= Tools and methods

== Protein Data Bank <pdb-section>
The Protein Data Bank (PDB) is the single global archive of three-dimensional
macromolecular structures. Established in 1971, its purpose is to serve as a
central repository for macromolecular data, ensuring their accessibility to all
@pdb-history.

// TODO: add references
Since 2003, it is managed by the wwPDB consortium @wwpdb, consisting of: 
- Research Collaboratory for Structural Bioinformatics (RCSB)
- Macromolecular Structure Database (MSD) at the European Bioinformatics Institute (EBI) 
- Protein Data Bank Japan (PDBj) at the Institute for Protein Research in Osaka University
- Electron Microscopy Data Bank (EMDB)

As of now (May 2024), it stores over two hundred and nineteen thousand
structures @pdb-entry-stats. Eighty-four percent of this data was obtained using
X-ray crystallography, nine percent using electron microscopy, and around six
percent by nuclear magnetic resonance @pdb-stats-summary.

// TODO: bridge to need for deposition validation

=== OneDep deposition system
Describe what the pipeline is made of and its outputs.
Describe how the standalone pipeline is used.

== MolProbity <molprobity>

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
  
#bibliography("bibliography.yaml") 