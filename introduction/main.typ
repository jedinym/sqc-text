#set page(margin: 1.75in, numbering: "1")
#set par(leading: 0.55em, justify: true)
#set text(font: "New Computer Modern")
#show raw: set text(font: "New Computer Modern Mono")
#show heading: set block(above: 1.4em, below: 1em, )


#align(center, text(17pt)[
  *Replication of PDB validation server functionality in the MetaCentrum*
  #linebreak()
  *environment*
])


#linebreak()
#linebreak()
#linebreak()

#align(center)[
  #text(size: 17pt)[
    Bachelor's Thesis 
  ]
]

#linebreak()
#linebreak()

#align(center)[
  #text(size: 17pt)[
    Martin Jediný
  ]
]

#align(bottom + center)[
  Brno, Spring 2024
]

#pagebreak()

= Introduction

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
models @pdb) is validated using community-developed tools @pdb-validation[p.
1917]. Based on the results of these tools, the validation pipeline generates a
report @pdb-validation that can be used for further refining the coordinate
models.

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
#bibliography("bibliography.yaml", style: "ieee")